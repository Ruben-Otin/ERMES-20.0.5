/*****************************************************************************************************
 *
 *  ERMESPETScSolver.cpp  -  ERMES 20.0 / PETSc direct interface
 *
 *  Solves the ERMES linear system A Xo = B with PETSc, reading and writing the ERMES binary files
 *  directly (no conversion to PETSc format). Files in the problem folder:
 *
 *    Matrix_A_cmplx.bin      input   complex<double> (real, imag) coefficients a_ij
 *    Matrix_A_int.bin        input   int32 pairs (i, j), 1-based
 *    Matrix_A_aux_cmplx.bin  input   auxiliary matrix coefficients c_mn   (Hermitic formats only)
 *    Matrix_A_aux_int.bin    input   auxiliary matrix indices (m, n)      (Hermitic formats only)
 *    Vector_B.bin            input   complex<double> vector B
 *    Vector_Xo.bin           output  complex<double> solution Xo, read back by ERMES
 *
 *  Matrix storage formats, detected automatically from the files:
 *
 *    Format               Detected when                     Assembled matrix
 *    -------------------  --------------------------------  --------------------------------------
 *    Symmetric            no aux files, one triangle        A = U + (U - D)^T
 *    Full-matrix          no aux files, both triangles      A = M
 *    Hermitic-Symmetric   aux files, aux has one triangle   A = U + (U - D)^H + R + (R - D_R)^T
 *    Hermitic-Full        aux files, aux has both triangles A = U + (U - D)^H + R
 *
 *    U: stored diagonal + upper diagonal of Matrix_A (Hermitian volumetric part in the Hermitic
 *    formats), R: stored auxiliary matrix (Robin boundary part), D, D_R: their diagonals,
 *    ^T transpose, ^H conjugate transpose. Diagonals are used as stored.
 *
 *  Parallel algorithm (no process ever holds the whole matrix):
 *
 *    1) Each MPI process reads a contiguous 1/P slice of the entries of each matrix file (MPI-IO).
 *    2) The transposed / conjugate-transposed off-diagonal entries are added on the fly.
 *    3) All entries go to PETSc with MatSetPreallocationCOO / MatSetValuesCOO, which sends each
 *       entry to the process owning its row and adds repeated entries (main + aux are summed).
 *    4) Each process reads its part of B, the system is solved with KSP, and each process writes
 *       its part of Xo.
 *
 *  Options (any PETSc option, e.g. -ksp_type, -pc_type, -mat_mumps_icntl_*, can also be used):
 *
 *    -ermes_folder <dir>       Folder with the ERMES files (default: current folder).
 *                              The old option -my_folder is also accepted.
 *    -ermes_monitor_every <n>  Iterative solvers: print the residual every n iterations
 *                              (solver residual and true relative residual ||b-Ax||/||b||)
 *
 *  Requirements: PETSc >= 3.19 configured with --with-scalar-type=complex (double precision).
 *
 *****************************************************************************************************/

static char help[] = "Solves an ERMES linear system with PETSc, reading/writing ERMES binary files directly.\n"
                     "  -ermes_folder <dir>        folder with the ERMES files (default: current folder)\n"
                     "  -ermes_monitor_every <n>   print the residual every n iterations\n\n";

#include <petscksp.h>
#if PETSC_VERSION_LT(3, 19, 0)
  #define PETSC_SUCCESS ((PetscErrorCode)0)
#endif
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>

#if !defined(PETSC_USE_COMPLEX)
  #error "ERMESPETScSolver requires PETSc configured with --with-scalar-type=complex"
#endif
#if !defined(PETSC_USE_REAL_DOUBLE)
  #error "ERMESPETScSolver requires double precision PETSc"
#endif

static_assert(sizeof(PetscScalar) == 2 * sizeof(double), "Unexpected PetscScalar layout");

/* Maximum number of items per MPI-IO call (MPI counts are 'int') */
static const long long kMaxChunk = 1LL << 26;

/* Bytes per entry in the ERMES files */
static const int kCmplxBytes = 16; /* complex<double> */
static const int kIndexBytes = 8;  /* 2 x int32       */

/*====================================================================================================
   File helpers
  ====================================================================================================*/

/* Joins folder and file name */
static std::string JoinPath(const char *dir, const char *name)
{
  std::string d(dir);
  if (!d.empty() && d.back() != '/') d += '/';
  return d + name;
}

/* True on all processes if the file exists (checked on process 0) */
static PetscErrorCode FileExists(MPI_Comm comm, const std::string &path, PetscBool *exists)
{
  PetscMPIInt rank;
  int         e = 0;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  if (rank == 0) e = (access(path.c_str(), F_OK) == 0);
  PetscCallMPI(MPI_Bcast(&e, 1, MPI_INT, 0, comm));
  *exists = e ? PETSC_TRUE : PETSC_FALSE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Collective open of a file for reading; returns its size in bytes */
static PetscErrorCode OpenForRead(MPI_Comm comm, const std::string &path, MPI_File *fh, MPI_Offset *bytes)
{
  PetscFunctionBeginUser;
  int err = MPI_File_open(comm, path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, fh);
  PetscCheck(err == MPI_SUCCESS, comm, PETSC_ERR_FILE_OPEN, "Cannot open ERMES file %s", path.c_str());
  PetscCallMPI(MPI_File_get_size(*fh, bytes));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Collective read of 'count' items starting at byte 'offset'. Large reads are split in chunks, and
   all processes do the same number of calls (a process with less data reads 0 items). */
static PetscErrorCode ReadAtAll(MPI_File fh, MPI_Offset offset, void *buf, long long count, MPI_Datatype type, int typeBytes, MPI_Comm comm)
{
  long long nRoundsLoc, nRounds, done = 0;
  char     *p = static_cast<char *>(buf);

  PetscFunctionBeginUser;
  nRoundsLoc = (count + kMaxChunk - 1) / kMaxChunk;
  PetscCallMPI(MPI_Allreduce(&nRoundsLoc, &nRounds, 1, MPI_LONG_LONG, MPI_MAX, comm));
  for (long long r = 0; r < nRounds; ++r) {
    int        n   = (int)std::min(kMaxChunk, count - done);
    int        got = 0;
    MPI_Status st;
    PetscCallMPI(MPI_File_read_at_all(fh, offset + (MPI_Offset)done * typeBytes, p + done * typeBytes, n, type, &st));
    PetscCallMPI(MPI_Get_count(&st, type, &got));
    PetscCheck(got == n, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Short read from ERMES file: expected %d items, got %d", n, got);
    done += n;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Collective write, same chunking logic as ReadAtAll */
static PetscErrorCode WriteAtAll(MPI_File fh, MPI_Offset offset, const void *buf, long long count, MPI_Datatype type, int typeBytes, MPI_Comm comm)
{
  long long   nRoundsLoc, nRounds, done = 0;
  const char *p = static_cast<const char *>(buf);

  PetscFunctionBeginUser;
  nRoundsLoc = (count + kMaxChunk - 1) / kMaxChunk;
  PetscCallMPI(MPI_Allreduce(&nRoundsLoc, &nRounds, 1, MPI_LONG_LONG, MPI_MAX, comm));
  for (long long r = 0; r < nRounds; ++r) {
    int        n = (int)std::min(kMaxChunk, count - done);
    MPI_Status st;
    PetscCallMPI(MPI_File_write_at_all(fh, offset + (MPI_Offset)done * typeBytes, p + done * typeBytes, n, type, &st));
    done += n;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*====================================================================================================
   Matrix parts (main and auxiliary matrix files)
  ====================================================================================================*/

/* How the stored triangle of a part is completed */
enum Expansion { EXPAND_NONE, EXPAND_TRANSPOSE, EXPAND_CONJUGATE };

/* One matrix stored in a pair of ERMES files (indices + values) */
struct MatrixPart {
  std::string          pathAi, pathAv;          /* index and value files                     */
  long long            nnz        = 0;          /* entries in the files                      */
  long long            kStart     = 0;          /* first entry read by this process          */
  long long            nLoc       = 0;          /* entries read by this process              */
  std::vector<int32_t> ij;                      /* local (row, col) pairs, until converted   */
  long long            triLoc[2]  = {0, 0};     /* local strictly upper / lower entries      */
  long long            triGlob[2] = {0, 0};     /* global strictly upper / lower entries     */
  Expansion            expand     = EXPAND_NONE;/* how the stored triangle is completed      */
  std::vector<char>    isOff;                   /* local entry is off-diagonal               */
  long long            cooOffset  = 0;          /* first position of this part in COO arrays */

  bool      OneTriangle() const { return (triGlob[0] == 0) != (triGlob[1] == 0); }
  long long NumAdded() const { return expand == EXPAND_NONE ? 0 : triLoc[0] + triLoc[1]; }
  long long NumCoo() const { return nLoc + NumAdded(); }
};

/* Reads this process's slice of the indices of a part and counts upper / lower entries */
static PetscErrorCode ReadPartIndices(MPI_Comm comm, MatrixPart &p)
{
  PetscMPIInt rank, size;
  MPI_File    fhAi, fhAv;
  MPI_Offset  bytesAi, bytesAv;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  PetscCallMPI(MPI_Comm_size(comm, &size));

  /* Sizes: nnz index pairs and nnz complex values */
  PetscCall(OpenForRead(comm, p.pathAi, &fhAi, &bytesAi));
  PetscCall(OpenForRead(comm, p.pathAv, &fhAv, &bytesAv));
  PetscCallMPI(MPI_File_close(&fhAv));
  PetscCheck(bytesAi % kIndexBytes == 0, comm, PETSC_ERR_FILE_UNEXPECTED, "%s: size is not a multiple of 8 bytes", p.pathAi.c_str());
  p.nnz = (long long)(bytesAi / kIndexBytes);
  PetscCheck((long long)bytesAv == p.nnz * kCmplxBytes, comm, PETSC_ERR_FILE_UNEXPECTED, "Inconsistent files: %s has %lld entries, %s has %lld", p.pathAi.c_str(), p.nnz, p.pathAv.c_str(), (long long)(bytesAv / kCmplxBytes));

  /* Contiguous slice of the entries for this process */
  p.kStart = p.nnz * rank / size;
  p.nLoc   = p.nnz * (rank + 1) / size - p.kStart;

  p.ij.resize((size_t)(2 * p.nLoc));
  PetscCall(ReadAtAll(fhAi, (MPI_Offset)p.kStart * kIndexBytes, p.ij.data(), 2 * p.nLoc, MPI_INT32_T, 4, comm));
  PetscCallMPI(MPI_File_close(&fhAi));

  /* Stored triangles (used to detect the format) */
  for (long long k = 0; k < p.nLoc; ++k) {
    if (p.ij[2 * k] < p.ij[2 * k + 1]) ++p.triLoc[0];
    else if (p.ij[2 * k] > p.ij[2 * k + 1]) ++p.triLoc[1];
  }
  PetscCallMPI(MPI_Allreduce(p.triLoc, p.triGlob, 2, MPI_LONG_LONG, MPI_SUM, comm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Fills the COO index arrays of a part: stored entries first, then the added transposed entries.
   Counts the indices out of range in *nBad. */
static void FillPartIndices(MatrixPart &p, long long N, PetscInt *ci, PetscInt *cj, long long *nBad)
{
  long long m = p.nLoc;

  if (p.expand != EXPAND_NONE) p.isOff.resize((size_t)p.nLoc);
  for (long long k = 0; k < p.nLoc; ++k) {
    const long long i = (long long)p.ij[2 * k] - 1; /* ERMES indices are 1-based */
    const long long j = (long long)p.ij[2 * k + 1] - 1;
    if (i < 0 || i >= N || j < 0 || j >= N) ++(*nBad);
    ci[k] = (PetscInt)i;
    cj[k] = (PetscInt)j;
    if (p.expand != EXPAND_NONE) {
      p.isOff[k] = (char)(i != j);
      if (i != j) { /* transposed position */
        ci[m] = (PetscInt)j;
        cj[m] = (PetscInt)i;
        ++m;
      }
    }
  }
  std::vector<int32_t>().swap(p.ij); /* free */
}

/* Reads this process's slice of the values of a part, and adds the transposed (or conjugate-
   transposed) values in the same order as FillPartIndices */
static PetscErrorCode ReadPartValues(MPI_Comm comm, const MatrixPart &p, PetscScalar *v)
{
  MPI_File   fh;
  MPI_Offset bytes;

  PetscFunctionBeginUser;
  PetscCall(OpenForRead(comm, p.pathAv, &fh, &bytes));
  PetscCall(ReadAtAll(fh, (MPI_Offset)p.kStart * kCmplxBytes, v, 2 * p.nLoc, MPI_DOUBLE, 8, comm));
  PetscCallMPI(MPI_File_close(&fh));

  if (p.expand != EXPAND_NONE) {
    long long m = p.nLoc;
    for (long long k = 0; k < p.nLoc; ++k)
      if (p.isOff[k]) v[m++] = (p.expand == EXPAND_CONJUGATE) ? PetscConj(v[k]) : v[k];
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*====================================================================================================
   Matrix assembly
  ====================================================================================================*/

/* Reads the ERMES matrix files in parallel, detects the storage format and assembles the
   distributed matrix A. Returns the detected format name and the stored entries of each part. */
static PetscErrorCode AssembleMatrix(MPI_Comm comm, const char *folder, long long N, Mat *A, std::string *format, long long nnzFile[2])
{
  MatrixPart mainPart, auxPart;
  PetscBool  hasAux = PETSC_FALSE;

  PetscFunctionBeginUser;
  mainPart.pathAi = JoinPath(folder, "Matrix_A_int.bin");
  mainPart.pathAv = JoinPath(folder, "Matrix_A_cmplx.bin");
  auxPart.pathAi  = JoinPath(folder, "Matrix_A_aux_int.bin");
  auxPart.pathAv  = JoinPath(folder, "Matrix_A_aux_cmplx.bin");

  /* ---- Read the indices and detect the format ---- */
  PetscCall(FileExists(comm, auxPart.pathAi, &hasAux));
  PetscCall(ReadPartIndices(comm, mainPart));
  if (hasAux) PetscCall(ReadPartIndices(comm, auxPart));

  if (hasAux) { /* Hermitian volumetric part + Robin part */
    mainPart.expand = EXPAND_CONJUGATE;
    auxPart.expand  = auxPart.OneTriangle() ? EXPAND_TRANSPOSE : EXPAND_NONE;
    *format         = auxPart.OneTriangle() ? "Hermitic-Symmetric" : "Hermitic-Full";
  } else {
    mainPart.expand = mainPart.OneTriangle() ? EXPAND_TRANSPOSE : EXPAND_NONE;
    *format         = mainPart.OneTriangle() ? "Symmetric" : "Full-matrix";
  }
  nnzFile[0] = mainPart.nnz;
  nnzFile[1] = hasAux ? auxPart.nnz : 0;

  /* ---- COO index arrays for all parts ---- */
  auxPart.cooOffset    = mainPart.NumCoo();
  const long long nCoo = mainPart.NumCoo() + (hasAux ? auxPart.NumCoo() : 0);

  std::vector<PetscInt> ci((size_t)nCoo), cj((size_t)nCoo);
  long long             nBad = 0, nBadGlob = 0;
  FillPartIndices(mainPart, N, ci.data(), cj.data(), &nBad);
  if (hasAux) FillPartIndices(auxPart, N, ci.data() + auxPart.cooOffset, cj.data() + auxPart.cooOffset, &nBad);
  PetscCallMPI(MPI_Allreduce(&nBad, &nBadGlob, 1, MPI_LONG_LONG, MPI_SUM, comm));
  PetscCheck(nBadGlob == 0, comm, PETSC_ERR_FILE_UNEXPECTED, "%lld matrix indices out of range (check that Vector_B.bin matches the matrix files)", nBadGlob);

  /* ---- Distributed matrix and nonzero pattern ---- */
  PetscCall(MatCreate(comm, A));
  PetscCall(MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, (PetscInt)N, (PetscInt)N));
  PetscCall(MatSetType(*A, MATAIJ));
  PetscCall(MatSetFromOptions(*A));
  PetscCall(MatSetPreallocationCOO(*A, (PetscCount)nCoo, ci.data(), cj.data()));
  std::vector<PetscInt>().swap(ci);
  std::vector<PetscInt>().swap(cj);

  /* ---- Values (entries at the same position, e.g. main + aux, are added) ---- */
  std::vector<PetscScalar> val((size_t)nCoo);
  PetscCall(ReadPartValues(comm, mainPart, val.data()));
  if (hasAux) PetscCall(ReadPartValues(comm, auxPart, val.data() + auxPart.cooOffset));
  PetscCall(MatSetValuesCOO(*A, val.data(), INSERT_VALUES));

  /* A complex symmetric matrix can use -pc_type cholesky (LDL^T) */
  if (*format == "Symmetric") {
    PetscCall(MatSetOption(*A, MAT_SYMMETRIC, PETSC_TRUE));
    PetscCall(MatSetOption(*A, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE));
    PetscCall(MatSetOption(*A, MAT_SYMMETRY_ETERNAL, PETSC_TRUE));
    PetscCall(MatSetOption(*A, MAT_HERMITIAN, PETSC_FALSE));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*====================================================================================================
   Vector input / output in ERMES format
  ====================================================================================================*/

/* Each process reads its own part of the vector */
static PetscErrorCode ReadVector(MPI_Comm comm, const std::string &path, Vec v)
{
  MPI_File     fh;
  MPI_Offset   bytes;
  PetscInt     rs, re;
  PetscScalar *pv;

  PetscFunctionBeginUser;
  PetscCall(OpenForRead(comm, path, &fh, &bytes));
  PetscCall(VecGetOwnershipRange(v, &rs, &re));
  PetscCall(VecGetArray(v, &pv));
  PetscCall(ReadAtAll(fh, (MPI_Offset)rs * kCmplxBytes, pv, 2LL * (re - rs), MPI_DOUBLE, 8, comm));
  PetscCall(VecRestoreArray(v, &pv));
  PetscCallMPI(MPI_File_close(&fh));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Each process writes its own part of the vector */
static PetscErrorCode WriteVector(MPI_Comm comm, const std::string &path, Vec v)
{
  MPI_File           fh;
  PetscInt           N, rs, re;
  const PetscScalar *pv;

  PetscFunctionBeginUser;
  PetscCall(VecGetSize(v, &N));
  PetscCall(VecGetOwnershipRange(v, &rs, &re));
  int err = MPI_File_open(comm, path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  PetscCheck(err == MPI_SUCCESS, comm, PETSC_ERR_FILE_OPEN, "Cannot create %s", path.c_str());
  PetscCallMPI(MPI_File_set_size(fh, (MPI_Offset)N * kCmplxBytes));
  PetscCall(VecGetArrayRead(v, &pv));
  PetscCall(WriteAtAll(fh, (MPI_Offset)rs * kCmplxBytes, pv, 2LL * (re - rs), MPI_DOUBLE, 8, comm));
  PetscCall(VecRestoreArrayRead(v, &pv));
  PetscCallMPI(MPI_File_close(&fh));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*====================================================================================================
   Residual monitor every n iterations (option -ermes_monitor_every <n>)
  ====================================================================================================*/

typedef struct {
  PetscInt  every; /* print every 'every' iterations        */
  PetscReal normB; /* ||b||, for the relative true residual */
} MonitorCtx;

/* Called by KSP at each iteration; prints only every mc->every iterations. The true residual
   b - A x is computed only when printing, so the extra cost is negligible for large n. */
static PetscErrorCode MonitorEvery(KSP ksp, PetscInt it, PetscReal rnorm, void *ctx)
{
  MonitorCtx *mc = (MonitorCtx *)ctx;
  Vec         res;
  PetscReal   tnorm;

  PetscFunctionBeginUser;
  if (it % mc->every) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(KSPBuildResidual(ksp, NULL, NULL, &res));
  PetscCall(VecNorm(res, NORM_2, &tnorm));
  PetscCall(VecDestroy(&res));
  PetscCall(PetscPrintf(PetscObjectComm((PetscObject)ksp), "  Iteration %8" PetscInt_FMT " : residual %e, ||b-Ax||/||b|| %e\n", it, (double)rnorm, (double)(mc->normB > 0 ? tnorm / mc->normB : tnorm)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/*====================================================================================================
   Main
  ====================================================================================================*/

int main(int argc, char **argv)
{
  MPI_Comm           comm;
  PetscMPIInt        rank, size;
  char               folder[PETSC_MAX_PATH_LEN] = ".";
  PetscBool          flg;
  PetscInt           its = 0;
  MonitorCtx         monitor = {0, 0.0};
  std::string        format;
  long long          N, nnzFile[2] = {0, 0};
  Mat                A   = NULL;
  Vec                b   = NULL, x = NULL, r = NULL;
  KSP                ksp = NULL;
  KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
  PetscReal          normB = 0.0, normR = 0.0;
  double             t0, t1, t2, t3, t4;

  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  comm = PETSC_COMM_WORLD;
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  PetscCallMPI(MPI_Comm_size(comm, &size));
  t0 = MPI_Wtime();

  /* ---------------------------------------- Options ---------------------------------------- */
  PetscCall(PetscOptionsGetString(NULL, NULL, "-ermes_folder", folder, sizeof(folder), &flg));
  if (!flg) PetscCall(PetscOptionsGetString(NULL, NULL, "-my_folder", folder, sizeof(folder), NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-ermes_monitor_every", &monitor.every, NULL));

  const std::string pathB = JoinPath(folder, "Vector_B.bin");
  const std::string pathX = JoinPath(folder, "Vector_Xo.bin");

  /* Remove any old solution, so that a failed run never leaves a previous Vector_Xo.bin behind */
  if (rank == 0) (void)remove(pathX.c_str());
  PetscCallMPI(MPI_Barrier(comm));

  /* ------------------------------ System size (from Vector_B.bin) ------------------------------ */
  {
    MPI_File   fh;
    MPI_Offset bytes;
    PetscCall(OpenForRead(comm, pathB, &fh, &bytes));
    PetscCallMPI(MPI_File_close(&fh));
    PetscCheck(bytes > 0 && bytes % kCmplxBytes == 0, comm, PETSC_ERR_FILE_UNEXPECTED, "%s: size is not a multiple of 16 bytes", pathB.c_str());
    N = (long long)(bytes / kCmplxBytes);
    PetscCheck(N <= (long long)std::numeric_limits<PetscInt>::max(), comm, PETSC_ERR_SUP, "System too large for 32-bit PetscInt: configure PETSc with --with-64-bit-indices");
  }

  PetscCall(PetscPrintf(comm, "------------------------------------------------------------\n"));
  PetscCall(PetscPrintf(comm, "Reading ERMES linear system...\n"));

  /* ------------------------------------ Matrix and vectors ------------------------------------ */
  PetscCall(AssembleMatrix(comm, folder, N, &A, &format, nnzFile));
  PetscCall(MatCreateVecs(A, &x, &b));
  PetscCall(ReadVector(comm, pathB, b));
  t1 = MPI_Wtime();

  PetscCall(PetscPrintf(comm, "- Unknowns      : %lld\n", N));
  PetscCall(PetscPrintf(comm, "- Matrix format : %s\n", format.c_str()));
  if (nnzFile[1] > 0) PetscCall(PetscPrintf(comm, "- Stored entries: %lld (Matrix_A) + %lld (Matrix_A_aux)\n", nnzFile[0], nnzFile[1]));
  else PetscCall(PetscPrintf(comm, "- Stored entries: %lld\n", nnzFile[0]));
  PetscCall(PetscPrintf(comm, "- Processes     : %d\n", (int)size));
  PetscCall(PetscPrintf(comm, "Reading ended.\n"));
  PetscCall(PetscPrintf(comm, "------------------------------------------------------------\n"));
  PetscCall(PetscPrintf(comm, "Solving linear system...\n"));

  /* ------------------------------------------ Solve ------------------------------------------ */
  PetscCall(KSPCreate(comm, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetFromOptions(ksp));
  if (monitor.every > 0) {
    PetscCall(VecNorm(b, NORM_2, &monitor.normB));
    PetscCall(KSPMonitorSet(ksp, MonitorEvery, &monitor, NULL));
  }
  PetscCall(KSPSetUp(ksp)); /* preconditioner setup / factorization */
  t2 = MPI_Wtime();
  PetscCall(KSPSolve(ksp, b, x));
  t3 = MPI_Wtime();

  PetscCall(KSPGetConvergedReason(ksp, &reason));
  PetscCall(KSPGetIterationNumber(ksp, &its));

  /* True relative residual ||b - A x|| / ||b|| */
  PetscCall(VecDuplicate(b, &r));
  PetscCall(MatMult(A, x, r));
  PetscCall(VecAYPX(r, -1.0, b));
  PetscCall(VecNorm(r, NORM_2, &normR));
  PetscCall(VecNorm(b, NORM_2, &normB));
  PetscCall(VecDestroy(&r));

  /* ------------------------------------- Write solution -------------------------------------- */
  PetscCall(WriteVector(comm, pathX, x));
  t4 = MPI_Wtime();

  PetscCall(PetscPrintf(comm, "Linear system solved.\n"));
  PetscCall(PetscPrintf(comm, "- Converged reason: %s (%" PetscInt_FMT " iterations)\n", KSPConvergedReasons[reason], its));
  PetscCall(PetscPrintf(comm, "- ||b-Ax||/||b||  : %e\n", (double)(normB > 0 ? normR / normB : normR)));
  PetscCall(PetscPrintf(comm, "- Time (s)        : read %.1f, setup %.1f, solve %.1f, write %.1f, total %.1f\n", t1 - t0, t2 - t1, t3 - t2, t4 - t3, t4 - t0));
  if (reason < 0) PetscCall(PetscPrintf(comm, "WARNING: the PETSc solver did not converge.\n"));
  PetscCall(PetscPrintf(comm, "Solution written in %s\n", pathX.c_str()));
  PetscCall(PetscPrintf(comm, "------------------------------------------------------------\n"));

  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
  return 0;
}
