# ERMES 20.0.5 — corrections and improvements with respect to ERMES 20.0.4

All changes keep the input format, the GiD interface and the results of correct
models unchanged, except where a result was wrong (marked **[results change]**).
Every modified spot carries a comment in the code explaining what changed and why.

## 1. Bug fixes

### Gaussian beam source (`ERMES/GenericRobin_{1st,2nd,3sb,1eg,2eg}_FullWave.cpp`)
- **[results change]** Amplitude normalization was a fourth root; correct is
  `sqrt( zR_pol zR_per / (|q_pol||q_per|) )`. The old code overdrove ports placed away
  from the waist by `sqrt(w/w0)`. No effect when the port is at the waist plane.
  Manual equation (3.11) corrected in the same way. Note: the old fourth root coincides with a
  2D beam ( slab models ) when both radii are equal; for 2D models set the radius of the invariant
  direction very large ( e.g. 1e3 ), which gives the exact 2D beam with 20.0.5 ( verified: |E| = 0.99
  and Poynting 1.33e-3 W/m^2 = 1/(2 eta0) at the waist, against 1.15 / 1.79e-3 with 20.0.4 ).
- Polarization vector is made orthogonal to the propagation vector (Gram–Schmidt)
  instead of silently producing a distorted beam when they are not perpendicular.
- Zero-length vectors and non-positive waist radii no longer produce NaN.
- **[results change, accuracy]** Nodal elements (1st, 2nd, 3sb) now evaluate the beam at
  the 12th-order Gauss points of each face, like the edge elements already did, instead
  of interpolating nodal values. Resolves the beam phase better on coarse port meshes.
  Curl and divergence terms are applied in one pass.
- **[results change, accuracy]** The same applies to the plane wave sources of the nodal Robin
  elements ( single plane wave and multiple plane waves, 1st, 2nd and 3sb ): the field
  E0 exp( jk k.r ) is now evaluated at the Gauss points instead of interpolating nodal values.
  Checked and unchanged: rectangular waveguide TE10 and coaxial TEM ports ( nodal elements
  already used Gauss points ), far field ( no source ), J sources and imported fluxes ( data
  given at the nodes, interpolated identically in nodal and edge elements ).

### Cold plasma (`ERMES/ColdPlasma.cpp/.h`)
- **[results change]** Hermitian storage (`Herm-Full`, `Herm-Symm`) with active damping
  produced a wrong matrix (the lower triangle is rebuilt as the conjugate, which is only
  valid for a lossless tensor). ERMES now switches to `Full_matrix` with a warning.
- B = 0 at a node gave a NaN tensor (division by |B| in the rotation). Now identity
  (the tensor is isotropic there).
- **[results change when grids differ]** The 1D B-field file positions were read and
  discarded: B was interpolated on the density grid, wrong if the grids differed and an
  out-of-bounds read if the B file was shorter. B now has its own grid (`mLo_B`).
- Both 1D grids are checked to be strictly ascending (required by the interpolation).
- Profile readers no longer drop lines that start with blanks, a sign or a decimal point;
  malformed data lines and out-of-range node ids stop with a clear error.
- Imported K-tensor file: header lines no longer consume a tensor nor write to index -1;
  a short binary file is an error (previously a silent zero tensor); bounds-checked lookup.
- All input files are checked for successful opening.
- **[results change, ~1e-4 relative]** CODATA 2018 constants from the new
  `ERMES/Utils/PhysicalConstants.h` (previously mixed values, e.g. me = 9.1095e-31).
  Ion masses subtract the removed electrons (neutral masses were ~0.05 % heavy for H+,
  shifting ion-cyclotron resonances by ~1 mm at R ~ 2 m).
- Species table built once in `Load_Data` instead of 15 heap-allocated species per node
  evaluation; absent species are skipped. Quasi-neutrality warning if sum(Z n_i)/n_e != 1.
- Exact resonance check fixed (it tested `w ± Wc` without the collision term).

### Iterative solvers (`linear_solvers/`)
- A NaN/Inf residual (solver breakdown) was reported as **"CONVERGENCE ACHIEVED"**
  (`NaN > tol` is false). Now reported as breakdown and ERMES stops with a failure status.
- The `_cp_` solvers ran to the iteration limit on NaN; they now stop.
- Warning when the true residual exceeds 10x the tolerance while the recursive estimate
  says converged (residual gap).
- `IniSolver(const char*)`, virtual destructor in `ComplexSolver`.

### Error handling and exit status
- All fatal errors returned exit status 0 (`exit(EXIT_SUCCESS)`), so SLURM `afterok`
  chains and scripts treated failed runs as successful. They now return `EXIT_FAILURE`
  and also print to stderr (`Modeler`, `ColdPlasma`, `Kernel` input reading, solver
  breakdown).
- `main()`: usage message instead of an uncaught exception when no argument is given;
  any escaping exception ends with a message and a failure status (previously abort).

### Portability / undefined behaviour
- `includes/compare_pdof.h`: comparator made `const` — required since C++17; the old code
  **did not compile with GCC >= 11** in its default mode.
- `includes/vector.h`: missing `#include <ios>`.
- `ERMES/Utils/*.h`: wrong relative include paths (`../includes`) that only worked through
  the `-I../external_libraries/` search path.
- `includes/node.h`: `pDof()` could fall off the end of a non-void function.
- `includes/csr_space.h`: `GetColumn` declared non-void without return.
- `sources/input_c_parser.cpp`: virtual destructors for `Statement` and `ValueStatement`
  (objects deleted through base pointers).

## 2. Performance

### Sparse matrix-vector products (`linear_solvers/complex_matrix_operations.*`)
- Every product used to allocate a full-size vector per thread, then add them serially
  in an `omp critical` section. Now:
  - per-thread scatter buffers are allocated once and reused,
  - each thread's buffer only spans the rows it can touch (upper-triangular storage),
    ~half the memory of the old implementation,
  - the reduction is parallel,
  - full-storage products write rows directly with no buffers at all.
- Dot products, norms and vector updates (all serial before) are OpenMP-parallel above
  20 000 entries, removing the serial fraction of each solver iteration.
- One templated kernel replaces ten hand-copied variants; semantics are identical
  (verified against dense references by `tests/test_matrix_operations.cpp`).
- Note: the sandbox where this was prepared has a single core, so the parallel speed-up
  could not be measured; single-thread speed is unchanged. Please benchmark on a node.

### Build
- `-O3` instead of `-O2`; optional `NATIVE=1` (`-march=native`).

## 3. Build system (`unix/`)
- New hand-written `Makefile` (old one kept as `Makefile_codeblocks`):
  explicit `-std=c++11`, warnings enabled for ERMES code (no more `-w -fpermissive`),
  header dependency tracking, automatic source discovery, `DEBUG=1`, `SANITIZE=1`,
  `NATIVE=1`, and `make test`.
- With the new flags the whole ERMES code compiles without warnings ( -O0 and -O3 ); the
  static link only reports the usual glibc `dlopen` note.
- `unix/README.txt` rewritten.
- `windows/.vs/` (50 MB Visual Studio cache) removed from the distribution. The Visual
  Studio project needs no change (only a new header was added).

## 4. New: internal direct solver (`linear_solvers/direct_complex_solver.*`)
- Sparse LU for complex matrices in every ERMES storage (full-wave symmetric,
  plasma Full, Herm-Symm, Herm-Full): approximate minimum degree ordering
  (AMD-type degree, aggressive absorption) + left-looking LU (Gilbert–Peierls) with
  threshold partial pivoting + up to 3 steps of iterative refinement.
- Selected with `Direct_solver` in GiD (the option has to be added to the problemtype,
  which is not part of this source tree) or with the environment variable
  `ERMES_SOLVER=Direct`, which overrides the GiD choice.
- Prints n, nnz(L+U), memory, fill ratio, off-diagonal pivots and the true residual;
  singular matrices stop ERMES with a failure status.
- Serial; for small and medium problems (one core, 3D 7-point grid: n = 27 000 in 45 s,
  280 MB). Large problems: iterative solvers or PETSc/MUMPS via the external solver.
- Symbolic analysis before the factorization: predicted nnz(L+U), memory and rough time are
  printed, and ERMES stops with an error if the factors would exceed 80% of the physical memory.
  Progress ( % of predicted work, elapsed and remaining time ) is printed every 5%.
- Tested by `tests/test_direct_solver.cpp` (all storages, zero diagonals needing
  pivoting, n = 1 … 729, singular detection; relative residuals ~1e-15).
- Added to the Visual Studio project.

### Multifrontal version ( 20.0.5 )
- The direct solver now uses a supernodal multifrontal factorization
  ( `linear_solvers/multifrontal_complex_solver.*` ): LDL^T for complex symmetric matrices,
  LU otherwise, relaxed supernodes, cache-blocked dense kernels parallel with OpenMP, static
  pivoting with iterative refinement and automatic fallback to the threshold-pivoting LU.
- One core: 4-7x faster than the first version on 3D FEM-like matrices; additionally
  parallel. Verified against all storage formats by `tests/test_direct_solver.cpp` and with
  1-8 threads.

### Direct solver optimization ( 20.0.5 )
- Blocked factorization of each front ( 64-column blocks ): almost all the work now goes through
  the dense kernel, including the Schur complement.
- Dense kernel with packed real / imaginary blocks and 4-column register blocking, so the inner
  loops vectorize; runtime CPU dispatch ( GCC ): generic, AVX2/FMA and AVX-512 versions in the
  same portable executable.
- Two-level parallelism: independent subtrees of the elimination tree are factorized in parallel
  ( one thread each ), then the large fronts at the top use all threads in the dense kernels.
- One core, 3D test grids: n = 27 000 in 4.6 s ( LDL^T ) / 8.5 s ( LU ), versus 45 s with the first
  direct solver ( ~10x ). Parallel speed-up not measured ( single-core test machine ).
- `Element` base class with a virtual destructor ( all elements are deleted polymorphically ).

### Direct solver robustness ( 20.0.5 )
- Fix: an indefinite matrix ( 165 856 unknowns, 32 518 perturbed pivots ) ended in NaN. Causes: no
  scaling, too small static perturbation ( 1e-10 ) and plain iterative refinement; the fallback to
  the threshold LU was also never reached ( pivot counter cleared before the check ).
- Now: Ruiz equilibration, pivoting inside the fronts ( partial pivoting for LU, symmetric 1x1 for
  LDL^T ), perturbation 1e-8, NaN/Inf detection during the factorization, GMRES(30) preconditioned
  by the factorization instead of iterative refinement, and automatic refactorization as LU if
  LDL^T fails. Saddle-point stress tests with 11% perturbed pivots converge in 2 GMRES iterations
  to 1e-14.
- External BLAS: called only outside OpenMP parallel regions ( top fronts, using the BLAS's own
  threads ); inside parallel subtrees ERMES's kernel is used. Concurrent calls to a single-threaded
  BLAS build ( not thread-safe ) were the likely cause of the NaN seen with OpenBLAS + METIS.
  The BLAS is used for the top fronts only if it is multithreaded ( openblas_get_parallel /
  bli_info_get_enable_threading ) or ERMES runs with one thread; otherwise ERMES's parallel kernel
  is faster. ERMES_BLAS=OFF forces ERMES's kernel. The .info file shows the kernel used.

### METIS bundled ( 20.0.5 )
- METIS 5 and GKlib sources ( KarypisLab, Apache-2.0 ) bundled in `external_libraries/metis` and
  `external_libraries/GKlib`, compiled automatically by the Makefile ( default, `METIS=0` to
  disable, `METIS_DIR` for an installed METIS ) and by the Visual Studio project ( all
  configurations, separate object folders because both contain a graph.c ). The vcpkg-based
  `Release_Libs` configuration was removed: no installation is needed anywhere.

### Linux executable with METIS ( 20.0.5 )
- The distributed static Linux executable now includes METIS 5 and GKlib ( KarypisLab sources,
  built statically, 32-bit indices ). Test, n = 119 164 FEM-like matrix: METIS 2.4e11 operations
  against 1.2e12 with minimum degree ( 5x less work, 2x less memory ).
- execs/licenses/: METIS / GKlib licenses and THIRD_PARTY_NOTICES.txt for the executables.

### Ordering speed ( 20.0.5 )
- The built-in minimum degree ordering ( used without METIS, e.g. the distributed Linux executable )
  was very slow for large FEM matrices ( 161 s for n = 165 856; it seemed to hang on a cluster ).
  Now: graph compression of indistinguishable unknowns ( Ax, Ay, Az, V of the same node ), weighted
  approximate minimum degree on the compressed graph and O(1) bucket degree lists. 30-40x faster
  ( n = 119 164, 4 unknowns per node, 27-point coupling: 0.8 s ) with similar fill.

### External libraries ( 20.0.5 )
- **Boost removed**: `boost::shared_ptr` replaced by `std::shared_ptr` ( C++11 ) everywhere
  ( Legacy folders excepted ). The 364 MB `external_libraries/boost` folder is no longer needed;
  missing standard includes that came through Boost added. `std::bind2nd` ( removed in C++17,
  rejected by Visual C++ 2022 ) replaced by plain loops in `includes/vector.h`; the whole code
  compiles as C++11 and C++17.
- Visual Studio project: reference to `external_libraries/libc.lib` removed ( a 32-bit lcc-win32
  library, ignored by the x64 linker ); the file can be deleted.
- Optional, open-source libraries for the direct solver ( `make BLAS=openblas METIS=1` ):
  optimized BLAS ( `cblas_zgemm`, OpenBLAS / BLIS ) for the dense kernels and METIS nested
  dissection ordering. Without them ERMES uses its own kernel and minimum degree ordering.
  The four combinations were tested ( BLAS with the reference BLAS, METIS with a stand-in
  library ); the real OpenBLAS / METIS were not available on the test machine.
- ERMES sets OpenBLAS / BLIS to one thread itself ( no environment variable needed ).
- Visual Studio: new configuration `Release_Libs | x64` ( METIS through vcpkg, static library and
  static runtime ); `Release | x64` is unchanged and needs no libraries. OpenBLAS is not used on
  Windows ( vcpkg's Visual C++ build has generic kernels only: 122 s vs 57 s with ERMES's kernel ).
- Measured ( Windows, n = 165 856, 16 threads ): previous direct solver 582 s; with METIS and the
  robust multifrontal solver 62 s, relative residual 8e-13, 1 perturbed pivot, 1 GMRES iteration.
- `INSTALL_OPTIONAL_LIBRARIES.txt`: step-by-step portable builds for Windows, Linux and clusters.

## 5. Compatibility fixes with the GiD problemtype and utilities
- External solver: any old `Vector_Xo.bin` is deleted before the external solver is called,
  and a missing or incomplete solution now stops ERMES with an error ( previously it wrote a
  zero or stale solution as if valid ). The solver exit status is printed. Also applies to the
  second ERMES call that only reads the solution.
- Old `Matrix_A_aux_*.bin` files are deleted when the current storage has no auxiliary matrix,
  so external solvers that detect the format from the files ( ERMESPETScSolver, NumPy scripts )
  cannot mix files from different runs.
- Failure to create the system matrix files is fatal ( was silently ignored ).
- Cold plasma: file paths in the plasma data file may contain spaces.
- Periodic boundary conditions: a face node without partner on the opposite face stops ERMES
  with an error ( it was paired with an uninitialized node id ).
- Variables reported as possibly uninitialized by the -O3 build initialized ( 2nd order edge
  face bases, BiCG residual when no iteration runs, edge node loop ).
- `make STATIC=1`: fully static executable for distribution with the GiD problemtype.

## 6. Version 20.0.5
- Version number 20.0.5 in folder, executable, manual and problemtype names, README files and
  scripts ( historical references to 20.0.4 kept where they describe the previous version ).
- ERMES prints its version at start-up ( first line of the `.info` file:
  `ERMES 20.0.5 analysis started...` ) and with `ERMES --version`.
- Manual: Direct solver, exit status, external-solver checks and Hermitian/damping storage
  documented; Gaussian beam equation (3.11) corrected ( square root ); duplicated LaTeX label
  fixed. PDF rebuilt.

## 7. Documentation
- Manual ( chapters 1-3 ) and GiD help updated: bundled METIS, no Boost / cbp2make, `make` and
  Release x64 builds, static Linux executable, final direct solver description. PDF rebuilt.
- `RELEASE_NOTES_20.0.5.txt` ( user summary, in the problemtype, source and Utilities folders ).
- Solver header comments updated to the final direct solver design.
- `README.txt`, `unix/README.txt`, `windows/README.txt` completed; new
  `linear_solvers/README.txt` (storages, solvers, how to choose, limits).

## 8. Not changed (recommendations)
- The per-element-order code duplication (5 copies of each boundary/volume element) is
  the main maintainability issue; merging it needs regression cases to validate.
- Suggested regression suite: plane wave, Gaussian-beam power conservation, TE10 port,
  cold-plasma slab vs dispersion relation.
- `Kij_Tol` drops element-matrix entries by absolute value before summation.
- Gaussian beam assumes vacuum (k0) at the port; a refractive-index input would need
  new GUI fields.
- `sources/update.cpp` does not compile and is excluded from the build (legacy).
- `Legacy_*` folders are kept untouched.
