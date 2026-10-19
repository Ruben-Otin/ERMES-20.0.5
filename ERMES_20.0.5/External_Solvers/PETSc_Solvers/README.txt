************************************************************************************
* ERMES 20.0.5 - ERMES-PETSc interface
************************************************************************************

ERMESPETScSolver reads the ERMES linear system A Xo = B directly and in parallel,
solves it with PETSc, and writes the solution Vector_Xo.bin for ERMES. No Python and
no intermediate files are needed. Files read from the problem folder:

- Matrix_A_cmplx.bin, Matrix_A_int.bin         : system matrix (coefficients, indices).
- Matrix_A_aux_cmplx.bin, Matrix_A_aux_int.bin : auxiliary matrix (Hermitic formats).
- Vector_B.bin                                 : right-hand side.

MATRIX FORMATS (detected automatically):

- Symmetric          : Matrix_A: diagonal+upper diagonal of a complex symmetric matrix.
- Full-matrix        : Matrix_A: entire matrix.
- Hermitic-Symmetric : Matrix_A: diagonal+upper diagonal of the Hermitian volumetric
                       contribution; Matrix_A_aux: diagonal + upper diagonal of the
                       complex symmetric Robin boundary contribution.
- Hermitic-Full      : Matrix_A: as above; Matrix_A_aux: entire Robin contribution,
                       added as stored.

Matrix_A_aux files present -> Hermitic format; the number of stored triangles selects 
the variant. The detected format is printed in the "*.info" file.

FILES:

- ERMESPETScSolver.cpp : Solver source code.
- makefile             : Builds the solver.
- ERMES2PETSc.sh       : Script called by ERMES (all settings at the top).

====================================================================================
= INSTALLATION
====================================================================================

1) Download PETSc from its git repository to the local folder "petsc":

   >> git clone -b release https://gitlab.com/petsc/petsc.git petsc

2) Set PETSC_DIR, PETSC_ARCH and LD_LIBRARY_PATH (they must also be set before
   compiling the solver). Use the same PETSC_DIR and PETSC_ARCH in section 1 of
   ERMES2PETSc.sh:

   >> export PETSC_DIR=$HOME/petsc
   >> export PETSC_ARCH=arch-complex-M1
   >> export LD_LIBRARY_PATH=$PETSC_DIR/$PETSC_ARCH/lib:$LD_LIBRARY_PATH

3) Configure PETSc with "./configure" inside the folder "petsc" (see the PETSc manual
   for the meaning of each option). Several configurations can be installed on the
   same machine, each with its own PETSC_ARCH. For example:

   >> ./configure PETSC_ARCH=arch-complex-M1
      --with-cc=gcc
      --with-cxx=g++
      --with-fc=gfortran
      --with-debugging=0
      --with-scalar-type=complex
      --with-64-bit-indices=1
      --download-mpich
      --download-hwloc
      --download-openblas
      --download-scalapack
      --download-mumps
      --download-metis
      --download-parmetis
      --download-ptscotch
      --download-cmake
      --download-bison
      --download-make

   Then build and check PETSc with the "make" commands printed at the end of configure.
   On a cluster, see "CLUSTER NOTES" below before configuring.

4) Compile the solver in this folder (PETSC_DIR and PETSC_ARCH set as in step 2):

   >> make ERMESPETScSolver

5) Give execution permissions to the script and the solver:

   >> chmod +x *.sh ERMESPETScSolver

====================================================================================
= USE FROM ERMES
====================================================================================

Select "External solver" on "Solver type" and type on "External solvers settings":

- Linux  : bash     /path_to/PETSc_Solvers/ERMES2PETSc.sh
- Windows: wsl bash /path_to/PETSc_Solvers/ERMES2PETSc.sh

ERMES runs the script from the problem folder. Inside a SLURM job, the solver runs
on all the nodes and tasks of the job. All the solver output, including errors, is
written to the ERMES "*.info" file.

The script can also be run directly (without ERMES) from a folder containing the
ERMES matrix and vector files:

>> cd /path_to/problem_folder
>> bash /path_to/PETSc_Solvers/ERMES2PETSc.sh

====================================================================================
= SCRIPT SETTINGS (ERMES2PETSc.sh)
====================================================================================

The settings are at the top of the script, in two sections. Nothing needs to be
changed below them.

1) SECTION 1 - RUN SETTINGS (change as needed for each problem):

   PETSC_DIR,     PETSc installation used by the solver. Several PETSc configurations
   PETSC_ARCH     can be installed on the same machine (see INSTALLATION, step 3);
                  select the one to use here.

   SolverType     Leave exactly ONE SolverType line active (without #):
                  - MUMPS LU: direct solver, works for every matrix format (default).
                  - MUMPS Cholesky: about half the memory and time of LU, but ONLY
                    valid for the "Symmetric" matrix format.
                  - LGMRES + SOR: iterative solver, low memory, but may converge
                    slowly or not at all.

   SolverOptions  Extra PETSc options added to the solver (may be empty). Useful ones:
                  -memory_view              memory used by the solver (default;
                                            negligible cost)
                  -log_view                 timing of every solver stage (small cost)
                  -mat_mumps_icntl_4 2      MUMPS info, including memory estimates
                  -mat_mumps_icntl_14 50    50% more MUMPS workspace (use it if MUMPS
                                            stops with INFOG(1)=-9)
                  -mat_mumps_icntl_28 2 -mat_mumps_icntl_29 2
                                            parallel ordering with ParMETIS: faster
                                            analysis and less memory on process 0 for
                                            very large problems (may increase fill-in)
                  -ermes_monitor_every <n>  residual every n iterations (iterative)
                  -ksp_monitor_true_residual
                                            residual at every iteration (iterative).
                                            Slows long solves (one extra matrix-vector
                                            product and one line per iteration):
                                            prefer -ermes_monitor_every.

2) SECTION 2 - CLUSTER SETTINGS (set once for each machine):

   NumParallTasks Number of MPI processes. Inside a SLURM job, all the tasks of the
                  job are used automatically; the default value (10) is used only
                  when running outside SLURM.

   RanksPerNode   MPI processes per node in SLURM jobs. Leave empty to use all the
                  tasks of the job. Set a lower value (e.g. 16 on 64-core nodes) if
                  the factorization runs out of memory: fewer processes per node
                  gives each one more memory. Total = RanksPerNode x number of nodes.

   Launcher       "petsc": PETSc's own mpiexec. Required when PETSc was configured
                  with --download-mpich (as in this README). Do not use srun or the
                  cluster's mpirun with this build.
                  "srun": only if PETSc was configured with the cluster's own MPI.

   BindTo         "core" (default) pins each MPI process to its own core, so that
                  processes do not migrate or share cores during the factorization.
                  Leave empty to disable. Only used with the "petsc" launcher.
                  Normally faster; to confirm on a new cluster, compare the solve
                  time with and without it once.

   MaxOutputLines Maximum number of solver output lines written to the "*.info" file
                  (default 0 = no limit). A normal run prints well under 100 lines.
                  A limit (e.g. 5000) stops an MPI failure, which can print endless
                  error backtraces from every process, from filling the disk with a
                  huge "*.info" file. If a limit is set, keep it well above the
                  expected output when using -ksp_monitor_true_residual on long
                  iterative solves, or the solver will be stopped at the limit.

   SolverFullPath Solver executable (default: same folder as the script).

3) AUTOMATIC SETTINGS (no action needed):

   - One thread per MPI process (OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1). PETSc,
     MUMPS and OpenBLAS configured as in this README (without --with-openmp) gain
     nothing from extra threads, and stray threads would compete for the cores.
   - LD_LIBRARY_PATH is set from PETSC_DIR and PETSC_ARCH, so the PETSc libraries
     are found without setting it before calling ERMES.
   - With the "petsc" launcher, cluster MPI settings that conflict with PETSc's own
     MPICH are removed (SLURM PMIx variables, Intel MPI PMI library), PMI version 1
     is set explicitly, and SLURM is prevented from pinning all the processes of a
     node to a single core. The thread and PMI settings are also passed explicitly
     to every MPI process. The same script therefore works on clusters with
     different SLURM and MPI setups.
   - Before starting, the script checks the settings, the PETSc installation, the
     solver executable and the ERMES input files, and stops with a one-line error
     if anything is missing or wrong.
   - The run information (solver, options, processes, binding, SLURM job and nodes,
     start and finish times) is written at the top of the "*.info" file.
   - Any old Vector_Xo.bin is deleted before solving, and any partial one after a
     failure, so ERMES never reads an old or incomplete solution.

====================================================================================
= CLUSTER NOTES
====================================================================================

Clusters often limit the amount of time and the number of processes that can be
run on login nodes and, typically, do not provide internet access from compute 
nodes. To overcome these limitations, either:
   
- Configure and build PETSc on a login node, limiting the parallel build to the 
  login-node rules (e.g. at most 4 cores), and run it inside tmux or screen, so 
  a dropped connection does not stop it:
	 
  >> tmux new -s petsc
  >> cd $HOME/petsc
  >> nice -19 ./configure PETSC_ARCH=arch-complex-M1 
     --with-cc=gcc --with-cxx=g++ --with-fc=gfortran 
     --with-debugging=0 --with-scalar-type=complex 
     --with-64-bit-indices=1 --download-mpich --download-hwloc 
     --download-openblas --download-scalapack --download-mumps 
     --download-metis --download-parmetis --download-ptscotch 
     --download-cmake --download-bison --download-make 
     --with-make-np=4

- Or download the packages on the login node first and build on a compute node:
   
  1) Creates a list of package URLs:
	 
     >> mkdir -p $HOME/petsc-pkgs
     >> cd $HOME/petsc
     >> ./configure PETSC_ARCH=arch-complex-M1
        --with-cc=gcc --with-cxx=g++ --with-fc=gfortran 
        --with-debugging=0 --with-scalar-type=complex --with-64-bit-indices=1 
        --download-mpich --download-hwloc --download-openblas --download-scalapack 
        --download-mumps --download-metis --download-parmetis --download-ptscotch 
        --download-cmake --download-bison --download-make 
        --with-packages-download-dir=$HOME/petsc-pkgs | tee pkg-list.txt
	 
  2) Download packages into $HOME/petsc-pkgs:	
	 
     >> cd $HOME/petsc-pkgs
     >> grep "\['" $HOME/petsc-3.26/pkg-list.txt | 
        grep -oE "https?://[^']+\.(tar\.gz|tgz|tar\.bz2|zip)" | 
        awk -F/ '!seen[$3 FS $4]++' | xargs -n1 wget -nc
     >> ls

  3) Configure using the downloaded packages:
	 
     >> ./configure PETSC_ARCH=arch-complex-M1
        --with-cc=gcc --with-cxx=g++ --with-fc=gfortran 
        --with-debugging=0 --with-scalar-type=complex --with-64-bit-indices=1 
        --download-mpich --download-hwloc --download-openblas --download-scalapack 
        --download-mumps --download-metis --download-parmetis --download-ptscotch 
        --download-cmake --download-bison --download-make 
        --with-packages-download-dir=$HOME/petsc-pkgs

An interrupted configure can be restarted with exactly the same options and the 
same PETSC_ARCH: packages already downloaded and built are reused.

Cluster modules can interfere with PETSc's MPICH (e.g. Intel MPI loaded by 
default). When running "make check" or PETSc programs by hand, put PETSc's bin 
folder first in PATH, so its mpiexec is used:

>> export PATH=$PETSC_DIR/$PETSC_ARCH/bin:$PATH

ERMES2PETSc.sh already handles this, always calling PETSc's mpiexec by its full path.

====================================================================================
= TROUBLESHOOTING
====================================================================================

-) "Found both env vars PMI_SIZE and PMIX_NAMESPACE" or "PMI_Init returned 14":
    The cluster's MPI environment is clashing with PETSc's MPICH. ERMES2PETSc.sh
    fixes this automatically. If it appears when running PETSc by hand, use:
    >> export SLURM_MPI_TYPE=none MPIR_CVAR_PMI_VERSION=1
    >> unset I_MPI_PMI_LIBRARY $(compgen -e | grep '^PMIX_')

-) "solver output exceeded ... lines":
    The output reached MaxOutputLines (only when a limit is set). Check the first 
    lines of the "*.info" file for the actual error.

-) "MUMPS stops with INFOG(1)=-9":
    Not enough MUMPS workspace: add -mat_mumps_icntl_14 50 (or higher) to 
    SolverOptions.

-) "Out of memory during the factorization"
    Use more nodes, set a lower RanksPerNode, or use Cholesky if the matrix format
    is "Symmetric".

-) "Quick MPI test on a cluster (inside an interactive job)":
    >> cd $PETSC_DIR/src/snes/tutorials && make ex19
    >> $PETSC_DIR/$PETSC_ARCH/bin/mpiexec -n 4 ./ex19

For more information, visit the PETSc website at https://petsc.org.
