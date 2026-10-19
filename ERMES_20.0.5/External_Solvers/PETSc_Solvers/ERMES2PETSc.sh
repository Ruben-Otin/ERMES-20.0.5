#!/bin/bash
###################################################################################################
#
#  ERMES2PETSc.sh  -  ERMES 20.0 / PETSc interface
#
#  Solves the ERMES linear system with ERMESPETScSolver, which reads the ERMES files directly in
#  parallel and writes the solution Vector_Xo.bin for ERMES.
#
#  Usage from ERMES:
#    1) Select "External solver" on ERMES "Solver type".
#    2) Type on ERMES "External solvers settings":
#       >> bash     /path_to/PETSc_Solvers/ERMES2PETSc.sh (Linux)
#       >> wsl bash /path_to/PETSc_Solvers/ERMES2PETSc.sh (Windows)
#
#  ERMES runs this script from the problem folder. The matrix format (Symmetric, Full-matrix,
#  Hermitic-Symmetric or Hermitic-Full) is detected automatically. Inside a SLURM job, the solver
#  runs on all the nodes and tasks of the job. All the solver output, including errors, goes to
#  the ERMES "*.info" file.
#
#  Edit only sections 1 and 2 below.
#
###################################################################################################

#==================================================================================================
#= 1. RUN SETTINGS  -  change as needed for each problem
#==================================================================================================

# ---- PETSc installation -------------------------------------------------------------------------

export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-complex-M1

# ---- Solver type: leave exactly ONE SolverType line active (without #) --------------------------
# 1) Direct solver, MUMPS LU. Works for every ERMES matrix format (recommended default).
# 2) Direct solver, MUMPS Cholesky. About half the memory and time of LU, but ONLY valid for the
#    "Symmetric" matrix format (not for Full-matrix or Hermitic formats).
# 3) Iterative solver, LGMRES + SOR. Low memory, but may converge slowly or not at all.

SolverType="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
#SolverType="-ksp_type preonly -pc_type cholesky -pc_factor_mat_solver_type mumps"
#SolverType="-ksp_type lgmres -pc_type sor -ksp_rtol 1e-6 -ksp_max_it 1000000"

# ---- Solver options: added to the solver type above (leave empty for none) ----------------------
#   -memory_view                 memory used by the solver (negligible cost)
#   -log_view                    timing of every solver stage (small cost)
#   -mat_mumps_icntl_4 2         MUMPS information, including memory estimates
#   -mat_mumps_icntl_14 50       50% more MUMPS workspace (if MUMPS stops with INFOG(1)=-9)
#   -mat_mumps_icntl_28 2 -mat_mumps_icntl_29 2
#                                parallel ordering with ParMETIS: faster analysis and less memory
#                                on process 0 for very large problems (may increase fill-in)
#   -ermes_monitor_every 1000    residual every 1000 iterations (iterative solver)
#   -ksp_monitor_true_residual   residual at every iteration (iterative solver). Slows long
#                                solves: one extra matrix-vector product and one line per iteration

SolverOptions="-memory_view"

#==================================================================================================
#= 2. CLUSTER SETTINGS  -  set for each machine
#==================================================================================================

# Number of MPI processes. Inside a SLURM job, all the tasks of the job are used; the number after
# ":-" is used only when running outside SLURM.
NumParallTasks="${SLURM_NTASKS:-10}"

# MPI processes per node (SLURM jobs). Leave empty to use all the tasks of the job. Set a lower
# value (e.g. "16" on 64-core nodes) if the factorization runs out of memory: fewer processes per
# node gives each one more memory. Total processes = RanksPerNode x number of nodes.
RanksPerNode=""

# MPI launcher: "petsc" = PETSc's own mpiexec (required when PETSc was configured with
# --download-mpich, as in the README). "srun" only if PETSc was configured with the cluster MPI.
Launcher="petsc"

# Process binding ("petsc" launcher): "core" pins each MPI process to its own core, so processes
# do not migrate or share cores during the factorization. Leave empty to disable.
BindTo="core"

# Maximum number of solver output lines written to the "*.info" file (0 = no limit). A normal run
# prints well under 100 lines. The limit stops an MPI failure, which can print endless error
# backtraces from every process, from filling the disk with a huge "*.info" file.
MaxOutputLines=0

# Solver executable (default: in the same folder as this script)
SolverFullPath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/ERMESPETScSolver"

#==================================================================================================
#= 3. SCRIPT  -  no changes needed below this line
#==================================================================================================

# Send error messages to the standard output, so they also appear in the ERMES *.info file
exec 2>&1

export LD_LIBRARY_PATH="$PETSC_DIR/$PETSC_ARCH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
FolderPath="$(pwd)"

# Stop with a clear message, instead of a cryptic MPI or PETSc error
Fail() {
    echo "ERROR: $1"
    echo "ERROR: PETSc solver not run. No solution written."
    exit 1
}

# ---- Checks before starting MPI ----
[ -n "$SolverType" ]                || Fail "no SolverType selected in section 1."
[ "$NumParallTasks" -gt 0 ] 2>/dev/null || Fail "NumParallTasks must be a positive number (now \"$NumParallTasks\")."
[ -z "$RanksPerNode" ] || [ "$RanksPerNode" -gt 0 ] 2>/dev/null || Fail "RanksPerNode must be empty or a positive number (now \"$RanksPerNode\")."
[ -d "$PETSC_DIR/$PETSC_ARCH/lib" ] || Fail "PETSc not found in $PETSC_DIR/$PETSC_ARCH (check PETSC_DIR and PETSC_ARCH)."
[ -x "$SolverFullPath" ]            || Fail "solver $SolverFullPath not found or not executable (compile it with make and run chmod +x)."
for f in Matrix_A_cmplx.bin Matrix_A_int.bin Vector_B.bin; do
    [ -s "$FolderPath/$f" ]         || Fail "ERMES file $FolderPath/$f not found or empty."
done

# ---- Threads ----
# One thread per MPI process. PETSc, MUMPS and OpenBLAS built as in the README (without
# --with-openmp) gain nothing from extra threads, and stray threads would compete for the cores.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ---- Processes per node ----
PerNodeFlag=""
if [ -n "$RanksPerNode" ]; then
    NumParallTasks=$(( RanksPerNode * ${SLURM_JOB_NUM_NODES:-1} ))
    if [ "$Launcher" = "srun" ]; then
        PerNodeFlag="--ntasks-per-node=$RanksPerNode"
    else
        PerNodeFlag="-ppn $RanksPerNode"
    fi
fi

# ---- MPI launcher and its environment ----
LaunchFlags=""
case "$Launcher" in
    srun)
        MpiExec="srun"
        ;;
    petsc)
        MpiExec="$PETSC_DIR/$PETSC_ARCH/bin/mpiexec"
        [ -x "$MpiExec" ] || Fail "PETSc mpiexec not found in $PETSC_DIR/$PETSC_ARCH/bin."

        # PETSc's MPICH talks to its own mpiexec through PMI version 1. Clusters whose SLURM sets
        # up PMIx (e.g. CUMULUS) make MPICH abort ("Found both env vars PMI_SIZE and
        # PMIX_NAMESPACE"), and Intel MPI settings (e.g. CSD3) do not apply to it. Remove them and
        # set PMI version 1 explicitly.
        export SLURM_MPI_TYPE=none
        PmixVars=$(compgen -e | grep '^PMIX_')
        [ -n "$PmixVars" ] && unset $PmixVars
        unset I_MPI_PMI_LIBRARY
        export MPIR_CVAR_PMI_VERSION=1

        # Inside SLURM, PETSc's mpiexec starts one helper per node with srun. Stop SLURM pinning
        # that helper (and so all the MPI processes it starts on the node) to a single core.
        export SLURM_CPU_BIND=none

        # Binding, and the environment passed explicitly to every process (PETSc's mpiexec does
        # not always forward it under SLURM)
        [ -n "$BindTo" ] && LaunchFlags="-bind-to $BindTo"
        LaunchFlags="$LaunchFlags -genv OMP_NUM_THREADS 1 -genv OPENBLAS_NUM_THREADS 1 -genv MPIR_CVAR_PMI_VERSION 1"
        ;;
    *)
        Fail "unknown Launcher \"$Launcher\" in section 2 (use petsc or srun)."
        ;;
esac

# ---- Run information ----
echo "------------------------------------------------------------"
echo "- PETSc  solver : $SolverType"
echo "- Solver options: ${SolverOptions:-none}"
echo "- Processes     : $NumParallTasks${RanksPerNode:+ ($RanksPerNode per node)}, 1 thread each"
[ "$Launcher" = "petsc" ] && echo "- Binding       : ${BindTo:-none}"
[ -n "$SLURM_JOB_ID" ]    && echo "- SLURM job     : $SLURM_JOB_ID on ${SLURM_JOB_NUM_NODES:-1} node(s): $SLURM_JOB_NODELIST"
echo "- Started       : $(date '+%Y-%m-%d %H:%M:%S')"

# Remove any old solution, so that a failed run is never taken as valid by ERMES
rm -f "$FolderPath/Vector_Xo.bin"

# ---- Solve ----
RunSolver() {
    "$MpiExec" -n "$NumParallTasks" $PerNodeFlag $LaunchFlags "$SolverFullPath" \
        -ermes_folder "$FolderPath" $SolverType $SolverOptions
}
if [ "$MaxOutputLines" -gt 0 ] 2>/dev/null; then
    # awk passes the output through line by line (flushed at once, so the "*.info" file shows the
    # progress live) and stops after MaxOutputLines lines; the solver then stops at its next write
    RunSolver 2>&1 | awk -v n="$MaxOutputLines" '
        { print; fflush() }
        NR >= n { print "ERROR: solver output exceeded " n " lines; solver stopped (see MaxOutputLines)."; fflush(); exit }'
    Status=${PIPESTATUS[0]}
else
    RunSolver
    Status=$?
fi

echo "- Finished      : $(date '+%Y-%m-%d %H:%M:%S')"

# ---- Result: a missing solution means the solver failed ----
if [ "$Status" -ne 0 ] || [ ! -s "$FolderPath/Vector_Xo.bin" ]; then
    echo "ERROR: PETSc solver failed (exit code $Status). No solution written."
    rm -f "$FolderPath/Vector_Xo.bin"
    exit 1
fi
exit 0
