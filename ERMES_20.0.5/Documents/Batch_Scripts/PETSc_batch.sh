###################################################################################################
# - Example of a PETSc call from Linux bash.
###################################################################################################

# Set GiD problem path
MY_PATH=$HOME/Simulations/ProblemName.gid/

# Set PETSc files paths
MPI_EXE=$PETSC_DIR/lib/petsc/bin/petscmpiexec
PETSC_EXE=$MY_PATH/PETSc/PETScSolver
SOLVER_INFO=$MY_PATH/PETSc/SolverInfo.info

# PETSc solver type and options
PETSC_SOLVER='-ksp_type richardson -pc_type lu -pc_factor_mat_solver_type mumps'
SOLVER_PARAMETERS='-ksp_max_it 1 -ksp_monitor_true_residual -memory_view -log_view'

# Run PETSc solver 
$MPI_EXE -n 36 $PETSC_EXE $PETSC_SOLVER $SOLVER_PARAMETERS -my_folder $MY_PATH > $SOLVER_INFO