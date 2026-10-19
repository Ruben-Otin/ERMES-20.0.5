###################################################################################################
# - Example of batch script to run ERMES on Linux bash.
###################################################################################################

# Set paths
ERMES_EXE=$HOME/GiDx64/gid-15.0.4/problemtypes/ERMES_20.0.5/ERMES.gid/execs/unix/ERMES_20.0.5
PROBLEM_NAME=GiDProblemName
PROBLEM_PATH=$HOME/Projects/$PROBLEM_NAME.gid
SOLVER_PATH=$HOME/Projects/PETSc

# Set working directory
cd $PROBLEM_PATH

# Execute ERMES ( non-zero exit status on error; stderr also goes to the .info file )
$ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info 2>&1 || { echo "ERMES failed, see $PROBLEM_NAME.info"; exit 1; }

# Solve the system with PETSc ( ERMES "External_solver" with empty settings writes the system and stops )
bash $SOLVER_PATH/ERMES2PETSc.sh >> $PROBLEM_PATH/$PROBLEM_NAME.info || { echo "PETSc solver failed, see $PROBLEM_NAME.info"; exit 1; }