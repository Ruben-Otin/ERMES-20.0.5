###################################################################################################
# - Example of batch script to run ERMES on Linux bash using GiD offscreen.
# - See file "GiD_offscreen.bch" for an example of how to run GiD in batch mode.
###################################################################################################

# Set paths
PROBLEM_NAME=GiDProblemName
PROBLEM_PATH=$HOME/Projects/$PROBLEM_NAME.gid

GID_PATH=$HOME/GiDx64/gid-15.0.4
ERMES_EXE=$GID_PATH/problemtypes/ERMES_20.0.5/ERMES.gid/execs/unix/ERMES_20.0.5
MACROS_PATH=$PROBLEM_PATH/GiDMacros
SOLVER_PATH=$HOME/Projects/PETSc

cd $PROBLEM_PATH

# Mesh with GiD offscreen
$GID_PATH/gid_offscreen -b $MACROS_PATH/Mesh.bch -n

# Delete previous files
rm -f $PROBLEM_NAME*dat Matrix_* Vector_* *.flavia.res *.info

# Calculate static currents after setting ERMES with GiD offscreen
$GID_PATH/gid_offscreen -b $MACROS_PATH/Calculate_J.bch -n
$ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info 2>&1 || { echo "ERMES failed, see $PROBLEM_NAME.info"; exit 1; }

# Move results files to a different folder
mv $PROBLEM_NAME.post.res Currents/$PROBLEM_NAME-J.post.res
mv $PROBLEM_NAME.post.msh Currents/$PROBLEM_NAME-J.post.msh 

# Delete previous files
rm -f $PROBLEM_NAME*dat Matrix_* Vector_* *.flavia.res *.info

# Calculate fields
$GID_PATH/gid_offscreen -b $MACROS_PATH/Calculate_E.bch -n
$ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info 2>&1 || { echo "ERMES failed, see $PROBLEM_NAME.info"; exit 1; }

# Solve with external solver
bash $SOLVER_PATH/ERMES2PETSc.sh >> $PROBLEM_PATH/$PROBLEM_NAME.info || { echo "PETSc solver failed, see $PROBLEM_NAME.info"; exit 1; }