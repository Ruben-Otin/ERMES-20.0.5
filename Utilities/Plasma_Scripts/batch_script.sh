#!/bin/bash

# SLURM job: builds the model with problem_generator.py and runs ERMES

#SBATCH -J MAST-U

#SBATCH -p cpu

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --mem=330GB 

#SBATCH --output=%j_output.txt 
#SBATCH --error=%j_error.txt   

#SBATCH --time=36:00:00

###################################################################################################

source ~/.bashrc

export SLURM_NTASKS=256

GID_PATH=/lustre/home/rotin/GiDx64/gid-15.0.4

ERMES_EXE=$GID_PATH/problemtypes/ERMES_20.0.5/ERMES.gid/execs/unix/ERMES_20.0.5

###################################################################################################

WORK_PATH=$HOME/projects/MAST-U

cd $WORK_PATH

# problem_generator.py: geometry + mesh (GiD) and plasma files only;
# ERMES is run by this script, below (RUN_ERMES=0)
export GID_EXE=$GID_PATH/gid_offscreen
export ERMES_EXE
export RUN_ERMES=0

python3 problem_generator.py || { echo "problem_generator.py failed"; exit 1; }

###################################################################################################

# Must match Problem_Name in problem_generator.py ('MASTU-2D-t<tor>_p<pol>')
PROBLEM_NAME=MASTU-2D-t20_p20

PROBLEM_PATH=$WORK_PATH/$PROBLEM_NAME.gid

cd $PROBLEM_PATH

# ERMES returns a non-zero exit status on error; stderr also goes to the .info file
$ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info 2>&1 || { echo "ERMES failed, see $PROBLEM_NAME.info"; exit 1; }

###################################################################################################
