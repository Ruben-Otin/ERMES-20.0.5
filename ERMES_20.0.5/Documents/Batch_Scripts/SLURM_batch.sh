#!/bin/bash
###################################################################################################
# - Example of SLURM batch file for running ERMES on a Linux bash HPC cluster.
###################################################################################################

#SBATCH -J JOB_NAME

#SBATCH -A ACCOUNT_NAME
#SBATCH -p MACHINE_NAME

#SBATCH --nodes=5
#SBATCH --ntasks=380
#SBATCH --exclusive

#SBATCH --time=24:00:00

#SBATCH --mail-type=ALL
#SBATCH --mail-user=user.mail@mail.com

###################################################################################################

# Load required system modules
source ~/util_modules.sh

export OMP_NUM_THREADS=1

# Set PETSc
export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-complex
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib

# Set ERMES executable path
ERMES_EXE=$HOME/GiDx64/gid-15.0.4/problemtypes/ERMES_20.0.5/ERMES.gid/execs/unix/ERMES_20.0.5

###################################################################################################

# Set active project path
WORK_PATH=$HOME/projects/active_project

PROBLEM_NAME=GiDProjectName

PROBLEM_PATH=$WORK_PATH/$PROBLEM_NAME.gid

cd $PROBLEM_PATH

###################################################################################################

# Run first case
CASE=M1

RESULTS_FOLDER=$WORK_PATH/Results

mkdir $RESULTS_FOLDER

FREQN=("f001" "f01" "f1")
IMAGS=("0.06283185307179587" "0.6283185307179586" "6.283185307179586")

for i in ${!FREQN[@]}; 
do
  PFREQ="ProblemFrequency = "
  PFREQ+=${IMAGS[$i]}
  PFREQ+=";"
  
  echo $PFREQ > $PROBLEM_PATH/$PROBLEM_NAME-16.dat
  
  # ERMES returns a non-zero status on error: report it and continue with the next case
  $ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info 2>&1 || echo "ERMES failed for case $CASE ${FREQN[$i]} ( exit status $? )"

  mv -f $PROBLEM_NAME.info       $RESULTS_FOLDER/$PROBLEM_NAME-$CASE-${FREQN[$i]}.info
  mv -f $PROBLEM_NAME.flavia.res $RESULTS_FOLDER/$PROBLEM_NAME-$CASE-${FREQN[$i]}.flavia.res  
  mv -f Fields/Volume-1-B.dat    $RESULTS_FOLDER/$PROBLEM_NAME-$CASE-Vol1-B-${FREQN[$i]}.dat
  
  rm -rf Integrals/ Fields/
done

###################################################################################################

# Run second case
CASE=M2

INFOS=$WORK_PATH/Infos
FLAVS=$WORK_PATH/Flavias
MODBS=$WORK_PATH/MOD-B

mkdir $INFOS
mkdir $FLAVS
mkdir $MODBS

FREQN=( "m2e6i0"    "m2e5i0"   "m2e4i0" )
REALS=( "-2000000"  "-200000"  "-20000" )
IMAGS=( "0.0"       "0.0"      "0.0"    )

for i in ${!FREQN[@]}; 
do
  echo ${REALS[$i]} >  $PROBLEM_PATH/$PROBLEM_NAME-17.dat
  echo ${IMAGS[$i]} >> $PROBLEM_PATH/$PROBLEM_NAME-17.dat
  # ERMES returns a non-zero status on error: report it and continue with the next case
  $ERMES_EXE $PROBLEM_NAME > $PROBLEM_PATH/$PROBLEM_NAME.info 2>&1 || echo "ERMES failed for case $CASE ${FREQN[$i]} ( exit status $? )"

  mv -f $PROBLEM_NAME.info       $INFOS/$PROBLEM_NAME-$CASE-${FREQN[$i]}.info
  mv -f $PROBLEM_NAME.flavia.res $FLAVS/$PROBLEM_NAME-$CASE-${FREQN[$i]}.flavia.res
  mv -f Fields/Surface-1-B.dat   $MODBS/$PROBLEM_NAME-$CASE-B-${FREQN[$i]}.dat
  
  rm -rf Integrals/ Fields/
done

###################################################################################################