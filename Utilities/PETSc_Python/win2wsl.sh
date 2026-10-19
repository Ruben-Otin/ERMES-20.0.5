#!/bin/bash

###################################################################################################
#
# - This script solves ERMES linear systems with PETSc from Windows WSL.
#
# - Check "README.txt" file provided in "/External_Solvers/PETSc" folder to install PETSc.
#
# - To call this script from ERMES, follow these steps: 
#
#   1) Select "External solver" on ERMES "Solver type".
#
#   2) Type "wsl ../PETSc/win2wsl.sh" on "External solvers settings".
# 
#   3) PETSc solver info will be saved on the ERMES "*.info" file.
#
# - To save the solver info in a different file use one of these two options: 
#
#   1) "wsl     ../PETSc/win2wsl.sh     >> SolverInfo.dat" in ERMES GUI. 
#
#   2) "python3 ../PETSc/ERMES2PETSc.py >> SolverInfo.dat" in this script.
#
# - Note: if end_of_line error appears, use in WSL: "dos2unix win2wsl.sh".  
#
###################################################################################################

export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-complex-M1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib

python3 ../PETSc/ERMES2PETSc.py