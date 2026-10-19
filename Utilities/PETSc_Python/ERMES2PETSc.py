######################################################################################################################################
# 
# - This script is the interface of ERMES with PETSc:
#
#   1) Converts ERMES linear system into PETSc format, 
#
#   2) Calls PETSc, and 
#
#   3) Converts back PETSc solution into ERMES format.
# 
# - Check "README.txt" file provided in "/External_Solvers/PETSc" folder to install PETSc.
#
# - To call this script from ERMES, follow these steps:
#    
#   1) Select "External solver" on ERMES "Solver type".
#
#   2) Type "python3 ../PETSc/ERMES2PETSc.py" on ERMES "External solvers settings".
# 
#   3) PETSc solver info will be saved on the ERMES "*.info" file.
#
# - To save the solver info in a different file use:
#
#   "python3 ../PETSc/ERMES2PETSc.py >> SolverInfo.dat" on ERMES "External solvers settings".
#
# - Before using it, costumize the values on the "Interface with PETSc" section, e.g.:
#
#   SolverType     = '-ksp_type richardson -pc_type lu -pc_factor_mat_solver_type mumps';
#   SolverOptions  = '-ksp_max_it 1 -ksp_monitor_true_residual -memory_view -log_view';
#   MatrixType     = 'Full_Matrix'; 
#   NumParallTasks = '16';
#   FolderPath     = os.getcwd() + '/';
#   SolverFullPath = FolderPath + '../PETSc/Solver/PETScSolver';
#
#   SolverType     = '-ksp_type lgmres -pc_type sor';
#   SolverOptions  = '-ksp_rtol 1e-4 -ksp_max_it 1000000 -ksp_converged_reason -log_view -memory_view -ksp_monitor_true_residual'; 
#   MatrixType     = 'Symmetric_Matrix'; 
#   NumParallTasks = '$SLURM_NTASKS';
#   FolderPath     = r'/mnt/d/Projects/Problem.gid/';
#   SolverFullPath = r'/mnt/d/Solvers/PETSc/Solver/PETScSolver';
#
# - Check PETSc documentation for more options and solvers.
#
######################################################################################################################################

# Importing required python modules
import Solver.PetscBinaryIO as PBIO
import numpy as np
import scipy.sparse
import subprocess as sp
import os

#############################################    Interface with PETSc    #############################################################

# Solver type and options ( check PETSc documentation )
SolverType    = '-ksp_type richardson -pc_type lu -pc_factor_mat_solver_type mumps';
SolverOptions = '-ksp_max_it 1 -ksp_monitor_true_residual -memory_view -log_view';

# Matrix format ( 'Symmetric_Matrix' | 'Full_Matrix' )
MatrixType = 'Full_Matrix'; 

# Number of parallel tasks ( '$SLURM_NTASKS' | 'int' e.g. '10' )
NumParallTasks = '10';

# Path to ERMES matrixes and vector files 
FolderPath = os.getcwd() + '/';

# Full path to the PETSc solver executable (PETScSolver) 
SolverFullPath = FolderPath + '../PETSc/Solver/PETScSolver';

######################################################################################################################################

# Displaying solver settings
print( '------------------------------------------------------------' );
print( '- PETSc  solver : ' + SolverType      );
print( '- Solver options: ' + SolverOptions   );
print( '- Matrix type   : ' + MatrixType      );
print( '- Processors    : ' + NumParallTasks  );
print( '                                    ' );
print( 'Converting ERMES linear system...   ' );

# Reading and converting B_vector from ERMES to PETSc format
B_vector = np.array( np.memmap( FolderPath + 'Vector_B.bin', dtype=np.dtype( ( np.float64, 2 ) ) ) );
B_vector = B_vector[:,0] + 1j * B_vector[:,1];

PBIO.PetscBinaryIO().writeBinaryFile( FolderPath + 'Vector_B_PETSc.dat', [ B_vector.view( PBIO.Vec ), ] );
del B_vector;

# Reading and converting A_matrix from ERMES to PETSc format
A_matrix = np.array( np.memmap( FolderPath + 'Matrix_A_cmplx.bin', dtype=np.dtype( ( np.float64, 2 ) ) ) );
A_indexs = np.array( np.memmap( FolderPath + 'Matrix_A_int.bin'  , dtype=np.dtype( ( np.int32  , 2 ) ) ) );

A_matrix = A_matrix[:,0] + 1j * A_matrix[:,1];
A_indexs = ( A_indexs[:,0]-1, A_indexs[:,1]-1 );
A_matrix = scipy.sparse.csr_matrix( ( A_matrix, ( A_indexs ) ) ); 
del A_indexs;

if ( MatrixType == 'Symmetric_Matrix' ): 
    A_matrix = A_matrix + A_matrix.transpose() - scipy.sparse.diags( A_matrix.diagonal() ); 

PBIO.PetscBinaryIO().writeBinaryFile( FolderPath + 'Matrix_A_PETSc.dat', [ A_matrix , ] );
del A_matrix;

print( 'Conversion ended.' );
print( '                 ' );
print( '------------------------------------------------------------' );
print( '                 ' );

######################################################################################################################################

# Calling PETSc
command = ( '$PETSC_DIR/lib/petsc/bin/petscmpiexec' + 
            ' -n '         + NumParallTasks         + 
            ' '            + SolverFullPath         + 
            ' '            + SolverType             + 
            ' '            + SolverOptions          + 
            ' -my_folder ' + FolderPath             );

process = sp.Popen( command, shell = True, stdout = sp.PIPE, stderr = sp.PIPE, encoding = 'utf8' );

while ( True ):
    output = process.stdout.readline();

    if ( ( output == '' ) and ( process.poll() is not None ) ):
        break;

    if ( output ):
        print( output.strip() );

######################################################################################################################################
   
# Reading PETSc solution X_vector and conversion to ERMES format.
X_vector = np.array( np.fromfile( FolderPath + 'Vector_X_PETSc.dat', dtype=np.dtype( ( np.float64, 2 ) ) ) );
X_vector = np.array( X_vector ).byteswap();

NewFile = open( FolderPath + 'Vector_Xo.bin', 'wb' );
X_vector.tofile( NewFile );
NewFile.close();

print( '------------------------------------------------------------' );
print( '                                                            ' );

######################################################################################################################################