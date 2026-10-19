##################################################################################################################
# 
# - Example of use of Python SciPy BiCG ( iterative solver, Jacobi preconditioner ) as an external solver of ERMES.
# 
# - Before using this example, install Python 3 with NumPy and SciPy and add "python" to the system PATH.
# 
# - In the ERMES GiD interface select: "Solver" -> "Solver type" = "External_solver" and write in
#   "External_solvers_settings":
#
#   1) Solver_executable = python
#   2) Solver_parameters = {python_script_path}\BiCG.py ( or "/" in Linux )
#
# - The matrix storage format ( symmetric, full, Hermitic-symmetric, Hermitic-full ) is detected automatically.
# - Works with old and new SciPy versions ( the "tol" argument was renamed "rtol" in SciPy 1.12 and removed in 1.14 ).
#
##################################################################################################################

# Required imports
import numpy as np
import scipy.sparse.linalg 
import scipy.sparse
import sys
import os

# GiD problem path ( ERMES runs the external solver from the problem folder )
Folder = os.getcwd();

def FilePath( name ):
    return os.path.join( Folder, name );

def ReadPart( IntFile, CmplxFile, N ):
    # ERMES format: (row, col) int32 pairs ( 1-based ) and (real, imag) float64 pairs
    Idx = np.fromfile( FilePath( IntFile   ), dtype=np.dtype( ( np.int32  , 2 ) ) );
    Val = np.fromfile( FilePath( CmplxFile ), dtype=np.dtype( ( np.float64, 2 ) ) );
    Val = Val[:,0] + 1j*Val[:,1];
    I = Idx[:,0] - 1; J = Idx[:,1] - 1;
    OneTriangle = bool( np.all( I <= J ) or np.all( I >= J ) );
    M = scipy.sparse.csr_matrix( ( Val, ( I, J ) ), shape=( N, N ) );
    return M, OneTriangle;

# Right hand side
B_vector = np.fromfile( FilePath( 'Vector_B.bin' ), dtype=np.dtype( ( np.float64, 2 ) ) );
B_vector = B_vector[:,0] + 1j*B_vector[:,1];
N = B_vector.size;

# System matrix. The storage format is detected from the files, as in ERMESPETScSolver
# ( no manual editing needed ):
#   Symmetric           : no aux files, one triangle    A = U + (U - D)^T
#   Full-matrix         : no aux files, both triangles  A = M
#   Hermitic-Symmetric  : aux files, aux one triangle   A = U + (U - D)^H + R + (R - D_R)^T
#   Hermitic-Full       : aux files, aux both triangles A = U + (U - D)^H + R
U, MainOneTriangle = ReadPart( 'Matrix_A_int.bin', 'Matrix_A_cmplx.bin', N );
D = scipy.sparse.diags( U.diagonal() );

if os.path.isfile( FilePath( 'Matrix_A_aux_int.bin' ) ):
    R, AuxOneTriangle = ReadPart( 'Matrix_A_aux_int.bin', 'Matrix_A_aux_cmplx.bin', N );
    A_matrix = U + ( U - D ).conj().transpose();
    if AuxOneTriangle:
        Format = 'Hermitic-Symmetric';
        A_matrix = A_matrix + R + ( R - scipy.sparse.diags( R.diagonal() ) ).transpose();
    else:
        Format = 'Hermitic-Full';
        A_matrix = A_matrix + R;
elif MainOneTriangle:
    Format   = 'Symmetric';
    A_matrix = U + ( U - D ).transpose();
else:
    Format   = 'Full-matrix';
    A_matrix = U;

A_matrix = scipy.sparse.csc_matrix( A_matrix );

# Relative tolerance and maximum number of iterations
Tolerance     = 1e-6;
MaxIterations = 100000;

# Jacobi preconditioner
Diag = A_matrix.diagonal();
Diag[ np.abs( Diag ) == 0 ] = 1.0;
Mid  = scipy.sparse.diags( 1.0 / Diag );

# BiCG solver
try:
    X_vector, info = scipy.sparse.linalg.bicg( A_matrix, B_vector, rtol=Tolerance, maxiter=MaxIterations, M=Mid );
except TypeError:
    X_vector, info = scipy.sparse.linalg.bicg( A_matrix, B_vector, tol=Tolerance, maxiter=MaxIterations, M=Mid );

Residual = np.linalg.norm( A_matrix @ X_vector - B_vector ) / np.linalg.norm( B_vector );

# Showing solver info
print( '------------------------------------------------------------'               );
print( '- Python BiConjugate Gradient solver ( matrix format: %s, N = %d )' % ( Format, N ) );
print( '- Residual: %r' % Residual                                                  );
print( '------------------------------------------------------------'               );

if ( info != 0 ) or not np.isfinite( Residual ):
    # No solution file is written: ERMES stops with an error instead of using a wrong solution
    sys.exit( "- Iterative solver did not converge ( info = %d ). No solution written." % info );

# Writing solution vector on file ( complex128, little endian, as read by ERMES )
np.asarray( X_vector, dtype='<c16' ).tofile( FilePath( 'Vector_Xo.bin' ) );
