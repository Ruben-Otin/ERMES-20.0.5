#################################################################################################################
# 
# - Example of use of an external solver with ERMES 8.0
# 
# - Before using this example write the python.exe path inside GiD: 
#   In GiD's "Solving Parameters" window:
#   - Solver path      =  F:\Python\WinPython-64bit-3.4.4.4Qt5\python-3.4.4.amd64\python.exe
#   - Input parameters = "F:/GiD/GiD 13.0.3/problemtypes/ERMES-ICRH-v8.0/ExternalSolvers/DirectSolver.py"
#
#################################################################################################################

import numpy as np
import scipy.sparse.linalg 
import scipy.sparse

# GiD problem path
import os
Folder = os.getcwd()

# Reading matrix  
b_vector = np.array(np.fromfile(Folder + r'\Vector_b.bin'      ,dtype=np.dtype((np.float64,2))));
A_matrix = np.array(np.fromfile(Folder + r'\Matrix_A_cmplx.bin',dtype=np.dtype((np.float64,2))));
A_indexs = np.array(np.fromfile(Folder + r'\Matrix_A_int.bin'  ,dtype=np.dtype((int       ,2))));

# Adapting formats
b_vector = b_vector[:,0] + 1j*b_vector[:,1]
A_matrix = A_matrix[:,0] + 1j*A_matrix[:,1]

A_matrix = scipy.sparse.csr_matrix((A_matrix, (A_indexs[:,0]-1, A_indexs[:,1]-1))) 

# Filling matrix lower diagonal (for symmetric matrices)
A_matrix = A_matrix + A_matrix.transpose() - scipy.sparse.diags(A_matrix.diagonal()) 
 
# Filling matrix lower diagonal (for Hermitian matrices)
#A = A + np.conjugate(A.transpose()) - scipy.sparse.diags(A.diagonal())  

# Solving system
x_complex = scipy.sparse.linalg.spsolve(A_matrix, b_vector) 

# Showing sover info
print('Python direct solver for symmetric matrices ends')
print('Problem solved: ' + Folder)
print('Residual: %r' % (np.linalg.norm(A_matrix * x_complex - b_vector)/np.linalg.norm(b_vector)))

# Writing solution vector on file
newFile = open(Folder + r'\Vector_Xo.dat', "wb")
newFile.write(x_complex)
newFile.close()



