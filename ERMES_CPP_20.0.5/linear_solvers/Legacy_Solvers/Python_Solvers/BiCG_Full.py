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

# Diagonal preconditioner
invD=[]
for i in range(0,len(A_matrix.diagonal())): 
    if (abs(A_matrix[i,i]) > 0) : invD.append(1.0/A_matrix[i,i])
    else                        : invD.append(1.0       )  
Mid = scipy.sparse.diags(np.array(invD))  

# BiCG solver
x_complex, info = scipy.sparse.linalg.bicg(A_matrix, b_vector, tol=1e-5, M=Mid)

if info != 0: raise RuntimeError("iterative solver did not converge")

# Showing sover info
print('Python BiCG solver ends')
print('Problem solved: ' + Folder)
print('Residual: %r' % (np.linalg.norm(A_matrix * x_complex - b_vector)/np.linalg.norm(b_vector)))

# Writing solution vector on file
newFile = open(Folder + r'\Vector_Xo.dat', "wb")
newFile.write(x_complex)
newFile.close()


