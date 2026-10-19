###############################################################################
# - This script calculates the element stiffness matrix for each tetrahedral 
#   element in the mesh for the RME_1st element type.
#
# - The matrix is saved in the binary file "Matrix_K_IVEM.bin".
###############################################################################

import IVEM_Utils as iu
import numpy      as np
import re

from pathlib import Path

###############################################################################
# Material properties
###############################################################################

mat = iu.Material();

mat.freq  = 3e8;                   

mat.sgm_r = 0.0; 
mat.sgm_i = 0.0;      

mat.epr_r = 1.0; 
mat.epr_i = 0.0;  
    
mat.mur_r = 1.0; 
mat.mur_i = 0.0;      

###############################################################################
# K matrix binary file name and path
###############################################################################

IVEM_K_FileName = '../Matrix_K_IVEM.bin'

###############################################################################
# Calculate K matrix and save it in binary file.
###############################################################################

# Volume elements list file
EVolFile = next( Path("..").glob("*-4.dat"), None );

# Volume elements list file format
patve = re.compile(r"VE\(\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)\s*\);");
  
# List of GiD mesh nodes:  No[ID] = p(x,y,z);
nIdl, nXl, nYl, nZl = iu.Read_Nodes_File();

# Element nodes Ids and Cartesian coordinates vectors
nId = np.zeros( 4, dtype=np.int64   );
nXe = np.zeros( 4, dtype=np.float64 );
nYe = np.zeros( 4, dtype=np.float64 ); 
nZe = np.zeros( 4, dtype=np.float64 );

# Open IVEM_K_FileName and write K matrix in binary
Kfile = open( IVEM_K_FileName, 'wb' );

# Read volume element list and calculate K matrix for each element
with open( EVolFile, "r" ) as f:
    
    for line in f:       
        m = patve.search( line );
        if m:
            nId[0] = int( m.group(1) );
            nId[1] = int( m.group(2) );
            nId[2] = int( m.group(3) );
            nId[3] = int( m.group(4) );
            eMatId = int( m.group(5) );
            for i in range( 4 ):    
                nXe[i] = nXl[ nId[i]-1 ]; 
                nYe[i] = nYl[ nId[i]-1 ]; 
                nZe[i] = nZl[ nId[i]-1 ];    
             
            # Calculate element stiffness matrix
            cK = iu.Get_K_Matrix_RME_1st( mat, nXe, nYe, nZe );
            
            # Write K matrix in binary files 
            np.array( cK, dtype='<c16' ).tofile( Kfile );
         
Kfile.close();

###############################################################################