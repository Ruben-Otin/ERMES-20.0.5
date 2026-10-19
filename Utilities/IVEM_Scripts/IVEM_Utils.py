###############################################################################
# - This script defines functions and objects used in IVEM_K_Matrix.py.
###############################################################################

import numpy as np
import re

from pathlib import Path

###############################################################################

class Material:  
    freq  = 0.0;                   # Frequency
    sgm_r = 0.0; sgm_i = 0.0;      # Conductivity 
    epr_r = 1.0; epr_i = 0.0;      # Relative permittivity
    mur_r = 1.0; mur_i = 0.0;      # Relative permeability
    
###############################################################################
    
def Read_Nodes_File():
    """
    - Reads nodes Ids and nodes coordinates from *-1.dat file.
    - Returns the vector lists: nIdl, nXl, nYl, nZl.
    """
    NodesFile = next( Path("..").glob("*-1.dat"), None );
    
    nIdl, nXl, nYl, nZl = [], [], [], [];
    
    pat = re.compile(r"No\[(\d+)\]\s*=\s*p\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^)]+)\s*\)\s*;");
    
    with open( NodesFile, 'r' ) as f:
        for line in f:
            m = pat.search( line );
            if m:
                nIdl.append( int  ( m.group(1) ) );
                nXl .append( float( m.group(2) ) );
                nYl .append( float( m.group(3) ) );
                nZl .append( float( m.group(4) ) );
                
    nIdl = np.asarray( nIdl, dtype=np.int64   );
    nXl  = np.asarray( nXl , dtype=np.float64 );
    nYl  = np.asarray( nYl , dtype=np.float64 );
    nZl  = np.asarray( nZl , dtype=np.float64 );

    return nIdl, nXl, nYl, nZl;
    
###############################################################################    

def Calculate_DN_Volume( nX, nY, nZ ):
    """
    - Returns derivatives and volume of a tetrahedron with nodal
      Cartesian coordinates nX, nY, nZ.
    """
    def X(i): return nX[i-1];
    def Y(i): return nY[i-1];
    def Z(i): return nZ[i-1];
    
    DN = np.zeros( ( 3, 4 ), dtype=np.float64 );
    
    Volume = ( X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - 
               X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) -
               X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + 
               X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) +
               X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - 
               X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) -
               X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + 
               X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3) );
    
    Volume = abs( Volume / 6.0 );
    
    Cte = 1.0 / ( 6.0 * Volume ); 

    DN[ 0 ][ 0 ] = Cte * ( Y(3)*Z(2) + Y(4)*Z(3) + Y(2)*Z(4) - 
                           Y(3)*Z(4) - Y(2)*Z(3) - Y(4)*Z(2) );
    
    DN[ 0 ][ 1 ] = Cte * ( Y(3)*Z(4) + Y(4)*Z(1) + Y(1)*Z(3) - 
                           Y(3)*Z(1) - Y(1)*Z(4) - Y(4)*Z(3) );
    
    DN[ 0 ][ 2 ] = Cte * ( Y(2)*Z(1) + Y(1)*Z(4) + Y(4)*Z(2) - 
                           Y(2)*Z(4) - Y(4)*Z(1) - Y(1)*Z(2) );
    
    DN[ 0 ][ 3 ] = Cte * ( Y(2)*Z(3) + Y(3)*Z(1) + Y(1)*Z(2) - 
                           Y(2)*Z(1) - Y(3)*Z(2) - Y(1)*Z(3) );
    
    DN[ 1 ][ 0 ] = Cte * ( X(3)*Z(4) + X(4)*Z(2) + X(2)*Z(3) - 
                           X(3)*Z(2) - X(2)*Z(4) - X(4)*Z(3) ); 
    
    DN[ 1 ][ 1 ] = Cte * ( X(3)*Z(1) + X(1)*Z(4) + X(4)*Z(3) - 
                           X(1)*Z(3) - X(3)*Z(4) - X(4)*Z(1) );
    
    DN[ 1 ][ 2 ] = Cte * ( X(2)*Z(4) + X(4)*Z(1) + X(1)*Z(2) - 
                           X(2)*Z(1) - X(4)*Z(2) - X(1)*Z(4) );
    
    DN[ 1 ][ 3 ] = Cte * ( X(1)*Z(3) + X(3)*Z(2) + X(2)*Z(1) - 
                           X(3)*Z(1) - X(1)*Z(2) - X(2)*Z(3) );
    
    DN[ 2 ][ 0 ] = Cte * ( X(4)*Y(3) + X(2)*Y(4) + X(3)*Y(2) - 
                           X(2)*Y(3) - X(4)*Y(2) - X(3)*Y(4) );
    
    DN[ 2 ][ 1 ] = Cte * ( X(3)*Y(4) + X(4)*Y(1) + X(1)*Y(3) - 
                           X(3)*Y(1) - X(4)*Y(3) - X(1)*Y(4) ); 
    
    DN[ 2 ][ 2 ] = Cte * ( X(2)*Y(1) + X(4)*Y(2) + X(1)*Y(4) - 
                           X(4)*Y(1) - X(2)*Y(4) - X(1)*Y(2) );
    
    DN[ 2 ][ 3 ] = Cte * ( X(2)*Y(3) + X(3)*Y(1) + X(1)*Y(2) - 
                           X(2)*Y(1) - X(3)*Y(2) - X(1)*Y(3) );
    
    return DN, Volume;
    
###############################################################################

def Get_K_Matrix_RME_1st( mat, nXe, nYe, nZe ):
    """
    - Returns RME_1st element stiffness matrix from a tetrahedron with nodal 
      Cartesian coordinates nX, nY, nZ and material properties mat.
    """
    wf = 2.0 * np.pi * mat.freq;
    w2 = wf * wf;
    
    eo = 8.8541878176e-12;
    mo = np.pi * 4.0e-7;    

    sg_real = mat.sgm_r;
    sg_imag = mat.sgm_i;
    
    ep_real = mat.epr_r * eo;
    ep_imag = mat.epr_i * eo;
    
    mu_real = mat.mur_r * mo;
    mu_imag = mat.mur_i * mo; 
    
    w2cEps = w2*ep_real - wf*sg_imag + 1j * ( w2*ep_imag + wf*sg_real );
    
    DN, mVolume = Calculate_DN_Volume( nXe, nYe, nZe ); 
    
    cMu       = mu_real + 1j * mu_imag;
    w2cEpsVol = mVolume * w2cEps;
    cMuInvVol = mVolume / cMu;
    
    NiNj_Matrix = np.zeros( ( 4, 4 ), dtype=np.complex128 );

    for i in range( 4 ):
        for j in range( 4 ):
            if( i!=j ): NiNj_Matrix[ i ][ j ] = w2cEpsVol / 20.0;
            else      : NiNj_Matrix[ i ][ j ] = w2cEpsVol / 10.0;	
				
    StiffMatrix = np.zeros( ( 12, 12 ), dtype=np.complex128 );
    
    for i in range( 4 ):
        for j in range( 4 ):

            eK_ij = cMuInvVol * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + 
                                  DN[ 1 ][ i ] * DN[ 1 ][ j ] + 
                                  DN[ 2 ][ i ] * DN[ 2 ][ j ] );
            
            eK_ij -= NiNj_Matrix[ i ][ j ]
    
            StiffMatrix[ i     ][ j     ] = eK_ij;
            StiffMatrix[ i + 4 ][ j + 4 ] = eK_ij;
            StiffMatrix[ i + 8 ][ j + 8 ] = eK_ij;
    
            eK_ij = cMuInvVol * ( DN[ 0 ][ i ] * DN[ 1 ][ j ] - 
                                  DN[ 1 ][ i ] * DN[ 0 ][ j ] );
    
            StiffMatrix[ i     ][ j + 4 ] =  eK_ij;
            StiffMatrix[ i + 4 ][ j     ] = -eK_ij;
    
            eK_ij = cMuInvVol * ( DN[ 0 ][ i ] * DN[ 2 ][ j ] - 
                                  DN[ 2 ][ i ] * DN[ 0 ][ j ] );             
    
            StiffMatrix[ i     ][ j + 8 ] =  eK_ij;
            StiffMatrix[ i + 8 ][ j     ] = -eK_ij;
    
            eK_ij = cMuInvVol * ( DN[ 1 ][ i ] * DN[ 2 ][ j ] - 
                                  DN[ 2 ][ i ] * DN[ 1 ][ j ] );
        
            StiffMatrix[ i + 4 ][ j + 8 ] =  eK_ij;
            StiffMatrix[ i + 8 ][ j + 4 ] = -eK_ij;
 
    return StiffMatrix;

###############################################################################