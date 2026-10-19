###############################################################################
# - This script computes the cold plasma permittivity tensor K ( relative, in
#   the local frame with z along B ) and writes it to the binary file
#   "Tensor_K_CPGE.bin" in the format expected by ERMES ( "Tensor" geometry
#   extrapolation ):
#
#     * one 3x3 complex tensor ( row major, little endian complex128 ) per data
#       line of the 3D electron density file, in the same order as that file;
#     * ERMES rotates K to the global frame with the B field of each node.
#
# - Conventions identical to ERMES ( ColdPlasma.cpp ):
#     * time dependence exp(-i w t), Stix tensor  K = [[ S, -iD, 0 ],
#                                                      [ iD,  S, 0 ],
#                                                      [ 0 ,  0, P ]]
#     * CODATA 2018 constants; ion mass = atomic mass - Z * electron mass.
#
# - Changes in ERMES 20.0.5 ( with respect to ERMES 20.0.4 ):
#     * K_xy = -iD ( was +iD: opposite handedness to the tensor computed inside
#       ERMES, i.e. results as if B were reversed );
#     * text / comment lines of the plasma files are skipped, with the same rule
#       as ERMES ( previous version failed on the "#" headers of the samples );
#     * the B field is matched to the density by node Id ( previously by line
#       position, which required both files to list the nodes in the same order );
#     * constants and masses updated as in ERMES.
###############################################################################

import numpy as np

###############################################################################
# Plasma parameters
###############################################################################

f   = 28e9     ;       # Problem frequency
Wf  = 2*np.pi*f;       # Problem frequency (angular)

eo  = 8.8541878128e-12;   # Vacuum electric permittivity  (CODATA 2018)
Qe  = 1.602176634e-19;    # Elementary charge
Me  = 9.1093837015e-31;   # Electron mass
uKg = 1.66053906660e-27;  # Atomic mass unit (kg)

def IonMass( AtomicMass_u, Z ):
    return AtomicMass_u * uKg - Z * Me;

MH1 = IonMass( 1.00782503223, 1 );  # H+
MH2 = IonMass( 2.01410177812, 1 );  # D+
MH3 = IonMass( 3.01604928199, 1 );  # T+

###############################################################################
# Plasma species
###############################################################################
   
# Number of species in the plasma   
NumSpecies = 4;
    
# Charge, mass, and relative density
Qs = np.zeros( NumSpecies, dtype=np.dtype( np.float64 ) );  
Ms = np.zeros( NumSpecies, dtype=np.dtype( np.float64 ) ); 
Nr = np.zeros( NumSpecies, dtype=np.dtype( np.float64 ) ); 

Qs[ 0 ] =-Qe; Ms[ 0 ] = Me ; Nr[ 0 ] = 1.0; # Electron
Qs[ 1 ] = Qe; Ms[ 1 ] = MH1; Nr[ 1 ] = 0.0; # Hidrogen
Qs[ 2 ] = Qe; Ms[ 2 ] = MH2; Nr[ 2 ] = 1.0; # Deuterium
Qs[ 3 ] = Qe; Ms[ 3 ] = MH3; Nr[ 3 ] = 0.0; # Tritium

###############################################################################
# Damping ( collision frequency in rad/s, w -> w + i*nu )
###############################################################################

eDamping = 1e9; 

Wd = np.zeros( NumSpecies, dtype=np.dtype( np.complex128 ) ); 

Wd[ 0 ] = Wf + 1j * eDamping;
Wd[ 1 ] = Wf;
Wd[ 2 ] = Wf;
Wd[ 3 ] = Wf;

###############################################################################
# Stix tensor components
###############################################################################

def Get_Stix_Tensor( ne, modB ):
    
    R = 1.0 + 1j * 0.0; 
    L = 1.0 + 1j * 0.0; 
    P = 1.0 + 1j * 0.0;
    
    for i in range( NumSpecies ):

        if Nr[i] == 0.0: continue;
        
        W2_ps   = Qs[i] * Qs[i] * ne * Nr[i] / ( eo * Ms[i] );
        Omega_s = Qs[i] * modB / Ms[i];
        
        R -= W2_ps / ( Wf * ( Wd[i] + Omega_s ) );
        L -= W2_ps / ( Wf * ( Wd[i] - Omega_s ) );
        P -= W2_ps / ( Wf * ( Wd[i]           ) );
    
    S = 0.5 * ( R + L );
    D = 0.5 * ( R - L );
    
    K = np.zeros( (3,3), dtype=np.dtype( np.complex128 ) );
    
    K[0][0] =  S   ; K[0][1] =-1j*D; K[0][2] = 0.0;
    K[1][0] =+1j*D ; K[1][1] =  S  ; K[1][2] = 0.0;
    K[2][0] = 0.0  ; K[2][1] = 0.0 ; K[2][2] = P  ;
    
    return K;

###############################################################################
# Plasma files
###############################################################################

ne_FileName = r'./Plasma/ne_3D.dat'
Be_FileName = r'./Plasma/Be_3D.dat'

def Is_Data_Line( line ):
    # Same rule as ERMES: the first non-blank character starts a number
    s = line.lstrip();
    if not s: return False;
    if s[0].isdigit(): return True;
    return ( s[0] in '+-.' ) and len( s ) > 1 and ( s[1].isdigit() or s[1] == '.' );

# Electron density ( node Id, ne ), in file order
ne_Ids = [];
ne     = [];
with open( ne_FileName, "r" ) as f:
    for line in f:
        if not Is_Data_Line( line ): continue;
        lsplit = line.split();
        ne_Ids.append( int  ( lsplit[0] ) );
        ne    .append( float( lsplit[1] ) );

# |B| per node Id
modB = {};
with open( Be_FileName, "r" ) as f:
    for line in f:
        if not Is_Data_Line( line ): continue;
        lsplit = line.split();
        Bx = float( lsplit[1] );
        By = float( lsplit[2] );
        Bz = float( lsplit[3] );
        modB[ int( lsplit[0] ) ] = np.sqrt( Bx**2 + By**2 + Bz**2 );

missing = [ n for n in ne_Ids if n not in modB ];
if missing:
    raise SystemExit( 'Error: %d nodes of %s have no B field in %s ( first: %d )' % ( len( missing ), ne_FileName, Be_FileName, missing[0] ) );

###############################################################################
# Write Tensor_K_CPGE.bin file ( same order as the data lines of ne file )
###############################################################################

Tensor_K_FileName = r'./Tensor_K_CPGE.bin'

with open( Tensor_K_FileName, 'wb' ) as Kfile:
    for nId, n_e in zip( ne_Ids, ne ):
        cK = Get_Stix_Tensor( n_e, modB[ nId ] );
        np.array( cK, dtype='<c16' ).tofile( Kfile );

print( 'Written %s ( %d tensors )' % ( Tensor_K_FileName, len( ne_Ids ) ) );

###############################################################################
