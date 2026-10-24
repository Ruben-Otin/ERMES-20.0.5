###############################################################################
# permittivity_tensor.py
#
# - Computes the cold plasma permittivity tensor K ( relative, in the local
#   frame with z along B ) at every node of the ERMES plasma files and writes
#   the binary file TENSOR_FILE in the format expected by ERMES ( "Tensor"
#   geometry extrapolation ):
#
#     * one 3x3 complex tensor ( row major, little endian complex128 ) per data
#       line of the electron density file, in the same order as that file;
#     * ERMES rotates K to the global frame with the B field of each node.
#
# - Input: the ne and Be files written by eqdsk_to_ermes_turb.py,
#   eqdsk_to_ermes_simp.py or torbeam_to_ermes.py ( "node ne" and
#   "node Bx By Bz" ). B is matched to ne by node Id.
#
# - Conventions identical to ERMES ( ColdPlasma.cpp ):
#     * time dependence exp(-i w t), Stix tensor  K = [[ S, -iD, 0 ],
#                                                      [ iD,  S, 0 ],
#                                                      [ 0 ,  0, P ]]
#     * collisions: w -> w + i*nu in the electron terms ( DAMPING = nu );
#     * CODATA 2018 constants; ion mass = atomic mass - Z * electron mass.
#
# - Changes in ERMES 20.0.5 ( with respect to ERMES 20.0.4 ):
#     * K_xy = -iD ( was +iD: opposite handedness to the tensor computed inside
#       ERMES, i.e. results as if B were reversed );
#     * text / comment lines of the plasma files are skipped, with the same rule
#       as ERMES;
#     * the B field is matched to the density by node Id ( previously by line
#       position );
#     * constants and masses updated as in ERMES.
###############################################################################

import numpy as np

###############################################################################
# USER PARAMETERS
###############################################################################

FREQUENCY = 28e9;      # Problem frequency [Hz]
DAMPING   = 1e9;       # Electron collision frequency [rad/s] ( CPlasma.Damping
                       # in problem_generator.py )

# Relative densities of the ion species ( CPlasma.H, D, T in problem_generator.py )
N_H = 0.0;
N_D = 1.0;
N_T = 0.0;

NE_FILE     = r'./Plasma/ne_2D-t20_p20.dat';   # node  ne [m^-3]
BE_FILE     = r'./Plasma/Be_2D-t20_p20.dat';   # node  Bx By Bz [T]
TENSOR_FILE = r'./Tensor_K_CPGE.bin';          # output

###############################################################################
# Constants ( CODATA 2018 )
###############################################################################

eo  = 8.8541878128e-12;   # Vacuum electric permittivity
Qe  = 1.602176634e-19;    # Elementary charge
Me  = 9.1093837015e-31;   # Electron mass
uKg = 1.66053906660e-27;  # Atomic mass unit (kg)

def IonMass( AtomicMass_u, Z ):
    return AtomicMass_u * uKg - Z * Me;

MH1 = IonMass( 1.00782503223, 1 );  # H+
MH2 = IonMass( 2.01410177812, 1 );  # D+
MH3 = IonMass( 3.01604928199, 1 );  # T+

Wf = 2 * np.pi * FREQUENCY;         # Angular frequency

###############################################################################
# Plasma species: charge, mass, relative density, (complex) frequency
###############################################################################

Qs = np.array( [ -Qe, Qe , Qe , Qe  ] );
Ms = np.array( [  Me, MH1, MH2, MH3 ] );
Nr = np.array( [ 1.0, N_H, N_D, N_T ] );
Wd = np.array( [ Wf + 1j * DAMPING, Wf, Wf, Wf ] );   # damping on electrons only

###############################################################################
# Stix tensor
###############################################################################

def Get_Stix_Tensor( ne, modB ):
    """
    - Cold plasma tensor K ( local frame, z along B ) for electron density
      ne [m^-3] and |B| [T].
    """
    R = 1.0 + 0j;
    L = 1.0 + 0j;
    P = 1.0 + 0j;
    
    for i in range( len( Qs ) ):

        if Nr[i] == 0.0: continue;
        
        W2_ps   = Qs[i] * Qs[i] * ne * Nr[i] / ( eo * Ms[i] );   # plasma frequency^2
        Omega_s = Qs[i] * modB / Ms[i];                         # signed cyclotron frequency
        
        R -= W2_ps / ( Wf * ( Wd[i] + Omega_s ) );
        L -= W2_ps / ( Wf * ( Wd[i] - Omega_s ) );
        P -= W2_ps / ( Wf * ( Wd[i]           ) );
    
    S = 0.5 * ( R + L );
    D = 0.5 * ( R - L );
    
    return np.array( [ [  S    , -1j * D, 0.0 ],
                       [ 1j * D,  S     , 0.0 ],
                       [ 0.0   ,  0.0   , P   ] ], dtype=np.complex128 );

###############################################################################
# Plasma files
###############################################################################

def Is_Data_Line( line ):
    """ Same rule as ERMES: the first non-blank character starts a number """
    s = line.lstrip();
    if not s: return False;
    if s[0].isdigit(): return True;
    return ( s[0] in '+-.' ) and len( s ) > 1 and ( s[1].isdigit() or s[1] == '.' );

def main():
    
    # Electron density ( node Id, ne ), in file order
    ne_Ids = [];
    ne     = [];
    with open( NE_FILE, "r" ) as f:
        for line in f:
            if not Is_Data_Line( line ): continue;
            lsplit = line.split();
            ne_Ids.append( int  ( lsplit[0] ) );
            ne    .append( float( lsplit[1] ) );
    
    # |B| per node Id
    modB = {};
    with open( BE_FILE, "r" ) as f:
        for line in f:
            if not Is_Data_Line( line ): continue;
            lsplit = line.split();
            Bx, By, Bz = ( float( v ) for v in lsplit[1:4] );
            modB[ int( lsplit[0] ) ] = np.sqrt( Bx**2 + By**2 + Bz**2 );
    
    missing = [ n for n in ne_Ids if n not in modB ];
    if missing:
        raise SystemExit( 'Error: %d nodes of %s have no B field in %s ( first: %d )'
                          % ( len( missing ), NE_FILE, BE_FILE, missing[0] ) );
    
    # Write the tensors ( same order as the data lines of the ne file )
    with open( TENSOR_FILE, 'wb' ) as Kfile:
        for nId, n_e in zip( ne_Ids, ne ):
            Get_Stix_Tensor( n_e, modB[ nId ] ).astype( '<c16' ).tofile( Kfile );
    
    print( 'Written %s ( %d tensors )' % ( TENSOR_FILE, len( ne_Ids ) ) );

if __name__ == "__main__":
    main();

###############################################################################
