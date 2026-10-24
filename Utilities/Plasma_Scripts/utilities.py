###############################################################################
# utilities.py
#
# Helper functions for problem_generator.py and eqdsk_to_ermes_simp.py:
#
#   - geometry       : Rot3DAxis_rad/_dgr, Calculate_WaveVector,
#                      LineSphere_Intersect, Angle_Between_Vectors,
#                      Reflected_Vector, LineLine_Intersect,
#                      Calculate_GiD_Points
#   - GiD / ERMES    : GaussBeam, ColdPlasma, ProblemSettings (data classes),
#                      Write_GiD_bch_File_2D
#   - polarisation   : Calculate_O_Mode
#   - mesh / fields  : Read_GiD_Nodes, Transform_FRZT_2_FXYZ
#
# Frame: z is the tokamak axis, R = sqrt(x^2 + y^2), phi = atan2(y, x),
# (R, phi, Z) right-handed. All lengths in metres.
#
# The plasma-file functions of older versions (Read_Equilibrium_Files,
# Generate_Plasma_Files, Get_Point_Equilibrium_Data,
# Get_Pcrtc_Equilibrium_Data) are now in eqdsk_to_ermes_turb.py, with the
# same signatures.
###############################################################################

import re
import numpy as np

from scipy.constants import e as QE, m_e as ME

###############################################################################

def Rot3DAxis_rad(vector, axis, rad):
    """
    - Rotates a vector around an axis an angle given in radians.
    - Returns a vector rotated rad radians around the axis.
    """
    vector = np.asarray(vector, dtype=float);
    axis   = np.asarray(axis  , dtype=float);
    
    norm = np.linalg.norm(axis);
    axis = axis/norm;
    
    R = np.zeros((3,3));

    R[0][0] = axis[0]*axis[0] * (1.0-np.cos(rad)) +           np.cos(rad);
    R[0][1] = axis[0]*axis[1] * (1.0-np.cos(rad)) - axis[2] * np.sin(rad);
    R[0][2] = axis[0]*axis[2] * (1.0-np.cos(rad)) + axis[1] * np.sin(rad);
    
    R[1][0] = axis[1]*axis[0] * (1.0-np.cos(rad)) + axis[2] * np.sin(rad);
    R[1][1] = axis[1]*axis[1] * (1.0-np.cos(rad)) +           np.cos(rad);
    R[1][2] = axis[1]*axis[2] * (1.0-np.cos(rad)) - axis[0] * np.sin(rad);
    
    R[2][0] = axis[2]*axis[0] * (1.0-np.cos(rad)) - axis[1] * np.sin(rad);
    R[2][1] = axis[2]*axis[1] * (1.0-np.cos(rad)) + axis[0] * np.sin(rad);
    R[2][2] = axis[2]*axis[2] * (1.0-np.cos(rad)) +           np.cos(rad);

    return R.dot(vector);

###############################################################################

def Rot3DAxis_dgr(vector, axis, dgr):
    """
    - Rotates a vector around an axis an angle given in degrees.
    - Returns a vector rotated dgr degress around the axis.
    """
    return Rot3DAxis_rad(vector, axis, dgr*np.pi/180);

###############################################################################

def Calculate_WaveVector(tor_dgr, pol_dgr, MirrorCoord, PlasmaCentr):
    """
    - Calculates the normalized wave vector of the Gaussian beam.
    - Starts from the direction mirror -> plasma centre, rotates it tor_dgr
      degrees about the tokamak axis (z) and then pol_dgr degrees about the
      horizontal axis k0 x z_hat (pol_dgr > 0 tilts the beam towards +z).
    """
    MirrorCoord = np.asarray(MirrorCoord, dtype=float);
    PlasmaCentr = np.asarray(PlasmaCentr, dtype=float);    
    
    # Rotation around Z-axis (tokamak axis)
    vector = PlasmaCentr - MirrorCoord; 
    axis   = [0, 0, 1]; 
    angle  = tor_dgr;
    k0     = Rot3DAxis_dgr(vector, axis, angle);

    # Rotation around the horizontal axis k0 x z_hat (pol > 0 tilts the beam
    # towards +z, as in the previous y-axis version of this function)
    vector = k0; 
    axis   = [k0[1], -k0[0], 0]; 
    angle  = pol_dgr;
    k0     = Rot3DAxis_dgr(vector, axis, angle);

    # Normalization
    k0 = k0 / np.linalg.norm(k0);
    
    return k0;

###############################################################################

def LineSphere_Intersect(r0, rV, eR, eC):
    """
    - Computes the closest intersection point between a ray and a sphere.
      Ray   : r(t) = r0 + rV * t,  t >= 0
      Sphere: ||r - eC||^2 = eR^2
    - Returns: closest intersection point, or None if no intersection
    """
    rV = np.asarray(rV, dtype=float);
    r0 = np.asarray(r0, dtype=float);
    eC = np.asarray(eC, dtype=float);

    OC = r0 - eC;

    a = np.dot(rV, rV);
    b = 2.0 * np.dot(rV, OC);
    c = np.dot(OC, OC) - eR*eR;

    discrt = b*b - 4*a*c;

    if ( discrt < 0 ) : return None;

    sqrt_d = np.sqrt(discrt);

    t1 = (-b - sqrt_d) / (2 * a);
    t2 = (-b + sqrt_d) / (2 * a);
    
    t_values = [t for t in (t1, t2) if t >= 0];
    
    if ( not t_values ) : return None;
    
    t_closest = min(t_values);
    
    P_closest = np.array( r0 + rV * t_closest );
    
    return P_closest;
    
###############################################################################

def Angle_Between_Vectors(v1, v2, deg=True):
    """
    - Computes the angle between two vectors.
    - Returns the angle in degrees (if deg==True ).
    - Returns the angle in radians (if deg==False).
    """
    v1 = np.asarray(v1, dtype=float);
    v2 = np.asarray(v2, dtype=float);

    dot = np.dot(v1, v2);
    norm_product = np.linalg.norm(v1) * np.linalg.norm(v2);

    # Avoid numerical issues
    cos_theta = np.clip(dot / norm_product, -1.0, 1.0);
    
    angle_rad = np.arccos( cos_theta );
    
    if( deg == True ): angle = np.degrees( angle_rad );
    else             : angle = angle_rad;

    return angle;

###############################################################################

def Reflected_Vector(k0, p0, pC):
    """
    Input:
    - k0: incident vector
    - p0: point on sphere
    - pC: sphere center point
      
    Output:
    - kr: reflected vector normalized
    """
    k0 = np.asarray(k0, dtype=float);
    p0 = np.asarray(p0, dtype=float);
    pC = np.asarray(pC, dtype=float);
    
    # Exterior normal at p0
    n0 = p0 - pC; 
    n0 = n0 / np.linalg.norm(n0);

    # Incident angle
    angle = Angle_Between_Vectors(-k0, n0);
    
    # Rotation axis
    t0 = np.cross(-k0, n0);

    # Reflected vector
    kr = Rot3DAxis_dgr(-k0, t0, 2*angle);
    
    # Normalization
    kr = kr / np.linalg.norm(kr);

    return kr;

###############################################################################

def LineLine_Intersect(p1, v1, p2, v2, tol=1e-9):
    """
    Input:
    - p1, v1 : point and direction of line 1
    - p2, v2 : point and direction of line 2
    
    Output:
    - None  -> no intersection 
    - Point -> unique intersection point
    """
    p1 = np.array(p1, dtype=float);
    v1 = np.array(v1, dtype=float);
    p2 = np.array(p2, dtype=float);
    v2 = np.array(v2, dtype=float);

    # Check directions
    if ( np.linalg.norm(np.cross(v1, v2)) < tol ): return None;
    
    # Solve: [p1 + t*v1 = p2 + s*v2] == [t*v1 - s*v2 = p2 - p1]
    A = np.column_stack((v1, -v2));
    b = p2 - p1;

    # Least-squares solution
    t_s, residuals, _, _ = np.linalg.lstsq(A, b, rcond=None);
    t, s = t_s;

    # Check consistency
    Pint_1 = p1 + t * v1;
    Pint_2 = p2 + s * v2;

    if ( np.linalg.norm(Pint_1 - Pint_2) < tol ):
        return Pint_1;
    else:
        return None;

###############################################################################

def Calculate_GiD_Points(P0, k0,  Rpl, Cpl, Dpl,
                         WidthR, WidthL, Length, 
                         PdeepI, PdeepR, Vdispl):
    """
    Input:
    - P0 : mirror location
    - k0 : beam wave vector 
    
    - Rpl : plasma radius
    - Cpl : plasma center
    - Dpl : distance of the plasma to beam boundary
      (if Dpl == 0.0 is calculated automatically)
    
    - WidthR : right length from the beam center
    - WidthL : left  length from the beam center
    - Length : length of the receiving boundary
    
    - PdeepI : length inside the plasma (incident beam)
    - PdeepR : length inside the plasma (reflected beam)
    - Vdispl : extrusion distance
    
    Output:
    - P1,..., P5 : corner points of the GiD geometry
    - VV         : extrusion point
    """
    P0_sph = LineSphere_Intersect(P0, k0, Rpl, Cpl);
    
    if ( P0_sph is None ):
        raise ValueError("Calculate_GiD_Points: the beam does not hit the "
                         "plasma sphere (check MirrorCoord, angles, PRadius)");
    
    n0 = P0_sph - Cpl; 
    n0 = n0 / np.linalg.norm(n0);
    
    t0 = np.cross(k0, n0);
    t0 = t0 / np.linalg.norm(t0);
    
    kpi = np.cross(k0, t0);
    kpi = kpi / np.linalg.norm(kpi);
    
    P0r = P0 + kpi * WidthR;
    P0l = P0 - kpi * WidthL;   
    
    P0r_sph = LineSphere_Intersect(P0r, k0, Rpl, Cpl);
    
    if ( Dpl == 0.0 ):
        k0rr = Reflected_Vector(k0, P0r_sph, Cpl);
        k0rr = k0rr / np.linalg.norm(k0rr);
        P1   = LineLine_Intersect(P0r_sph, k0rr, P0l, k0);
    else:
        P0l_sph = LineSphere_Intersect(P0l, k0, Rpl, Cpl);
        P1      = P0l_sph - k0 * Dpl;
    
    P2 = LineLine_Intersect(P1, kpi, P0r, k0);
    P3 = P0r_sph + k0 * PdeepI;
    
    k0r = Reflected_Vector(k0, P0_sph, Cpl);
    k0r = k0r / np.linalg.norm(k0r);    
    
    t0r = np.cross(k0r, n0);
    t0r = t0r / np.linalg.norm(t0r);
   
    k1r = np.cross(k0r, t0r);
    k1r = k1r / np.linalg.norm(k1r);
 
    P5 = P1 + k1r * Length;
    
    P5_sph = LineSphere_Intersect(P5, -k0r, Rpl, Cpl);

    P4 = P5_sph - k0r * PdeepR;
    
    VV = P1 + Vdispl * t0;
     
    return P1, P2, P3, P4, P5, VV;

###############################################################################

class GaussBeam:
    """
    Gaussian beam data written by Write_GiD_bch_File_2D:
    - modEpar, phaEpar : amplitude and phase of E along Ep
    - modEper, phaEper : amplitude and phase of E along k0 x Ep
    - Rpar, Rper       : waist radii [m]
    - OO_*             : waist position [m]
    - Ep_*             : polarization vector
    - k0_*             : propagation direction (normalized)
    """
    modEpar = 0.0; 
    phaEpar = 0.0;
    
    modEper = 0.0; 
    phaEper = 0.0;
    
    Rpar = 0.0; 
    Rper = 0.0;
    
    OO_X = 0.0; 
    OO_Y = 0.0; 
    OO_Z = 0.0;
    
    Ep_X = 0.0;
    Ep_Y = 0.0;
    Ep_Z = 0.0;
    
    k0_X = 0.0;
    k0_Y = 0.0;
    k0_Z = 0.0;
    
###############################################################################
    
class ColdPlasma:
    """
    Cold plasma material data written by Write_GiD_bch_File_2D:
    - ne_File, Be_File : ERMES plasma files (node ne / node Bx By Bz)
    - Damping          : electron collision frequency [rad/s]
    - H, D, T          : relative densities of the ion species
    """
    ne_File = '../Electron_Density_Data.dat'; 
    Be_File = '../Bext_Field_Data.dat';
    
    Damping = 1e9;
    
    H = 0.0; 
    D = 0.0;
    T = 0.0;
    
###############################################################################
    
class ProblemSettings:
    """
    ERMES problem data written by Write_GiD_bch_File_2D:
    - ElemType  : element type
    - Frequency : frequency [Hz]
    - SolvName, SolvFile : external solver command and script
    """
    ElemType = 'RME_1st'; 
    
    Frequency = 28e9;
    
    SolvName = 'python3';
    SolvFile = '../PETSc/ERMES2PETSc.py'; 

###############################################################################

def Write_GiD_bch_File_2D(GiDFileName,
                          GeoFileName, 
                          MeshSize   , 
                          GBeam      , 
                          CPlasma    , 
                          PSettings  , 
                          P1, P2, P3, P4, P5, VV):
    """                         
    Writes a GiD batch (bch) file that creates the geometry, assigns the
    ERMES materials and conditions, meshes and writes the ERMES input files.
    
    Input:
    - GiDFileName : GiD project to create, <path>/<Problem_Name>.gid
    - GeoFileName : batch file to write (full path)
    - MeshSize    : mesh element size
    - GBeam       : Gaussian beam properties
    - CPlasma     : cold plasma properties
    - PSettings   : problem settings
    - P1,..., P5  : corner points of the GiD geometry
    - VV          : extrusion point
    """
    File = open( GeoFileName, 'w');

    ## Login
    #File.write('*****TCL GiD_Login login user@email.com password \n');
    #File.write('*****TCL GiD_Login start_session gid 15 \n\n');
    
    # Clean GiD enviorement
    File.write('Mescape Files New No escape escape escape escape \n');
    File.write('data defaults ProblemType UNKNOWN \n\n');
    
    # Load problem type
    File.write('Mescape Data Defaults ProblemType yes ERMES_20.0.5/ERMES \n');
    File.write('escape escape escape escape escape escape \n\n');
    
    # Save file
    File.write('Mescape Files SaveAs ');
    File.write('-alsoresults:0 -geoversion:current -- \n');
    File.write('{%s} \n\n' % (GiDFileName));  
        
    # Create geometry
    File.write('Mescape Geometry Create Line \n');
    File.write('%.6f,%.6f,%.6f \n' % (P1[0],P1[1],P1[2]));
    File.write('%.6f,%.6f,%.6f \n' % (P2[0],P2[1],P2[2]));
    File.write('%.6f,%.6f,%.6f \n' % (P3[0],P3[1],P3[2]));
    File.write('%.6f,%.6f,%.6f \n' % (P4[0],P4[1],P4[2]));
    File.write('%.6f,%.6f,%.6f \n' % (P5[0],P5[1],P5[2]));
    File.write('%.6f,%.6f,%.6f \n' % (P1[0],P1[1],P1[2]));
    File.write('escape escape \n\n');

    ## Add in Linux
    # File.write('Mescape Geometry Create Line Join 1 5 \n');
    # File.write('escape escape \n\n');

    File.write('Mescape Utilities Collapse model Yes \n');
    File.write('Mescape Utilities Repair Yes \n\n');
    
    File.write('Mescape Geometry Create NurbsSurface 1 2 3 4 5 ');
    File.write('escape escape escape \n');
    File.write('\'GetPointCoord Silent FJoin 1 escape \n'); 
    File.write('\'GetPointCoord Silent FJoin 1 escape \n\n'); 
    
    File.write('Mescape Utilities Copy Surfaces DoExtrude Volumes ');
    File.write('MaintainLayers \nTranslation FJoin 1 FNoJoin ');
    File.write('%.6f,%.6f,%.6f ' % (VV[0],VV[1],VV[2])); 
    File.write('1 escape \n\n');
    
    File.write('Mescape Geometry Create Contact ');
    File.write('SeparatedVolume 1 7 escape escape \n\n');
    
    File.write('Mescape Utilities Collapse model Yes \n');
    File.write('Mescape Utilities Repair Yes \n');
    File.write('\'Zoom Frame escape \n\n');
    
    File.write('Mescape Files Save \n\n');
        
    # Gaussian beam
    File.write('Mescape Data Materials NewMaterial \n');
    File.write('Robin_coefficients_Gauss_beam \n');
    File.write('Robin_coefficients_Gauss_beam Yes \n');
    File.write('Robin_coeffs_Gauss_beam \n');
    File.write('%.6f %.6f { } \n'      % (GBeam.modEpar,GBeam.phaEpar));
    File.write('%.6f %.6f { } \n'      % (GBeam.modEper,GBeam.phaEper));
    File.write('%.6f %.6f { } \n'      % (GBeam.Rpar,GBeam.Rper));
    File.write('%.6f %.6f %.6f { } \n' % (GBeam.OO_X,GBeam.OO_Y,GBeam.OO_Z));
    File.write('%.6f %.6f %.6f { } \n' % (GBeam.Ep_X,GBeam.Ep_Y,GBeam.Ep_Z));
    File.write('%.6f %.6f %.6f \n'     % (GBeam.k0_X,GBeam.k0_Y,GBeam.k0_Z));
    File.write('escape \n\n');
    
    # Cold plasma
    File.write('Mescape Data Materials NewMaterial \n');
    File.write('Cold_plasma Cold_plasma Yes cold_plasma \n');
    File.write('Full_3D 0.0 0.0 { } \n');
    File.write('0.0 0.0 0.0 { } \n');
    File.write('0.0 0.0 0.0 { } \n');
    File.write('0.0 0.0 0.0 { } \n');
    File.write('{%s} \n' % (CPlasma.ne_File));
    File.write('{%s} \n' % (CPlasma.Be_File));
    File.write('%.6f %.6f %.6f { } \n' % (CPlasma.H,CPlasma.D,CPlasma.T));
    File.write('0.0 0.0 { } \n');
    File.write('0.0 0.0 { } \n');
    File.write('0.0 0.0 0.0 0.0 { } \n');
    File.write('0.0 0.0 0.0 \n');  
    File.write('Full_matrix Off 1e-6 Off 1e6 { } \n');
    File.write('ComplexFreq_ele \n');
    File.write('%.6g \n' % (CPlasma.Damping));
    File.write('escape \n\n');
    
    # Problem settings
    File.write('Mescape Data ProblemData Cold_plasma 3D { } \n');
    File.write('{%s} \n' % (PSettings.ElemType));
    File.write('Off { } On Off Off { } Nodes Off Release m 1.0 1e-9 { } \n');
    File.write('Area_weighted On On { } \n');
    File.write('Manual Periodic 0.0 0.0 Off Off Off { } \n');
    File.write('Off Binary 1 %.6g { } \n' % (PSettings.Frequency));
    File.write('0 0.0 1.0 { } 0 1e9 3e9 0.1e9 { } \n');
    File.write('{%s} { } { } 8 1e-3 1 \n' % (PSettings.Solver));
    File.write('Diagonal 1e6 1e6 250 250 Nil_vector Final_step { } { } \n');
    File.write('{%s} {%s} \n' % (PSettings.SolvName,PSettings.SolvFile));
    File.write('1 0 0 0 1 0 0 {Time step} T T/32 \n');
    File.write('0 0 0 0 1 1 1 0 0 0 0 1 1 1 1 0 1 0 \n\n');
   
    File.write('Mescape Files Save \n\n');
    
    # Assign materials and boundary conditions
    File.write('Mescape Data Materials ');  
    File.write('AssignMaterial Cold_plasma Volumes 1 escape \n');
    File.write('Mescape escape escape escape escape \n\n');
    
    File.write('Mescape Data Conditions AssignCond ');
    File.write('Far_field_condition_Cold_plasma \n');
    File.write('Change FW 3 4 5 6 escape \n');
    File.write('Mescape escape escape escape escape \n\n');    
    
    File.write('Mescape Data Conditions AssignCond ');
    File.write('Robin_condition_Full_wave \n');
    File.write('Change Robin_coefficients_Gauss_beam 2 escape \n');
    File.write('Mescape escape escape escape escape \n\n');  
    
    File.write('Mescape Data Conditions AssignCond ');
    File.write('Field_surface_integral \n');
    File.write('Change Surface_1 2 escape \n');
    File.write('Mescape escape escape escape escape \n\n');   

    File.write('Mescape Data Conditions AssignCond ');
    File.write('Field_surface_integral \n');
    File.write('Change Surface_2 6 escape \n');
    File.write('Mescape escape escape escape escape \n\n');              
    
    File.write('Mescape Files Save \n\n');
    
    # Mesh
    if ( MeshSize > 0 ):
        File.write('Mescape Meshing Generate ');
        File.write('%.6f ' % (MeshSize));
        File.write('MeshingParametersFrom=Model escape \n\n');
    
        File.write('Mescape Files Save \n');
        File.write('Mescape Files WriteCalcFile \n');
        File.write('Mescape Files Save \n\n');
    
    File.close();

    ## Log-out
    # File.write('*****TCL GiD_Login logout \n\n');
    
    return 0;
    
###############################################################################

def Calculate_O_Mode(k0, Be, f):
    """
    O-mode polarization of the launched beam, evaluated with the B field at
    the critical density (cold plasma, X = 1).

    Input:
    - k0 : beam wave vector (normalized)
    - Be : external B field [T] (Cartesian)
    - f  : beam frequency [Hz]

    Output:
    - cEx   : complex E component along vEp (perpendicular to k0 and B);
              used as GaussBeam Epar
    - cEh   : complex E component in the plane of k0 and B (transverse to
              k0); used as GaussBeam Eper
    - vEp   : polarization vector vEp = (B x k0)/|B x k0|
    - theta : angle between k0 and B [rad]

    Notes:
    - The amplitude of cEx is fixed (0.5827); cEx and cEh are not normalized.
    - cEh diverges as theta -> 90 degrees (cos(theta) -> 0).
    """
    k0 = np.asarray(k0, dtype=float);
    Be = np.asarray(Be, dtype=float);
    
    # Vertical axis
    x0 = np.cross(Be, k0);
    x0 = x0 / np.linalg.norm(x0);
    
    vEp = x0;
    
    # Angle between k0 and Be
    theta = Angle_Between_Vectors(k0, Be, deg=False);
            
    # Y parameter
    w0 = 2 * np.pi * f  ;
    
    modBe = np.linalg.norm(Be);
    we    = QE * modBe / ME   ;
 
    Y = we / w0;
    
    # O-mode vertical component
    cEx = -1j * 0.5827;
    
    # O-mode horizontal component
    cEz =  1j * cEx * 0.5 * ( 
        Y * np.sin(theta)**2 + 
        np.sqrt( Y**2 * np.sin(theta)**4 + 4.0  * np.cos(theta)**2 ) 
    );
    
    cEy = -cEz * np.sin(theta) / np.cos(theta);
    
    cEh = np.sqrt( cEy*np.conj(cEy) + cEz*np.conj(cEz) );
    
    if ( abs(theta) > np.pi/2.0 ): cEh = -cEh;
        
    return cEx, cEh, vEp, theta;
    
###############################################################################

def Read_GiD_Nodes(Mesh_File):
    """
    - Reads the nodes of an ERMES/GiD input file <Problem_Name>-1.dat,
      lines "No[ID] = p(x,y,z);".
    - Returns: node ids (int64) and coordinates x, y, z (float), in file order.
    """
    pat = re.compile(
        r"No\[(\d+)\]\s*=\s*p\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^)]+)\s*\)\s*;"
    );
    
    Id, Xs, Ys, Zs = [], [], [], [];
    
    with open( Mesh_File, 'r' ) as f:
        for line in f:
            m = pat.search( line );
            if m:
                Id.append( int  ( m.group(1) ) );
                Xs.append( float( m.group(2) ) );
                Ys.append( float( m.group(3) ) );
                Zs.append( float( m.group(4) ) );
    
    if not Id:
        raise RuntimeError("No nodes 'No[ID] = p(x,y,z);' found in " + Mesh_File);
    
    return ( np.asarray(Id, dtype=np.int64), np.asarray(Xs, float),
             np.asarray(Ys, float), np.asarray(Zs, float) );

###############################################################################

def Transform_FRZT_2_FXYZ(nX, nY, nZ, Fr, Fz, Ft):
    """
    - Coordinate system transformation (R, phi, Z) -> (x, y, z)
    - The tokamak axis is the z axis: R = sqrt(x^2 + y^2), phi = atan2(y, x),
      with (R, phi, Z) right-handed (phi counter-clockwise seen from +z),
      as in the EQDSK (B_t = F/R along +phi) and Scotty.
    - Input : node coordinates nX, nY, nZ and the components Fr, Fz, Ft
              (note the order: R, Z, toroidal).
    - Output: Cartesian components Fx, Fy, Fz (0 on the axis R = 0).
    """
    nX = np.asarray(nX, float); nY = np.asarray(nY, float)
    Fr = np.asarray(Fr, float); Fz = np.asarray(Fz, float); Ft = np.asarray(Ft, float)
    modr = np.sqrt(nX**2 + nY**2)
    with np.errstate(invalid="ignore", divide="ignore"):
        cphi = np.where(modr > 0, nX / modr, 1.0)
        sphi = np.where(modr > 0, nY / modr, 0.0)
    cFx = np.where(modr > 0, Fr * cphi - Ft * sphi, 0.0)
    cFy = np.where(modr > 0, Fr * sphi + Ft * cphi, 0.0)
    cFz = Fz.copy()
    return cFx, cFy, cFz;

###############################################################################
