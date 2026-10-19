###############################################################################
# - This script defines functions used in problem_generator.py
###############################################################################

import re
import os
import numpy as np

from scipy.interpolate import UnivariateSpline, RectBivariateSpline

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
    """
    MirrorCoord = np.asarray(MirrorCoord, dtype=float);
    PlasmaCentr = np.asarray(PlasmaCentr, dtype=float);    
    
    # Rotation around Y-axis 
    vector = PlasmaCentr - MirrorCoord; 
    axis   = [0, 1, 0]; 
    angle  = tor_dgr;
    k0     = Rot3DAxis_dgr(vector, axis, angle);

    # Rotation around axis perpendicular to r0 
    vector = k0; 
    axis   = [-k0[2], 0, k0[0]]; 
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
    ne_File = '..\Electron_Density_Data.dat'; 
    Be_File = '..\Bext_Field_Data.dat';
    
    Damping = 1e9;
    
    H = 0.0; 
    D = 0.0;
    T = 0.0;
    
###############################################################################
    
class ProblemSettings:
    ElemType = 'RME_1st'; 
    
    Freqency = 28e9;
    
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
    Writes a bch file to create the geometry in GiD
    
    Input:
    - GeoFileName : GiD file name (full path)
    - GiDFileName : batch file name (full path)
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
    File.write('Off Binary 1 %.6g { } \n' % (PSettings.Freqency));
    File.write('0 0.0 1.0 { } 0 1e9 3e9 0.1e9 { } \n');
    File.write('External_solver { } { } 8 1e-3 1 \n');
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
    Input:
    - k0 : beam wave vector 
    - Be : external B field
    - f  : beam frequency

    Output:
    - cEx   : E perpendicular component 
    - cEh   : E horizontal component
    - vEp   : E polarization vector
    - theta : angle between k0 and Bext in rad
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
    qe = 1.602176634e-19; 
    me = 9.109383560e-31;
    w0 = 2 * np.pi * f  ;
    
    modBe = np.linalg.norm(Be);
    we    = qe * modBe / me   ;
 
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

def find_nearest_index(array, value):
    """
    - Index search
    """  
    array = np.asarray(array);
    idx = (np.abs(array - value)).argmin();
    return idx;

###############################################################################

def preprocess_line(line):
    """
    - Line processing
    """  
    ## Add spaces before the negative sign, 
    ## if not preceded by 'e' or whitespace
    return re.sub(r'(?<![\sEe])([+-])(?=\d)', r' \1', line);

###############################################################################

def Transform_FRZT_2_FXYZ(nX, nY, nZ, Fr, Fz, Ft):
    """
    - Coordinate system transformation
    """   
    # Array size
    asize = len( nX );
    
    # B field in cartesian coordinates
    cFx = np.zeros( asize ); 
    cFy = np.zeros( asize );
    cFz = np.zeros( asize );

    # Cylindrical unitary vectors
    rhat = np.zeros( 2 ); 
    that = np.zeros( 2 );
    
    for i in range( asize ):
        
        modr = np.sqrt( nX[i]**2 + nZ[i]**2 );
    
        if modr > 0: 
            rhat[0] = nX[i] / modr; 
            rhat[1] = nZ[i] / modr;
            
            that[0] = -rhat[1]; 
            that[1] =  rhat[0]; 
            
            cFx[i] = Fr[i] * rhat[0] + Ft[i] * that[0];
            cFy[i] = Fz[i];
            cFz[i] = Fr[i] * rhat[1] + Ft[i] * that[1];
            
        else:
            cFx[i] = 0.0;
            cFy[i] = Fz[i];
            cFz[i] = 0.0;
            
    return cFx, cFy, cFz;

###############################################################################

def Read_Equilibrium_Files(EQDSK_File, Norm_ne_File, Norm_rp_File):
    """
    - Reads EQDSK equilibrium files
    """   
    with open( EQDSK_File, 'r' ) as fid:
        ## Read the file header (first 51 characters) and discard it
        fid.read(51)
        
        ## Read the simulation grid info (3 float numbers)
        temp0  = [float(val) for val in fid.readline().split()]
        nw, nh = temp0[1], temp0[2]
        
        ## Read 0D parameters (20 float numbers)
        temp0 = fid.readline()
        temp0 = [float(val) for val in preprocess_line(temp0).split()]
        
        temp1 = fid.readline()
        temp1 = [float(val) for val in preprocess_line(temp1).split()]
        
        temp2 = fid.readline()
        temp2 = [float(val) for val in preprocess_line(temp2).split()]
           
        temp3 = fid.readline()
        temp3 = [float(val) for val in preprocess_line(temp3).split()]
        
        parameters_0D = {
            'rdim'     : temp0[0],
            'zdim'     : temp0[1],
            'rcentr'   : temp0[2],
            'rleft'    : temp0[3],
            'zmid'     : temp0[4],
            
            'rmaxis'   : temp1[0],
            'zmaxis'   : temp1[1],
            'simag'    : temp1[2],
            'sibry'    : temp1[3],
            'bcentr'   : temp1[4],
            
            'current'  : temp2[0],
            'simag2'   : temp2[1],
            'xdum1'    : temp2[2],
            'rmagaxis2': temp2[3],
            'xdum2'    : temp2[4],
            
            'zmaxis2'  : temp3[0],
            'xdum3'    : temp3[1],
            'sibry2'   : temp3[2],
            'xdum4'    : temp3[3],
            'xdum5'    : temp3[4],  
        }
         
        templine = fid.readline()
        templine = [float(val) for val in preprocess_line(templine).split()]
        
        ## Read 1D profile (nw float numbers)
        fpol = []
        
        indextrack = 0
        while len(fpol) < nw:
            fpol.append(templine[indextrack])
            indextrack += 1
            if indextrack == len(templine):
                templine = fid.readline()
                templine = (
                [float(val) for val in preprocess_line(templine).split()] 
                )
                indextrack = 0   
                
        pres = []
        
        while len(pres) < nw:
            pres.append(templine[indextrack])
            indextrack += 1
            if indextrack == len(templine):
                templine = fid.readline()
                templine = (
                [float(val) for val in preprocess_line(templine).split()] 
                )
                indextrack = 0
                
        ffprim = []
        
        while len(ffprim) < nw:
            ffprim.append(templine[indextrack])
            indextrack += 1
            if indextrack == len(templine):
                templine = fid.readline()
                templine = ( 
                [float(val) for val in preprocess_line(templine).split()] 
                )
                indextrack = 0
                
        pprime = []
        
        while len(pprime) < nw:
            pprime.append(templine[indextrack])
            indextrack += 1
            if indextrack == len(templine):
                templine = fid.readline()
                templine = ( 
                [float(val) for val in preprocess_line(templine).split()] 
                )
                indextrack = 0
                
        psirz = []
        
        while len(psirz) < nh:
            sub = []
            while len(sub)< nw:
                sub.append(templine[indextrack])
                indextrack += 1
                if indextrack == len(templine):
                    templine = fid.readline()
                    templine = ( 
                    [float(val) for val in preprocess_line(templine).split()] 
                    )
                    indextrack = 0
                if len(sub) == nw:
                    psirz.append(sub)
                    
        qpsi = []
        
        while len(qpsi) < nw:
            qpsi.append(templine[indextrack])
            indextrack += 1
            if indextrack == len(templine):
                templine = fid.readline()
                templine = ( 
                [float(val) for val in preprocess_line(templine).split()] 
                )
                indextrack = 0
    
        nbbs, limitr = int(templine[0]),int(templine[1])
    
        rbbs = []
        zbbs = []
        
        oddevenline = 1 
        templine = fid.readline()
        templine = [float(val) for val in preprocess_line(templine).split()]
        while len(rbbs)< nbbs or len(zbbs) < nbbs:
            if oddevenline%2 == 1 and len(templine) == 5:
                rbbs.append(templine[0])
                zbbs.append(templine[1])
                rbbs.append(templine[2])
                zbbs.append(templine[3])
                rbbs.append(templine[4])
                templine = fid.readline()
                templine = ( 
                [float(val) for val in preprocess_line(templine).split()] 
                )
                oddevenline += 1
            elif oddevenline%2 == 0 and len(templine) == 5:
                zbbs.append(templine[0])
                rbbs.append(templine[1])
                zbbs.append(templine[2])
                rbbs.append(templine[3])
                zbbs.append(templine[4])
                templine = fid.readline()
                templine = ( 
                [float(val) for val in preprocess_line(templine).split()] 
                )
                oddevenline += 1
            elif oddevenline%2 == 1: #odd number, clear up the lines
                for i in range(0,len(templine)):
                    if i%2 == 0:
                        rbbs.append(templine[i])
                    else:
                        zbbs.append(templine[i])
            elif oddevenline%2 == 0: #odd number, clear up the lines
                for i in range(0,len(templine)):
                    if i%2 == 0:
                        zbbs.append(templine[i])
                    else:
                        rbbs.append(templine[i])   
    
        rlim = []
        zlim = []
        
        while len(rlim)< limitr or len(zlim) < limitr:
            templine = fid.readline()
            templine = ( 
            [float(val) for val in preprocess_line(templine).split()] 
            )
            if oddevenline%2 == 1:
                for i in range(0,len(templine)):
                    if i%2 == 0:
                        if len(rlim) != limitr:
                            rlim.append(templine[i])
                    else:
                        if len(zlim) != limitr:
                            zlim.append(templine[i])
                oddevenline += 1
            else:
                for i in range(0,len(templine)):
                    if i%2 == 0:
                        if len(zlim) != limitr:
                            zlim.append(templine[i])
                    else:
                        if len(rlim) != limitr:
                            rlim.append(templine[i])
                oddevenline += 1
            
        Rgrid = np.linspace(parameters_0D['rleft'], parameters_0D['rleft'] + 
                            parameters_0D['rdim' ], 
                            int(nw));
        
        Zgrid = np.linspace(parameters_0D['zmid'] - parameters_0D['zdim'] / 2,
                            parameters_0D['zmid'] + parameters_0D['zdim'] / 2, 
                            int(nh));
    
    ## Flux on the RZ plane
    psirz2d = np.transpose( np.array( psirz ) );
    
    ## Normalize such that 0 at magnetic axis, 1 at LCFS.
    psirz2dnorm = np.divide( psirz2d-parameters_0D['sibry2'], 
                             parameters_0D['simag2'] );
    
    ## Normalize such that minimum is 0. Divide by abs of min value
    max_psi     =  np.max( psirz2dnorm ); ### rescaling.
    psirz2dnorm =  psirz2dnorm - max_psi; ### shifting..
    psirz2dnorm = -psirz2dnorm / max_psi;
    
    ## Getting B_R, B_T, and B_Z
    B_R = np.zeros_like( psirz2dnorm );
    B_T = np.zeros_like( psirz2dnorm );
    B_Z = np.zeros_like( psirz2dnorm );
    
    interp_order     = 3;
    interp_smoothing = 0;
    
    Rin_new = [];
    
    for i in range( len( Rgrid ) ):
        Rin_join = [Rgrid[i]]*len(Rgrid);
        Rin_new = Rin_new + Rin_join;
    
    interp_psirz2d = RectBivariateSpline(
        Rgrid,
        Zgrid,
        psirz2d,
        kx=interp_order,
        ky=interp_order,
        s=0
        ) 
    
    ## From 0 to 1(LCFS) to max.
    npsi = len( fpol );
    pol_flux_array_EQDSK_norm = np.linspace( 0, 1, npsi ); 
    
    interp_rBphi = UnivariateSpline(
        pol_flux_array_EQDSK_norm, 
        fpol,                      
        w    = None,
        bbox = [None, None],
        k    = interp_order,
        s    = interp_smoothing,
        ext  = 0,
        check_finite = False,
    )
    
    ## First derivative
    dpolflux_dR = interp_psirz2d( Rgrid, Zgrid, dx=1, dy=0 );
    dpolflux_dZ = interp_psirz2d( Rgrid, Zgrid, dx=0, dy=1 );
    
    for i in range(int(nw)):
        B_R[i,:] = -1*dpolflux_dZ[i,:] / Rgrid[i];
        B_Z[i,:] =    dpolflux_dR[i,:] / Rgrid[i];
    
    for i in range( int( nw ) ):     #R
        for j in range( int( nh ) ): #Z
            B_T[i,j] = interp_rBphi(psirz2dnorm[i,j]) / Rgrid[i];
    
    ## Parse GiD mesh nodes and write ERMES ne and Be files
    ne_data   =   np.loadtxt( Norm_ne_File, dtype = float );
    ro_data   = ( np.loadtxt( Norm_rp_File, dtype = float ) )**2;
    
    ne_interp = UnivariateSpline   ( ro_data, ne_data, k=5, s=0 );
    Ps_spline = RectBivariateSpline( Rgrid, Zgrid, psirz2dnorm, 
                                     kx=5 , ky=5 , s=0 );
    
    Br_spline = RectBivariateSpline( Rgrid, Zgrid, B_R, kx=5, ky=5, s=0 );
    Bt_spline = RectBivariateSpline( Rgrid, Zgrid, B_T, kx=5, ky=5, s=0 );
    Bz_spline = RectBivariateSpline( Rgrid, Zgrid, B_Z, kx=5, ky=5, s=0 );
    
    return ne_interp, Ps_spline, Br_spline, Bt_spline, Bz_spline;
    
###############################################################################

def Generate_Plasma_Files(EQDSK_File    , 
                          Norm_ne_File  , Norm_rp_File,
                          Problem_Folder, Problem_Name, 
                          Plasma_Folder , Plasma_Name ):
    """
    - Generates ERMES plasma files from equilibrium data files
    """
    ## Compose file names
    ERMES_Folder = Problem_Folder + Problem_Name + '.gid/';
    
    ne_3D_File = Plasma_Folder + 'ne_' + Plasma_Name + '.dat';
    Be_3D_File = Plasma_Folder + 'Be_' + Plasma_Name + '.dat';
        
    ## Read equilibrium file
    ne_interp, Ps_spline, Br_spline, Bt_spline, Bz_spline = (
        Read_Equilibrium_Files(EQDSK_File, Norm_ne_File, Norm_rp_File) 
    );
    
    ## Parse GiD mesh nodes:  No[ID] = p(x,y,z);
    Id, Xs, Ys, Zs = [], [], [], [];
    
    pat = re.compile(
      r"No\[(\d+)\]\s*=\s*p\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^)]+)\s*\)\s*;"
    );
    
    # Nodes list
    ERMES_1_dat_File = ERMES_Folder + Problem_Name + '-1.dat';
    
    # Check if the file exists
    if os.path.exists(ERMES_1_dat_File):
       print( 'Reading nodes list file:\n' + ERMES_1_dat_File );
    else:
       print( '\nNodes list file:    \n' + 
               ERMES_1_dat_File + '  \n' + 
              'must be created !!!!!.\n' ); 
       return;       
    
    # Open and read file
    with open( ERMES_1_dat_File, 'r' ) as f:
        for line in f:
            m = pat.search( line );
            if m:
                Id.append( int  ( m.group(1) ) );
                Xs.append( float( m.group(2) ) );
                Ys.append( float( m.group(3) ) );
                Zs.append( float( m.group(4) ) );
                
    print( 'File read.' );
                
    nId = np.asarray( Id, dtype=np.int64 );
    nX  = np.asarray( Xs, float );
    nY  = np.asarray( Ys, float );
    nZ  = np.asarray( Zs, float );
    
    ## Map (nX, nY, nZ) -> (R, phi, Z) 
    Rnod = np.sqrt( nX**2 + nZ**2 );
    Znod = nY;
    
    ## Evaluate Psi, Br, Bt, Bz at node (R,Z)
    Ps = Ps_spline.ev( Rnod, Znod );
    
    Br = Br_spline.ev( Rnod, Znod );
    Bt = Bt_spline.ev( Rnod, Znod );
    Bz = Bz_spline.ev( Rnod, Znod );
    
    ne = ne_interp( Ps );
    
    ## Transform to cartesian coordinates
    cBx, cBy, cBz = Transform_FRZT_2_FXYZ( nX, nY, nZ, Br, Bz, Bt );
    
    ## Write ouputs ( sorted by NodeID )
    order      = np.argsort( nId );
    nid_sorted = nId[ order ];
    
    np.savetxt( 
        ne_3D_File, 
        np.column_stack( [ nid_sorted, ne[order] ] ),
        fmt=["%d", "%.8e"]
    );
    
    np.savetxt(
        Be_3D_File,
        np.column_stack( [ nid_sorted, 
                           cBx[order], 
                           cBy[order], 
                           cBz[order] ] ),
        fmt=["%d", "%.8e", "%.8e", "%.8e"]
    );
    
###############################################################################

def Get_Point_Equilibrium_Data(P0, EQDSK_File, Norm_ne_File, Norm_rp_File):
    """
    - Gives ne and Bext at P0
    """                   
    ## Read equilibrium file
    ne_interp, Ps_spline, Br_spline, Bt_spline, Bz_spline = (
    Read_Equilibrium_Files(EQDSK_File, Norm_ne_File, Norm_rp_File) );
    
    ## Get P0 point
    nX  = np.asarray([ P0[0] ], float);
    nY  = np.asarray([ P0[1] ], float);
    nZ  = np.asarray([ P0[2] ], float);
    
    ## Map (nX, nY, nZ) -> (R, phi, Z) 
    Rnod = np.sqrt( nX**2 + nZ**2 );
    Znod = nY;
    
    ## Evaluate Psi, Br, Bt, Bz at node (R,Z)
    Ps = Ps_spline.ev( Rnod, Znod );
    
    Br = Br_spline.ev( Rnod, Znod );
    Bt = Bt_spline.ev( Rnod, Znod );
    Bz = Bz_spline.ev( Rnod, Znod );
    
    ne = ne_interp( Ps );
    
    cBx, cBy, cBz = Transform_FRZT_2_FXYZ( nX, nY, nZ, Br, Bz, Bt );
    
    Bext = np.array( [ cBx[0], cBy[0], cBz[0] ] );
    
    return ne, Bext;
    
###############################################################################

def Get_Pcrtc_Equilibrium_Data(P0, k0, Freq, Rpl, Cpl,
                               EQDSK_File, Norm_ne_File, Norm_rp_File,
                               step=1e-3, tol=1e-3, max_iter=100000):
    """
    Returns:
        nc  : critical density
        Bcr : magnetic field at location where ne = nc
    """
    # Critical density
    wo = 2.0*np.pi* Freq;
    qe = 1.602176634e-19; 
    me = 9.109383560e-31;
    eo = 8.854187817e-12;
 
    nc = eo * me * ( wo / qe )**2;

    # Starting point on plasma boundary 
    Pn = LineSphere_Intersect(P0, k0, Rpl, Cpl);

    ne_prev, B_prev = Get_Point_Equilibrium_Data(
        Pn, EQDSK_File, Norm_ne_File, Norm_rp_File
    );

    # If already at critical density (within tolerance)
    if np.isclose(ne_prev, nc, rtol=tol):
        return nc, B_prev;

    # Determine marching direction
    direction = k0 if ne_prev < nc else -k0;

    # March until bracketing nc
    for ni in range(max_iter):

        Pn_next = Pn + direction * step;

        ne_next, B_next = Get_Point_Equilibrium_Data(
            Pn_next, EQDSK_File, Norm_ne_File, Norm_rp_File
        );

        # Check if nc is bracketed between previous and next step
        if ( (ne_prev - nc) * (ne_next - nc) <= 0 ):

            # Linear interpolation for better accuracy
            t = (nc - ne_prev) / (ne_next - ne_prev);

            Pcr = Pn + t * (Pn_next - Pn);
            
            ni, Bcr = Get_Point_Equilibrium_Data(
                Pcr, EQDSK_File, Norm_ne_File, Norm_rp_File
            );

            return ni[0], Bcr;

        # Advance
        Pn      = Pn_next;
        ne_prev = ne_next;

    raise RuntimeError("Critical density not found within max_iter");
    
###############################################################################