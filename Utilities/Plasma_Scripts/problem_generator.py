###############################################################################
# - This script generates a GiD geometry, apply boundary conditions,
#   mesh and calculates the EM fields with ERMES.
#
# - The functions 'gu.*' used in this script are defined in the file
#   'utilities.py'.
###############################################################################

import numpy     as np
import utilities as gu
import subprocess
import os

###############################################################################

# Mirror XYZ coordinates
MirrorCoord = np.array([-1818.389, 0.000, 284.600]) * 1e-3;

# Beam waist (distance behind mirror and size)
WaistDist = (0.5 * (597.95 + 578.22) - 200) * 1e-3;
WaistSize = (0.5 * ( 47.18 +  47.06)      ) * 1e-3;

# Frequency
Freq = 28e9;
WavL = 3e8/Freq;

# Plasma radius and center
PRadius = 1.42;
PCenter = np.array([0,0,0]);

# GiD mesh size
MeshSize = 0.003;

# Toroidal and poloidal angles
tor_dgr = 20;
pol_dgr = 20;

###############################################################################
# Boundary size:
#
# - WidthR : right length from the beam center
# - WidthL : left  length from the beam center
#
# - PdeepI : length inside the plasma (incident beam)
# - PdeepR : length inside the plasma (reflected beam)
#
# - Length : length of the receiving boundary
# - Distpl : distance of the plasma to beam boundary
#   (if Distpl == 0.0 is calculated automatically)
#
# - Vdispl : extrusion distance
###############################################################################

WidthR = 0.250-0.000; 
WidthL = 0.250-0.050; 

PdeepI = 0.110-0.010; 
PdeepR = 0.110-0.010;

Length = 0.670-0.020; 
Distpl = 0.520-0.000;

Vdispl = MeshSize; 

###############################################################################

# Absolute path where GiD problem files are saved
Absolute_Path = os.getcwd() + r'/';

# Case name Id
NameId = 't' + str(tor_dgr) + '_p' + str(pol_dgr);

# GiD problem name and folder 
Problem_Name = 'MASTU-2D-' + NameId;

# BCH files
BCH_Folder = Absolute_Path + r'BatchFiles/';
BCH_Path   = BCH_Folder + Problem_Name + ".bch";

# Plasma files 
Plasma_NameId = '2D-' + NameId; 
Plasma_Folder = Absolute_Path + r'Plasma/';

# Equilibrium files
EQDSK_File   = r'./EQDSK_Samples/mast-u-sample.eqdsk';
Norm_ne_File = r'./EQDSK_Samples/rho_ne.txt';
Norm_rp_File = r'./EQDSK_Samples/rho_flux.txt';

# Create folders
if not os.path.exists(BCH_Folder   ): os.makedirs(BCH_Folder   );
if not os.path.exists(Plasma_Folder): os.makedirs(Plasma_Folder);

###############################################################################

# Gaussian beam wave vector
k0 = gu.Calculate_WaveVector(tor_dgr, pol_dgr, MirrorCoord, PCenter);

# Waist position
WaistCoord = MirrorCoord - k0 * WaistDist;

# E polarization
ncrt, Bcrt = gu.Get_Pcrtc_Equilibrium_Data(
    MirrorCoord, k0, Freq, PRadius, PCenter,
    EQDSK_File, Norm_ne_File, Norm_rp_File,
    step=1e-4, tol=1e-4, max_iter=100000
);

cEpar, cEper, vEp, theta = gu.Calculate_O_Mode(k0, Bcrt, Freq);

# Gauss beam properties
GBeam = gu.GaussBeam();

GBeam.modEpar = +abs(cEpar);
GBeam.phaEpar = np.angle(cEpar,deg=False);

# -:Left-polarized, +:Right-polarized
GBeam.modEper = -abs(cEper);
GBeam.phaEper = np.angle(cEper,deg=False);

GBeam.Rpar = WaistSize;
GBeam.Rper = WaistSize;

GBeam.OO_X = WaistCoord[0];
GBeam.OO_Y = WaistCoord[1];
GBeam.OO_Z = WaistCoord[2];

GBeam.Ep_X = vEp[0];
GBeam.Ep_Y = vEp[1];
GBeam.Ep_Z = vEp[2];

GBeam.k0_X = k0[0];
GBeam.k0_Y = k0[1];
GBeam.k0_Z = k0[2];

###############################################################################

# Cold plasma properties
CPlasma = gu.ColdPlasma();

CPlasma.ne_File = Plasma_Folder + 'ne_' + Plasma_NameId + '.dat';
CPlasma.Be_File = Plasma_Folder + 'Be_' + Plasma_NameId + '.dat';

CPlasma.Damping = 1e9;

CPlasma.H = 0.0;
CPlasma.D = 1.0;
CPlasma.T = 0.0;

###############################################################################

# Problem settings
PSettings = gu.ProblemSettings();

PSettings.ElemType = 'RME_1st'; 

PSettings.Freqency = Freq;

PSettings.SolvName = 'wsl';
PSettings.SolvFile = '../PETSc/win2wsl.sh >> SolverInfo.info'; 

###############################################################################

print(); 
print("--------------------------");
print("Problem settings          ");
print("--------------------------");  
print("Name    = %s"   % (Problem_Name)); 
print("Freq    = %.6g" % (Freq)); 
print("Element = %s"   % (PSettings.ElemType)); 
print("Damping = %.6g" % (CPlasma.Damping));
print();
print("B_crtic = %f,%f,%f"     % (Bcrt[0], Bcrt[1], Bcrt[2]));
print("n_crtic = %.6g"         % (ncrt));
print("ang_Bk0 = %.6g degrees" % (np.degrees(theta)));
print();
print("k0      = %f,%f,%f"            % (k0[0], k0[1], k0[2]));
print("Epol    = %f,%f,%f"            % (vEp[0], vEp[1], vEp[2]));
print("|Epar|  = %.6g,  Phase = %.6g" % (GBeam.modEpar, GBeam.phaEpar));
print("|Eper|  = %.6g,  Phase = %.6g" % (GBeam.modEper, GBeam.phaEper));
print(); 

###############################################################################

print("Creating GiD bch file..."); 

# Calculate GiD geometry points
P1,P2,P3,P4,P5,VV = gu.Calculate_GiD_Points( 
    MirrorCoord, k0, PRadius, PCenter, Distpl, 
    WidthR, WidthL, Length, PdeepI, PdeepR, Vdispl                                                
);

# Write GiD batch geometry generator file
gu.Write_GiD_bch_File_2D( 
    Absolute_Path + Problem_Name + '.gid',
    BCH_Path, MeshSize, GBeam, CPlasma, 
    PSettings, P1, P2, P3, P4, P5, VV 
);

print("GiD bch file created.\n"); 

###############################################################################

print("Executing GiD...");

# Windows
GiD_Path = r'D:/MySoftware/GiD 15.0.3/gid_offscreen.bat';
command  = [GiD_Path, "-b", BCH_Path, "-n"]; print(command);
results  = subprocess.run(command, capture_output=True, text=True);

## Linux
# GiD_Path = r"$HOME/GiDx64/gid-15.0.4/gid_offscreen";
# ERR_File = Absolute_Path + r"BatchFiles/" + Problem_Name + ".err";
# LOG_File = Absolute_Path + r"BatchFiles/" + Problem_Name + ".log";
# command  = [GiD_Path , "-offscreen", 
#             "-error" , ERR_File    , 
#             "-output", LOG_File    ,
#             "-b"     , BCH_Path    ,
#             "-n"
# ]; 
# print(command);
# results = subprocess.run(command, check=True);

print("GiD finished.\n");

###############################################################################

print("Creating plasma files..."); 

gu.Generate_Plasma_Files( 
    EQDSK_File, Norm_ne_File, Norm_rp_File, Absolute_Path, 
    Problem_Name, Plasma_Folder, Plasma_NameId 
);

print("Plasma files created.\n"); 

###############################################################################

print("Executing ERMES..."); 

# Problemtype folder
PType_Folder = r'D:/MySoftware/GiD 15.0.3/problemtypes/ERMES_20.0.5'

# Windows
ERMES_Exec = PType_Folder + r'/ERMES.gid/execs/win/ERMES_20.0.5.exe';

## Linux
# ERMES_Exec = PType_Folder +  r"/ERMES.gid/execs/unix/ERMES_20.0.5";

PROBL_Path = Absolute_Path + Problem_Name + r".gid/";
PROBL_File = PROBL_Path + Problem_Name;

with open(PROBL_File + r".info", "w") as f:
    subprocess.run(
        [ERMES_Exec, PROBL_File],
        stdout = f,
        stderr = subprocess.STDOUT,   # ERMES error messages also go to the .info file
        cwd    = PROBL_Path,
        check  = True
    );   # check=True: raises an exception if ERMES stops with an error ( non-zero exit status )
      
print("ERMES finished."); 

###############################################################################