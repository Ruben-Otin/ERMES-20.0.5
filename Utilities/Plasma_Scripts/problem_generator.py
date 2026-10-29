###############################################################################
# problem_generator.py
#
# - Generates a 2D (one-element-thick) ERMES model of a Gaussian beam
#   launched into a tokamak plasma (MAST-U DBS benchmark):
#     1. beam direction (mirror + toroidal/poloidal angles) and O- or X-mode
#        polarisation (B at the critical density),
#     2. GiD batch file (geometry, materials, conditions, mesh) -> run GiD,
#     3. ERMES plasma files ne_<id>.dat, Be_<id>.dat at the mesh nodes,
#     4. run ERMES (optional).
#
# - 'gu.*'  : utilities.py           (geometry, GiD batch file, polarisation)
# - 'e2e.*' : eqdsk_to_ermes_turb.py  (EQDSK equilibrium, plasma files)
#
# - Frame: z = tokamak axis, lengths in metres.
#
# - Environment variables (optional, used by batch_script.sh):
#     RUN_GID=0 / RUN_ERMES=0 : skip the GiD / ERMES step
#     GID_EXE, ERMES_EXE      : executables (override the defaults below)
###############################################################################
# From bash:
# FREQ=28e9
# TOR=20
# POL=20
# EPOL="X"               # "O" or "X"
# MESH=0.003
# ETYP="RME_1st"
# SOLV="Direct_solver"   # optional: "Direct_solver" or "External_solver"
#
# python3 problem_generator_v2.py $FREQ $TOR $POL $EPOL $MESH $ETYP $SOLV
#
# Without arguments (e.g. from Spyder) the defaults below are used.
###############################################################################

import numpy               as np
import utilities           as gu
import eqdsk_to_ermes_turb as e2e   # EQDSK -> ERMES plasma files (+ turbulence)
import subprocess
import platform
import sys
import os

###############################################################################

# Mirror XYZ coordinates (z = tokamak axis). Same physical position as the
# previous y-axis version (-1818.389, 0, 284.600) mm, rotated (x, y, z) -> (x, -z, y)
MirrorCoord = np.array([-1818.389, -284.600, 0.000]) * 1e-3;

# Beam waist (distance behind mirror and size)
WaistDist = (580 - 200) * 1e-3;
WaistSize = (47.10    ) * 1e-3;

# Plasma radius and center
PRadius = 1.42;
PCenter = np.array([0,0,0]);

###############################################################################
# Case parameters (defaults, overridden by the command-line arguments)
###############################################################################

Freq     = 28e9;            # Frequency
tor_dgr  = 20;              # Toroidal angle
pol_dgr  = 20;              # Poloidal angle
Epol     = 'O';             # Polarization: 'O', 'X'
MeshSize = 0.003;           # GiD mesh size
ElemType = 'RME_1st';       # Element type
Solver   = 'Direct_solver'; # 'Direct_solver' or 'External_solver'

if len(sys.argv) > 1:
    if len(sys.argv) < 7:
        sys.exit("Usage: python3 problem_generator_v2.py "
                 "FREQ TOR POL EPOL MESH ETYP [SOLV]");
    Freq     = float(sys.argv[1]);
    tor_dgr  = float(sys.argv[2]);
    pol_dgr  = float(sys.argv[3]);
    Epol     = sys.argv[4];
    MeshSize = float(sys.argv[5]);
    ElemType = sys.argv[6];
    if len(sys.argv) > 7: Solver = sys.argv[7];

Epol = Epol.strip().upper();

if Epol not in ('O', 'X'):
    sys.exit("Polarization must be 'O' or 'X' (got '%s')" % Epol);

WavL = 3e8/Freq;

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
NameId = 't%g_p%g_e%s' % (tor_dgr, pol_dgr, Epol);

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

# Plasma options (eqdsk_to_ermes_turb.py), used for BOTH the launch polarisation
# (B at the critical density) and the ERMES plasma files, so they stay consistent:
#  - toroidal_convention: "right-handed" ((R,phi,Z) right-handed, as the EQDSK
#    and Scotty) or "legacy" (old utilities.py frame, left-handed: mirrored
#    field-line helicity; only to reproduce old runs)
#  - turbulence: True adds the synthetic turbulent density (parameters turb_*,
#    see torbeam_to_ermes.py), e.g. turb_Lr = 0.015
#  - blobs: True adds field-aligned blobs/holes (parameters blob_*, hole_*)
Plasma_Options = dict(
    toroidal_convention = "right-handed",
    turbulence          = False,
);

# Create folders
os.makedirs(BCH_Folder   , exist_ok=True);
os.makedirs(Plasma_Folder, exist_ok=True);

###############################################################################
# Executables and run steps
###############################################################################

Windows = ( platform.system() == "Windows" );

if Windows:
    GiD_Exec   = r'D:/MySoftware/GiD 15.0.3/gid_offscreen.bat';
    ERMES_Exec = r'D:/MySoftware/GiD 15.0.3/problemtypes/ERMES_20.0.5/ERMES.gid/execs/win/ERMES_20.0.5.exe';
else:
    GiD_Exec   = os.path.expanduser(r'~/GiDx64/gid-15.0.4/gid_offscreen');
    ERMES_Exec = os.path.expanduser(r'~/GiDx64/gid-15.0.4/problemtypes/ERMES_20.0.5/ERMES.gid/execs/unix/ERMES_20.0.5');

GiD_Exec   = os.environ.get("GID_EXE"  , GiD_Exec  );
ERMES_Exec = os.environ.get("ERMES_EXE", ERMES_Exec);

RUN_GID   = os.environ.get("RUN_GID"  , "1") != "0";
RUN_ERMES = os.environ.get("RUN_ERMES", "1") != "0";

###############################################################################

# Gaussian beam wave vector
k0 = gu.Calculate_WaveVector(tor_dgr, pol_dgr, MirrorCoord, PCenter);

# Waist position
WaistCoord = MirrorCoord - k0 * WaistDist;

# E polarization
ncrt, Bcrt = e2e.Get_Pcrtc_Equilibrium_Data(
    MirrorCoord, k0, Freq, PRadius, PCenter,
    EQDSK_File, Norm_ne_File, Norm_rp_File,
    step=1e-4, tol=1e-4, max_iter=100000, **Plasma_Options
);

if ( Epol == 'O' ):
    cEpar, cEper, vEp, theta = gu.Calculate_O_Mode(k0, Bcrt, Freq);
else:
    cEpar, cEper, vEp, theta = gu.Calculate_X_Mode(k0, Bcrt, Freq);

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

PSettings.ElemType = ElemType; 

PSettings.Frequency = Freq;

# Solver type: 'Direct_solver' or 'External_solver'
PSettings.Solver = Solver;

# If PSettings.Solver ='External_solver'
PSettings.SolvName = 'bash';
PSettings.SolvFile = '../PETSc/ERMES2PETSc.sh'; 

###############################################################################

print(); 
print("--------------------------");
print("Problem settings          ");
print("--------------------------");  
print("Name    = %s"   % (Problem_Name)); 
print("Freq    = %.6g" % (Freq)); 
print("Mode    = %s"   % (Epol)); 
print("Element = %s"   % (PSettings.ElemType)); 
print("Solver  = %s"   % (PSettings.Solver)); 
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

if RUN_GID:
    print("Executing GiD...");
    
    if Windows:
        command = [GiD_Exec, "-b", BCH_Path, "-n"];
    else:
        command = [GiD_Exec, "-offscreen",
                   "-error" , BCH_Folder + Problem_Name + ".err",
                   "-output", BCH_Folder + Problem_Name + ".log",
                   "-b"     , BCH_Path, "-n"];
    print(command);
    
    results = subprocess.run(command, capture_output=True, text=True);
    
    # Only a warning: a missing mesh is detected below (Generate_Plasma_Files)
    if results.returncode != 0:
        print(results.stdout); print(results.stderr);
        print("WARNING: GiD exit status %d" % results.returncode);
    
    print("GiD finished.\n");

###############################################################################

print("Creating plasma files..."); 

# Raises FileNotFoundError if GiD did not write <Problem_Name>.gid/<Problem_Name>-1.dat
e2e.Generate_Plasma_Files( 
    EQDSK_File, Norm_ne_File, Norm_rp_File, Absolute_Path, 
    Problem_Name, Plasma_Folder, Plasma_NameId, **Plasma_Options 
);

print("Plasma files created.\n"); 

###############################################################################

if RUN_ERMES:
    print("Executing ERMES..."); 
    
    PROBL_Path = Absolute_Path + Problem_Name + r".gid/";
    PROBL_File = PROBL_Path + Problem_Name;
    
    # ERMES output and error messages go to the .info file; check=True raises
    # an exception if ERMES stops with an error (non-zero exit status)
    with open(PROBL_File + r".info", "w") as f:
        subprocess.run(
            [ERMES_Exec, PROBL_File],
            stdout = f,
            stderr = subprocess.STDOUT,
            cwd    = PROBL_Path,
            check  = True
        );
    
    print("ERMES finished."); 

###############################################################################
