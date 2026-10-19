###############################################################################
# - This script reads and visualize an EQDSK equilibrium file. 
###############################################################################

import numpy as np
import re
import matplotlib.pyplot as plt

from scipy.constants   import mu_0
from scipy.interpolate import UnivariateSpline, RectBivariateSpline

###############################################################################
## Files names
###############################################################################

EQDSK_eqfile = './EQDSK_Samples/mast-u-sample.eqdsk'
Norm_ne_file = "./EQDSK_Samples/rho_ne.txt"
Norm_te_file = "./EQDSK_Samples/rho_te.txt"
Norm_rp_file = "./EQDSK_Samples/rho_flux.txt"

###############################################################################
## Read eqdsk file
###############################################################################

def find_nearest_index(array, value):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return idx

def preprocess_line(line):
    ## Add spaces before the negative sign if not preceded by 'e' or whitespace
    return re.sub(r'(?<![\sEe])([+-])(?=\d)', r' \1', line)

## Open the file
with open( EQDSK_eqfile, 'r' ) as fid:
    
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
            templine = [float(val) for val in preprocess_line(templine).split()]
            indextrack = 0   
            
    pres = []
    while len(pres) < nw:
        pres.append(templine[indextrack])
        indextrack += 1
        if indextrack == len(templine):
            templine = fid.readline()
            templine = [float(val) for val in preprocess_line(templine).split()]
            indextrack = 0
            
    ffprim = []
    while len(ffprim) < nw:
        ffprim.append(templine[indextrack])
        indextrack += 1
        if indextrack == len(templine):
            templine = fid.readline()
            templine = [float(val) for val in preprocess_line(templine).split()]
            indextrack = 0
            
    pprime = []
    while len(pprime) < nw:
        pprime.append(templine[indextrack])
        indextrack += 1
        if indextrack == len(templine):
            templine = fid.readline()
            templine = [float(val) for val in preprocess_line(templine).split()]
            indextrack = 0
            
    psirz = []
    while len(psirz) < nh:
        sub = []
        while len(sub)< nw:
            sub.append(templine[indextrack])
            indextrack += 1
            if indextrack == len(templine):
                templine = fid.readline()
                templine = [float(val) for val in preprocess_line(templine).split()]
                indextrack = 0
            if len(sub) == nw:
                psirz.append(sub)
                
    qpsi = []
    while len(qpsi) < nw:
        qpsi.append(templine[indextrack])
        indextrack += 1
        if indextrack == len(templine):
            templine = fid.readline()
            templine = [float(val) for val in preprocess_line(templine).split()]
            indextrack = 0

    nbbs, limitr = int(templine[0]),int(templine[1])

    oddevenline = 1 
    rbbs = []
    zbbs = []
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
            templine = [float(val) for val in preprocess_line(templine).split()]
            oddevenline += 1
        elif oddevenline%2 == 0 and len(templine) == 5:
            zbbs.append(templine[0])
            rbbs.append(templine[1])
            zbbs.append(templine[2])
            rbbs.append(templine[3])
            zbbs.append(templine[4])
            templine = fid.readline()
            templine = [float(val) for val in preprocess_line(templine).split()]
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
        templine = [float(val) for val in preprocess_line(templine).split()]
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

###############################################################################
## Plot flux countour
###############################################################################

plt.figure(figsize=(8,12));

sLvL = 32;

cs = plt.contour(Rgrid, Zgrid, psirz2d.T, levels=sLvL, colors='k');
plt.clabel(cs, inline=True, fontsize=8);

plt.contourf(Rgrid, Zgrid, psirz2d.T, levels=sLvL, cmap="turbo");
plt.colorbar();

plt.xlabel("R (m)", fontsize=14);
plt.ylabel("Z (m)", fontsize=14);
plt.title ("$\Psi$ (poloidal flux)",fontsize=14);

plt.grid();
plt.tight_layout();
plt.show();

###############################################################################
## Normalize such that 0 at magnetic axis, 1 at LCFS.
###############################################################################

psirz2dnorm = np.divide(psirz2d-parameters_0D['sibry2'],parameters_0D['simag2']);

## Normalize such that minimum is 0. Divide by abs of min value
min_psi_abs = np.abs(np.min(psirz2dnorm)); ### rescaling.
min_psi     = (np.min(psirz2dnorm));       ### rescaling.
max_psi     = (np.max(psirz2dnorm));       ### rescaling.
psirz2dnorm =  psirz2dnorm - max_psi;      ### shifting..
psirz2dnorm = -psirz2dnorm / max_psi;

###############################################################################
## Plot normalized flux
###############################################################################

plt.figure(figsize=(8,12));

sLvL = 32;

cs = plt.contour(Rgrid, Zgrid, psirz2dnorm.T, levels=sLvL, colors='k');
plt.clabel(cs, inline=True, fontsize=8);

plt.contourf(Rgrid, Zgrid, psirz2dnorm.T, levels=sLvL, cmap="turbo");
plt.colorbar()

plt.xlabel("R (m)", fontsize=14);
plt.ylabel("Z (m)", fontsize=14);
plt.title ("$\Psi_{N}$ (normalized poloidal flux)", fontsize=14);

plt.grid();
plt.tight_layout();
plt.show();

###############################################################################
## Getting B_R, B_T, and B_Z
###############################################################################

B_R = np.zeros_like(psirz2dnorm)
B_T = np.zeros_like(psirz2dnorm)
B_Z = np.zeros_like(psirz2dnorm)

interp_order     = 3
interp_smoothing = 0

## Flatten psi2d - joins up the rows tgt into one array
psi2d_flatten = psirz2d.flatten()

## Should get a 129*129 size
Zin_new = np.array(list(Zgrid)*len(Zgrid))

Rin_new = []
for i in range(len(Rgrid)):
    Rin_join = [Rgrid[i]]*len(Rgrid)
    Rin_new = Rin_new + Rin_join

interp_psirz2d = RectBivariateSpline(
    Rgrid,
    Zgrid,
    psirz2d,
    kx=interp_order,
    ky=interp_order,
    s=0
    ) 

## From 0 to 1(LCFS) to max.
npsi = len(fpol)
pol_flux_array_EQDSK_norm = np.linspace(0,1,npsi) 

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
dpolflux_dR = interp_psirz2d(Rgrid, Zgrid, dx=1, dy=0) 
dpolflux_dZ = interp_psirz2d(Rgrid, Zgrid, dx=0, dy=1)

for i in range(int(nw)):
    B_R[i,:] = -1*dpolflux_dZ[i,:] / Rgrid[i]
    B_Z[i,:] =    dpolflux_dR[i,:] / Rgrid[i]

for i in range(int(nw)): #R
    for j in range(int(nh)): #Z
        B_T[i,j] = interp_rBphi( psirz2dnorm[i,j] ) / Rgrid[i]

interp_B_T_grid = RectBivariateSpline(
    Rgrid,
    Zgrid,
    B_T,
    kx = interp_order,
    ky = interp_order,
    s  = 0
    ) 

###############################################################################
## Plot B_R, B_Z
###############################################################################

plt.figure(figsize=(8,12));

Bmag = np.sqrt( B_R**2 + B_Z**2 );
plt.pcolormesh( Rgrid, Zgrid, Bmag.T, shading='auto', cmap='plasma' );
plt.colorbar();

plt.streamplot( Rgrid, Zgrid, B_R.T, B_Z.T,color='white',density=1.3 );

plt.xlabel("R (m)", fontsize=14);
plt.ylabel("Z (m)", fontsize=14);
plt.title('|B| poloidal (T)', fontsize=14);

plt.show();

###############################################################################
## Getting J_R, J_T, and J_Z
###############################################################################

Fgrid = np.zeros_like( psirz2dnorm )

interp_psirz2dnorm = RectBivariateSpline(
    Rgrid,
    Zgrid,
    psirz2dnorm,
    kx = interp_order,
    ky = interp_order,
    s  = 0
    ) 

for i in range(int(nw)): 
    for j in range(int(nh)): 
        Fgrid[i,j] = interp_rBphi( interp_psirz2dnorm( Rgrid[i], Zgrid[j] ) ) / mu_0
        
interp_Fgrid = RectBivariateSpline(
    Rgrid,
    Zgrid,
    Fgrid,
    kx = interp_order,
    ky = interp_order,
    s  = 0
    ) 

## First derivatives
dF_dR = interp_Fgrid( Rgrid, Zgrid, dx=1, dy=0 ) 
dF_dZ = interp_Fgrid( Rgrid, Zgrid, dx=0, dy=1 )

J_R = np.zeros_like( psirz2dnorm )
J_Z = np.zeros_like( psirz2dnorm )

for i in range(int(nw)):
    J_R[i,:] = -1*dF_dZ[i,:] / Rgrid[i]
    J_Z[i,:] =    dF_dR[i,:] / Rgrid[i]
   
J_T = np.zeros_like( psirz2dnorm )

interp_pprime = UnivariateSpline(
    pol_flux_array_EQDSK_norm, 
    pprime,                    
    w    = None,
    bbox = [None, None],
    k    = interp_order,
    s    = interp_smoothing,
    ext  = 0,
    check_finite = False,
)

interp_ffprim = UnivariateSpline(
    pol_flux_array_EQDSK_norm, 
    ffprim,                    
    w    = None,
    bbox = [None, None],
    k    = interp_order,
    s    = interp_smoothing,
    ext  = 0,
    check_finite = False,
)

for i in range(int(nw)):
    for j in range(int(nh)):
        J_T[i,j] = ( interp_pprime( interp_psirz2dnorm( Rgrid[i], Zgrid[j] ) ) * ( Rgrid[i] ) + 
                     interp_ffprim( interp_psirz2dnorm( Rgrid[i], Zgrid[j] ) ) / ( Rgrid[i] * mu_0 ) );
        
###############################################################################
## Plot J_R, J_Z
###############################################################################

plt.figure( figsize=(8,12) );

Jmag = np.sqrt( J_R**2 + J_Z**2 );

plt.pcolormesh( Rgrid, Zgrid, Jmag.T, shading='auto', cmap='plasma' );
plt.colorbar();

plt.streamplot( Rgrid, Zgrid, J_R.T, J_Z.T, color='white', density=1.3 );

plt.xlabel("R (m)", fontsize=14);
plt.ylabel("Z (m)", fontsize=14);
plt.title('|J| poloidal (T)', fontsize=14);

plt.show();

###############################################################################
## Get electron density ne
###############################################################################

ne_dens          = np.loadtxt( Norm_ne_file, dtype = float );
psi              = np.loadtxt( Norm_rp_file, dtype = float )**2;

ne_interp        = UnivariateSpline( psi, ne_dens, s=0 );
pol_flux_interp  = RectBivariateSpline(Rgrid, Zgrid, psirz2dnorm);

ne_2D            = np.zeros( np.shape( psirz2dnorm ) );

for i in range(len(Rgrid)):
    for j in range(len(Zgrid)):
        ne_2D[i,j] = ne_interp(pol_flux_interp(Rgrid[i],Zgrid[j]))[0][0];
 
###############################################################################
## Plot ne
###############################################################################

plt.figure(figsize=(8,12));

plt.contourf(Rgrid, Zgrid, ne_2D.T, levels=np.linspace(0,ne_dens.max()))
plt.colorbar();

plt.xlabel("R (m)", fontsize=14);
plt.ylabel("Z (m)", fontsize=14);
plt.title('Electron density ($m^{-3}$)',fontsize=14); 
 
plt.show();

plt.plot( np.sqrt(psi), ne_dens ); plt.grid(); plt.show();

###############################################################################
## Get electron temperature Te
###############################################################################

te_dens   = np.loadtxt( Norm_te_file, dtype = float );
te_interp = UnivariateSpline( psi, te_dens, s=0 );
te_2D     = np.zeros( np.shape( psirz2dnorm ) );

for i in range(len(Rgrid)):
    for j in range(len(Zgrid)):
        te_2D[i,j] = te_interp(pol_flux_interp(Rgrid[i],Zgrid[j]))[0][0];
  
###############################################################################
## Plot Te
###############################################################################

plt.figure(figsize=(8,12));

plt.contourf(Rgrid, Zgrid, te_2D.T, levels=np.linspace(0,te_dens.max()))
plt.colorbar();

plt.xlabel("R (m)", fontsize=14);
plt.ylabel("Z (m)", fontsize=14);
plt.title('Electron temperature (eV)',fontsize=14);

plt.show();

plt.plot( np.sqrt(psi), te_dens ); plt.grid(); plt.show();

###############################################################################











