###############################################################################
# - This script reads an EQDSK equilibrium file and generates 
#   the ne and Be plasma files for ERMES. 
#
# - It requires an EQDSK plasma equilibrium file and the list 
#   of nodes of the problem.
###############################################################################

import numpy as np
import re

from scipy.interpolate import UnivariateSpline, RectBivariateSpline

###############################################################################
## Files names
###############################################################################

EQDSK_File = './EQDSK_Samples/mast-u-sample.eqdsk'

###############################################################################

Norm_ne_File = './EQDSK_Samples/rho_ne.txt'
Norm_rp_File = './EQDSK_Samples/rho_flux.txt'

###############################################################################

Problem_Name   = 'MASTU-2D-t20_p20';
Problem_Folder = './';
ERMES_Folder   = Problem_Folder + Problem_Name + '.gid/';

ExtrFormat3D   = True;

###############################################################################

Plasma_Name   = '2D-t20_p20'; 
Plasma_Folder = './';

ne_3D_File = Plasma_Folder + 'ne_' + Plasma_Name + '.dat';
Be_3D_File = Plasma_Folder + 'Be_' + Plasma_Name + '.dat';

ExtrFormat3D = True;

###############################################################################

print('EQDSK file: ' + EQDSK_File  );
print('ERMES file: ' + ERMES_Folder);
print('ne    file: ' + ne_3D_File  );
print('Be    file: ' + Be_3D_File  );

###############################################################################
## Read flux surfaces eqilibrium file
###############################################################################

def find_nearest_index(array, value):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return idx

def preprocess_line(line):
    ## Add spaces before the negative sign if not preceded by 'e' or whitespace
    return re.sub(r'(?<![\sEe])([+-])(?=\d)', r' \1', line)

## Open the file
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
## Normalize such that 0 at magnetic axis, 1 at LCFS.
###############################################################################

psirz2dnorm = np.divide(psirz2d-parameters_0D['sibry2'],parameters_0D['simag2']);

## Normalize such that minimum is 0. Divide by abs of min value
min_psi_abs = np.abs( np.min( psirz2dnorm ) ); ### rescaling.
min_psi     = np.min( psirz2dnorm );           ### rescaling.
max_psi     = np.max( psirz2dnorm );           ### rescaling.

psirz2dnorm =  psirz2dnorm - max_psi;          ### shifting..
psirz2dnorm = -psirz2dnorm / max_psi;

###############################################################################
## Getting B_R, B_T, and B_Z
###############################################################################

B_R = np.zeros_like( psirz2dnorm );
B_T = np.zeros_like( psirz2dnorm );
B_Z = np.zeros_like( psirz2dnorm );

interp_order     = 3;
interp_smoothing = 0;

## Flatten psi2d - joins up the rows tgt into one array
psi2d_flatten = psirz2d.flatten();

## Should get a 129*129 size
Zin_new = np.array( list(Zgrid)*len(Zgrid) );

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

interp_B_T_grid = RectBivariateSpline(
    Rgrid,
    Zgrid,
    B_T,
    kx = interp_order,
    ky = interp_order,
    s  = 0
    ) 

###############################################################################
## Transform Fr, Fz, Ft to cartesian Fx, Fy, Fz
###############################################################################

def Transform_FRZT_2_FXYZ( nX, nY, nZ, Fr, Fz, Ft ):
    
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
## Parse GiD mesh nodes and write ERMES ne and Be files
###############################################################################

ne_data   =   np.loadtxt( Norm_ne_File, dtype = float );
ro_data   = ( np.loadtxt( Norm_rp_File, dtype = float ) )**2;

ne_interp = UnivariateSpline   ( ro_data, ne_data, k=5, s=0 );
Ps_spline = RectBivariateSpline( Rgrid, Zgrid, psirz2dnorm, kx=5, ky=5, s=0 );

Br_spline = RectBivariateSpline( Rgrid, Zgrid, B_R, kx=5, ky=5, s=0 );
Bt_spline = RectBivariateSpline( Rgrid, Zgrid, B_T, kx=5, ky=5, s=0 );
Bz_spline = RectBivariateSpline( Rgrid, Zgrid, B_Z, kx=5, ky=5, s=0 );

## Parse GiD mesh nodes:  No[ID] = p(x,y,z);
Id, Xs, Ys, Zs = [], [], [], [];

pat = re.compile(
    r"No\[(\d+)\]\s*=\s*p\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^)]+)\s*\)\s*;"
);

with open( ERMES_Folder + Problem_Name + '-1.dat', 'r' ) as f:
    for line in f:
        m = pat.search( line );
        if m:
            Id.append( int  ( m.group(1) ) );
            Xs.append( float( m.group(2) ) );
            Ys.append( float( m.group(3) ) );
            Zs.append( float( m.group(4) ) );
            
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
if ExtrFormat3D:
    cBx, cBy, cBz = Transform_FRZT_2_FXYZ( nX, nY, nZ, Br, Bz, Bt );
else:
    cBx = Br; cBy = Bz; cBz = Bt;

## Write ouputs ( sorted by NodeID )
order      = np.argsort( nId );
nid_sorted = nId[ order ];

np.savetxt( 
    ne_3D_File, 
    np.column_stack( [ nid_sorted, ne[ order ] ] ),
    fmt=["%d", "%.8e"]
);

np.savetxt(
    Be_3D_File,
    np.column_stack( [ nid_sorted, cBx[ order ], cBy[ order ], cBz[ order ] ] ),
    fmt=["%d", "%.8e", "%.8e", "%.8e"]
);

###############################################################################