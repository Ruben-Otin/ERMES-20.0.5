###############################################################################
# - This script reads an EQDSK equilibrium file and generates 
#   imported J source files for ERMES tetrahedral elements. 
#
# - It requieres a EQDSK plasma equilibrium file and the list 
#   of elements of the problem.
###############################################################################

import numpy as np
import re
import os

from scipy.constants   import mu_0
from scipy.interpolate import UnivariateSpline, RectBivariateSpline

###############################################################################
## Input file names
###############################################################################

EQDSK_File = './EQDSK_Samples/mast-u-sample.eqdsk';

###############################################################################

Problem_Name   = 'MASTU-2D-t20_p20';
Problem_Folder = './';
ERMES_Folder   = Problem_Folder + Problem_Name + '.gid/';

###############################################################################
## Extrusion format (cylindrical or axysimmetric) and file format
############################################################################### 

ExtrFormat3D = True; 
VariablPhase = True;

###############################################################################
## Output file names
###############################################################################

Jsources_FileName = 'Jsource-1.dat'
Jsources_Folder   = ERMES_Folder    + 'Export_J_Sources/';
Je_3D_File        = Jsources_Folder + Jsources_FileName  ;

os.makedirs( Jsources_Folder, exist_ok=True );

###############################################################################
## Current parameters
###############################################################################

## Global current phase
Jphase = 0.0;

# J cut option
JCutOption = 1;

## Option 1 - by coord
RJ_max = 6.8;
ZJ_max = 7.5;
Ps_max = 1.0;

## Option 2 - by material Id
Mat_J_Id = 2;
PsMJ_max = 1.0;

###############################################################################

print( 'EQDSK file: ' + EQDSK_File   );
print( 'ERMES file: ' + ERMES_Folder );
print( 'Je    file: ' + Je_3D_File   );

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
## Getting J_R, J_T, and J_Z
###############################################################################

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
## Parse GiD mesh nodes:  No[ID] = p(x,y,z);
###############################################################################

Ps_spline = RectBivariateSpline( Rgrid, Zgrid, psirz2dnorm, kx=5, ky=5, s=0 );
Jr_spline = RectBivariateSpline( Rgrid, Zgrid, J_R, kx=5, ky=5, s=0 );
Jt_spline = RectBivariateSpline( Rgrid, Zgrid, J_T, kx=5, ky=5, s=0 );
Jz_spline = RectBivariateSpline( Rgrid, Zgrid, J_Z, kx=5, ky=5, s=0 );

Id, Xs, Ys, Zs = [], [], [], [];

patno = re.compile(r"No\[(\d+)\]\s*=\s*p\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^)]+)\s*\)\s*;");

with open( ERMES_Folder + Problem_Name + '-1.dat', "r" ) as f:
    for line in f:
        m = patno.search( line );
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
Jr = Jr_spline.ev( Rnod, Znod );
Jt = Jt_spline.ev( Rnod, Znod );
Jz = Jz_spline.ev( Rnod, Znod );

## Transform to cartesian coordinates
if ExtrFormat3D:    
    cJx, cJy, cJz = Transform_FRZT_2_FXYZ( nX, nY, nZ, Jr, Jz, Jt );
else:    
    cJx = Jr; cJy = Jz; cJz = Jt;

###############################################################################
## Parse GiD mesh thetrahedra:  VE(Id1,Id2,Id3,Id4,mId); 
###############################################################################

patve = re.compile(r"VE\(\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)\s*\);")

Jfile = open( Je_3D_File, 'w' );

if( VariablPhase == True ):
    Jfile.write( 'V' + '\n' );
else:
    Jfile.write( str( Jphase ) + '\n' );

with open( ERMES_Folder + Problem_Name + '-4.dat', "r" ) as f:
    
    for line in f:       
        m = patve.search( line );
        
        if m:
            nId = np.zeros( 4, dtype=np.int64 );
            
            nId[0] = int( m.group(1) );
            nId[1] = int( m.group(2) );
            nId[2] = int( m.group(3) );
            nId[3] = int( m.group(4) );
            eMatId = int( m.group(5) );
            
            avgJx = 0.0; avgJy = 0.0; avgJz = 0.0;
            
            ## Cut option 2
            if( JCutOption == 2 and eMatId == Mat_J_Id):
                for i in range( len(nId) ): 
                    vId = nId[i] - 1;
                    
                    if( Ps[vId] <= PsMJ_max ):
                        
                        avgJx += cJx[vId]; 
                        avgJy += cJy[vId]; 
                        avgJz += cJz[vId];
             
            ## Cut option 1
            elif( JCutOption == 1 ):
                for i in range( len(nId) ):
                    vId = nId[i] - 1;
                    
                    if( ( Ps[vId]        <= Ps_max                     ) and
                        ( abs( nY[vId] ) <  ZJ_max                     ) and
                        ( np.sqrt( nX[vId]**2 + nZ[vId]**2 ) <  RJ_max ) ):
                        
                        avgJx += cJx[vId]; 
                        avgJy += cJy[vId]; 
                        avgJz += cJz[vId];
            
            ## Do not write in file cero J elements
            if( avgJx == 0.0 and avgJy == 0.0 and avgJz == 0.0 ): continue;
            
            ## Averaging the value of each node
            avgJx *= 0.25; avgJy *= 0.25; avgJz *= 0.25;
            
            if( VariablPhase == True ):
                Jfile.write(
                    '%i %i %i %i %.8f %.8f %.8f %.8f %.8f %.8f\n'%(
                        nId[0], nId[1], nId[2], nId[3], 
                        avgJx, Jphase,
                        avgJy, Jphase,
                        avgJz, Jphase
                    )
                );
            else:    
                Jfile.write(
                    '%i %i %i %i %.8f %.8f %.8f\n'%(
                        nId[0], nId[1], nId[2], nId[3], 
                        avgJx, avgJy, avgJz
                    )   
                );
            
Jfile.close();