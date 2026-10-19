###############################################################################
# - Plane wave Robin boundary condition for EDG_1st and RME_1st elements.
# - This script generates the files that can be imported into ERMES as a
#   plane wave Robin boundary condiction.
###############################################################################

import numpy as np
import re
from pathlib import Path

###############################################################################

class PlaneWave:  
    freq = 0.0;                    # Frequency

    modE = 0.0; phaE = 0.0;        # Module and phase of E

    ex = 0.0; ey = 0.0; ez = 0.0;  # Polarization vector
    kx = 0.0; ky = 0.0; kz = 0.0;  # Wave vector
    
    sgm_r = 0.0; sgm_i = 0.0;      # Conductivity 
    epr_r = 1.0; epr_i = 0.0;      # Relative permittivity
    mur_r = 1.0; mur_i = 0.0;      # Relative permeability
    
###############################################################################

class Node:   
    ID = 0;   # Node ID

    X = 0.0;  # Node cartesian coordinate X
    Y = 0.0;  # Node cartesian coordinate Y
    Z = 0.0;  # Node cartesian coordinate Z
    
###############################################################################

def Get_Element_Nodes_List( NodesIds, nXs, nYs, nZs ):
    """
    - Returns a list of element nodes objects.
    """
    eNodes = [];
    
    for i in range( len( NodesIds ) - 1 ):   
        Ni    = Node();  
        Ni.ID = NodesIds[ i ];
        
        Ni.X = nXs[ Ni.ID - 1 ];
        Ni.Y = nYs[ Ni.ID - 1 ];
        Ni.Z = nZs[ Ni.ID - 1 ];
        
        eNodes.append( Ni );
        
    return eNodes;

###############################################################################

def Read_Nodes_File( FolderName = '..' ):
    """
    - Reads nodes Ids and nodes coordinates from *-1.dat file.
    - Returns vectors: nId, nXc, nYc, nZc.
    """
    NodesFile = next( Path(FolderName).glob("*-1.dat"), None );
    
    nId, nXc, nYc, nZc = [], [], [], [];
    
    pat = re.compile(r"No\[(\d+)\]\s*=\s*p\(\s*([^,]+)\s*,\s*([^,]+)\s*,\s*([^)]+)\s*\)\s*;");
    
    with open( NodesFile, 'r' ) as f:
        for line in f:
            m = pat.search( line );
            if m:
                nId.append( int  ( m.group(1) ) );
                nXc.append( float( m.group(2) ) );
                nYc.append( float( m.group(3) ) );
                nZc.append( float( m.group(4) ) );
                
    nId = np.asarray( nId, dtype=np.int64 );
    nXc = np.asarray( nXc, float );
    nYc = np.asarray( nYc, float );
    nZc = np.asarray( nZc, float );

    return nId, nXc, nYc, nZc;

###############################################################################

def Read_IRBC_Elemets_File( FolderName = '..' ):
    """
    - Reads IRBC surface elements from *-19.dat file.
    - Returns element list: [ Node1, Node2, Node3, SurfaceId ].
    """
    ElementsFile = next( Path(FolderName).glob("*-19.dat"), None );
    
    ElementList = [];
    
    with open( ElementsFile, 'r' ) as f:
        for line in f:
            iIds = list( map( int, line.split() ) );
            
            ElementList.append( [ iIds[0], iIds[1], iIds[2], iIds[3] ] );
  
    ElementList = np.asarray( ElementList, dtype=np.int64 );
    
    return ElementList;

###############################################################################

def GaussPoints2D_Order12():
    """
    - Returns 2D order 12 Gauss points coordinates and weights 
    """	
    numGaussPoints = 144;
    
    cX = np.zeros( numGaussPoints, dtype=np.dtype( np.float64 ) );
    cY = np.zeros( numGaussPoints, dtype=np.dtype( np.float64 ) );

    cX[  0] = 9.7843793683414970e-001; cY[  0] = 1.9879538455522941e-004;
    cX[  1] = 9.2894210126441101e-001; cY[  1] = 6.5513129222256026e-004;
    cX[  2] = 8.5455254376493583e-001; cY[  2] = 1.3409794217013301e-003;
    cX[  3] = 7.5959888952522703e-001; cY[  3] = 2.2164220017696101e-003;
    cX[  4] = 6.4960065027725500e-001; cY[  4] = 3.2305708846247349e-003;
    cX[  5] = 5.3095084931281777e-001; cY[  5] = 4.3244844228933522e-003;
    cX[  6] = 4.1054508120145772e-001; cY[  6] = 5.4345874213983974e-003;
    cX[  7] = 2.9538088426258025e-001; cY[  7] = 6.4963647959178127e-003;
    cX[  8] = 1.9215105452985404e-001; cY[  8] = 7.4481110894631360e-003;
    cX[  9] = 1.0685449088347654e-001; cY[  9] = 8.2345183567499140e-003;
    cX[ 10] = 4.4446463155407789e-002; cY[ 10] = 8.8099005813592913e-003;
    cX[ 11] = 8.5390549884274414e-003; cY[ 11] = 9.1409554975809375e-003;
    cX[ 12] = 9.7843793683414970e-001; cY[ 12] = 1.0337148873274264e-003;
    cX[ 13] = 9.2894210126441101e-001; cY[ 13] = 3.4066131436586238e-003;
    cX[ 14] = 8.5455254376493583e-001; cY[ 14] = 6.9729505788766223e-003;
    cX[ 15] = 7.5959888952522703e-001; cY[ 15] = 1.1525159021952914e-002;
    cX[ 16] = 6.4960065027725500e-001; cY[ 16] = 1.6798625508709152e-002;
    cX[ 17] = 5.3095084931281777e-001; cY[ 17] = 2.2486859732492814e-002;
    cX[ 18] = 4.1054508120145772e-001; cY[ 18] = 2.8259277430161614e-002;
    cX[ 19] = 2.9538088426258025e-001; cY[ 19] = 3.3780407015356881e-002;
    cX[ 20] = 1.9215105452985404e-001; cY[ 20] = 3.8729386664948143e-002;
    cX[ 21] = 1.0685449088347654e-001; cY[ 21] = 4.2818620937240690e-002;
    cX[ 22] = 4.4446463155407789e-002; cY[ 22] = 4.5810547398778052e-002;
    cX[ 23] = 8.5390549884274414e-003; cY[ 23] = 4.7531997804615700e-002;
    cX[ 24] = 9.7843793683414970e-001; cY[ 24] = 2.4806865366578204e-003;
    cX[ 25] = 9.2894210126441101e-001; cY[ 25] = 8.1751162382154730e-003;
    cX[ 26] = 8.5455254376493583e-001; cY[ 26] = 1.6733535362464594e-002;
    cX[ 27] = 7.5959888952522703e-001; cY[ 27] = 2.7657826320482415e-002;
    cX[ 28] = 6.4960065027725500e-001; cY[ 28] = 4.0312976667629133e-002;
    cX[ 29] = 5.3095084931281777e-001; cY[ 29] = 5.3963477622276650e-002;
    cX[ 30] = 4.1054508120145772e-001; cY[ 30] = 6.7816000249278960e-002;
    cX[ 31] = 2.9538088426258025e-001; cY[ 31] = 8.1065487121377039e-002;
    cX[ 32] = 1.9215105452985404e-001; cY[ 32] = 9.2941941003815823e-002;
    cX[ 33] = 1.0685449088347654e-001; cY[ 33] = 1.0275519660153920e-001;
    cX[ 34] = 4.4446463155407789e-002; cY[ 34] = 1.0993515674605736e-001;
    cX[ 35] = 8.5390549884274414e-003; cY[ 35] = 1.1406625604397524e-001;
    cX[ 36] = 9.7843793683414970e-001; cY[ 36] = 4.4491381685421369e-003;
    cX[ 37] = 9.2894210126441101e-001; cY[ 37] = 1.4662159507148631e-002;
    cX[ 38] = 8.5455254376493583e-001; cY[ 38] = 3.0011776891446976e-002;
    cX[ 39] = 7.5959888952522703e-001; cY[ 39] = 4.9604611031249100e-002;
    cX[ 40] = 6.4960065027725500e-001; cY[ 40] = 7.2301760230110687e-002;
    cX[ 41] = 5.3095084931281777e-001; cY[ 41] = 9.6784081522855489e-002;
    cX[ 42] = 4.1054508120145772e-001; cY[ 42] = 1.2162873087279910e-001;
    cX[ 43] = 2.9538088426258025e-001; cY[ 43] = 1.4539182906563652e-001;
    cX[ 44] = 1.9215105452985404e-001; cY[ 44] = 1.6669237772200932e-001;
    cX[ 45] = 1.0685449088347654e-001; cY[ 45] = 1.8429255791096372e-001;
    cX[ 46] = 4.4446463155407789e-002; cY[ 46] = 1.9716989418684217e-001;
    cX[ 47] = 8.5390549884274414e-003; cY[ 47] = 2.0457906551614960e-001;
    cX[ 48] = 9.7843793683414970e-001; cY[ 48] = 6.8154285750310689e-003;
    cX[ 49] = 9.2894210126441101e-001; cY[ 49] = 2.2460282664008218e-002;
    cX[ 50] = 8.5455254376493583e-001; cY[ 50] = 4.5973650191324156e-002;
    cX[ 51] = 7.5959888952522703e-001; cY[ 51] = 7.5987004824005083e-002;
    cX[ 52] = 6.4960065027725500e-001; cY[ 52] = 1.1075571583312008e-001;
    cX[ 53] = 5.3095084931281777e-001; cY[ 53] = 1.4825904924304636e-001;
    cX[ 54] = 4.1054508120145772e-001; cY[ 54] = 1.8631741621251197e-001;
    cX[ 55] = 2.9538088426258025e-001; cY[ 55] = 2.2271900508647624e-001;
    cX[ 56] = 1.9215105452985404e-001; cY[ 56] = 2.5534832844688155e-001;
    cX[ 57] = 1.0685449088347654e-001; cY[ 57] = 2.8230922883734993e-001;
    cX[ 58] = 4.4446463155407789e-002; cY[ 58] = 3.0203542350701656e-001;
    cX[ 59] = 8.5390549884274414e-003; cY[ 59] = 3.1338518970490681e-001;
    cX[ 60] = 9.7843793683414970e-001; cY[ 60] = 9.4308862505256327e-003;
    cX[ 61] = 9.2894210126441101e-001; cY[ 61] = 3.1079537937634207e-002;
    cX[ 62] = 8.5455254376493583e-001; cY[ 62] = 6.3616287765712168e-002;
    cX[ 63] = 7.5959888952522703e-001; cY[ 63] = 1.0514743000003748e-001;
    cX[ 64] = 6.4960065027725500e-001; cY[ 64] = 1.5325882240838168e-001;
    cX[ 65] = 5.3095084931281777e-001; cY[ 65] = 2.0515426339360834e-001;
    cX[ 66] = 4.1054508120145772e-001; cY[ 66] = 2.5781773507677480e-001;
    cX[ 67] = 2.9538088426258025e-001; cY[ 67] = 3.0818863108564265e-001;
    cX[ 68] = 1.9215105452985404e-001; cY[ 68] = 3.5333963423326181e-001;
    cX[ 69] = 1.0685449088347654e-001; cY[ 69] = 3.9064692635657489e-001;
    cX[ 70] = 4.4446463155407789e-002; cY[ 70] = 4.1794315520517711e-001;
    cX[ 71] = 8.5390549884274414e-003; cY[ 71] = 4.3364845573088551e-001;
    cX[ 72] = 9.7843793683414970e-001; cY[ 72] = 1.2131176915324662e-002;
    cX[ 73] = 9.2894210126441101e-001; cY[ 73] = 3.9978360797954825e-002;
    cX[ 74] = 8.5455254376493583e-001; cY[ 74] = 8.1831168469351970e-002;
    cX[ 75] = 7.5959888952522703e-001; cY[ 75] = 1.3525368047473552e-001;
    cX[ 76] = 6.4960065027725500e-001; cY[ 76] = 1.9714052731436327e-001;
    cX[ 77] = 5.3095084931281777e-001; cY[ 77] = 2.6389488729357380e-001;
    cX[ 78] = 4.1054508120145772e-001; cY[ 78] = 3.3163718372176743e-001;
    cX[ 79] = 2.9538088426258025e-001; cY[ 79] = 3.9643048465177699e-001;
    cX[ 80] = 1.9215105452985404e-001; cY[ 80] = 4.5450931123688404e-001;
    cX[ 81] = 1.0685449088347654e-001; cY[ 81] = 5.0249858275994841e-001;
    cX[ 82] = 4.4446463155407789e-002; cY[ 82] = 5.3761038163941499e-001;
    cX[ 83] = 8.5390549884274414e-003; cY[ 83] = 5.5781248928068683e-001;
    cX[ 84] = 9.7843793683414970e-001; cY[ 84] = 1.4746634590819232e-002;
    cX[ 85] = 9.2894210126441101e-001; cY[ 85] = 4.8597616071580835e-002;
    cX[ 86] = 8.5455254376493583e-001; cY[ 86] = 9.9473806043740023e-002;
    cX[ 87] = 7.5959888952522703e-001; cY[ 87] = 1.6441410565076797e-001;
    cX[ 88] = 6.4960065027725500e-001; cY[ 88] = 2.3964363388962495e-001;
    cX[ 89] = 5.3095084931281777e-001; cY[ 89] = 3.2079010144413594e-001;
    cX[ 90] = 4.1054508120145772e-001; cY[ 90] = 4.0313750258603037e-001;
    cX[ 91] = 2.9538088426258025e-001; cY[ 91] = 4.8190011065094357e-001;
    cX[ 92] = 1.9215105452985404e-001; cY[ 92] = 5.5250061702326447e-001;
    cX[ 93] = 1.0685449088347654e-001; cY[ 93] = 6.1083628027917358e-001;
    cX[ 94] = 4.4446463155407789e-002; cY[ 94] = 6.5351811333757581e-001;
    cX[ 95] = 8.5390549884274414e-003; cY[ 95] = 6.7807575530666586e-001;
    cX[ 96] = 9.7843793683414970e-001; cY[ 96] = 1.7112924997308165e-002;
    cX[ 97] = 9.2894210126441101e-001; cY[ 97] = 5.6395739228440420e-002;
    cX[ 98] = 8.5455254376493583e-001; cY[ 98] = 1.1543567934361722e-001;
    cX[ 99] = 7.5959888952522703e-001; cY[ 99] = 1.9079649944352398e-001;
    cX[100] = 6.4960065027725500e-001; cY[100] = 2.7809758949263436e-001;
    cX[101] = 5.3095084931281777e-001; cY[101] = 3.7226506916432683e-001;
    cX[102] = 4.1054508120145772e-001; cY[102] = 4.6782618792574326e-001;
    cX[103] = 2.9538088426258025e-001; cY[103] = 5.5922728667178334e-001;
    cX[104] = 1.9215105452985404e-001; cY[104] = 6.4115656774813679e-001;
    cX[105] = 1.0685449088347654e-001; cY[105] = 7.0885295120555991e-001;
    cX[106] = 4.4446463155407789e-002; cY[106] = 7.5838364265775016e-001;
    cX[107] = 8.5390549884274414e-003; cY[107] = 7.8688187949542310e-001;
    cX[108] = 9.7843793683414970e-001; cY[108] = 1.9081376629192476e-002;
    cX[109] = 9.2894210126441101e-001; cY[109] = 6.2882782497373571e-002;
    cX[110] = 8.5455254376493583e-001; cY[110] = 1.2871392087259956e-001;
    cX[111] = 7.5959888952522703e-001; cY[111] = 2.1274328415429061e-001;
    cX[112] = 6.4960065027725500e-001; cY[112] = 3.1008637305511583e-001;
    cX[113] = 5.3095084931281777e-001; cY[113] = 4.1508567306490557e-001;
    cX[114] = 4.1054508120145772e-001; cY[114] = 5.2163891854926325e-001;
    cX[115] = 2.9538088426258025e-001; cY[115] = 6.2355362861604269e-001;
    cX[116] = 1.9215105452985404e-001; cY[116] = 7.1490700446633004e-001;
    cX[117] = 1.0685449088347654e-001; cY[117] = 7.9039031251498426e-001;
    cX[118] = 4.4446463155407789e-002; cY[118] = 8.4561838009853485e-001;
    cX[119] = 8.5390549884274414e-003; cY[119] = 8.7739468896759731e-001;
    cX[120] = 9.7843793683414970e-001; cY[120] = 2.0528348278522873e-002;
    cX[121] = 9.2894210126441101e-001; cY[121] = 6.7651285591930427e-002;
    cX[122] = 8.5455254376493583e-001; cY[122] = 1.3847450565618755e-001;
    cX[123] = 7.5959888952522703e-001; cY[123] = 2.2887595145282014e-001;
    cX[124] = 6.4960065027725500e-001; cY[124] = 3.3360072421403586e-001;
    cX[125] = 5.3095084931281777e-001; cY[125] = 4.4656229095468941e-001;
    cX[126] = 4.1054508120145772e-001; cY[126] = 5.6119564136838074e-001;
    cX[127] = 2.9538088426258025e-001; cY[127] = 6.7083870872206286e-001;
    cX[128] = 1.9215105452985404e-001; cY[128] = 7.6911955880519789e-001;
    cX[129] = 1.0685449088347654e-001; cY[129] = 8.5032688817928281e-001;
    cX[130] = 4.4446463155407789e-002; cY[130] = 9.0974298944581422e-001;
    cX[131] = 8.5390549884274414e-003; cY[131] = 9.4392894720695697e-001;
    cX[132] = 9.7843793683414970e-001; cY[132] = 2.1363267781295069e-002;
    cX[133] = 9.2894210126441101e-001; cY[133] = 7.0402767443366487e-002;
    cX[134] = 8.5455254376493583e-001; cY[134] = 1.4410647681336283e-001;
    cX[135] = 7.5959888952522703e-001; cY[135] = 2.3818468847300342e-001;
    cX[136] = 6.4960065027725500e-001; cY[136] = 3.4716877883812025e-001;
    cX[137] = 5.3095084931281777e-001; cY[137] = 4.6472466626428888e-001;
    cX[138] = 4.1054508120145772e-001; cY[138] = 5.8402033137714393e-001;
    cX[139] = 2.9538088426258025e-001; cY[139] = 6.9812275094150189e-001;
    cX[140] = 1.9215105452985404e-001; cY[140] = 8.0040083438068288e-001;
    cX[141] = 1.0685449088347654e-001; cY[141] = 8.8491099075977353e-001;
    cX[142] = 4.4446463155407789e-002; cY[142] = 9.4674363626323288e-001;
    cX[143] = 8.5390549884274414e-003; cY[143] = 9.8231998951399158e-001;
    
    W = np.zeros( numGaussPoints, dtype=np.dtype( np.float64 ) );
    
    W[  0] = 1.8284741016797483e-005;
    W[  1] = 1.0484265146702545e-004;
    W[  2] = 2.9333645786585916e-004;
    W[  3] = 5.8678524201024478e-004;
    W[  4] = 9.5416018550717355e-004;
    W[  5] = 1.3353885606628456e-003;
    W[  6] = 1.6537158491122695e-003;
    W[  7] = 1.8326749408828461e-003;
    W[  8] = 1.8138115492054414e-003;
    W[  9] = 1.5710979332129132e-003;
    W[ 10] = 1.1187124356320090e-003;
    W[ 11] = 5.1102355005257921e-004;
    W[ 12] = 4.1448732115333276e-005;
    W[ 13] = 2.3766237492376082e-004;
    W[ 14] = 6.6494921916439883e-004;
    W[ 15] = 1.3301530649501910e-003;
    W[ 16] = 2.1629362914066790e-003;
    W[ 17] = 3.0271231443719506e-003;
    W[ 18] = 3.7487227826615782e-003;
    W[ 19] = 4.1543958762858026e-003;
    W[ 20] = 4.1116354309666658e-003;
    W[ 21] = 3.5614405093773315e-003;
    W[ 22] = 2.5359512620936526e-003;
    W[ 23] = 1.1584128105121938e-003;
    W[ 24] = 6.2044937122139279e-005;
    W[ 25] = 3.5575870131352181e-004;
    W[ 26] = 9.9536778055530619e-004;
    W[ 27] = 1.9911167137275374e-003;
    W[ 28] = 3.2377165561084968e-003;
    W[ 29] = 4.5313247370490769e-003;
    W[ 30] = 5.6114929810488949e-003;
    W[ 31] = 6.2187482649022569e-003;
    W[ 32] = 6.1547398138413608e-003;
    W[ 33] = 5.3311486550106611e-003;
    W[ 34] = 3.7960856357100482e-003;
    W[ 35] = 1.7340373594472681e-003;
    W[ 36] = 7.8745888534725199e-005;
    W[ 37] = 4.5152008106229280e-004;
    W[ 38] = 1.2632959905232383e-003;
    W[ 39] = 2.5270757304526204e-003;
    W[ 40] = 4.1092291951631086e-003;
    W[ 41] = 5.7510444721038433e-003;
    W[ 42] = 7.1219670982855090e-003;
    W[ 43] = 7.8926803766357703e-003;
    W[ 44] = 7.8114424451249456e-003;
    W[ 45] = 6.7661610636030739e-003;
    W[ 46] = 4.8178973209279021e-003;
    W[ 47] = 2.2007970183494595e-003;
    W[ 48] = 9.0499631523127739e-005;
    W[ 49] = 5.1891472331804152e-004;
    W[ 50] = 1.4518576623410843e-003;
    W[ 51] = 2.9042712793335185e-003;
    W[ 52] = 4.7225796155992011e-003;
    W[ 53] = 6.6094549859451240e-003;
    W[ 54] = 8.1850038155382790e-003;
    W[ 55] = 9.0707550464729902e-003;
    W[ 56] = 8.9773914054723739e-003;
    W[ 57] = 7.7760895771997783e-003;
    W[ 58] = 5.5370247307320059e-003;
    W[ 59] = 2.5292916611131273e-003;
    W[ 60] = 9.6567180156891032e-005;
    W[ 61] = 5.5370536574959005e-004;
    W[ 62] = 1.5491974727612535e-003;
    W[ 63] = 3.0989881741586064e-003;
    W[ 64] = 5.0392050096721534e-003;
    W[ 65] = 7.0525859567009850e-003;
    W[ 66] = 8.7337674721684726e-003;
    W[ 67] = 9.6789039025858871e-003;
    W[ 68] = 9.5792806954095588e-003;
    W[ 69] = 8.2974375749328075e-003;
    W[ 70] = 5.9082546050932168e-003;
    W[ 71] = 2.6988680439612046e-003;
    W[ 72] = 9.6567180156891005e-005;
    W[ 73] = 5.5370536574958994e-004;
    W[ 74] = 1.5491974727612531e-003;
    W[ 75] = 3.0989881741586055e-003;
    W[ 76] = 5.0392050096721517e-003;
    W[ 77] = 7.0525859567009833e-003;
    W[ 78] = 8.7337674721684692e-003;
    W[ 79] = 9.6789039025858836e-003;
    W[ 80] = 9.5792806954095554e-003;
    W[ 81] = 8.2974375749328040e-003;
    W[ 82] = 5.9082546050932151e-003;
    W[ 83] = 2.6988680439612037e-003;
    W[ 84] = 9.0499631523128146e-005;
    W[ 85] = 5.1891472331804379e-004;
    W[ 86] = 1.4518576623410908e-003;
    W[ 87] = 2.9042712793335315e-003;
    W[ 88] = 4.7225796155992220e-003;
    W[ 89] = 6.6094549859451535e-003;
    W[ 90] = 8.1850038155383154e-003;
    W[ 91] = 9.0707550464730318e-003;
    W[ 92] = 8.9773914054724138e-003;
    W[ 93] = 7.7760895771998139e-003;
    W[ 94] = 5.5370247307320311e-003;
    W[ 95] = 2.5292916611131385e-003;
    W[ 96] = 7.8745888534724928e-005;
    W[ 97] = 4.5152008106229122e-004;
    W[ 98] = 1.2632959905232337e-003;
    W[ 99] = 2.5270757304526113e-003;
    W[100] = 4.1092291951630939e-003;
    W[101] = 5.7510444721038225e-003;
    W[102] = 7.1219670982854838e-003;
    W[103] = 7.8926803766357426e-003;
    W[104] = 7.8114424451249179e-003;
    W[105] = 6.7661610636030496e-003;
    W[106] = 4.8178973209278848e-003;
    W[107] = 2.2007970183494517e-003;
    W[108] = 6.2044937122139130e-005;
    W[109] = 3.5575870131352100e-004;
    W[110] = 9.9536778055530403e-004;
    W[111] = 1.9911167137275326e-003;
    W[112] = 3.2377165561084894e-003;
    W[113] = 4.5313247370490665e-003;
    W[114] = 5.6114929810488819e-003;
    W[115] = 6.2187482649022430e-003;
    W[116] = 6.1547398138413470e-003;
    W[117] = 5.3311486550106490e-003;
    W[118] = 3.7960856357100396e-003;
    W[119] = 1.7340373594472642e-003;
    W[120] = 4.1448732115333527e-005;
    W[121] = 2.3766237492376226e-004;
    W[122] = 6.6494921916440274e-004;
    W[123] = 1.3301530649501988e-003;
    W[124] = 2.1629362914066920e-003;
    W[125] = 3.0271231443719684e-003;
    W[126] = 3.7487227826616008e-003;
    W[127] = 4.1543958762858268e-003;
    W[128] = 4.1116354309666901e-003;
    W[129] = 3.5614405093773527e-003;
    W[130] = 2.5359512620936678e-003;
    W[131] = 1.1584128105122007e-003;
    W[132] = 1.8284741016797375e-005;
    W[133] = 1.0484265146702483e-004;
    W[134] = 2.9333645786585737e-004;
    W[135] = 5.8678524201024132e-004;
    W[136] = 9.5416018550716780e-004;
    W[137] = 1.3353885606628376e-003;
    W[138] = 1.6537158491122595e-003;
    W[139] = 1.8326749408828350e-003;
    W[140] = 1.8138115492054305e-003;
    W[141] = 1.5710979332129039e-003;
    W[142] = 1.1187124356320023e-003;
    W[143] = 5.1102355005257607e-004;
    
    return cX, cY, W;

###############################################################################

def Get_CoPlanar_XY_Axis( Nodes ):
    """
    - Returns X, Y axis of the co-planar coordinate system 
    """	
    # Z axis of the co-planar coordinate system 
    Z_axis = Calculate_Exterior_Normal( Nodes );
    
    # X axis of the co-planar coordinate system 
    X_axis = np.zeros( 3, dtype=np.dtype(np.float64) );
    
    X_axis[ 0 ] = Nodes[ 1 ].X - Nodes[ 0 ].X; 
    X_axis[ 1 ] = Nodes[ 1 ].Y - Nodes[ 0 ].Y; 
    X_axis[ 2 ] = Nodes[ 1 ].Z - Nodes[ 0 ].Z; 
    
    X_axis = X_axis / np.linalg.norm( X_axis );  
    
    # Y axis of the co-planar coordinate system ( Y_axis = X_axis x Z_axis )
    Y_axis = np.zeros( 3, dtype=np.dtype(np.float64) );
    
    Y_axis[ 0 ] = X_axis[ 1 ] * Z_axis[ 2 ] - X_axis[ 2 ] * Z_axis[ 1 ];
    Y_axis[ 1 ] = X_axis[ 2 ] * Z_axis[ 0 ] - X_axis[ 0 ] * Z_axis[ 2 ];
    Y_axis[ 2 ] = X_axis[ 0 ] * Z_axis[ 1 ] - X_axis[ 1 ] * Z_axis[ 0 ];
    
    Y_axis = Y_axis / np.linalg.norm( Y_axis ); 
    
    return X_axis, Y_axis;

###############################################################################

def Get_Vertices_Surface_Coordinates( Nodes ):
    """
    - Returns co-planar coordinates of vertice nodes (pX, pY) 
    """	
    # XY axis of the co-planar coordinate system 
    X_axis, Y_axis = Get_CoPlanar_XY_Axis( Nodes );

    # Position vectors of vertice nodes
    R_2 = np.zeros( 3, dtype=np.dtype(np.float64) );

    R_2[ 0 ] = Nodes[1].X - Nodes[0].X; 
    R_2[ 1 ] = Nodes[1].Y - Nodes[0].Y; 
    R_2[ 2 ] = Nodes[1].Z - Nodes[0].Z; 
    
    R_3 = np.zeros( 3, dtype=np.dtype(np.float64) );
    
    R_3[ 0 ] = Nodes[2].X - Nodes[0].X; 
    R_3[ 1 ] = Nodes[2].Y - Nodes[0].Y;
    R_3[ 2 ] = Nodes[2].Z - Nodes[0].Z;
    
    # Co-planar coordinates of vertice nodes
    pX = np.zeros( 3, dtype=np.dtype(np.float64) );
    pY = np.zeros( 3, dtype=np.dtype(np.float64) );
    
    pX[ 0 ] = 0.0; 
    pY[ 0 ] = 0.0;
    
    pX[ 1 ] = R_2[0]*X_axis[0] + R_2[1]*X_axis[1] + R_2[2]*X_axis[2]; 
    pY[ 1 ] = R_2[0]*Y_axis[0] + R_2[1]*Y_axis[1] + R_2[2]*Y_axis[2];
    
    pX[ 2 ] = R_3[0]*X_axis[0] + R_3[1]*X_axis[1] + R_3[2]*X_axis[2]; 
    pY[ 2 ] = R_3[0]*Y_axis[0] + R_3[1]*Y_axis[1] + R_3[2]*Y_axis[2];
    
    return pX, pY;

###############################################################################

def Calculate_DN( Area, Nodes ):
    """
    - Returns Lagrange 1st order base derivatives - DN[ X, Y ][ 0, 1, 2 ]
    """	
    Cte = 1.0 / ( 2.0 * Area ); 
        
    # Plane coordinates
    pX, pY = Get_Vertices_Surface_Coordinates( Nodes );
    
    # Derivatives
    DN = np.zeros((2, 3), dtype=np.dtype(np.float64));

    # dNi/dx
    DN[ 0 ][ 0 ] = Cte * ( pY[ 1 ] - pY[ 2 ] );
    DN[ 0 ][ 1 ] = Cte * ( pY[ 2 ] - pY[ 0 ] );
    DN[ 0 ][ 2 ] = Cte * ( pY[ 0 ] - pY[ 1 ] );
                           
    # dNi/dy              
    DN[ 1 ][ 0 ] = Cte * ( pX[ 2 ] - pX[ 1 ] ); 
    DN[ 1 ][ 1 ] = Cte * ( pX[ 0 ] - pX[ 2 ] );
    DN[ 1 ][ 2 ] = Cte * ( pX[ 1 ] - pX[ 0 ] );
    
    return DN;

###############################################################################

def Calculate_Edge_Lengths( Nodes ):
    """
    - Returns edge lengths with sign
    """	
    # Edge coordinates vector
    R = np.zeros((3, 3), dtype=np.dtype(np.float64));
		
    # Edge: [0 1]
    R[ 0 ][ 0 ] = Nodes[ 1 ].X - Nodes[ 0 ].X;
    R[ 0 ][ 1 ] = Nodes[ 1 ].Y - Nodes[ 0 ].Y; 
    R[ 0 ][ 2 ] = Nodes[ 1 ].Z - Nodes[ 0 ].Z;

	# Edge: [1 2]
    R[ 1 ][ 0 ] = Nodes[ 2 ].X - Nodes[ 1 ].X;
    R[ 1 ][ 1 ] = Nodes[ 2 ].Y - Nodes[ 1 ].Y; 
    R[ 1 ][ 2 ] = Nodes[ 2 ].Z - Nodes[ 1 ].Z;

    # Edge: [2 0]
    R[ 2 ][ 0 ] = Nodes[ 0 ].X - Nodes[ 2 ].X;
    R[ 2 ][ 1 ] = Nodes[ 0 ].Y - Nodes[ 2 ].Y; 
    R[ 2 ][ 2 ] = Nodes[ 0 ].Z - Nodes[ 2 ].Z;

	# Edge lengths
    eLength = np.zeros( 3, dtype=np.dtype(np.float64) );
    
    for i in range( 3 ): eLength[ i ] = np.linalg.norm( R[ i ] );
        
    # Check directions
    if( Nodes[ 0 ].ID > Nodes[ 1 ].ID ): eLength[ 0 ] *= -1.0;
    if( Nodes[ 1 ].ID > Nodes[ 2 ].ID ): eLength[ 1 ] *= -1.0;
    if( Nodes[ 2 ].ID > Nodes[ 0 ].ID ): eLength[ 2 ] *= -1.0;
        
    return eLength;

###############################################################################

def Whitney2D_Ni_1st( cX, cY, Area, Nodes ):
    """
    - Returns 1st order Whitney base for 2D triangles
    """	
    NumEdges  = 3;
    NumPoints = len( cX );
    
    Nx = np.zeros((NumEdges, NumPoints), dtype=np.dtype(np.float64));
    Ny = np.zeros((NumEdges, NumPoints), dtype=np.dtype(np.float64));

    DN = Calculate_DN( Area, Nodes );

    dL1_dX = DN[ 0 ][ 0 ]; dL1_dY = DN[ 1 ][ 0 ];
    dL2_dX = DN[ 0 ][ 1 ]; dL2_dY = DN[ 1 ][ 1 ];
    dL3_dX = DN[ 0 ][ 2 ]; dL3_dY = DN[ 1 ][ 2 ];

    eLength = Calculate_Edge_Lengths( Nodes );

    for p in range( NumPoints ):
        # Natural coordinates
        L2 = cX[ p ];
        L3 = cY[ p ];
        L1 = 1.0 - L2 - L3;

		# Edge: [0 1]
        Nx[ 0 ][ p ] = ( L1 * dL2_dX - L2 * dL1_dX ) * eLength[ 0 ];
        Ny[ 0 ][ p ] = ( L1 * dL2_dY - L2 * dL1_dY ) * eLength[ 0 ];
		    
		# Edge: [1 2]
        Nx[ 1 ][ p ] = ( L2 * dL3_dX - L3 * dL2_dX ) * eLength[ 1 ];
        Ny[ 1 ][ p ] = ( L2 * dL3_dY - L3 * dL2_dY ) * eLength[ 1 ];
		
		# Edge: [2 0]
        Nx[ 2 ][ p ] = ( L3 * dL1_dX - L1 * dL3_dX ) * eLength[ 2 ];
        Ny[ 2 ][ p ] = ( L3 * dL1_dY - L1 * dL3_dY ) * eLength[ 2 ];
        
    return Nx, Ny;   

###############################################################################

def Calculate_jK( PWave ):
    """
    - Returns jk = 1j*w*sqrt(ep*mu).
    """
    wf = 2.0 * np.pi * PWave.freq;
    w2 = wf * wf;

    eo = 8.8541878176e-12;
    mo = np.pi * 4.0e-7;            
    
    sg_real = PWave.sgm_r;
    sg_imag = PWave.sgm_i;

    ep_real = PWave.epr_r * eo; 
    ep_imag = PWave.epr_i * eo; 

    mu_real = PWave.mur_r * mo; 
    mu_imag = PWave.mur_i * mo; 

    cMu   = mu_real + 1j * mu_imag;
    cEpw2 = ep_real * w2 - sg_imag * wf + 1j * ( ep_imag * w2 + sg_real * wf );
    
    return ( 1j * np.sqrt( cEpw2 * cMu ) );    

###############################################################################

def Calculate_jK_div_Mu( PWave ):
    """
    - Returns jk/mu = 1j*w*sqrt(ep/mu).
    """
    wf = 2.0 * np.pi * PWave.freq;
    w2 = wf * wf;

    eo = 8.8541878176e-12;
    mo = np.pi * 4.0e-7;            
    
    sg_real = PWave.sgm_r;
    sg_imag = PWave.sgm_i;

    ep_real = PWave.epr_r * eo; 
    ep_imag = PWave.epr_i * eo; 

    mu_real = PWave.mur_r * mo; 
    mu_imag = PWave.mur_i * mo; 

    cMu   = mu_real + 1j * mu_imag;
    cEpw2 = ep_real * w2 - sg_imag * wf + 1j * ( ep_imag * w2 + sg_real * wf );
    
    return ( 1j * np.sqrt( cEpw2 / cMu ) );    

###############################################################################

def Calculate_Exterior_Normal( Nodes ):
    """
    - Returns the exterior normal of a surface element:  
      n = v2 x v1 / ||v2 x v1||
    """
    v1 = np.zeros( 3, dtype=np.dtype( np.float64 ) );
    v2 = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    v2[ 0 ] = Nodes[ 2 ].X - Nodes[ 0 ].X;
    v2[ 1 ] = Nodes[ 2 ].Y - Nodes[ 0 ].Y;
    v2[ 2 ] = Nodes[ 2 ].Z - Nodes[ 0 ].Z;
		    
    v1[ 0 ] = Nodes[ 1 ].X - Nodes[ 0 ].X;
    v1[ 1 ] = Nodes[ 1 ].Y - Nodes[ 0 ].Y;
    v1[ 2 ] = Nodes[ 1 ].Z - Nodes[ 0 ].Z;

    Next = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Next[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
    Next[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
    Next[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

    Next = Next / np.linalg.norm( Next );
    
    return Next;

###############################################################################

def Calculate_Area( Nodes ):
    """
    - Returns the area of a surface element: 
      Area = 0.5 * ||v2 x v1||
    """
    v1 = np.zeros( 3, dtype=np.dtype( np.float64 ) );
    v2 = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    v2[ 0 ] = Nodes[ 2 ].X - Nodes[ 0 ].X;
    v2[ 1 ] = Nodes[ 2 ].Y - Nodes[ 0 ].Y;
    v2[ 2 ] = Nodes[ 2 ].Z - Nodes[ 0 ].Z;
		    
    v1[ 0 ] = Nodes[ 1 ].X - Nodes[ 0 ].X;
    v1[ 1 ] = Nodes[ 1 ].Y - Nodes[ 0 ].Y;
    v1[ 2 ] = Nodes[ 1 ].Z - Nodes[ 0 ].Z;

    # Area = 0.5 * ||v2 x v1||
    Next = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Next[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
    Next[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
    Next[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

    Area = 0.5 * np.linalg.norm( Next );
    
    return Area;

###############################################################################

def Calculate_NiNj_Matrix( Nodes ):
    """
    - Returns Surf_Integral( NiNj ) matrix
    """
    NiNj_Matrix = np.zeros( (3, 3), dtype=np.dtype( np.float64 ) );
    
    Area = Calculate_Area( Nodes );
    
    for i in range( 3 ):
    	for j in range( 3 ):
            if( i != j ):
                NiNj_Matrix[ i ][ j ] = Area / 12.0;
            else:
                NiNj_Matrix[ i ][ j ] = Area / 6.0;		         

    return NiNj_Matrix;

###############################################################################

def Get_U_Vector_PW_RME_1st( PWave, Nodes ):
    """
    - Returns the residual vector for RME_1st element type.
    """
    # Nodes and Degrees of Freedom (DoF) of the element
    NumNodes = 3;
    NumDOFs  = NumNodes*3;
    
    # Initializing U vector
    U_vector = np.zeros( NumDOFs, dtype=np.dtype( np.complex128 ) );
    
    # Module and phase of E
    modul_F = PWave.modE; 
    phase_F = PWave.phaE;
    
    if( PWave.modE == 0.0 ): return U_vector;

    # E polarization vector ( normalized )
    Fpol = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Fpol[0] = PWave.ex; 
    Fpol[1] = PWave.ey; 
    Fpol[2] = PWave.ez;

    Fpol = Fpol / np.linalg.norm( Fpol );

    # Complex E vector
    cF = np.zeros( 3, dtype=np.dtype( np.complex128 ) );

    cF[0] = Fpol[0]; 
    cF[1] = Fpol[1]; 
    cF[2] = Fpol[2];

    cF = cF*( modul_F*np.cos( phase_F ) + 1j*modul_F*np.sin( phase_F ) );
    
    # Wave propagation vector ( normalized )
    Kdir = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Kdir[0] = PWave.kx;
    Kdir[1] = PWave.ky;
    Kdir[2] = PWave.kz;

    Kdir = Kdir / np.linalg.norm( Kdir );
    
    # Complex value jK
    jK = Calculate_jK( PWave );
    
    # Complex E vector on element nodes
    cFn = np.zeros( (NumNodes, NumNodes), dtype=np.dtype( np.complex128 ) );	
    
    for i in range( NumNodes ):
        KdirR  = Kdir[0]*Nodes[i].X + Kdir[1]*Nodes[i].Y + Kdir[2]*Nodes[i].Z;
        cFn[i] = cF;
        cFn[i] = cFn[i]*np.exp( jK * KdirR );
	                                                                         
	# Element exterior normal
    Next = Calculate_Exterior_Normal( Nodes );

	# k*n
    KdirNext = Kdir[0]*Next[0] + Kdir[1]*Next[1] + Kdir[2]*Next[2];

	# Complex value jK/Mu
    jK_div_Mu = Calculate_jK_div_Mu( PWave );
    
    # NiNj surface integrals
    NiNj_Matrix = Calculate_NiNj_Matrix( Nodes );
    
	# Curl element coefficient matrix ( n x n x E - n x k x E )
    CMatrix = np.zeros( (3, 3), dtype=np.dtype( np.float64 ) );

    CMatrix[0][0] = ( Next[0]*Next[0] - 1.0 ) - ( Kdir[0]*Next[0] - KdirNext );
    CMatrix[0][1] = ( Next[0]*Next[1]       ) - ( Kdir[0]*Next[1]            );
    CMatrix[0][2] = ( Next[0]*Next[2]       ) - ( Kdir[0]*Next[2]            );
    										  
    CMatrix[1][0] = ( Next[1]*Next[0]       ) - ( Kdir[1]*Next[0]            );
    CMatrix[1][1] = ( Next[1]*Next[1] - 1.0 ) - ( Kdir[1]*Next[1] - KdirNext );
    CMatrix[1][2] = ( Next[1]*Next[2]       ) - ( Kdir[1]*Next[2]            );
    										  
    CMatrix[2][0] = ( Next[2]*Next[0]       ) - ( Kdir[2]*Next[0]            );
    CMatrix[2][1] = ( Next[2]*Next[1]       ) - ( Kdir[2]*Next[1]            );
    CMatrix[2][2] = ( Next[2]*Next[2] - 1.0 ) - ( Kdir[2]*Next[2] - KdirNext );
    
    # Residual vector
    for i in range( NumNodes ):
        intg_NiCx = 0.0 + 1j * 0.0;
        intg_NiCy = 0.0 + 1j * 0.0;
        intg_NiCz = 0.0 + 1j * 0.0;

        for j in range( NumNodes ):
            intg_NiCx += ( CMatrix[0][0]*cFn[j][0] + 
                           CMatrix[0][1]*cFn[j][1] + 
                           CMatrix[0][2]*cFn[j][2] ) * NiNj_Matrix[i][j];
            
            intg_NiCy += ( CMatrix[1][0]*cFn[j][0] + 
                           CMatrix[1][1]*cFn[j][1] + 
                           CMatrix[1][2]*cFn[j][2] ) * NiNj_Matrix[i][j];
            
            intg_NiCz += ( CMatrix[2][0]*cFn[j][0] + 
                           CMatrix[2][1]*cFn[j][1] + 
                           CMatrix[2][2]*cFn[j][2] ) * NiNj_Matrix[i][j];

        U_vector[i             ] += jK_div_Mu * intg_NiCx;
        U_vector[i +   NumNodes] += jK_div_Mu * intg_NiCy;
        U_vector[i + 2*NumNodes] += jK_div_Mu * intg_NiCz;
        
    # Div element coefficient matrix ( n . E - k . E )
    DMatrix = np.zeros( (3, 3), dtype=np.dtype( np.float64 ) );

    DMatrix[0][0] = ( Next[0]*Next[0] ) - ( Next[0]*Kdir[0] );
    DMatrix[0][1] = ( Next[0]*Next[1] ) - ( Next[0]*Kdir[1] );
    DMatrix[0][2] = ( Next[0]*Next[2] ) - ( Next[0]*Kdir[2] );
								        
    DMatrix[1][0] = ( Next[1]*Next[0] ) - ( Next[1]*Kdir[0] );
    DMatrix[1][1] = ( Next[1]*Next[1] ) - ( Next[1]*Kdir[1] );
    DMatrix[1][2] = ( Next[1]*Next[2] ) - ( Next[1]*Kdir[2] );
									         
    DMatrix[2][0] = ( Next[2]*Next[0] ) - ( Next[2]*Kdir[0] );
    DMatrix[2][1] = ( Next[2]*Next[1] ) - ( Next[2]*Kdir[1] );
    DMatrix[2][2] = ( Next[2]*Next[2] ) - ( Next[2]*Kdir[2] );

    # Residual vector
    for i in range( NumNodes ):
        intg_NiDx = 0.0 + 1j * 0.0;
        intg_NiDy = 0.0 + 1j * 0.0;
        intg_NiDz = 0.0 + 1j * 0.0;

        for j in range( NumNodes ):
            intg_NiDx += ( DMatrix[0][0]*cFn[j][0] + 
                           DMatrix[0][1]*cFn[j][1] + 
                           DMatrix[0][2]*cFn[j][2] ) * NiNj_Matrix[i][j];
            
            intg_NiDy += ( DMatrix[1][0]*cFn[j][0] + 
                           DMatrix[1][1]*cFn[j][1] + 
                           DMatrix[1][2]*cFn[j][2] ) * NiNj_Matrix[i][j];
            
            intg_NiDz += ( DMatrix[2][0]*cFn[j][0] + 
                           DMatrix[2][1]*cFn[j][1] + 
                           DMatrix[2][2]*cFn[j][2] ) * NiNj_Matrix[i][j];

        U_vector[i             ] -= jK_div_Mu * intg_NiDx;
        U_vector[i +   NumNodes] -= jK_div_Mu * intg_NiDy;
        U_vector[i + 2*NumNodes] -= jK_div_Mu * intg_NiDz;

    return U_vector;

###############################################################################

def Get_P_Matrix_PW_RME_1st( PWave, Nodes ):
    """
    - Returns the stiffness matrix for RME_1st element type.
    """
    # Nodes and Degrees of Freedom (DoF) of the element
    NumNodes = 3;
    NumDOFs  = NumNodes*3;
    
	# Setting P matrix
    P_matrix = np.zeros((NumDOFs, NumDOFs), dtype=np.dtype(np.complex128));
    
    # Robin coefficient
    cRCoeff = -Calculate_jK_div_Mu( PWave );

    # NiNj surface integrals matrix
    NiNj_Matrix = Calculate_NiNj_Matrix( Nodes );

    # Stiffness matrix
    for i in range( NumNodes ):
        for j in range( NumNodes ):
            P_matrix[i             ][j             ] = cRCoeff*NiNj_Matrix[i][j];
            P_matrix[i +   NumNodes][j +   NumNodes] = cRCoeff*NiNj_Matrix[i][j];
            P_matrix[i + 2*NumNodes][j + 2*NumNodes] = cRCoeff*NiNj_Matrix[i][j];
            
    return P_matrix;

###############################################################################

def Get_U_Vector_PW_EDG_1st( PWave, Nodes ):
    """
    - Returns the residual vector for EDG_1st element type.
    """
    NumEdges = 3;
    NumVerts = 3;
    NumNodes = 6;
    NumDOFs  = 18;
    
    # Initializing U vector
    U_vector = np.zeros( NumDOFs, dtype=np.dtype( np.complex128 ) );
    
    # Module and phase of E
    modul_F = PWave.modE; 
    phase_F = PWave.phaE;
    
    if( PWave.modE == 0.0 ): return U_vector;

    # E polarization vector ( normalized )
    Fpol = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Fpol[0] = PWave.ex; 
    Fpol[1] = PWave.ey; 
    Fpol[2] = PWave.ez;

    Fpol = Fpol / np.linalg.norm( Fpol );

    # Complex E vector
    cF = np.zeros( 3, dtype=np.dtype( np.complex128 ) );

    cF[0] = Fpol[0]; 
    cF[1] = Fpol[1]; 
    cF[2] = Fpol[2];

    cF = cF*( modul_F*np.cos( phase_F ) + 1j*modul_F*np.sin( phase_F ) );
    
    # Wave propagation vector ( normalized )
    Kdir = np.zeros( 3, dtype=np.dtype( np.float64 ) );

    Kdir[0] = PWave.kx;
    Kdir[1] = PWave.ky;
    Kdir[2] = PWave.kz;

    Kdir = Kdir / np.linalg.norm( Kdir );
    
	# Element exterior normal
    Next = Calculate_Exterior_Normal( Nodes );

	# k*n
    KdirNext = Kdir[0]*Next[0] + Kdir[1]*Next[1] + Kdir[2]*Next[2];    
    
	# Curl element coefficient matrix ( n x n x E - n x k x E )
    CMatrix = np.zeros( (3, 3), dtype=np.dtype( np.float64 ) );

    CMatrix[0][0] = ( Next[0]*Next[0] - 1.0 ) - ( Kdir[0]*Next[0] - KdirNext );
    CMatrix[0][1] = ( Next[0]*Next[1]       ) - ( Kdir[0]*Next[1]            );
    CMatrix[0][2] = ( Next[0]*Next[2]       ) - ( Kdir[0]*Next[2]            );
    										  
    CMatrix[1][0] = ( Next[1]*Next[0]       ) - ( Kdir[1]*Next[0]            );
    CMatrix[1][1] = ( Next[1]*Next[1] - 1.0 ) - ( Kdir[1]*Next[1] - KdirNext );
    CMatrix[1][2] = ( Next[1]*Next[2]       ) - ( Kdir[1]*Next[2]            );
    										  
    CMatrix[2][0] = ( Next[2]*Next[0]       ) - ( Kdir[2]*Next[0]            );
    CMatrix[2][1] = ( Next[2]*Next[1]       ) - ( Kdir[2]*Next[1]            );
    CMatrix[2][2] = ( Next[2]*Next[2] - 1.0 ) - ( Kdir[2]*Next[2] - KdirNext );    
    
    # Complex value jK
    jK = Calculate_jK( PWave );
    
    # Gauss points
    cX, cY, W = GaussPoints2D_Order12();
    
    # Jacobian
    Area = Calculate_Area( Nodes );

    detJ = 2.0 * Area;
    
    for gp in range( len(W) ): W[ gp ] *= detJ;     
    
    # XY axis of the co-planar coordinate system 
    X_axis, Y_axis = Get_CoPlanar_XY_Axis( Nodes );    
    
    # Input field on Gauss points 
    cFx = np.zeros( len(W), dtype=np.dtype( np.complex128 ) );
    cFy = np.zeros( len(W), dtype=np.dtype( np.complex128 ) );
    
    for gp in range( len(W) ):
        
        alpha = cX[gp];
        beta  = cY[gp];
        
        rX = ( alpha * ( Nodes[1].X - Nodes[0].X ) + 
               beta  * ( Nodes[2].X - Nodes[0].X ) + Nodes[0].X );
        
        rY = ( alpha * ( Nodes[1].Y - Nodes[0].Y ) + 
               beta  * ( Nodes[2].Y - Nodes[0].Y ) + Nodes[0].Y );
        
        rZ = ( alpha * ( Nodes[1].Z - Nodes[0].Z ) + 
               beta  * ( Nodes[2].Z - Nodes[0].Z ) + Nodes[0].Z );
        
        KdirR = Kdir[0] * rX + Kdir[1] * rY + Kdir[2] * rZ;
        
        gpcF = cF * np.exp( jK * KdirR );

        cFx3D = ( CMatrix[0][0] * gpcF[0] + 
                  CMatrix[0][1] * gpcF[1] + 
                  CMatrix[0][2] * gpcF[2] ); 
        
        cFy3D = ( CMatrix[1][0] * gpcF[0] + 
                  CMatrix[1][1] * gpcF[1] + 
                  CMatrix[1][2] * gpcF[2] ); 
        
        cFz3D = ( CMatrix[2][0] * gpcF[0] + 
                  CMatrix[2][1] * gpcF[1] + 
                  CMatrix[2][2] * gpcF[2] ); 

        cFx[gp] = cFx3D * X_axis[0] + cFy3D * X_axis[1] + cFz3D * X_axis[2]; 
        cFy[gp] = cFx3D * Y_axis[0] + cFy3D * Y_axis[1] + cFz3D * Y_axis[2];  
        
    # 1st order Whitney base functions
    Nx, Ny = Whitney2D_Ni_1st( cX, cY, Area, Nodes );
        
    # Complex value jK/Mu
    jK_div_Mu = Calculate_jK_div_Mu( PWave );
    
    # Residual vector
    for i in range( NumEdges ):
        
        sIntg_FNi = 0.0 + 1j * 0.0;
        
        for gp in range( len(W) ):
            sIntg_FNi += ( W[ gp ] * ( Nx[ i ][ gp ] * cFx[ gp ] + 
                                       Ny[ i ][ gp ] * cFy[ gp ] ) );
            
        U_vector[ i + NumNodes + NumVerts ] += jK_div_Mu * sIntg_FNi; 

    return U_vector;

###############################################################################

def Get_P_Matrix_PW_EDG_1st( PWave, Nodes ):
    """
    - Returns the stiffness matrix for EDG_1st element type.
    """
    NumEdges = 3;
    NumVerts = 3;
    NumNodes = 6;
    NumDOFs  = 18;

	# Setting P matrix
    P_matrix = np.zeros((NumDOFs, NumDOFs), dtype=np.dtype(np.complex128));
    
    # Robin coefficient
    cRCoeff = -Calculate_jK_div_Mu( PWave );
    
    # Gauss points
    cX, cY, W = GaussPoints2D_Order12();
    
    # Jacobian
    Area = Calculate_Area( Nodes );

    detJ = 2.0 * Area;
    
    for gp in range( len(W) ): W[ gp ] *= detJ; 

    # Edge matrix
    Nx, Ny = Whitney2D_Ni_1st( cX, cY, Area, Nodes );
    
    EDG_matrix = np.zeros((NumEdges, NumEdges), dtype=np.dtype(np.complex128));
    
    # Upper diagonal
    for i in range( NumEdges ):
        for j in range( i, NumEdges ):
            
            sIntg_NiNj = 0.0;
            
            for gp in range( len(W) ):
                sIntg_NiNj += W[ gp ] * ( Nx[ i ][ gp ] * Nx[ j ][ gp ] + 
                                          Ny[ i ][ gp ] * Ny[ j ][ gp ] );
            
            EDG_matrix[ i ][ j ] = cRCoeff * sIntg_NiNj;

    # Lower diagonal
    for i in range( NumEdges ):
        for j in range( i+1, NumEdges ):
            
            EDG_matrix[ j ][ i ] = EDG_matrix[ i ][ j ];

    # Edge element matrix inclusion into stiffness matrix  
    edg2dof = NumNodes + NumVerts;
    
    for i in range( NumEdges ):
        for j in range( NumEdges ):
            
            P_matrix[ i + edg2dof ][ j + edg2dof ] = EDG_matrix[ i ][ j ];
            
    return P_matrix;

###############################################################################
# - Defines a plane wave object and generates the files containing the 
#   Imported Robin Boundary Condition (IRBC) matrix and vector.
###############################################################################
if __name__ == "__main__":
    
    PWave = PlaneWave();

    # Frequency
    PWave.freq = 12e8;

    # Module and phase of E
    PWave.modE = 1.0; PWave.phaE = 0.0;

    # Polarization vector
    PWave.ex = 1.0; PWave.ey = 1.0; PWave.ez = 0.0;

    # Wave vector
    PWave.kx = 1.0; PWave.ky = -1.0; PWave.kz = 0.0;

    # Material properties ( jk = j*w*sqrt(ep*mu) )
    PWave.sgm_r = 0.0; PWave.sgm_i = 0.0;
    PWave.epr_r = 1.0; PWave.epr_i = 0.0;
    PWave.mur_r = 1.0; PWave.mur_i = 0.0;

    # Elemet type (RME_1st | EDG_1st)
    ElementType = 'EDG_1st'
    
    # GiD problem folder
    ProblemName = r'../test_problem.gid'

    # IRBC matrix and vector file names
    IRBC_U_FileName = ProblemName + r'/Vector_U_IRBC.bin'
    IRBC_P_FileName = ProblemName + r'/Matrix_P_IRBC.bin'

    # Read GiD mesh nodes:  No[ID] = p(x,y,z);
    nIds, nXs, nYs, nZs = Read_Nodes_File( ProblemName );

    # Read IRBC surface elements
    eList = Read_IRBC_Elemets_File( ProblemName );

    # Write in 'Vector_U_IRBC.bin' file
    Ufile = open( IRBC_U_FileName, 'wb' );

    # Write in 'Vector_P_IRBC.bin' file
    Pfile = open( IRBC_P_FileName, 'wb' );

    # Loop over IRBC surface elements
    for i in range(len( eList )):
        
        # Get element nodes list
        eNodes = Get_Element_Nodes_List( eList[i], nXs, nYs, nZs );
        
        # Get IRBC matrix and vector
        if  (ElementType == 'RME_1st'):
            cU = Get_U_Vector_PW_RME_1st( PWave, eNodes );
            cP = Get_P_Matrix_PW_RME_1st( PWave, eNodes );
        elif(ElementType == 'EDG_1st'):
            cU = Get_U_Vector_PW_EDG_1st( PWave, eNodes );
            cP = Get_P_Matrix_PW_EDG_1st( PWave, eNodes );
         
        # Write in binary files 
        np.array( cU, dtype='<c16' ).tofile( Ufile );
        np.array( cP, dtype='<c16' ).tofile( Pfile );
             
    Ufile.close();
    Pfile.close();
    
###############################################################################
