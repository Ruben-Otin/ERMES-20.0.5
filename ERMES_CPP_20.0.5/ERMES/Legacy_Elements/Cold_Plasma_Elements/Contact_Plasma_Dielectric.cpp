//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// - Trans(T)*StiffMatrix*T to eliminate one of the nodes in a contact pair for Dielectric-Plasma interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Modeler::ContactHO_TKT_Plasma( std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
{
	if ( mContactPairs.size() == 0 ) return;
	
	int local_i = 0;

    std::vector<int>::iterator hoit;

    for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
    {
        if ( mContactPairs.find(*hoit) != mContactPairs.end() )
        {
            int RNodeId = (*hoit);
            int LNodeId = mContactPairs[RNodeId].ctcNode;
            
            ///////////////////////////////////////////////////   Modify Id vector    ////////////////////////////////////////////////
            int step    = HONodesId.size();
			int step2   = 2*step;
            int dofSize = 3*step;

            Node::Pointer pNode = mpModel->GetNode( LNodeId - 1 );
 
            // Carefull!!!, this changes the content of (*hoit)
            HONodesId[local_i        ] = pNode->Id();

			IdVector [local_i        ] = pNode->pDofcEx()->EquationId();
			IdVector [local_i + step ] = pNode->pDofcEy()->EquationId();
			IdVector [local_i + step2] = pNode->pDofcEz()->EquationId();

            //////////////////////////////////////////////    Inv(T)*StiffnesMatrix*T    /////////////////////////////////////////////
            std::complex<double> dof_x, dof_y, dof_z;

            Matrix< std::complex<double> > ContactMatrix(3,3);

            Get_Contact_Matrix_Plasma( ContactMatrix, RNodeId );

            // StiffnesMatrix*T
            for (int i=0; i<dofSize; i++)
            { 
                dof_x = StiffMatrix[i][local_i        ];
                dof_y = StiffMatrix[i][local_i + step ];
                dof_z = StiffMatrix[i][local_i + step2];

                StiffMatrix[i][local_i        ] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                StiffMatrix[i][local_i + step ] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                StiffMatrix[i][local_i + step2] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
            }

            // Inv(T)*StiffnesMatrix    
            for (int i=0; i<dofSize; i++)
            { 
                dof_x = StiffMatrix[local_i        ][i];
                dof_y = StiffMatrix[local_i + step ][i];
                dof_z = StiffMatrix[local_i + step2][i];

				StiffMatrix[local_i        ][i] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                StiffMatrix[local_i + step ][i] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                StiffMatrix[local_i + step2][i] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
            }

        }//end if( mContactPairs.find(*hoit) != mContactPairs.end() )

        local_i++; 

    }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) 	
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// - Calculates the contact matrix for Dielectric-Plasma interfaces
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Modeler::Get_Contact_Matrix_Plasma( Matrix< std::complex<double> >& ContactMatrix, int RNodesId )
{
    ContactPairData& pContactPair = mContactPairs[RNodesId];

    Node::Pointer prNode = mpModel->GetNode( RNodesId             - 1 );
    Node::Pointer plNode = mpModel->GetNode( pContactPair.ctcNode - 1 );

    std::vector<double> n = pContactPair.ctcNormal;
    std::vector<double> t(3), b(3);
   
    TangencialCoordinates(n,t,b);
      
    // Material properties constants
    Properties::Pointer Properties;
    
    double eo   = 8.8541878176e-12;
    double freq = mProblemFrequency;
    
    std::complex<double> cZero( 0.0, 0.0 ); 
    std::complex<double> cUnit( 0.0, 1.0 );
    
    double sigma; 
    double epr_real, epr_imag;
    double mur_real, mur_imag;
    
    Matrix< std::complex<double> > TEn(3,3);

    // Material R
    Properties = mpModel->GetProperties( pContactPair.matR );
    
    sigma      = (*Properties)( IHL_ELECTRIC_CONDUCTIVITY  );
    epr_real   = (*Properties)( REAL_ELECTRIC_PERMITTIVITY );
    epr_imag   = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY );
    mur_real   = (*Properties)( REAL_MAGNETIC_PERMEABILITY );
    mur_imag   = (*Properties)( IMAG_MAGNETIC_PERMEABILITY );
    
    // Plasma permittivity tensor
    if ( (sigma == 0.0) && (epr_real == 1.0) && (epr_imag == 0.0) && (mur_real == 1.0) && (mur_imag == 0.0) )
    {
        mpColdPlasma->Get_PermittivityTensor_InNode( TEn, prNode, freq ); 
    }
    // IHL permittivity tensor 
    else
    {
        std::complex<double> cEps( eo*epr_real , eo*epr_imag + ( sigma / freq ) );      

        TEn[0][0] = cEps ; TEn[0][1] = cZero; TEn[0][2] = cZero; 
        TEn[1][0] = cZero; TEn[1][1] = cEps ; TEn[1][2] = cZero; 
        TEn[2][0] = cZero; TEn[2][1] = cZero; TEn[2][2] = cEps ;
    }
    // [TR] * {Er} = [TL] * {El}
    Matrix< std::complex<double> > TR(3,3);
    
    TR[0][0] = n[0]*TEn[0][0] + n[1]*TEn[1][0] + n[2]*TEn[2][0]; 
    TR[0][1] = n[0]*TEn[0][1] + n[1]*TEn[1][1] + n[2]*TEn[2][1]; 
    TR[0][2] = n[0]*TEn[0][2] + n[1]*TEn[1][2] + n[2]*TEn[2][2]; 
    
    TR[1][0] = t[0]; TR[1][1] = t[1]; TR[1][2] = t[2]; 
    TR[2][0] = b[0]; TR[2][1] = b[1]; TR[2][2] = b[2];
    
    // Material L
    Properties = mpModel->GetProperties( pContactPair.matL );
    
    sigma      = (*Properties)( IHL_ELECTRIC_CONDUCTIVITY  );
    epr_real   = (*Properties)( REAL_ELECTRIC_PERMITTIVITY );
    epr_imag   = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY );
    mur_real   = (*Properties)( REAL_MAGNETIC_PERMEABILITY );
    mur_imag   = (*Properties)( IMAG_MAGNETIC_PERMEABILITY );
    
    // Plasma permittivity tensor
    if ( (sigma == 0.0) && (epr_real == 1.0) && (epr_imag == 0.0) && (mur_real == 1.0) && (mur_imag == 0.0) )
    {
        mpColdPlasma->Get_PermittivityTensor_InNode( TEn, plNode, freq ); 
    }
    // IHL permittivity tensor 
    else
    {
        std::complex<double> cEps( eo*epr_real , eo*epr_imag + ( sigma / freq ) );  

        TEn[0][0] = cEps ; TEn[0][1] = cZero; TEn[0][2] = cZero; 
        TEn[1][0] = cZero; TEn[1][1] = cEps ; TEn[1][2] = cZero; 
        TEn[2][0] = cZero; TEn[2][1] = cZero; TEn[2][2] = cEps ;
    }
    // [TR] * {Er} = [TL] * {El}
    Matrix< std::complex<double> > TL(3,3);
    
    TL[0][0] = n[0]*TEn[0][0] + n[1]*TEn[1][0] + n[2]*TEn[2][0]; 
    TL[0][1] = n[0]*TEn[0][1] + n[1]*TEn[1][1] + n[2]*TEn[2][1]; 
    TL[0][2] = n[0]*TEn[0][2] + n[1]*TEn[1][2] + n[2]*TEn[2][2]; 
    
    TL[1][0] = t[0]; TL[1][1] = t[1]; TL[1][2] = t[2]; 
    TL[2][0] = b[0]; TL[2][1] = b[1]; TL[2][2] = b[2];
    
    // {Er} = ( [invTR] * [TL] ) * {El}
    Matrix< std::complex<double> > invTR ( 3, 3 );
    
    std::complex<double> detTR = Determinant( TR, 3 );
    
    invTR[0][0] = ( TR[1][1] * TR[2][2] - TR[1][2] * TR[2][1] ) / detTR; 
    invTR[0][1] = ( TR[0][2] * TR[2][1] - TR[0][1] * TR[2][2] ) / detTR; 
    invTR[0][2] = ( TR[0][1] * TR[1][2] - TR[0][2] * TR[1][1] ) / detTR; 
    
    invTR[1][0] = ( TR[1][2] * TR[2][0] - TR[1][0] * TR[2][2] ) / detTR; 
    invTR[1][1] = ( TR[0][0] * TR[2][2] - TR[0][2] * TR[2][0] ) / detTR; 
    invTR[1][2] = ( TR[0][2] * TR[1][0] - TR[0][0] * TR[1][2] ) / detTR; 
    
    invTR[2][0] = ( TR[1][0] * TR[2][1] - TR[1][1] * TR[2][0] ) / detTR;  
    invTR[2][1] = ( TR[0][1] * TR[2][0] - TR[0][0] * TR[2][1] ) / detTR;  
    invTR[2][2] = ( TR[0][0] * TR[1][1] - TR[0][1] * TR[1][0] ) / detTR;  
    
    // Contact matrix = [invTR] * [TL] 
    ContactMatrix.Resize(3,3);
    ContactMatrix = invTR * TL;    
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// - Inv(T)*b in contact nodes for Dielectric-Plasma interfaces
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Modeler::ContactHO_TR_Plasma( std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
{
    if ( mContactPairs.size() == 0 ) return;
	
	int local_i = 0;

    std::vector<int>::iterator hoit;

    for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
    {
        if ( mContactPairs.find(*hoit) != mContactPairs.end() )
        {
            int RNodeId = (*hoit);
            int LNodeId = mContactPairs[RNodeId].ctcNode;
 
            ///////////////////////////////////////////////   Modify Id vector    ///////////////////////////////////////////
            int step  = HONodesId.size();
			int step2 = 2*step;

            Node::Pointer pNode = mpModel->GetNode( LNodeId - 1 );

            // Carefull!!!, this changes the content of (*hoit)
			HONodesId[local_i        ] = pNode->Id();

			IdVector [local_i        ] = pNode->pDofcEx()->EquationId();
			IdVector [local_i + step ] = pNode->pDofcEy()->EquationId();
			IdVector [local_i + step2] = pNode->pDofcEz()->EquationId();

            ///////////////////////////////////////////////////  Inv(T)*b  //////////////////////////////////////////////////
            std::complex<double> dof_x, dof_y, dof_z;

            Matrix< std::complex<double> > ContactMatrix;

            Get_Contact_Matrix_Plasma( ContactMatrix, RNodeId );
           
			dof_x = ResVector[local_i        ];
            dof_y = ResVector[local_i + step ];
            dof_z = ResVector[local_i + step2];

            ResVector[local_i        ] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
            ResVector[local_i + step ] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
            ResVector[local_i + step2] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
        }

        local_i++; 
    }  
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// - Returns the determinant of a square matrix M of size n (for complex matrices)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::complex<double> Modeler::Determinant( Matrix< std::complex<double> >& M, int n ) 
{   
    std::complex<double> det( 0.0, 0.0 );

    Matrix< std::complex<double> > subM(n,n);

    if (n == 1)
    {
        return M[0][0];
    }

    if (n == 2)
    {
        return ( (M[0][0] * M[1][1]) - (M[1][0] * M[0][1]) );
    }

    for ( int x = 0; x < n; x++ ) 
    {
        int subi = 0; 
        for ( int i = 1; i < n; i++ ) 
        {
            int subj = 0;
            for ( int j = 0; j < n; j++ ) 
            {
                if (j == x) continue;

                subM[subi][subj] = M[i][j];
                subj++;
            }
            subi++;
        }
        det = det + ( pow( -1, x ) * M[0][x] * Determinant( subM, n-1 ) );
    }
   
    return det;
}