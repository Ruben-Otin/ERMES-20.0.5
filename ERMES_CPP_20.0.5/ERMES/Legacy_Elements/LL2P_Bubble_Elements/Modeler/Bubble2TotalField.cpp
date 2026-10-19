    //**************************************************************************************************************************************************************
    //* - Set PBC with debugging info
    //**************************************************************************************************************************************************************
	void Modeler::SetPBC()
	{
        // Setting PBC for Front-Back surfaces
		if ( mPBCFrontElements.size() > 0 ) Set_FrontBack_PBC();

        //// Check if Front-Back PBC is compatible with LL2P 3sb elements
        //if ( ( mElementOrder == 0 ) && ( ( mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) > 0 ) )
        //{
        //    Send_Error_Msg( "PBC for LL2P 3sb elements is only valid when Front-Back surface meshes are equal", "Please, re-mesh or change element type.", 1 );            
        //}
	
        // Setting PBC for Right-Left surfaces
        if ( mPBCRightElements.size() > 0 ) Set_RightLeft_PBC();

        //// Check if Right-Left PBC is compatible with LL2P 3sb elements
        //if ( ( mElementOrder == 0 ) && ( ( mPBC_NodeEdgePairs_Right.size() + mPBC_NodeElementPairs_Right.size() ) > 0 ) )
        //{
        //    Send_Error_Msg( "PBC for LL2P 3sb elements is only valid when Right-Left surface meshes are equal", "Please, re-mesh or change element type.", 1 );                         
        //}

        // Front-Back info
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator it;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            if ( mPBC_NodeNodePairs_Front.find((*it)->Id()) != mPBC_NodeNodePairs_Front.end() )
            {
                int node1 = mPBC_NodeNodePairs_Front[(*it)->Id()];
                std::printf("N: %4i - ( %+f, %+f, %+f ) - B: %4i - ( %+f, %+f, %+f )\n", 
                (*it)->Id(),          (*it)->X(),          (*it)->Y(),          (*it)->Z(), 
                node1      , nodes[node1-1]->X(), nodes[node1-1]->Y(), nodes[node1-1]->Z() );
            }
			else if ( mPBC_NodeEdgePairs_Front.find((*it)->Id()) != mPBC_NodeEdgePairs_Front.end() )
			{
                std::printf("E: %4i - ( %+f, %+f, %+f ) - B: %4i, %4i - %+f\n", (*it)->Id(), (*it)->X(), (*it)->Y(), (*it)->Z(), 
                mPBC_NodeEdgePairs_Front[(*it)->Id()][0], mPBC_NodeEdgePairs_Front[(*it)->Id()][1], 
                mPBC_NodeEdgeCoord_Front[(*it)->Id()]  );
			}
			else if ( mPBC_NodeElementPairs_Front.find((*it)->Id()) != mPBC_NodeElementPairs_Front.end() )
			{
                std::printf("F: %4i - ( %+f, %+f, %+f ) - B: %4i, %4i, %4i - %+f, %+f\n", (*it)->Id(), (*it)->X(), (*it)->Y(), (*it)->Z(),
                mPBC_NodeElementPairs_Front[(*it)->Id()][0], mPBC_NodeElementPairs_Front[(*it)->Id()][1], mPBC_NodeElementPairs_Front[(*it)->Id()][2],
                mPBC_NodeElementCoord_Front[(*it)->Id()][0], mPBC_NodeElementCoord_Front[(*it)->Id()][1] );
			}
            else 
            {
                std::printf("*: %4i - ( %+f, %+f, %+f ) - Not assigned\n", (*it)->Id(), (*it)->X(), (*it)->Y(), (*it)->Z() );
            }
        }

        // Right-Left info
        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            if ( mPBC_NodeNodePairs_Right.find((*it)->Id()) != mPBC_NodeNodePairs_Right.end() )
            {
                int node1 = mPBC_NodeNodePairs_Right[(*it)->Id()];
                std::printf("N: %4i - ( %+f, %+f, %+f ) - L: %4i - ( %+f, %+f, %+f )\n", 
                (*it)->Id(),          (*it)->X(),          (*it)->Y(),          (*it)->Z(), 
                node1      , nodes[node1-1]->X(), nodes[node1-1]->Y(), nodes[node1-1]->Z() );
               
            }
			else if ( mPBC_NodeEdgePairs_Right.find((*it)->Id()) != mPBC_NodeEdgePairs_Right.end() )
			{
                std::printf("E: %4i - ( %+f, %+f, %+f ) - L: %4i, %4i - %+f\n", (*it)->Id(), (*it)->X(), (*it)->Y(), (*it)->Z(), 
                mPBC_NodeEdgePairs_Right[(*it)->Id()][0], mPBC_NodeEdgePairs_Right[(*it)->Id()][1],
                mPBC_NodeEdgeCoord_Right[(*it)->Id()]  );				
			}
			else if ( mPBC_NodeElementPairs_Right.find((*it)->Id()) != mPBC_NodeElementPairs_Right.end() )
			{
                std::printf("F: %4i - ( %+f, %+f, %+f ) - L: %4i, %4i, %4i - %+f, %+f\n", (*it)->Id(), (*it)->X(), (*it)->Y(), (*it)->Z(),
                mPBC_NodeElementPairs_Right[(*it)->Id()][0], mPBC_NodeElementPairs_Right[(*it)->Id()][1], mPBC_NodeElementPairs_Right[(*it)->Id()][2],
				mPBC_NodeElementCoord_Right[(*it)->Id()][0], mPBC_NodeElementCoord_Right[(*it)->Id()][1] );
			}
            else 
            {
                std::printf("*: %4i - ( %+f, %+f, %+f ) - Not assigned\n", (*it)->Id(), (*it)->X(), (*it)->Y(), (*it)->Z() );
            }
        }
	}
	
	//**************************************************************************************************************************************************************
    //* - Trans(T)*StiffMatrix*T for PBC nodes - Transforming E_bubble -> E_total
    //**************************************************************************************************************************************************************
    void Modeler::PBC_Transf2TotalE( std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix )
    {
        if ( ( mElementOrder != 0 ) || ( ( mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) == 0 ) )
        {
            return;
        }

        int mNumNodes = 17;
        int mNumDofs  = mNumNodes * 3;

		// P1 basis functions on nodal points
		std::vector<double> nX( mNumNodes ), nY( mNumNodes ), nZ( mNumNodes );

		// Vertices 0,1,2,3
		nX[0] = 0.00; nY[0] = 0.00; nZ[0] = 0.00;
		nX[1] = 1.00; nY[1] = 0.00; nZ[1] = 0.00;
		nX[2] = 0.00; nY[2] = 1.00; nZ[2] = 0.00;
		nX[3] = 0.00; nY[3] = 0.00; nZ[3] = 1.00;
		
		// Face [0 1 2]
		nX[4] = 0.25; nY[4] = 0.25; nZ[4] = 0.00;
		nX[5] = 0.50; nY[5] = 0.25; nZ[5] = 0.00;
		nX[6] = 0.25; nY[6] = 0.50; nZ[6] = 0.00;

		// Face [0 1 3]
		nX[7] = 0.25; nY[7] = 0.00; nZ[7] = 0.25;
		nX[8] = 0.50; nY[8] = 0.00; nZ[8] = 0.25;
		nX[9] = 0.25; nY[9] = 0.00; nZ[9] = 0.50;

		// Face [1 2 3]
		nX[10] = 0.50; nY[10] = 0.25; nZ[10] = 0.25;
		nX[11] = 0.25; nY[11] = 0.50; nZ[11] = 0.25;
		nX[12] = 0.25; nY[12] = 0.25; nZ[12] = 0.50;

		// Face [0 2 3]
		nX[13] = 0.00; nY[13] = 0.25; nZ[13] = 0.25;
		nX[14] = 0.00; nY[14] = 0.50; nZ[14] = 0.25;
		nX[15] = 0.00; nY[15] = 0.25; nZ[15] = 0.50;

		// Volume [0 1 2 3]
		nX[16] = 0.25; nY[16] = 0.25; nZ[16] = 0.25; 

		// P1 basis on nodes
		Matrix<double> Nn;

		Lagrange3D_Ni_1st( Nn, nX, nY, nZ );

		Matrix< std::complex<double> > NE( mNumDofs, mNumDofs );

		for( int i=0; i<mNumDofs; i++ ) 
		{
			for( int j=0; j<mNumDofs; j++ ) 
			{
			    if ( i == j ) NE[i][i] = std::complex<double>( 1.00, 0.00 );
				else          NE[i][j] = std::complex<double>( 0.00, 0.00 );
			}
		}

		for( int i=4; i<mNumNodes; i++ ) 
		{
			if ( ( mPBC_NodeEdgePairs_Front.find   ( HONodesId[i] ) != mPBC_NodeEdgePairs_Front.end()    ) || 
                 ( mPBC_NodeElementPairs_Front.find( HONodesId[i] ) != mPBC_NodeElementPairs_Front.end() )  )
            {
                for( int j=0; j<4; j++ ) 
			    {
			        NE[ i               ][ j               ] = std::complex<double>( -Nn[j][i], 0.00 );
				    NE[ i + mNumNodes   ][ j + mNumNodes   ] = std::complex<double>( -Nn[j][i], 0.00 );
				    NE[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = std::complex<double>( -Nn[j][i], 0.00 );
		        }
            }
		}

		Matrix< std::complex<double> > Stiff_x_NE = StiffMatrix * NE;
		Matrix< std::complex<double> > NE_tran  = NE.Transpose();
		StiffMatrix = ( NE_tran * Stiff_x_NE );
    }
 
    //**************************************************************************************************************************************************************
    //* - Trans(T)*b for PBC nodes - Transforming E_bubble -> E_total
    //**************************************************************************************************************************************************************
	void Modeler::PBC_Transf2TotalR( std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector )
	{
        if ( ( mElementOrder != 0 ) || ( ( mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) == 0 ) )
        {
            return;
        }

        int mNumNodes = 17;
        int mNumDofs  = mNumNodes * 3;

		// P1 basis functions on nodal points
		std::vector<double> nX( mNumNodes ), nY( mNumNodes ), nZ( mNumNodes );

		// Vertices 0,1,2,3
		nX[0] = 0.00; nY[0] = 0.00; nZ[0] = 0.00;
		nX[1] = 1.00; nY[1] = 0.00; nZ[1] = 0.00;
		nX[2] = 0.00; nY[2] = 1.00; nZ[2] = 0.00;
		nX[3] = 0.00; nY[3] = 0.00; nZ[3] = 1.00;
		
		// Face [0 1 2]
		nX[4] = 0.25; nY[4] = 0.25; nZ[4] = 0.00;
		nX[5] = 0.50; nY[5] = 0.25; nZ[5] = 0.00;
		nX[6] = 0.25; nY[6] = 0.50; nZ[6] = 0.00;

		// Face [0 1 3]
		nX[7] = 0.25; nY[7] = 0.00; nZ[7] = 0.25;
		nX[8] = 0.50; nY[8] = 0.00; nZ[8] = 0.25;
		nX[9] = 0.25; nY[9] = 0.00; nZ[9] = 0.50;

		// Face [1 2 3]
		nX[10] = 0.50; nY[10] = 0.25; nZ[10] = 0.25;
		nX[11] = 0.25; nY[11] = 0.50; nZ[11] = 0.25;
		nX[12] = 0.25; nY[12] = 0.25; nZ[12] = 0.50;

		// Face [0 2 3]
		nX[13] = 0.00; nY[13] = 0.25; nZ[13] = 0.25;
		nX[14] = 0.00; nY[14] = 0.50; nZ[14] = 0.25;
		nX[15] = 0.00; nY[15] = 0.25; nZ[15] = 0.50;

		// Volume [0 1 2 3]
		nX[16] = 0.25; nY[16] = 0.25; nZ[16] = 0.25; 

		// P1 basis on nodes
		Matrix<double> Nnn;

		Lagrange3D_Ni_1st( Nnn, nX, nY, nZ );
		Matrix< std::complex<double> > NE( mNumDofs, mNumDofs );

		for( int i=0; i<mNumDofs; i++ ) 
		{
			for( int j=0; j<mNumDofs; j++ ) 
			{
			    if ( i == j ) NE[i][i] = std::complex<double>( 1.00, 0.00 );
				else          NE[i][j] = std::complex<double>( 0.00, 0.00 );
			}
		}

		for( int i=4; i<mNumNodes; i++ ) 
		{
			if ( ( mPBC_NodeEdgePairs_Front.find   ( HONodesId[i] ) != mPBC_NodeEdgePairs_Front.end()    ) || 
                 ( mPBC_NodeElementPairs_Front.find( HONodesId[i] ) != mPBC_NodeElementPairs_Front.end() )  )
            {
                for( int j=0; j<4; j++ ) 
			    {
			        NE[ i               ][ j               ] = std::complex<double>( -Nnn[j][i], 0.00 );
				    NE[ i + mNumNodes   ][ j + mNumNodes   ] = std::complex<double>( -Nnn[j][i], 0.00 );
				    NE[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = std::complex<double>( -Nnn[j][i], 0.00 );
		        }
            }
		}

		Matrix< std::complex<double> > NE_tran  = NE.Transpose();
		Vector< std::complex<double> > NE_tran_x_ResV = NE_tran * ResVector;
		ResVector = NE_tran_x_ResV;
    }