
#include "VolumeElement_2eg_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{ 
    //******************************************************************************************************************************************
    // - Global indexes of the DOFs
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( mNumDofs );

		if( mPotentialsOn )
		{
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                EquationId[ i               ] = mNodes[ i ]->pDofcAx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcAy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcAz()->EquationId();
            }	

            for( int i=0; i<4; i++ ) 
            { 
                EquationId[ i + 3*mNumNodes ] = mNodes[ i ]->pDofcVs()->EquationId();
            }	
		}											   
		else										   
		{											   
            for( int i=0; i<mNumNodes; i++ ) 		   
            { 										   
                EquationId[ i               ] = mNodes[ i ]->pDofcEx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcEy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcEz()->EquationId();
            }		
		}
    }

    //******************************************************************************************************************************************
    // - Volume of the element
    //******************************************************************************************************************************************
    double VolumeElement_2eg_FullWave::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //******************************************************************************************************************************************
    // - Lagrangian 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Calculate_DN( double DN[ 3 ][ 4 ] )
    {
        double Cte = 1.0 / ( 6.0 * mVolume ); 

        // dNi/dx
        DN[ 0 ][ 0 ] = Cte * ( Y(3)*Z(2) + Y(4)*Z(3) + Y(2)*Z(4) - Y(3)*Z(4) - Y(2)*Z(3) - Y(4)*Z(2) );
        DN[ 0 ][ 1 ] = Cte * ( Y(3)*Z(4) + Y(4)*Z(1) + Y(1)*Z(3) - Y(3)*Z(1) - Y(1)*Z(4) - Y(4)*Z(3) );
        DN[ 0 ][ 2 ] = Cte * ( Y(2)*Z(1) + Y(1)*Z(4) + Y(4)*Z(2) - Y(2)*Z(4) - Y(4)*Z(1) - Y(1)*Z(2) );
        DN[ 0 ][ 3 ] = Cte * ( Y(2)*Z(3) + Y(3)*Z(1) + Y(1)*Z(2) - Y(2)*Z(1) - Y(3)*Z(2) - Y(1)*Z(3) );
        
        // dNi/dy
        DN[ 1 ][ 0 ] = Cte * ( X(3)*Z(4) + X(4)*Z(2) + X(2)*Z(3) - X(3)*Z(2) - X(2)*Z(4) - X(4)*Z(3) ); 
        DN[ 1 ][ 1 ] = Cte * ( X(3)*Z(1) + X(1)*Z(4) + X(4)*Z(3) - X(1)*Z(3) - X(3)*Z(4) - X(4)*Z(1) );
        DN[ 1 ][ 2 ] = Cte * ( X(2)*Z(4) + X(4)*Z(1) + X(1)*Z(2) - X(2)*Z(1) - X(4)*Z(2) - X(1)*Z(4) );
        DN[ 1 ][ 3 ] = Cte * ( X(1)*Z(3) + X(3)*Z(2) + X(2)*Z(1) - X(3)*Z(1) - X(1)*Z(2) - X(2)*Z(3) );
        
        // dNi/dz
        DN[ 2 ][ 0 ] = Cte * ( X(4)*Y(3) + X(2)*Y(4) + X(3)*Y(2) - X(2)*Y(3) - X(4)*Y(2) - X(3)*Y(4) );
        DN[ 2 ][ 1 ] = Cte * ( X(3)*Y(4) + X(4)*Y(1) + X(1)*Y(3) - X(3)*Y(1) - X(4)*Y(3) - X(1)*Y(4) ); 
        DN[ 2 ][ 2 ] = Cte * ( X(2)*Y(1) + X(4)*Y(2) + X(1)*Y(4) - X(4)*Y(1) - X(2)*Y(4) - X(1)*Y(2) ); 
        DN[ 2 ][ 3 ] = Cte * ( X(2)*Y(3) + X(3)*Y(1) + X(1)*Y(2) - X(2)*Y(1) - X(3)*Y(2) - X(1)*Y(3) );    
    }

    //******************************************************************************************************************************************
    // - Element edge lengths with sign
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Calculate_Edge_Lengths( Vector<double>& eLength )
    {
		// Edge coordinates vector
		Vector< Vector<double> > R( mNumEdges ); 
		
		for( int i=0; i<mNumEdges; i++ ) 
		{
			R[ i ].resize( 3 );
		}

		// Edge: [0 1]
		R[ 0 ][ 0 ] = X(2) - X(1); 
		R[ 0 ][ 1 ] = Y(2) - Y(1); 
		R[ 0 ][ 2 ] = Z(2) - Z(1);

		// Edge: [1 2]
		R[ 1 ][ 0 ] = X(3) - X(2); 
		R[ 1 ][ 1 ] = Y(3) - Y(2); 
		R[ 1 ][ 2 ] = Z(3) - Z(2);

		// Edge: [2 0]
		R[ 2 ][ 0 ] = X(1) - X(3); 
		R[ 2 ][ 1 ] = Y(1) - Y(3); 
		R[ 2 ][ 2 ] = Z(1) - Z(3);

		// Edge: [0 3]
		R[ 3 ][ 0 ] = X(4) - X(1); 
		R[ 3 ][ 1 ] = Y(4) - Y(1); 
		R[ 3 ][ 2 ] = Z(4) - Z(1);

		// Edge: [1 3]
		R[ 4 ][ 0 ] = X(4) - X(2); 
		R[ 4 ][ 1 ] = Y(4) - Y(2); 
		R[ 4 ][ 2 ] = Z(4) - Z(2);

		// Edge: [2 3]
		R[ 5 ][ 0 ] = X(4) - X(3); 
		R[ 5 ][ 1 ] = Y(4) - Y(3); 
		R[ 5 ][ 2 ] = Z(4) - Z(3);

		// Edge lengths
		for( int i=0; i<mNumEdges; i++ )
		{
			eLength[ i ] = sqrt( R[ i ][ 0 ]*R[ i ][ 0 ] + R[ i ][ 1 ]*R[ i ][ 1 ] + R[ i ][ 2 ]*R[ i ][ 2 ] );
		}
		
		// Check directions
		if( mNodes[ 0 ]->Id() > mNodes[ 1 ]->Id() ) eLength[ 0 ] *= -1.0;
		if( mNodes[ 1 ]->Id() > mNodes[ 2 ]->Id() ) eLength[ 1 ] *= -1.0;
		if( mNodes[ 2 ]->Id() > mNodes[ 0 ]->Id() ) eLength[ 2 ] *= -1.0;
		if( mNodes[ 0 ]->Id() > mNodes[ 3 ]->Id() ) eLength[ 3 ] *= -1.0;
		if( mNodes[ 1 ]->Id() > mNodes[ 3 ]->Id() ) eLength[ 4 ] *= -1.0;
		if( mNodes[ 2 ]->Id() > mNodes[ 3 ]->Id() ) eLength[ 5 ] *= -1.0;
	}

    //******************************************************************************************************************************************
    // - 1st order Whitney base for 3D tetrahedrons
    //******************************************************************************************************************************************
	void VolumeElement_2eg_FullWave::Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz,
		                                               Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ )
	{
		double L1, L2, L3, L4;
		
		int numPoints = cX.size();

		Nx.Resize( mNumEdges, numPoints );
		Ny.Resize( mNumEdges, numPoints );
		Nz.Resize( mNumEdges, numPoints );
		
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		double dL1_dX = DN[ 0 ][ 0 ], dL1_dY = DN[ 1 ][ 0 ], dL1_dZ = DN[ 2 ][ 0 ];
		double dL2_dX = DN[ 0 ][ 1 ], dL2_dY = DN[ 1 ][ 1 ], dL2_dZ = DN[ 2 ][ 1 ];
		double dL3_dX = DN[ 0 ][ 2 ], dL3_dY = DN[ 1 ][ 2 ], dL3_dZ = DN[ 2 ][ 2 ];
		double dL4_dX = DN[ 0 ][ 3 ], dL4_dY = DN[ 1 ][ 3 ], dL4_dZ = DN[ 2 ][ 3 ];

		Vector<double> eLength( mNumEdges ); Calculate_Edge_Lengths( eLength );

		for( int p=0; p<numPoints; p++ ) 
        {
            L2 = cX[ p ];
            L3 = cY[ p ];
			L4 = cZ[ p ];
            L1 = 1.0 - L2 - L3 - L4;

			// Edge: [0 1]
			Nx[ 0 ][ p ] = ( L1 * dL2_dX - L2 * dL1_dX ) * eLength[ 0 ];
            Ny[ 0 ][ p ] = ( L1 * dL2_dY - L2 * dL1_dY ) * eLength[ 0 ];
		    Nz[ 0 ][ p ] = ( L1 * dL2_dZ - L2 * dL1_dZ ) * eLength[ 0 ];

			// Edge: [1 2]
			Nx[ 1 ][ p ] = ( L2 * dL3_dX - L3 * dL2_dX ) * eLength[ 1 ];
            Ny[ 1 ][ p ] = ( L2 * dL3_dY - L3 * dL2_dY ) * eLength[ 1 ];
		    Nz[ 1 ][ p ] = ( L2 * dL3_dZ - L3 * dL2_dZ ) * eLength[ 1 ];

			// Edge: [2 0]
			Nx[ 2 ][ p ] = ( L3 * dL1_dX - L1 * dL3_dX ) * eLength[ 2 ];
            Ny[ 2 ][ p ] = ( L3 * dL1_dY - L1 * dL3_dY ) * eLength[ 2 ];
		    Nz[ 2 ][ p ] = ( L3 * dL1_dZ - L1 * dL3_dZ ) * eLength[ 2 ];

			// Edge: [0 3]
			Nx[ 3 ][ p ] = ( L1 * dL4_dX - L4 * dL1_dX ) * eLength[ 3 ];
            Ny[ 3 ][ p ] = ( L1 * dL4_dY - L4 * dL1_dY ) * eLength[ 3 ];
		    Nz[ 3 ][ p ] = ( L1 * dL4_dZ - L4 * dL1_dZ ) * eLength[ 3 ];

			// Edge: [1 3]
			Nx[ 4 ][ p ] = ( L2 * dL4_dX - L4 * dL2_dX ) * eLength[ 4 ];
            Ny[ 4 ][ p ] = ( L2 * dL4_dY - L4 * dL2_dY ) * eLength[ 4 ];
		    Nz[ 4 ][ p ] = ( L2 * dL4_dZ - L4 * dL2_dZ ) * eLength[ 4 ];

			// Edge: [2 3]
			Nx[ 5 ][ p ] = ( L3 * dL4_dX - L4 * dL3_dX ) * eLength[ 5 ];
            Ny[ 5 ][ p ] = ( L3 * dL4_dY - L4 * dL3_dY ) * eLength[ 5 ];
		    Nz[ 5 ][ p ] = ( L3 * dL4_dZ - L4 * dL3_dZ ) * eLength[ 5 ];
        }
	}

    //******************************************************************************************************************************************
    // - Edge IDs and opposite node IDs for a Whitney face base
    //******************************************************************************************************************************************
	void VolumeElement_2eg_FullWave::Face_Base( int edgeBase[ 2 ], int nodeBase[ 2 ], int nodeID[ 3 ], int edgeID[ 3 ], int oppnID[ 3 ] )
	{
		if( ( nodeID[ 0 ] < nodeID[ 1 ] ) && ( nodeID[ 0 ] < nodeID[ 2 ] ) )
		{
			if( nodeID[ 1 ] < nodeID[ 2 ] )
			{
			    edgeBase[ 0 ] = edgeID[ 0 ]; nodeBase[ 0 ] = oppnID[ 0 ]; 
				edgeBase[ 1 ] = edgeID[ 2 ]; nodeBase[ 1 ] = oppnID[ 2 ]; 
			}
			else
			{
			    edgeBase[ 0 ] = edgeID[ 2 ]; nodeBase[ 0 ] = oppnID[ 2 ]; 
				edgeBase[ 1 ] = edgeID[ 0 ]; nodeBase[ 1 ] = oppnID[ 0 ]; 
			}
		}
		else if( ( nodeID[ 1 ] < nodeID[ 0 ] ) && ( nodeID[ 1 ] < nodeID[ 2 ] ) )
		{
			if( nodeID[ 0 ] < nodeID[ 2 ] )
			{
			    edgeBase[ 0 ] = edgeID[ 0 ]; nodeBase[ 0 ] = oppnID[ 0 ]; 
				edgeBase[ 1 ] = edgeID[ 1 ]; nodeBase[ 1 ] = oppnID[ 1 ]; 
			}
			else
			{
			    edgeBase[ 0 ] = edgeID[ 1 ]; nodeBase[ 0 ] = oppnID[ 1 ]; 
				edgeBase[ 1 ] = edgeID[ 0 ]; nodeBase[ 1 ] = oppnID[ 0 ]; 	
			}
		}
		else if( ( nodeID[ 2 ] < nodeID[ 0 ] ) && ( nodeID[ 2 ] < nodeID[ 1 ] ) )
		{
			if( nodeID[ 0 ] < nodeID[ 1 ] )
			{
			    edgeBase[ 0 ] = edgeID[ 2 ]; nodeBase[ 0 ] = oppnID[ 2 ]; 
				edgeBase[ 1 ] = edgeID[ 1 ]; nodeBase[ 1 ] = oppnID[ 1 ]; 	
			}
			else
			{
			    edgeBase[ 0 ] = edgeID[ 1 ]; nodeBase[ 0 ] = oppnID[ 1 ]; 
				edgeBase[ 1 ] = edgeID[ 2 ]; nodeBase[ 1 ] = oppnID[ 2 ]; 	
			}
		}		
	}

    //******************************************************************************************************************************************
    // - Whitney face bases ordering - Fy:[0-3] - Fz:[4-7]
    //******************************************************************************************************************************************
	void VolumeElement_2eg_FullWave::Order_Face_Bases( int n[ 8 ], int e[ 8 ] )
	{
		// Local face base IDs 
		int edgeBase[ 2 ];
		int nodeBase[ 2 ];

		// Global face base IDs
		int nodeID[ 3 ];
		int edgeID[ 3 ];
		int oppnID[ 3 ];
        
		// Face: [0 1 2]
		nodeID[ 0 ] = mNodes[ 0 ]->Id(); 
		nodeID[ 1 ] = mNodes[ 1 ]->Id(); 
		nodeID[ 2 ] = mNodes[ 2 ]->Id();

		edgeID[ 0 ] = 0; edgeID[ 1 ] = 1; edgeID[ 2 ] = 2;
		oppnID[ 0 ] = 2; oppnID[ 1 ] = 0; oppnID[ 2 ] = 1;

		Face_Base( edgeBase, nodeBase, nodeID, edgeID, oppnID );

		// Fy - Face: [0 1 2]
		e[ 0 ] = edgeBase[ 0 ]; 
		n[ 0 ] = nodeBase[ 0 ]; 

		// Fz - Face: [0 1 2]
		e[ 4 ] = edgeBase[ 1 ]; 
		n[ 4 ] = nodeBase[ 1 ]; 

		// Face: [0 1 3]
		nodeID[ 0 ] = mNodes[ 0 ]->Id(); 
		nodeID[ 1 ] = mNodes[ 1 ]->Id(); 
		nodeID[ 2 ] = mNodes[ 3 ]->Id();

		edgeID[ 0 ] = 0; edgeID[ 1 ] = 4; edgeID[ 2 ] = 3;
		oppnID[ 0 ] = 3; oppnID[ 1 ] = 0; oppnID[ 2 ] = 1;

		Face_Base( edgeBase, nodeBase, nodeID, edgeID, oppnID );

		// Fy - Face: [0 1 3]
		e[ 1 ] = edgeBase[ 0 ]; 
		n[ 1 ] = nodeBase[ 0 ]; 

		// Fz - Face: [0 1 3]
		e[ 5 ] = edgeBase[ 1 ]; 
		n[ 5 ] = nodeBase[ 1 ]; 

		// Face: [1 2 3]
		nodeID[ 0 ] = mNodes[ 1 ]->Id(); 
		nodeID[ 1 ] = mNodes[ 2 ]->Id(); 
		nodeID[ 2 ] = mNodes[ 3 ]->Id();

		edgeID[ 0 ] = 1; edgeID[ 1 ] = 5; edgeID[ 2 ] = 4;
		oppnID[ 0 ] = 3; oppnID[ 1 ] = 1; oppnID[ 2 ] = 2;

		Face_Base( edgeBase, nodeBase, nodeID, edgeID, oppnID );

		// Fy - Face: [1 2 3] 
		e[ 2 ] = edgeBase[ 0 ]; 
		n[ 2 ] = nodeBase[ 0 ]; 

		// Fz - Face: [1 2 3] 
		e[ 6 ] = edgeBase[ 1 ]; 
		n[ 6 ] = nodeBase[ 1 ]; 	

		// Face: [0 2 3]
		nodeID[ 0 ] = mNodes[ 0 ]->Id(); 
		nodeID[ 1 ] = mNodes[ 2 ]->Id(); 
		nodeID[ 2 ] = mNodes[ 3 ]->Id();

		edgeID[ 0 ] = 2; edgeID[ 1 ] = 5; edgeID[ 2 ] = 3;
		oppnID[ 0 ] = 3; oppnID[ 1 ] = 0; oppnID[ 2 ] = 2;

		Face_Base( edgeBase, nodeBase, nodeID, edgeID, oppnID );

		// Fy - Face: [0 2 3] 
		e[ 3 ] = edgeBase[ 0 ]; 
		n[ 3 ] = nodeBase[ 0 ]; 

		// Fz - Face: [0 2 3] 
		e[ 7 ] = edgeBase[ 1 ]; 
		n[ 7 ] = nodeBase[ 1 ]; 
	}

    //******************************************************************************************************************************************
    // - 2nd order Whitney base for 3D tetrahedrons
    //******************************************************************************************************************************************
	void VolumeElement_2eg_FullWave::Whitney3D_Ni_2nd( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz,
		                                               Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ )
	{
		int numPoints = cX.size();

		Nx.Resize( mNumBases, numPoints );
		Ny.Resize( mNumBases, numPoints );
		Nz.Resize( mNumBases, numPoints );

		int n[ 8 ], e[ 8 ]; Order_Face_Bases( n, e );  

		Matrix<double> L; Lagrange3D_Ni_1st( L, cX, cY, cZ );

	    Matrix<double> Wx, Wy, Wz; Whitney3D_Ni_1st( Wx, Wy, Wz, cX, cY, cZ );

		for( int p=0; p<numPoints; p++ ) 
        {
			// Edge: [0 1]
			Nx[ 0 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wx[ 0 ][ p ];
            Ny[ 0 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wy[ 0 ][ p ];
		    Nz[ 0 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wz[ 0 ][ p ];
														   
			Nx[ 1 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wx[ 0 ][ p ];
            Ny[ 1 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wy[ 0 ][ p ];
		    Nz[ 1 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wz[ 0 ][ p ];

			// Edge: [1 2]
			Nx[ 2 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wx[ 1 ][ p ];
            Ny[ 2 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wy[ 1 ][ p ];
		    Nz[ 2 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wz[ 1 ][ p ];
														   
			Nx[ 3 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wx[ 1 ][ p ];
            Ny[ 3 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wy[ 1 ][ p ];
		    Nz[ 3 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wz[ 1 ][ p ];

			// Edge: [2 0]
			Nx[ 4 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wx[ 2 ][ p ];
            Ny[ 4 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wy[ 2 ][ p ];
		    Nz[ 4 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wz[ 2 ][ p ];
														   
			Nx[ 5 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wx[ 2 ][ p ];
            Ny[ 5 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wy[ 2 ][ p ];
		    Nz[ 5 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wz[ 2 ][ p ];
			
			// Edge: [0 3]
			Nx[ 6 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wx[ 3 ][ p ];
            Ny[ 6 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wy[ 3 ][ p ];
		    Nz[ 6 ][ p ] = ( 3.0 * L[ 0 ][ p ] - 1.0 ) * Wz[ 3 ][ p ];
														   
			Nx[ 7 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wx[ 3 ][ p ];
            Ny[ 7 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wy[ 3 ][ p ];
		    Nz[ 7 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wz[ 3 ][ p ];			

			// Edge: [1 3]
			Nx[ 8 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wx[ 4 ][ p ];
            Ny[ 8 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wy[ 4 ][ p ];
		    Nz[ 8 ][ p ] = ( 3.0 * L[ 1 ][ p ] - 1.0 ) * Wz[ 4 ][ p ];
													   
			Nx[ 9 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wx[ 4 ][ p ];
            Ny[ 9 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wy[ 4 ][ p ];
		    Nz[ 9 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wz[ 4 ][ p ];				

			// Edge: [2 3]
			Nx[ 10 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wx[ 5 ][ p ];
            Ny[ 10 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wy[ 5 ][ p ];
		    Nz[ 10 ][ p ] = ( 3.0 * L[ 2 ][ p ] - 1.0 ) * Wz[ 5 ][ p ];
													   
			Nx[ 11 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wx[ 5 ][ p ];
            Ny[ 11 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wy[ 5 ][ p ];
		    Nz[ 11 ][ p ] = ( 3.0 * L[ 3 ][ p ] - 1.0 ) * Wz[ 5 ][ p ];			

			// Face: [0 1 2] - Fy
			Nx[ 12 ][ p ] = 4.5 * L[ n[ 0 ] ][ p ] * Wx[ e[ 0 ] ][ p ];
            Ny[ 12 ][ p ] = 4.5 * L[ n[ 0 ] ][ p ] * Wy[ e[ 0 ] ][ p ];
		    Nz[ 12 ][ p ] = 4.5 * L[ n[ 0 ] ][ p ] * Wz[ e[ 0 ] ][ p ];
					
			// Face: [0 1 3] - Fy
			Nx[ 13 ][ p ] = 4.5 * L[ n[ 1 ] ][ p ] * Wx[ e[ 1 ] ][ p ];
            Ny[ 13 ][ p ] = 4.5 * L[ n[ 1 ] ][ p ] * Wy[ e[ 1 ] ][ p ];
		    Nz[ 13 ][ p ] = 4.5 * L[ n[ 1 ] ][ p ] * Wz[ e[ 1 ] ][ p ];			

			// Face: [1 2 3] - Fy
			Nx[ 14 ][ p ] = 4.5 * L[ n[ 2 ] ][ p ] * Wx[ e[ 2 ] ][ p ];
            Ny[ 14 ][ p ] = 4.5 * L[ n[ 2 ] ][ p ] * Wy[ e[ 2 ] ][ p ];
		    Nz[ 14 ][ p ] = 4.5 * L[ n[ 2 ] ][ p ] * Wz[ e[ 2 ] ][ p ];
				
			// Face: [0 2 3] - Fy
			Nx[ 15 ][ p ] = 4.5 * L[ n[ 3 ] ][ p ] * Wx[ e[ 3 ] ][ p ];
            Ny[ 15 ][ p ] = 4.5 * L[ n[ 3 ] ][ p ] * Wy[ e[ 3 ] ][ p ];
		    Nz[ 15 ][ p ] = 4.5 * L[ n[ 3 ] ][ p ] * Wz[ e[ 3 ] ][ p ];		

			// Face: [0 1 2] - Fz
			Nx[ 16 ][ p ] = 4.5 * L[ n[ 4 ] ][ p ] * Wx[ e[ 4 ] ][ p ];
            Ny[ 16 ][ p ] = 4.5 * L[ n[ 4 ] ][ p ] * Wy[ e[ 4 ] ][ p ];
		    Nz[ 16 ][ p ] = 4.5 * L[ n[ 4 ] ][ p ] * Wz[ e[ 4 ] ][ p ];
									  
			// Face: [0 1 3] - Fz
			Nx[ 17 ][ p ] = 4.5 * L[ n[ 5 ] ][ p ] * Wx[ e[ 5 ] ][ p ];
            Ny[ 17 ][ p ] = 4.5 * L[ n[ 5 ] ][ p ] * Wy[ e[ 5 ] ][ p ];
		    Nz[ 17 ][ p ] = 4.5 * L[ n[ 5 ] ][ p ] * Wz[ e[ 5 ] ][ p ];		

			// Face: [1 2 3] - Fz
			Nx[ 18 ][ p ] = 4.5 * L[ n[ 6 ] ][ p ] * Wx[ e[ 6 ] ][ p ];
            Ny[ 18 ][ p ] = 4.5 * L[ n[ 6 ] ][ p ] * Wy[ e[ 6 ] ][ p ];
		    Nz[ 18 ][ p ] = 4.5 * L[ n[ 6 ] ][ p ] * Wz[ e[ 6 ] ][ p ];
									
			// Face: [0 2 3] - Fz
			Nx[ 19 ][ p ] = 4.5 * L[ n[ 7 ] ][ p ] * Wx[ e[ 7 ] ][ p ];
            Ny[ 19 ][ p ] = 4.5 * L[ n[ 7 ] ][ p ] * Wy[ e[ 7 ] ][ p ];
		    Nz[ 19 ][ p ] = 4.5 * L[ n[ 7 ] ][ p ] * Wz[ e[ 7 ] ][ p ];		
        }
	}

    //******************************************************************************************************************************************
    // - Derivatives of the 1st order Whitney base 
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Derivatives_Whitney_1st( Vector<double>& dNx_dX, Vector<double>& dNx_dY, Vector<double>& dNx_dZ,
		                                                      Vector<double>& dNy_dX, Vector<double>& dNy_dY, Vector<double>& dNy_dZ,
		                                                      Vector<double>& dNz_dX, Vector<double>& dNz_dY, Vector<double>& dNz_dZ )
    {
		// Li derivatives
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		double dL1_dX = DN[ 0 ][ 0 ], dL1_dY = DN[ 1 ][ 0 ], dL1_dZ = DN[ 2 ][ 0 ];
		double dL2_dX = DN[ 0 ][ 1 ], dL2_dY = DN[ 1 ][ 1 ], dL2_dZ = DN[ 2 ][ 1 ];
		double dL3_dX = DN[ 0 ][ 2 ], dL3_dY = DN[ 1 ][ 2 ], dL3_dZ = DN[ 2 ][ 2 ];
		double dL4_dX = DN[ 0 ][ 3 ], dL4_dY = DN[ 1 ][ 3 ], dL4_dZ = DN[ 2 ][ 3 ];

		// Element edges lengths
		Vector<double> eLength( mNumEdges ); Calculate_Edge_Lengths( eLength );

		// Whitney base derivatives
		dNx_dX.resize( mNumEdges ); 
		dNx_dY.resize( mNumEdges ); 
		dNx_dZ.resize( mNumEdges ); 
				   			
		dNy_dX.resize( mNumEdges ); 
		dNy_dY.resize( mNumEdges ); 
		dNy_dZ.resize( mNumEdges );
				  			
		dNz_dX.resize( mNumEdges ); 
		dNz_dY.resize( mNumEdges ); 
		dNz_dZ.resize( mNumEdges );

        // Edge: [0 1]
		dNx_dX[ 0 ] = ( dL1_dX * dL2_dX - dL2_dX * dL1_dX ) * eLength[ 0 ];
		dNx_dY[ 0 ] = ( dL1_dY * dL2_dX - dL2_dY * dL1_dX ) * eLength[ 0 ];
		dNx_dZ[ 0 ] = ( dL1_dZ * dL2_dX - dL2_dZ * dL1_dX ) * eLength[ 0 ];
			    													   
  		dNy_dX[ 0 ] = ( dL1_dX * dL2_dY - dL2_dX * dL1_dY ) * eLength[ 0 ];
		dNy_dY[ 0 ] = ( dL1_dY * dL2_dY - dL2_dY * dL1_dY ) * eLength[ 0 ];
		dNy_dZ[ 0 ] = ( dL1_dZ * dL2_dY - dL2_dZ * dL1_dY ) * eLength[ 0 ];
			    													   
  		dNz_dX[ 0 ] = ( dL1_dX * dL2_dZ - dL2_dX * dL1_dZ ) * eLength[ 0 ];
		dNz_dY[ 0 ] = ( dL1_dY * dL2_dZ - dL2_dY * dL1_dZ ) * eLength[ 0 ];
		dNz_dZ[ 0 ] = ( dL1_dZ * dL2_dZ - dL2_dZ * dL1_dZ ) * eLength[ 0 ];

		// Edge: [1 2]
		dNx_dX[ 1 ] = ( dL2_dX * dL3_dX - dL3_dX * dL2_dX ) * eLength[ 1 ];
		dNx_dY[ 1 ] = ( dL2_dY * dL3_dX - dL3_dY * dL2_dX ) * eLength[ 1 ];
		dNx_dZ[ 1 ] = ( dL2_dZ * dL3_dX - dL3_dZ * dL2_dX ) * eLength[ 1 ];
			    	 												   
 		dNy_dX[ 1 ] = ( dL2_dX * dL3_dY - dL3_dX * dL2_dY ) * eLength[ 1 ];
		dNy_dY[ 1 ] = ( dL2_dY * dL3_dY - dL3_dY * dL2_dY ) * eLength[ 1 ];
		dNy_dZ[ 1 ] = ( dL2_dZ * dL3_dY - dL3_dZ * dL2_dY ) * eLength[ 1 ];
			    	 										   		   
 		dNz_dX[ 1 ] = ( dL2_dX * dL3_dZ - dL3_dX * dL2_dZ ) * eLength[ 1 ];
		dNz_dY[ 1 ] = ( dL2_dY * dL3_dZ - dL3_dY * dL2_dZ ) * eLength[ 1 ];
		dNz_dZ[ 1 ] = ( dL2_dZ * dL3_dZ - dL3_dZ * dL2_dZ ) * eLength[ 1 ];

		// Edge: [2 0]
		dNx_dX[ 2 ] = ( dL3_dX * dL1_dX - dL1_dX * dL3_dX ) * eLength[ 2 ];
		dNx_dY[ 2 ] = ( dL3_dY * dL1_dX - dL1_dY * dL3_dX ) * eLength[ 2 ];
		dNx_dZ[ 2 ] = ( dL3_dZ * dL1_dX - dL1_dZ * dL3_dX ) * eLength[ 2 ];
			    	 												     
		dNy_dX[ 2 ] = ( dL3_dX * dL1_dY - dL1_dX * dL3_dY ) * eLength[ 2 ];
		dNy_dY[ 2 ] = ( dL3_dY * dL1_dY - dL1_dY * dL3_dY ) * eLength[ 2 ];
		dNy_dZ[ 2 ] = ( dL3_dZ * dL1_dY - dL1_dZ * dL3_dY ) * eLength[ 2 ];
			    	 										   		   
 		dNz_dX[ 2 ] = ( dL3_dX * dL1_dZ - dL1_dX * dL3_dZ ) * eLength[ 2 ];
		dNz_dY[ 2 ] = ( dL3_dY * dL1_dZ - dL1_dY * dL3_dZ ) * eLength[ 2 ];
		dNz_dZ[ 2 ] = ( dL3_dZ * dL1_dZ - dL1_dZ * dL3_dZ ) * eLength[ 2 ];

		// Edge: [0 3]
		dNx_dX[ 3 ] = ( dL1_dX * dL4_dX - dL4_dX * dL1_dX ) * eLength[ 3 ];
		dNx_dY[ 3 ] = ( dL1_dY * dL4_dX - dL4_dY * dL1_dX ) * eLength[ 3 ];
		dNx_dZ[ 3 ] = ( dL1_dZ * dL4_dX - dL4_dZ * dL1_dX ) * eLength[ 3 ];
			  	  												   	   
 		dNy_dX[ 3 ] = ( dL1_dX * dL4_dY - dL4_dX * dL1_dY ) * eLength[ 3 ];
		dNy_dY[ 3 ] = ( dL1_dY * dL4_dY - dL4_dY * dL1_dY ) * eLength[ 3 ];
		dNy_dZ[ 3 ] = ( dL1_dZ * dL4_dY - dL4_dZ * dL1_dY ) * eLength[ 3 ];
			    	 												     
		dNz_dX[ 3 ] = ( dL1_dX * dL4_dZ - dL4_dX * dL1_dZ ) * eLength[ 3 ];
		dNz_dY[ 3 ] = ( dL1_dY * dL4_dZ - dL4_dY * dL1_dZ ) * eLength[ 3 ];
		dNz_dZ[ 3 ] = ( dL1_dZ * dL4_dZ - dL4_dZ * dL1_dZ ) * eLength[ 3 ];

		// Edge: [1 3]
		dNx_dX[ 4 ] = ( dL2_dX * dL4_dX - dL4_dX * dL2_dX ) * eLength[ 4 ];
		dNx_dY[ 4 ] = ( dL2_dY * dL4_dX - dL4_dY * dL2_dX ) * eLength[ 4 ];
		dNx_dZ[ 4 ] = ( dL2_dZ * dL4_dX - dL4_dZ * dL2_dX ) * eLength[ 4 ];
			  	  												   	   
 		dNy_dX[ 4 ] = ( dL2_dX * dL4_dY - dL4_dX * dL2_dY ) * eLength[ 4 ];
		dNy_dY[ 4 ] = ( dL2_dY * dL4_dY - dL4_dY * dL2_dY ) * eLength[ 4 ];
		dNy_dZ[ 4 ] = ( dL2_dZ * dL4_dY - dL4_dZ * dL2_dY ) * eLength[ 4 ];
			    	 												     
		dNz_dX[ 4 ] = ( dL2_dX * dL4_dZ - dL4_dX * dL2_dZ ) * eLength[ 4 ];
		dNz_dY[ 4 ] = ( dL2_dY * dL4_dZ - dL4_dY * dL2_dZ ) * eLength[ 4 ];
		dNz_dZ[ 4 ] = ( dL2_dZ * dL4_dZ - dL4_dZ * dL2_dZ ) * eLength[ 4 ];

		// Edge: [2 3]
		dNx_dX[ 5 ] = ( dL3_dX * dL4_dX - dL4_dX * dL3_dX ) * eLength[ 5 ];
		dNx_dY[ 5 ] = ( dL3_dY * dL4_dX - dL4_dY * dL3_dX ) * eLength[ 5 ];
		dNx_dZ[ 5 ] = ( dL3_dZ * dL4_dX - dL4_dZ * dL3_dX ) * eLength[ 5 ];
			  	  												   	   
 		dNy_dX[ 5 ] = ( dL3_dX * dL4_dY - dL4_dX * dL3_dY ) * eLength[ 5 ];
		dNy_dY[ 5 ] = ( dL3_dY * dL4_dY - dL4_dY * dL3_dY ) * eLength[ 5 ];
		dNy_dZ[ 5 ] = ( dL3_dZ * dL4_dY - dL4_dZ * dL3_dY ) * eLength[ 5 ];
			  	  												   	   
 		dNz_dX[ 5 ] = ( dL3_dX * dL4_dZ - dL4_dX * dL3_dZ ) * eLength[ 5 ];
		dNz_dY[ 5 ] = ( dL3_dY * dL4_dZ - dL4_dY * dL3_dZ ) * eLength[ 5 ];
		dNz_dZ[ 5 ] = ( dL3_dZ * dL4_dZ - dL4_dZ * dL3_dZ ) * eLength[ 5 ];		
	}

    //******************************************************************************************************************************************
    // - Derivatives of the 2nd order Whitney base on cX, cY, cZ points
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Derivatives_On_Points( Matrix<double>& dNx_dX, Matrix<double>& dNx_dY, Matrix<double>& dNx_dZ,
		                                                    Matrix<double>& dNy_dX, Matrix<double>& dNy_dY, Matrix<double>& dNy_dZ,
		                                                    Matrix<double>& dNz_dX, Matrix<double>& dNz_dY, Matrix<double>& dNz_dZ,
		                                                    Vector<double>&     cX, Vector<double>&     cY, Vector<double>&     cZ )
    {
		// Number of points
		int numPoints = cX.size();

		// Whitney 2nd order base derivatives
		dNx_dX.Resize( mNumBases, numPoints ); 
		dNx_dY.Resize( mNumBases, numPoints ); 
		dNx_dZ.Resize( mNumBases, numPoints ); 
				   		
		dNy_dX.Resize( mNumBases, numPoints ); 
		dNy_dY.Resize( mNumBases, numPoints ); 
		dNy_dZ.Resize( mNumBases, numPoints );
				  		
		dNz_dX.Resize( mNumBases, numPoints ); 
		dNz_dY.Resize( mNumBases, numPoints ); 
		dNz_dZ.Resize( mNumBases, numPoints );

		// Nodal 1st order base
		Matrix<double> L; Lagrange3D_Ni_1st( L, cX, cY, cZ );

		// Edge 1st order base
	    Matrix<double> Wx, Wy, Wz; Whitney3D_Ni_1st( Wx, Wy, Wz, cX, cY, cZ );

		// Li derivatives
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );
		
		double dL_dX[ 4 ], dL_dY[ 4 ], dL_dZ[ 4 ];

		dL_dX[ 0 ] = DN[ 0 ][ 0 ], dL_dY[ 0 ] = DN[ 1 ][ 0 ], dL_dZ[ 0 ] = DN[ 2 ][ 0 ];
		dL_dX[ 1 ] = DN[ 0 ][ 1 ], dL_dY[ 1 ] = DN[ 1 ][ 1 ], dL_dZ[ 1 ] = DN[ 2 ][ 1 ];
		dL_dX[ 2 ] = DN[ 0 ][ 2 ], dL_dY[ 2 ] = DN[ 1 ][ 2 ], dL_dZ[ 2 ] = DN[ 2 ][ 2 ];
		dL_dX[ 3 ] = DN[ 0 ][ 3 ], dL_dY[ 3 ] = DN[ 1 ][ 3 ], dL_dZ[ 3 ] = DN[ 2 ][ 3 ];

		// Wi derivatives
		Vector<double> dWx_dX, dWx_dY, dWx_dZ;
		Vector<double> dWy_dX, dWy_dY, dWy_dZ;
		Vector<double> dWz_dX, dWz_dY, dWz_dZ;

		Derivatives_Whitney_1st( dWx_dX, dWx_dY, dWx_dZ, dWy_dX, dWy_dY, dWy_dZ, dWz_dX, dWz_dY, dWz_dZ );
		
		// Face base order
		int n[ 8 ], e[ 8 ]; Order_Face_Bases( n, e );  

		for( int p=0; p<numPoints; p++ ) 
        {
			// Edge: [0 1]
			dNx_dX[ 0 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wx[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dX[ 0 ];
			dNx_dY[ 0 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wx[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dY[ 0 ];
			dNx_dZ[ 0 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wx[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dZ[ 0 ];

			dNy_dX[ 0 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wy[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dX[ 0 ];
			dNy_dY[ 0 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wy[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dY[ 0 ];
			dNy_dZ[ 0 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wy[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dZ[ 0 ];

			dNz_dX[ 0 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wz[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dX[ 0 ];
			dNz_dY[ 0 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wz[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dY[ 0 ];
			dNz_dZ[ 0 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wz[ 0 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dZ[ 0 ];

			dNx_dX[ 1 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wx[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dX[ 0 ];
			dNx_dY[ 1 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wx[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dY[ 0 ];
			dNx_dZ[ 1 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wx[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dZ[ 0 ];
					
			dNy_dX[ 1 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wy[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dX[ 0 ];
			dNy_dY[ 1 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wy[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dY[ 0 ];
			dNy_dZ[ 1 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wy[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dZ[ 0 ];
					
			dNz_dX[ 1 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wz[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dX[ 0 ];
			dNz_dY[ 1 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wz[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dY[ 0 ];
			dNz_dZ[ 1 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wz[ 0 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dZ[ 0 ];

			// Edge: [1 2]
			dNx_dX[ 2 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wx[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dX[ 1 ];
			dNx_dY[ 2 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wx[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dY[ 1 ];
			dNx_dZ[ 2 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wx[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dZ[ 1 ];
											  
			dNy_dX[ 2 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wy[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dX[ 1 ];
			dNy_dY[ 2 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wy[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dY[ 1 ];
			dNy_dZ[ 2 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wy[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dZ[ 1 ];
											 
			dNz_dX[ 2 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wz[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dX[ 1 ];
			dNz_dY[ 2 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wz[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dY[ 1 ];
			dNz_dZ[ 2 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wz[ 1 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dZ[ 1 ];

			dNx_dX[ 3 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wx[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dX[ 1 ];
			dNx_dY[ 3 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wx[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dY[ 1 ];
			dNx_dZ[ 3 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wx[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dZ[ 1 ];
										
			dNy_dX[ 3 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wy[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dX[ 1 ];
			dNy_dY[ 3 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wy[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dY[ 1 ];
			dNy_dZ[ 3 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wy[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dZ[ 1 ];
									
			dNz_dX[ 3 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wz[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dX[ 1 ];
			dNz_dY[ 3 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wz[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dY[ 1 ];
			dNz_dZ[ 3 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wz[ 1 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dZ[ 1 ];		

			// Edge: [2 0]
			dNx_dX[ 4 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wx[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dX[ 2 ];
			dNx_dY[ 4 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wx[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dY[ 2 ];
			dNx_dZ[ 4 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wx[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dZ[ 2 ];
										  
			dNy_dX[ 4 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wy[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dX[ 2 ];
			dNy_dY[ 4 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wy[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dY[ 2 ];
			dNy_dZ[ 4 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wy[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dZ[ 2 ];
											 			 
			dNz_dX[ 4 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wz[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dX[ 2 ];
			dNz_dY[ 4 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wz[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dY[ 2 ];
			dNz_dZ[ 4 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wz[ 2 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dZ[ 2 ];
														 
			dNx_dX[ 5 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wx[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dX[ 2 ];
			dNx_dY[ 5 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wx[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dY[ 2 ];
			dNx_dZ[ 5 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wx[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dZ[ 2 ];
											  
			dNy_dX[ 5 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wy[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dX[ 2 ];
			dNy_dY[ 5 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wy[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dY[ 2 ];
			dNy_dZ[ 5 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wy[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dZ[ 2 ];
												   
			dNz_dX[ 5 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wz[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dX[ 2 ];
			dNz_dY[ 5 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wz[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dY[ 2 ];
			dNz_dZ[ 5 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wz[ 2 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dZ[ 2 ];	

			// Edge: [0 3]
			dNx_dX[ 6 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wx[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dX[ 3 ];
			dNx_dY[ 6 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wx[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dY[ 3 ];
			dNx_dZ[ 6 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wx[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWx_dZ[ 3 ];
										  	 
			dNy_dX[ 6 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wy[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dX[ 3 ];
			dNy_dY[ 6 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wy[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dY[ 3 ];
			dNy_dZ[ 6 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wy[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWy_dZ[ 3 ];
											 
			dNz_dX[ 6 ][ p ] = ( 3.0 * dL_dX[ 0 ] ) * Wz[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dX[ 3 ];
			dNz_dY[ 6 ][ p ] = ( 3.0 * dL_dY[ 0 ] ) * Wz[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dY[ 3 ];
			dNz_dZ[ 6 ][ p ] = ( 3.0 * dL_dZ[ 0 ] ) * Wz[ 3 ][ p ] + ( 3.0 * L[ 0 ][ p ] - 1.0 ) * dWz_dZ[ 3 ];
										
			dNx_dX[ 7 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wx[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dX[ 3 ];
			dNx_dY[ 7 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wx[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dY[ 3 ];
			dNx_dZ[ 7 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wx[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dZ[ 3 ];
											
			dNy_dX[ 7 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wy[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dX[ 3 ];
			dNy_dY[ 7 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wy[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dY[ 3 ];
			dNy_dZ[ 7 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wy[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dZ[ 3 ];
											
			dNz_dX[ 7 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wz[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dX[ 3 ];
			dNz_dY[ 7 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wz[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dY[ 3 ];
			dNz_dZ[ 7 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wz[ 3 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dZ[ 3 ];	

			// Edge: [1 3]
			dNx_dX[ 8 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wx[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dX[ 4 ];
			dNx_dY[ 8 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wx[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dY[ 4 ];
			dNx_dZ[ 8 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wx[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWx_dZ[ 4 ];
									  	 	
			dNy_dX[ 8 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wy[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dX[ 4 ];
			dNy_dY[ 8 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wy[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dY[ 4 ];
			dNy_dZ[ 8 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wy[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWy_dZ[ 4 ];
										
			dNz_dX[ 8 ][ p ] = ( 3.0 * dL_dX[ 1 ] ) * Wz[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dX[ 4 ];
			dNz_dY[ 8 ][ p ] = ( 3.0 * dL_dY[ 1 ] ) * Wz[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dY[ 4 ];
			dNz_dZ[ 8 ][ p ] = ( 3.0 * dL_dZ[ 1 ] ) * Wz[ 4 ][ p ] + ( 3.0 * L[ 1 ][ p ] - 1.0 ) * dWz_dZ[ 4 ];
						
			dNx_dX[ 9 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wx[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dX[ 4 ];
			dNx_dY[ 9 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wx[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dY[ 4 ];
			dNx_dZ[ 9 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wx[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dZ[ 4 ];
														
			dNy_dX[ 9 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wy[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dX[ 4 ];
			dNy_dY[ 9 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wy[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dY[ 4 ];
			dNy_dZ[ 9 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wy[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dZ[ 4 ];
													
			dNz_dX[ 9 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wz[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dX[ 4 ];
			dNz_dY[ 9 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wz[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dY[ 4 ];
			dNz_dZ[ 9 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wz[ 4 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dZ[ 4 ];	

			// Edge: [2 3]
			dNx_dX[ 10 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wx[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dX[ 5 ];
			dNx_dY[ 10 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wx[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dY[ 5 ];
			dNx_dZ[ 10 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wx[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWx_dZ[ 5 ];
			
			dNy_dX[ 10 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wy[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dX[ 5 ];
			dNy_dY[ 10 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wy[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dY[ 5 ];
			dNy_dZ[ 10 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wy[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWy_dZ[ 5 ];
			
			dNz_dX[ 10 ][ p ] = ( 3.0 * dL_dX[ 2 ] ) * Wz[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dX[ 5 ];
			dNz_dY[ 10 ][ p ] = ( 3.0 * dL_dY[ 2 ] ) * Wz[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dY[ 5 ];
			dNz_dZ[ 10 ][ p ] = ( 3.0 * dL_dZ[ 2 ] ) * Wz[ 5 ][ p ] + ( 3.0 * L[ 2 ][ p ] - 1.0 ) * dWz_dZ[ 5 ];
					
			dNx_dX[ 11 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wx[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dX[ 5 ];
			dNx_dY[ 11 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wx[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dY[ 5 ];
			dNx_dZ[ 11 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wx[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWx_dZ[ 5 ];
			
			dNy_dX[ 11 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wy[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dX[ 5 ];
			dNy_dY[ 11 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wy[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dY[ 5 ];
			dNy_dZ[ 11 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wy[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWy_dZ[ 5 ];
			
			dNz_dX[ 11 ][ p ] = ( 3.0 * dL_dX[ 3 ] ) * Wz[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dX[ 5 ];
			dNz_dY[ 11 ][ p ] = ( 3.0 * dL_dY[ 3 ] ) * Wz[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dY[ 5 ];
			dNz_dZ[ 11 ][ p ] = ( 3.0 * dL_dZ[ 3 ] ) * Wz[ 5 ][ p ] + ( 3.0 * L[ 3 ][ p ] - 1.0 ) * dWz_dZ[ 5 ];	

			// Face: [0 1 2] - Fy
			dNx_dX[ 12 ][ p ] = 4.5 * ( dL_dX[ n[ 0 ] ] * Wx[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWx_dX[ e[ 0 ] ] );
			dNx_dY[ 12 ][ p ] = 4.5 * ( dL_dY[ n[ 0 ] ] * Wx[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWx_dY[ e[ 0 ] ] );
			dNx_dZ[ 12 ][ p ] = 4.5 * ( dL_dZ[ n[ 0 ] ] * Wx[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWx_dZ[ e[ 0 ] ] );

			dNy_dX[ 12 ][ p ] = 4.5 * ( dL_dX[ n[ 0 ] ] * Wy[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWy_dX[ e[ 0 ] ] );
			dNy_dY[ 12 ][ p ] = 4.5 * ( dL_dY[ n[ 0 ] ] * Wy[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWy_dY[ e[ 0 ] ] );
			dNy_dZ[ 12 ][ p ] = 4.5 * ( dL_dZ[ n[ 0 ] ] * Wy[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWy_dZ[ e[ 0 ] ] );
					
			dNz_dX[ 12 ][ p ] = 4.5 * ( dL_dX[ n[ 0 ] ] * Wz[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWz_dX[ e[ 0 ] ] );
			dNz_dY[ 12 ][ p ] = 4.5 * ( dL_dY[ n[ 0 ] ] * Wz[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWz_dY[ e[ 0 ] ] );
			dNz_dZ[ 12 ][ p ] = 4.5 * ( dL_dZ[ n[ 0 ] ] * Wz[ e[ 0 ] ][ p ] + L[ n[ 0 ] ][ p ] * dWz_dZ[ e[ 0 ] ] );

			// Face: [0 1 3] - Fy
			dNx_dX[ 13 ][ p ] = 4.5 * ( dL_dX[ n[ 1 ] ] * Wx[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWx_dX[ e[ 1 ] ] );
			dNx_dY[ 13 ][ p ] = 4.5 * ( dL_dY[ n[ 1 ] ] * Wx[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWx_dY[ e[ 1 ] ] );
			dNx_dZ[ 13 ][ p ] = 4.5 * ( dL_dZ[ n[ 1 ] ] * Wx[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWx_dZ[ e[ 1 ] ] );
			
			dNy_dX[ 13 ][ p ] = 4.5 * ( dL_dX[ n[ 1 ] ] * Wy[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWy_dX[ e[ 1 ] ] );
			dNy_dY[ 13 ][ p ] = 4.5 * ( dL_dY[ n[ 1 ] ] * Wy[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWy_dY[ e[ 1 ] ] );
			dNy_dZ[ 13 ][ p ] = 4.5 * ( dL_dZ[ n[ 1 ] ] * Wy[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWy_dZ[ e[ 1 ] ] );
			
			dNz_dX[ 13 ][ p ] = 4.5 * ( dL_dX[ n[ 1 ] ] * Wz[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWz_dX[ e[ 1 ] ] );
			dNz_dY[ 13 ][ p ] = 4.5 * ( dL_dY[ n[ 1 ] ] * Wz[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWz_dY[ e[ 1 ] ] );
			dNz_dZ[ 13 ][ p ] = 4.5 * ( dL_dZ[ n[ 1 ] ] * Wz[ e[ 1 ] ][ p ] + L[ n[ 1 ] ][ p ] * dWz_dZ[ e[ 1 ] ] );

			// Face: [1 2 3] - Fy
			dNx_dX[ 14 ][ p ] = 4.5 * ( dL_dX[ n[ 2 ] ] * Wx[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWx_dX[ e[ 2 ] ] );
			dNx_dY[ 14 ][ p ] = 4.5 * ( dL_dY[ n[ 2 ] ] * Wx[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWx_dY[ e[ 2 ] ] );
			dNx_dZ[ 14 ][ p ] = 4.5 * ( dL_dZ[ n[ 2 ] ] * Wx[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWx_dZ[ e[ 2 ] ] );
			
			dNy_dX[ 14 ][ p ] = 4.5 * ( dL_dX[ n[ 2 ] ] * Wy[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWy_dX[ e[ 2 ] ] );
			dNy_dY[ 14 ][ p ] = 4.5 * ( dL_dY[ n[ 2 ] ] * Wy[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWy_dY[ e[ 2 ] ] );
			dNy_dZ[ 14 ][ p ] = 4.5 * ( dL_dZ[ n[ 2 ] ] * Wy[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWy_dZ[ e[ 2 ] ] );
			
			dNz_dX[ 14 ][ p ] = 4.5 * ( dL_dX[ n[ 2 ] ] * Wz[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWz_dX[ e[ 2 ] ] );
			dNz_dY[ 14 ][ p ] = 4.5 * ( dL_dY[ n[ 2 ] ] * Wz[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWz_dY[ e[ 2 ] ] );
			dNz_dZ[ 14 ][ p ] = 4.5 * ( dL_dZ[ n[ 2 ] ] * Wz[ e[ 2 ] ][ p ] + L[ n[ 2 ] ][ p ] * dWz_dZ[ e[ 2 ] ] );

			// Face: [0 2 3] - Fy
			dNx_dX[ 15 ][ p ] = 4.5 * ( dL_dX[ n[ 3 ] ] * Wx[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWx_dX[ e[ 3 ] ] );
			dNx_dY[ 15 ][ p ] = 4.5 * ( dL_dY[ n[ 3 ] ] * Wx[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWx_dY[ e[ 3 ] ] );
			dNx_dZ[ 15 ][ p ] = 4.5 * ( dL_dZ[ n[ 3 ] ] * Wx[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWx_dZ[ e[ 3 ] ] );
			
			dNy_dX[ 15 ][ p ] = 4.5 * ( dL_dX[ n[ 3 ] ] * Wy[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWy_dX[ e[ 3 ] ] );
			dNy_dY[ 15 ][ p ] = 4.5 * ( dL_dY[ n[ 3 ] ] * Wy[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWy_dY[ e[ 3 ] ] );
			dNy_dZ[ 15 ][ p ] = 4.5 * ( dL_dZ[ n[ 3 ] ] * Wy[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWy_dZ[ e[ 3 ] ] );
			
			dNz_dX[ 15 ][ p ] = 4.5 * ( dL_dX[ n[ 3 ] ] * Wz[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWz_dX[ e[ 3 ] ] );
			dNz_dY[ 15 ][ p ] = 4.5 * ( dL_dY[ n[ 3 ] ] * Wz[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWz_dY[ e[ 3 ] ] );
			dNz_dZ[ 15 ][ p ] = 4.5 * ( dL_dZ[ n[ 3 ] ] * Wz[ e[ 3 ] ][ p ] + L[ n[ 3 ] ][ p ] * dWz_dZ[ e[ 3 ] ] );			

			// Face: [0 1 2] - Fz
			dNx_dX[ 16 ][ p ] = 4.5 * ( dL_dX[ n[ 4 ] ] * Wx[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWx_dX[ e[ 4 ] ] );
			dNx_dY[ 16 ][ p ] = 4.5 * ( dL_dY[ n[ 4 ] ] * Wx[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWx_dY[ e[ 4 ] ] );
			dNx_dZ[ 16 ][ p ] = 4.5 * ( dL_dZ[ n[ 4 ] ] * Wx[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWx_dZ[ e[ 4 ] ] );
			
			dNy_dX[ 16 ][ p ] = 4.5 * ( dL_dX[ n[ 4 ] ] * Wy[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWy_dX[ e[ 4 ] ] );
			dNy_dY[ 16 ][ p ] = 4.5 * ( dL_dY[ n[ 4 ] ] * Wy[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWy_dY[ e[ 4 ] ] );
			dNy_dZ[ 16 ][ p ] = 4.5 * ( dL_dZ[ n[ 4 ] ] * Wy[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWy_dZ[ e[ 4 ] ] );
			
			dNz_dX[ 16 ][ p ] = 4.5 * ( dL_dX[ n[ 4 ] ] * Wz[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWz_dX[ e[ 4 ] ] );
			dNz_dY[ 16 ][ p ] = 4.5 * ( dL_dY[ n[ 4 ] ] * Wz[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWz_dY[ e[ 4 ] ] );
			dNz_dZ[ 16 ][ p ] = 4.5 * ( dL_dZ[ n[ 4 ] ] * Wz[ e[ 4 ] ][ p ] + L[ n[ 4 ] ][ p ] * dWz_dZ[ e[ 4 ] ] );
			
			// Face: [0 1 3] - Fz
			dNx_dX[ 17 ][ p ] = 4.5 * ( dL_dX[ n[ 5 ] ] * Wx[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWx_dX[ e[ 5 ] ] );
			dNx_dY[ 17 ][ p ] = 4.5 * ( dL_dY[ n[ 5 ] ] * Wx[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWx_dY[ e[ 5 ] ] );
			dNx_dZ[ 17 ][ p ] = 4.5 * ( dL_dZ[ n[ 5 ] ] * Wx[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWx_dZ[ e[ 5 ] ] );
			
			dNy_dX[ 17 ][ p ] = 4.5 * ( dL_dX[ n[ 5 ] ] * Wy[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWy_dX[ e[ 5 ] ] );
			dNy_dY[ 17 ][ p ] = 4.5 * ( dL_dY[ n[ 5 ] ] * Wy[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWy_dY[ e[ 5 ] ] );
			dNy_dZ[ 17 ][ p ] = 4.5 * ( dL_dZ[ n[ 5 ] ] * Wy[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWy_dZ[ e[ 5 ] ] );
			
			dNz_dX[ 17 ][ p ] = 4.5 * ( dL_dX[ n[ 5 ] ] * Wz[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWz_dX[ e[ 5 ] ] );
			dNz_dY[ 17 ][ p ] = 4.5 * ( dL_dY[ n[ 5 ] ] * Wz[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWz_dY[ e[ 5 ] ] );
			dNz_dZ[ 17 ][ p ] = 4.5 * ( dL_dZ[ n[ 5 ] ] * Wz[ e[ 5 ] ][ p ] + L[ n[ 5 ] ][ p ] * dWz_dZ[ e[ 5 ] ] );

			// Face: [1 2 3] - Fz
			dNx_dX[ 18 ][ p ] = 4.5 * ( dL_dX[ n[ 6 ] ] * Wx[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWx_dX[ e[ 6 ] ] );
			dNx_dY[ 18 ][ p ] = 4.5 * ( dL_dY[ n[ 6 ] ] * Wx[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWx_dY[ e[ 6 ] ] );
			dNx_dZ[ 18 ][ p ] = 4.5 * ( dL_dZ[ n[ 6 ] ] * Wx[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWx_dZ[ e[ 6 ] ] );
			
			dNy_dX[ 18 ][ p ] = 4.5 * ( dL_dX[ n[ 6 ] ] * Wy[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWy_dX[ e[ 6 ] ] );
			dNy_dY[ 18 ][ p ] = 4.5 * ( dL_dY[ n[ 6 ] ] * Wy[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWy_dY[ e[ 6 ] ] );
			dNy_dZ[ 18 ][ p ] = 4.5 * ( dL_dZ[ n[ 6 ] ] * Wy[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWy_dZ[ e[ 6 ] ] );
			
			dNz_dX[ 18 ][ p ] = 4.5 * ( dL_dX[ n[ 6 ] ] * Wz[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWz_dX[ e[ 6 ] ] );
			dNz_dY[ 18 ][ p ] = 4.5 * ( dL_dY[ n[ 6 ] ] * Wz[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWz_dY[ e[ 6 ] ] );
			dNz_dZ[ 18 ][ p ] = 4.5 * ( dL_dZ[ n[ 6 ] ] * Wz[ e[ 6 ] ][ p ] + L[ n[ 6 ] ][ p ] * dWz_dZ[ e[ 6 ] ] );

			// Face: [0 2 3] - Fz
			dNx_dX[ 19 ][ p ] = 4.5 * ( dL_dX[ n[ 7 ] ] * Wx[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWx_dX[ e[ 7 ] ] );
			dNx_dY[ 19 ][ p ] = 4.5 * ( dL_dY[ n[ 7 ] ] * Wx[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWx_dY[ e[ 7 ] ] );
			dNx_dZ[ 19 ][ p ] = 4.5 * ( dL_dZ[ n[ 7 ] ] * Wx[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWx_dZ[ e[ 7 ] ] );
			
			dNy_dX[ 19 ][ p ] = 4.5 * ( dL_dX[ n[ 7 ] ] * Wy[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWy_dX[ e[ 7 ] ] );
			dNy_dY[ 19 ][ p ] = 4.5 * ( dL_dY[ n[ 7 ] ] * Wy[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWy_dY[ e[ 7 ] ] );
			dNy_dZ[ 19 ][ p ] = 4.5 * ( dL_dZ[ n[ 7 ] ] * Wy[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWy_dZ[ e[ 7 ] ] );
			
			dNz_dX[ 19 ][ p ] = 4.5 * ( dL_dX[ n[ 7 ] ] * Wz[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWz_dX[ e[ 7 ] ] );
			dNz_dY[ 19 ][ p ] = 4.5 * ( dL_dY[ n[ 7 ] ] * Wz[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWz_dY[ e[ 7 ] ] );
			dNz_dZ[ 19 ][ p ] = 4.5 * ( dL_dZ[ n[ 7 ] ] * Wz[ e[ 7 ] ][ p ] + L[ n[ 7 ] ][ p ] * dWz_dZ[ e[ 7 ] ] );
        }
	}

    //******************************************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
                                                                  Matrix< std::complex<double> >& EleStiffMatrix, 
                                                                  Vector< std::complex<double> >& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		if( mPotentialsOn ) 
		{
            for( int j=0; j<mNumVerts; j++ ) 
            { 
                if( FixedValue.find( mNodes[ j ]->Id() ) != FixedValue.end() )
                {    
                    std::complex<double> cVoltage = FixedValue[ mNodes[ j ]->Id() ];
		    
                    for( int i=0; i<mNumDofs; i++ ) 
                    {
                        ResidualVector[ i ] -= EleStiffMatrix[ i ][ j + 3*mNumNodes ] * cVoltage;
                    }
                }
            }
		}
    }

    //******************************************************************************************************************************************
    // - Stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );	

		Get_CurlCurl_NiNj_Matrix( StiffMatrix ); 

		if( mStabilizerOn )
		{
		    Add_Stabilization_Matrix( StiffMatrix );
		}
		
		if( mPotentialsOn ) 
		{
		    Add_Potentials_AV_Matrix( StiffMatrix ); 
		}
	}

    //******************************************************************************************************************************************
    // - Cur-Curl + NiNj contribution to the element stiffness matrix 
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Get_CurlCurl_NiNj_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{
        // Element material properties
		double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;

		double pi = 3.14159265358979;
        double eo = 8.8541878176e-12;
	    double mo = pi * 4.0e-7;    
		
		double sg_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
		double sg_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        double ep_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double ep_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
		double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

		std::complex<double> w2cEps;

		if( mComplxFreqOn )
		{
		    std::complex<double> jUnit( 0.0, 1.0 );
		    
		    std::complex<double> sLaplace  = mComplxFreq;
		    std::complex<double> sLaplace2 = sLaplace * sLaplace;
		    
            w2cEps = sLaplace * ( sg_real + jUnit * sg_imag ) - sLaplace2 * ( ep_real + jUnit * ep_imag );		 
		}
		else
		{
            w2cEps = std::complex<double>( w2 * ep_real - wf * sg_imag, w2 * ep_imag + wf * sg_real );
		}

		std::complex<double> cMuInv = 1.0 / std::complex<double>( mu_real, mu_imag );

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 2nd order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_2nd( Nx, Ny, Nz, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ;

		Derivatives_On_Points( dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ, cX, cY, cZ );

		// Edge base matrix
		Matrix< std::complex<double> > EdgeBaseMatrix( mNumBases, mNumBases );

		// Edge base matrix upper diagonal 
		for( int i=0; i<mNumBases; i++ )
		{
			for( int j=i; j<mNumBases; j++ )
			{
				double vIntg_Curli_Curlj = 0.0;
				double vIntg_Ni_Nj       = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Curli_Curlj += W[ gp ] * ( ( dNz_dY[ i ][ gp ] - dNy_dZ[ i ][ gp ] ) * ( dNz_dY[ j ][ gp ] - dNy_dZ[ j ][ gp ] ) + 
					                                 ( dNx_dZ[ i ][ gp ] - dNz_dX[ i ][ gp ] ) * ( dNx_dZ[ j ][ gp ] - dNz_dX[ j ][ gp ] ) +
					                                 ( dNy_dX[ i ][ gp ] - dNx_dY[ i ][ gp ] ) * ( dNy_dX[ j ][ gp ] - dNx_dY[ j ][ gp ] ) );
					
					vIntg_Ni_Nj += W[ gp ] * ( Nx[ i ][ gp ] * Nx[ j ][ gp ] + Ny[ i ][ gp ] * Ny[ j ][ gp ] + Nz[ i ][ gp ] * Nz[ j ][ gp ] );
				}

				EdgeBaseMatrix[ i ][ j ] = ( cMuInv * vIntg_Curli_Curlj ) - ( w2cEps * vIntg_Ni_Nj );
			}
		}

		// Edge base matrix lower diagonal  
        for( int i=0; i<mNumBases; i++ )
        {
            for( int j=i+1; j<mNumBases; j++ )
            {
                EdgeBaseMatrix[ j ][ i ] = EdgeBaseMatrix[ i ][ j ];
            }
        }      

		// Edge base matrix inclusion in stiffnes matrix - [ FYY ]
		int NUM_FYs = 2*mNumEdges + mNumFaces;

		int INI_I = 0;
		int INI_J = 0;

		int END_I = NUM_FYs; 
		int END_J = NUM_FYs;

		int POS_I = mNumNodes + mNumVerts - INI_I;
		int POS_J = mNumNodes + mNumVerts - INI_J;

		for( int i=INI_I; i<END_I; i++ )
        {
            for( int j=INI_J; j<END_J; j++ )
            {
                StiffMatrix[ i + POS_I ][ j + POS_J ] = EdgeBaseMatrix[ i ][ j ];
            }
        }    

		// Edge base matrix inclusion in stiffnes matrix - [ FZZ ]
		int NUM_FZs = mNumFaces;

		INI_I = NUM_FYs;
		INI_J = NUM_FYs;

		END_I = INI_I + NUM_FZs; 
		END_J = INI_J + NUM_FZs;

		POS_I = 2*mNumNodes + mNumVerts + 2*mNumEdges - INI_I;
		POS_J = 2*mNumNodes + mNumVerts + 2*mNumEdges - INI_J;

		for( int i=INI_I; i<END_I; i++ )
        {
            for( int j=INI_J; j<END_J; j++ )
            {
                StiffMatrix[ i + POS_I ][ j + POS_J ] = EdgeBaseMatrix[ i ][ j ];
            }
        }    

		// Edge base matrix inclusion in stiffness matrix - [ FYZ ] - [ FZY ]
		INI_I = 0;
		INI_J = NUM_FYs;

		END_I = NUM_FYs; 
		END_J = NUM_FYs + NUM_FZs;

		POS_I =   mNumNodes + mNumVerts               - INI_I; 
		POS_J = 2*mNumNodes + mNumVerts + 2*mNumEdges - INI_J;  

		for( int i=INI_I; i<END_I; i++ )
        {
            for( int j=INI_J; j<END_J; j++ )
            {
                StiffMatrix[ i + POS_I ][ j + POS_J ] = EdgeBaseMatrix[ i ][ j ];
				StiffMatrix[ j + POS_J ][ i + POS_I ] = EdgeBaseMatrix[ i ][ j ];
            }
        }    
	}

	//******************************************************************************************************************************************
    // - Edge element stabilization term with Lagrange multipliers
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Add_Stabilization_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
        // Stabilization constant
		std::complex<double> cStabCte;

		if( mPotentialsOn )
		{
		    cStabCte = 1.0;
		}
		else 
		{
            double eo      = 8.8541878176e-12;
            double freqc   = (*mProperties)( FREQUENCY );

		    double sg_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
		    double sg_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

            double ep_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		    double ep_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;
		    
		    std::complex<double> cEp( ep_real - ( sg_imag / freqc ), ep_imag + ( sg_real / freqc ) );

		    cStabCte = cEp / eo;
		}

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 2nd order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_2nd( Nx, Ny, Nz, cX, cY, cZ );

		// 1st order Lagrange base derivatives 
		double DN[ 3 ][ 4 ]; 
		
		Calculate_DN( DN );

		// [ Wi * grad( Nj ) ] - [ grad( Ni ) * Wj ] - [ FY ]
		int NUM_FYs = 2*mNumEdges + mNumFaces;

		int INI_I = 0;
		int INI_J = 0;

		int END_I = NUM_FYs; 
		int END_J = mNumVerts;

		int POS_I =   mNumNodes + mNumVerts - INI_I;
		int POS_J = 2*mNumNodes             - INI_J;

		for( int i=INI_I; i<END_I; i++ )
        {
            for( int j=INI_J; j<END_J; j++ )
            {
				double vIntg_Wi_DNj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Wi_DNj += W[ gp ] * ( Nx[ i ][ gp ] * DN[ 0 ][ j ] + Ny[ i ][ gp ] * DN[ 1 ][ j ] + Nz[ i ][ gp ] * DN[ 2 ][ j ] );
				}

				StiffMatrix[ i + POS_I ][ j + POS_J ] += cStabCte * vIntg_Wi_DNj;
				StiffMatrix[ j + POS_J ][ i + POS_I ] += cStabCte * vIntg_Wi_DNj;
            }
        }    

		// [ Wi * grad( Nj ) ] - [ grad( Ni ) * Wj ] - [ FZ ]
		int NUM_FZs = mNumFaces;

		INI_I = NUM_FYs;
		INI_J = 0;

		END_I = INI_I + NUM_FZs; 
		END_J = mNumVerts;

		POS_I = 2*mNumNodes + mNumVerts + 2*mNumEdges - INI_I;
		POS_J = 2*mNumNodes                           - INI_J;

		for( int i=INI_I; i<END_I; i++ )
        {
            for( int j=INI_J; j<END_J; j++ )
            {
				double vIntg_Wi_DNj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Wi_DNj += W[ gp ] * ( Nx[ i ][ gp ] * DN[ 0 ][ j ] + Ny[ i ][ gp ] * DN[ 1 ][ j ] + Nz[ i ][ gp ] * DN[ 2 ][ j ] );
				}

				StiffMatrix[ i + POS_I ][ j + POS_J ] += cStabCte * vIntg_Wi_DNj;
				StiffMatrix[ j + POS_J ][ i + POS_I ] += cStabCte * vIntg_Wi_DNj;
            }
        } 

		// [ DNi * DNj ] 
		for( int i=0; i<mNumVerts; i++ )
        {
            for( int j=0; j<mNumVerts; j++ )
            {
				double vIntg_DNi_DNj = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );

				StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] -= cStabCte * vIntg_DNi_DNj;
            }
        } 
	}

    //******************************************************************************************************************************************
    // - AV-potentials contribution [AV], [VA], [VV] to stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Add_Potentials_AV_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 		
		// Element material properties
		double wf = (*mProperties)( FREQUENCY );

        double w2 = wf * wf;
        double eo = 8.8541878176e-12;
	    
		double sgm_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
	    double sgm_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        double eps_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double eps_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		std::complex<double> nw2cEps;

		if( mComplxFreqOn )
		{
		    std::complex<double> jUnit( 0.0, 1.0 );
		    
		    std::complex<double> sLaplace  = mComplxFreq;
		    std::complex<double> sLaplace2 = sLaplace * sLaplace;
		    
            nw2cEps = sLaplace2 * ( eps_real + jUnit * eps_imag ) - sLaplace * ( sgm_real + jUnit * sgm_imag );       	 
		}
		else
		{
            nw2cEps = std::complex<double>( - w2 * eps_real + wf * sgm_imag, - w2 * eps_imag - wf * sgm_real );
		}

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 2nd order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_2nd( Nx, Ny, Nz, cX, cY, cZ );

		// 1st order Lagrange base derivatives 
		double DN[ 3 ][ 4 ]; 
		
		Calculate_DN( DN );

		// [ Wi * grad( Nj ) ] - [ grad( Ni ) * Wj ] - [ FY ]
		int NUM_FYs = 2*mNumEdges + mNumFaces;

		int INI_I = 0;
		int INI_J = 0;

		int END_I = NUM_FYs; 
		int END_J = mNumVerts;

		int POS_I =   mNumNodes + mNumVerts - INI_I;
		int POS_J = 3*mNumNodes             - INI_J;

		for( int i=INI_I; i<END_I; i++ )
        {
            for( int j=INI_J; j<END_J; j++ )
            {
				double vIntg_Wi_DNj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Wi_DNj += W[ gp ] * ( Nx[ i ][ gp ] * DN[ 0 ][ j ] + Ny[ i ][ gp ] * DN[ 1 ][ j ] + Nz[ i ][ gp ] * DN[ 2 ][ j ] );
				}

				StiffMatrix[ i + POS_I ][ j + POS_J ] += nw2cEps * vIntg_Wi_DNj;
				StiffMatrix[ j + POS_J ][ i + POS_I ] += nw2cEps * vIntg_Wi_DNj;
            }
        }   

		// [ Wi * grad( Nj ) ] - [ grad( Ni ) * Wj ] - [ FZ ]
		int NUM_FZs = mNumFaces;

		INI_I = NUM_FYs;
		INI_J = 0;

		END_I = INI_I + NUM_FZs; 
		END_J = mNumVerts;

		POS_I = 2*mNumNodes + mNumVerts + 2*mNumEdges - INI_I;
		POS_J = 3*mNumNodes                           - INI_J;

		for( int i=INI_I; i<END_I; i++ )
        {
            for( int j=INI_J; j<END_J; j++ )
            {
				double vIntg_Wi_DNj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Wi_DNj += W[ gp ] * ( Nx[ i ][ gp ] * DN[ 0 ][ j ] + Ny[ i ][ gp ] * DN[ 1 ][ j ] + Nz[ i ][ gp ] * DN[ 2 ][ j ] );
				}

				StiffMatrix[ i + POS_I ][ j + POS_J ] += nw2cEps * vIntg_Wi_DNj;
				StiffMatrix[ j + POS_J ][ i + POS_I ] += nw2cEps * vIntg_Wi_DNj;
            }
        }   

		// [ grad( Ni ) * grad( Nj ) ]
		for( int i=0; i<mNumVerts; i++ )
        {
            for( int j=0; j<mNumVerts; j++ ) 
        	{
				StiffMatrix[ i + 3*mNumNodes ][ j + 3*mNumNodes ] = nw2cEps * mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );
        	}
        }
	}

    //******************************************************************************************************************************************
    // - Rotational on nodes
    //******************************************************************************************************************************************
	void VolumeElement_2eg_FullWave::Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes )
	{
		// Nodes natural coordinates
		Vector<double> nX, nY, nZ;

		Get_Nodes_Natural_Coordinates( nX, nY, nZ );

        // Derivatives on nodes
		Matrix<double> dNx_dX, dNx_dY, dNx_dZ; 
		Matrix<double> dNy_dX, dNy_dY, dNy_dZ; 
		Matrix<double> dNz_dX, dNz_dY, dNz_dZ;

        Derivatives_On_Points( dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ, nX, nY, nZ );

        // Ft = [ At | Et ] 
		Vector< std::complex<double> > cFt( mNumBases );

		int NUM_FYs = 2*mNumEdges + mNumFaces;

		if( mPotentialsOn )
		{
		    for( int e=0; e<NUM_FYs; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }

			for( int e=NUM_FYs; e<mNumBases; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAz, *mNodes[ e ] );
		    }
		}
		else
		{
		    for( int e=0; e<NUM_FYs; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }

			for( int e=NUM_FYs; e<mNumBases; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEz, *mNodes[ e ] );
		    }
		}

		// Rotational on nodes
		rot_OnNodes.resize( mNumNodes );
        
		for( int n=0; n<mNumNodes; n++ )
		{
            Vector< std::complex<double> > cRot( 3, std::complex<double>( 0.0, 0.0 ) );
		 
		    for ( int e=0; e<mNumBases; e++ )
		    {
		        cRot[ 0 ] += cFt[ e ] * ( dNz_dY[ e ][ n ] - dNy_dZ[ e ][ n ] );
		        cRot[ 1 ] += cFt[ e ] * ( dNx_dZ[ e ][ n ] - dNz_dX[ e ][ n ] );
		        cRot[ 2 ] += cFt[ e ] * ( dNy_dX[ e ][ n ] - dNx_dY[ e ][ n ] );
		    }	

		    rot_OnNodes[ n ] = cRot;
		}
	}
	 
    //******************************************************************************************************************************************
    // - Rotational on Gauss points
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs )
	{
		// GiD internal Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );		

        // Derivatives on Gauss points
		Matrix<double> dNx_dX, dNx_dY, dNx_dZ; 
		Matrix<double> dNy_dX, dNy_dY, dNy_dZ; 
		Matrix<double> dNz_dX, dNz_dY, dNz_dZ;

        Derivatives_On_Points( dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ, cX, cY, cZ );

        // Ft = [ At | Et ] 
		Vector< std::complex<double> > cFt( mNumBases );

		int NUM_FYs = 2*mNumEdges + mNumFaces;

		if( mPotentialsOn )
		{
		    for( int e=0; e<NUM_FYs; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }

			for( int e=NUM_FYs; e<mNumBases; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAz, *mNodes[ e ] );
		    }
		}
		else
		{
		    for( int e=0; e<NUM_FYs; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }

			for( int e=NUM_FYs; e<mNumBases; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEz, *mNodes[ e ] );
		    }
		}

		// Rotational on Gauss points
		rot_OnGPs.resize( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
            Vector< std::complex<double> > cRot( 3, std::complex<double>( 0.0, 0.0 ) );
		 
		    for ( int e=0; e<mNumBases; e++ )
		    {
		        cRot[ 0 ] += cFt[ e ] * ( dNz_dY[ e ][ gp ] - dNy_dZ[ e ][ gp ] );
		        cRot[ 1 ] += cFt[ e ] * ( dNx_dZ[ e ][ gp ] - dNz_dX[ e ][ gp ] );
		        cRot[ 2 ] += cFt[ e ] * ( dNy_dX[ e ][ gp ] - dNx_dY[ e ][ gp ] );
		    }	

		    rot_OnGPs[ gp ] = cRot;
		}
	}

    //******************************************************************************************************************************************
    // - Electric field on nodes
    //******************************************************************************************************************************************
	void VolumeElement_2eg_FullWave::Calculate_E_field_OnNodes( cVector2Type& Ef_OnNodes )
	{
		// Nodes natural coordinates
		Vector<double> nX, nY, nZ;

		Get_Nodes_Natural_Coordinates( nX, nY, nZ );

        // 2nd order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_2nd( Nx, Ny, Nz, nX, nY, nZ );
	   
        // Ft = [ At | Et ] 
		Vector< std::complex<double> > cFt( mNumBases );

		int NUM_FYs = 2*mNumEdges + mNumFaces;

		if( mPotentialsOn )
		{
		    for( int e=0; e<NUM_FYs; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }

			for( int e=NUM_FYs; e<mNumBases; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAz, *mNodes[ e ] );
		    }
		}
		else
		{
		    for( int e=0; e<NUM_FYs; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }

			for( int e=NUM_FYs; e<mNumBases; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEz, *mNodes[ e ] );
		    }
		}

		// Field on nodes
        cVector2Type cF_OnNodes( mNumNodes );

		for( int n=0; n<mNumNodes; n++ )
		{
		    Vector< std::complex<double> > cF_OnNode( 3, std::complex<double>( 0.0, 0.0 ) );

			for( int e=0; e<mNumBases; e++ )
			{
		        cF_OnNode[ 0 ] += cFt[ e ] * Nx[ e ][ n ];
		        cF_OnNode[ 1 ] += cFt[ e ] * Ny[ e ][ n ];
		        cF_OnNode[ 2 ] += cFt[ e ] * Nz[ e ][ n ];
			}

			cF_OnNodes[ n ] = cF_OnNode;
		}

		// Electric field on nodes
		if( mPotentialsOn )
		{
            double DN[ 3 ][ 4 ]; Calculate_DN( DN );

			std::complex<double> jw;

			if( mComplxFreqOn )
			{
                jw = mComplxFreq;
			}
			else
			{
			    jw = std::complex<double>( 0.0, (*mProperties)( FREQUENCY ) );
			}
		    
			Vector< std::complex<double> > cGradVs( 3, std::complex<double>( 0.0, 0.0 ) ); 

			for( int i=0; i<mNumVerts; i++ )
            {
                cGradVs[ 0 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 0 ][ i ];
                cGradVs[ 1 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 1 ][ i ];
                cGradVs[ 2 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 2 ][ i ];
            }

            // E = jw * ( A + grad( V ) )
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        cF_OnNodes[ n ][ 0 ] = jw * ( cF_OnNodes[ n ][ 0 ] + cGradVs[ 0 ] );
		    	cF_OnNodes[ n ][ 1 ] = jw * ( cF_OnNodes[ n ][ 1 ] + cGradVs[ 1 ] );
		    	cF_OnNodes[ n ][ 2 ] = jw * ( cF_OnNodes[ n ][ 2 ] + cGradVs[ 2 ] );
		    }			
		}

		Ef_OnNodes = cF_OnNodes;
	}

    //******************************************************************************************************************************************
    // - Electric field on Gauss points
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs )
	{
		// GiD internal Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );

        // 2nd order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_2nd( Nx, Ny, Nz, cX, cY, cZ );
	   
        // Ft = [ At | Et ] 
		Vector< std::complex<double> > cFt( mNumBases );

		int NUM_FYs = 2*mNumEdges + mNumFaces;

		if( mPotentialsOn )
		{
		    for( int e=0; e<NUM_FYs; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }

			for( int e=NUM_FYs; e<mNumBases; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAz, *mNodes[ e ] );
		    }
		}
		else
		{
		    for( int e=0; e<NUM_FYs; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }

			for( int e=NUM_FYs; e<mNumBases; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEz, *mNodes[ e ] );
		    }
		}

		// Field on Gauss points
		cVector2Type cF_OnGPs( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
		    Vector< std::complex<double> > cF_OnGP( 3, std::complex<double>( 0.0, 0.0 ) );

			for( int e=0; e<mNumBases; e++ )
			{
		        cF_OnGP[ 0 ] += cFt[ e ] * Nx[ e ][ gp ];
		        cF_OnGP[ 1 ] += cFt[ e ] * Ny[ e ][ gp ];
		        cF_OnGP[ 2 ] += cFt[ e ] * Nz[ e ][ gp ];
			}

			cF_OnGPs[ gp ] = cF_OnGP;
		}

		// Electric field on Gauss points
		if( mPotentialsOn )
		{
            double DN[ 3 ][ 4 ]; Calculate_DN( DN );

			std::complex<double> jw;

			if( mComplxFreqOn )
			{
                jw = mComplxFreq;
			}
			else
			{
			    jw = std::complex<double>( 0.0, (*mProperties)( FREQUENCY ) );
			}
		    
			Vector< std::complex<double> > cGradVs( 3, std::complex<double>( 0.0, 0.0 ) ); 

			for( int i=0; i<mNumVerts; i++ )
            {
                cGradVs[ 0 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 0 ][ i ];
                cGradVs[ 1 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 1 ][ i ];
                cGradVs[ 2 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 2 ][ i ];
            }

            // E = jw * ( A + grad( V ) )
		    for( int gp=0; gp<numGaussPoints; gp++ )
		    {
		        cF_OnGPs[ gp ][ 0 ] = jw * ( cF_OnGPs[ gp ][ 0 ] + cGradVs[ 0 ] );
		    	cF_OnGPs[ gp ][ 1 ] = jw * ( cF_OnGPs[ gp ][ 1 ] + cGradVs[ 1 ] );
		    	cF_OnGPs[ gp ][ 2 ] = jw * ( cF_OnGPs[ gp ][ 2 ] + cGradVs[ 2 ] );
		    }			
		}

		Ef_OnGPs = cF_OnGPs;
	}

    //******************************************************************************************************************************************
    // - Natural coordinates of element nodes
    //******************************************************************************************************************************************
    void VolumeElement_2eg_FullWave::Get_Nodes_Natural_Coordinates( Vector<double>& nX, Vector<double>& nY, Vector<double>& nZ )
	{
		// Element nodes natural coordinates
		nX.resize( mNumNodes );
		nY.resize( mNumNodes );
		nZ.resize( mNumNodes );

		double a = 1.0 / 3.0;
		double b = 2.0 / 3.0;

		// Corners 
		nX[  0 ] = 0.0; nY[  0 ] = 0.0; nZ[  0 ] = 0.0;
		nX[  1 ] = 1.0; nY[  1 ] = 0.0; nZ[  1 ] = 0.0;
		nX[  2 ] = 0.0; nY[  2 ] = 1.0; nZ[  2 ] = 0.0;
		nX[  3 ] = 0.0; nY[  3 ] = 0.0; nZ[  3 ] = 1.0;
		
		// Edges 
		nX[  4 ] = a  ; nY[  4 ] = 0.0; nZ[  4 ] = 0.0; 
		nX[  5 ] = b  ; nY[  5 ] = 0.0; nZ[  5 ] = 0.0; 
			 					 					 		      
		nX[  6 ] = b  ; nY[  6 ] = a  ; nZ[  6 ] = 0.0; 
		nX[  7 ] = a  ; nY[  7 ] = b  ; nZ[  7 ] = 0.0;  
			 					 					 		      
		nX[  8 ] = 0.0; nY[  8 ] = b  ; nZ[  8 ] = 0.0; 
		nX[  9 ] = 0.0; nY[  9 ] = a  ; nZ[  9 ] = 0.0;  

		nX[ 10 ] = 0.0; nY[ 10 ] = 0.0; nZ[ 10 ] = a  ;
		nX[ 11 ] = 0.0; nY[ 11 ] = 0.0; nZ[ 11 ] = b  ; 

		nX[ 12 ] = b  ; nY[ 12 ] = 0.0; nZ[ 12 ] = a  ;
		nX[ 13 ] = a  ; nY[ 13 ] = 0.0; nZ[ 13 ] = b  ; 

		nX[ 14 ] = 0.0; nY[ 14 ] = b  ; nZ[ 14 ] = a  ;
		nX[ 15 ] = 0.0; nY[ 15 ] = a  ; nZ[ 15 ] = b  ;  

		// Faces
		nX[ 16 ] = a  ; nY[ 16 ] = a  ; nZ[ 16 ] = 0.0;
		nX[ 17 ] = a  ; nY[ 17 ] = 0.0; nZ[ 17 ] = a  ; 
		nX[ 18 ] = a  ; nY[ 18 ] = a  ; nZ[ 18 ] = a  ;
		nX[ 19 ] = 0.0; nY[ 19 ] = a  ; nZ[ 19 ] = a  ; 
	}

    //******************************************************************************************************************************************
    // - Inner GiD Gauss points 
    //******************************************************************************************************************************************
    int VolumeElement_2eg_FullWave::GetInnerGiDGaussPoints( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, int numResultsOnGPs )
	{
		if( numResultsOnGPs <= 1 ) 
		{
			GaussPoints3D_InnerGiD_0rd( cX, cY, cZ );
		}
		else    
		{
			GaussPoints3D_InnerGiD_1st( cX, cY, cZ );
		}
	
        return cX.size();
	}
} 