
#include "JSource_1eg_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
	//*****************************************************************************************************************
    //*******************************************  Cartesian mode  ****************************************************
    //*****************************************************************************************************************  
    // 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = Module of the cartesian current X component.
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Phase  of the cartesian current X component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Module of the cartesian current Y component.
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = Phase  of the cartesian current Y component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Module of the cartesian current Z component.
    // SINUSOIDAL_SURFACE_CURRENT[ 5 ] = Phase  of the cartesian current Z component.
    //  
    //*****************************************************************************************************************
    //*******************************************  Axysimmetric mode  *************************************************
    //*****************************************************************************************************************  
    //
    // COMPLEX_IBC[ 0 ] = Axisymmetric current module.
    // COMPLEX_IBC[ 1 ] = Axisymmetric current phase.
    //
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of curvature.    
    //
    // COMPLEX_NEUMANN_FLOW[ 3 ] = X component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Y component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 5 ] = Z component of the curvature axis vector.
    // 
    //*****************************************************************************************************************
    //*******************************************  Helicoidal mode  ***************************************************
    //*****************************************************************************************************************    
    // 
    // COMPLEX_IBC_2o[ 0 ] = Toroid major radius.
	// 
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of curvature.    
    //
    // COMPLEX_NEUMANN_FLOW[ 3 ] = X component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Y component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 5 ] = Z component of the curvature axis vector.
	// 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = Module of the poloidal current component.
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Phase  of the poloidal current component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Module of the toroidal current component.
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = Phase  of the toroidal current component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 5 ] = 0 
	//                                  
    //*****************************************************************************************************************
    //*******************************************  Plasma mode  *******************************************************
    //*****************************************************************************************************************  
	//  
    // COMPLEX_IBC_2o[ 0 ] = Toroid major radius.
    // 
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of curvature.    
    //
    // COMPLEX_NEUMANN_FLOW[ 3 ] = X component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Y component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 5 ] = Z component of the curvature axis vector.
	// 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = Module of the poloidal current component.
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Phase  of the poloidal current component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Module of the toroidal current component.
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = Phase  of the toroidal current component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Plasma mode index. 
	// SINUSOIDAL_SURFACE_CURRENT[ 5 ] = 1 
	// 
    //*****************************************************************************************************************
	
	//*****************************************************************************************************************
    //* - Global index of the nodes
    //*****************************************************************************************************************	
    void JSource_1eg_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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

    //*****************************************************************************************************************
    //* - Element volume
    //*****************************************************************************************************************
    double JSource_1eg_FullWave::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
        
        return fabs( Det / 6.0 );
    }

    //*****************************************************************************************************************
    // - 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //*****************************************************************************************************************
    void JSource_1eg_FullWave::Calculate_DN( double DN[ 3 ][ 4 ] )
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

    //*****************************************************************************************************************
    // - Element edge lengths with sign
    //*****************************************************************************************************************
    void JSource_1eg_FullWave::Calculate_Edge_Lengths( Vector<double>& eLength )
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

    //*****************************************************************************************************************
    // - 1st order Whitney base for 3D tetrahedrons
    //*****************************************************************************************************************
	void JSource_1eg_FullWave::Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz,
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

	//*****************************************************************************************************************
    //* - Cartesian sources
    //*****************************************************************************************************************
    void JSource_1eg_FullWave::Get_CartesianXYZ_Sources( Vector< std::complex<double> >& cJx, 
		                                                 Vector< std::complex<double> >& cJy, 
		                                                 Vector< std::complex<double> >& cJz ) 
    {
		// Cartesian source
        Vector<double> vJ = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT );

		if( ( vJ[ 0 ] != 0.0 ) || ( vJ[ 2 ] != 0.0 ) || ( vJ[ 4 ] != 0.0 ) ) 
		{
			std::complex<double> Element_cJx( vJ[ 0 ] * cos( vJ[ 1 ] ), vJ[ 0 ] * sin( vJ[ 1 ] ) );
			std::complex<double> Element_cJy( vJ[ 2 ] * cos( vJ[ 3 ] ), vJ[ 2 ] * sin( vJ[ 3 ] ) );
			std::complex<double> Element_cJz( vJ[ 4 ] * cos( vJ[ 5 ] ), vJ[ 4 ] * sin( vJ[ 5 ] ) );

			for( int i=0; i<mNumVerts; i++ )
			{
				cJx[ i ] += Element_cJx;
				cJy[ i ] += Element_cJy;
				cJz[ i ] += Element_cJz;
			}
		}

		if( mNodal_J.size() > 0 )
		{
			for( int i=0; i<mNumVerts; i++ )
			{
				cJx[ i ] += mNodal_J[ i ][ 0 ];
				cJy[ i ] += mNodal_J[ i ][ 1 ];
				cJz[ i ] += mNodal_J[ i ][ 2 ];
			}		
		}
	}

	//*****************************************************************************************************************
    //* - Axisymmetric sources
    //*****************************************************************************************************************
    void JSource_1eg_FullWave::Get_Axisymmetric_Sources( Vector< std::complex<double> >& cJx, 
		                                                 Vector< std::complex<double> >& cJy, 
		                                                 Vector< std::complex<double> >& cJz )  
    {
        // Axisymmetric source
		Vector<double> vJa = (*mProperties)( COMPLEX_IBC );

		if( vJa[ 0 ] == 0.0 )
		{
		    return;
		}

		// Complex axisymmetric current sources
		std::complex<double> cJa( vJa[ 0 ] * cos( vJa[ 1 ] ), vJa[ 0 ] * sin( vJa[ 1 ] ) );

        // Origin of coordinates
		double ox = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ];
		double oy = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ];
		double oz = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 2 ];

		// Rotation axis
		double ax = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ];
		double ay = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];
		double az = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ];

		// Transform sources to cartesian and sent to vertice nodes
		for( int i=0; i<mNumVerts; i++ )
		{
			// Position vector 
			double px = mNodes[ i ]->X() - ox;
			double py = mNodes[ i ]->Y() - oy;
			double pz = mNodes[ i ]->Z() - oz;
               
			// r = a x p 
			double rx = ay*pz - az*py;
			double ry = az*px - ax*pz;
			double rz = ax*py - ay*px;

			double norm = sqrt( rx*rx + ry*ry + rz*rz );

			if( norm > 0.0 )
			{
				rx = rx / norm;
				ry = ry / norm;
				rz = rz / norm;
			}
			else
			{
				rx = 0.0;
				ry = 0.0;
				rz = 0.0;
			}

			cJx[ i ] += cJa * rx;
			cJy[ i ] += cJa * ry;
			cJz[ i ] += cJa * rz;
		}
	}

	//*****************************************************************************************************************
    //* - Plasma modes J sources
    //*****************************************************************************************************************
    void JSource_1eg_FullWave::Get_PlasmaModeJs_Sources( Vector< std::complex<double> >& cJx, 
		                                                 Vector< std::complex<double> >& cJy, 
		                                                 Vector< std::complex<double> >& cJz ) 
    {
		// Torus major radius
		double Ro = (*mProperties)( COMPLEX_IBC_2o )[ 0 ];
		
		if( Ro <= 0.0 ) 
		{
		    return;
		}

		// Plasma mode J current sources 
        Vector<double> vJp = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT );

		if( ( vJp[ 0 ] == 0.0 ) && ( vJp[ 2 ] == 0.0 ) ) 
		{
		    return;
		}

		// Complex poloidal and toroidal current 
        std::complex<double> cJpol( vJp[ 0 ] * cos( vJp[ 1 ] ), vJp[ 0 ] * sin( vJp[ 1 ] ) );
		std::complex<double> cJtor( vJp[ 2 ] * cos( vJp[ 3 ] ), vJp[ 2 ] * sin( vJp[ 3 ] ) );

        // Origin of coordinates
		Vector<double> OO_A( 3 );

		OO_A[ 0 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ];
		OO_A[ 1 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ];
		OO_A[ 2 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 2 ];

	    // Rotation axis
		Vector<double> nAxi( 3 );

		nAxi[ 0 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ];
		nAxi[ 1 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];
		nAxi[ 2 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ];

	    nAxi *= ( 1.0 / sqrt( nAxi[ 0 ]*nAxi[ 0 ] + nAxi[ 1 ]*nAxi[ 1 ] + nAxi[ 2 ]*nAxi[ 2 ] ) );

		// Phi reference axis
		Vector<double> nRef( 3 );

		if( sqrt( nAxi[ 0 ]*nAxi[ 0 ] + nAxi[ 1 ]*nAxi[ 1 ] ) > 0.0 )
        {
            nRef[ 0 ] = nAxi[ 1 ];
		    nRef[ 1 ] =-nAxi[ 0 ];
		    nRef[ 2 ] = 0.0;
        }
        else
        {
            nRef[ 0 ] = nAxi[ 2 ];
            nRef[ 1 ] = 0.0;
            nRef[ 2 ] = 0.0;
        }

		nRef *= ( 1.0 / sqrt( nRef[ 0 ]*nRef[ 0 ] + nRef[ 1 ]*nRef[ 1 ] + nRef[ 2 ]*nRef[ 2 ] ) );

		// Phi reference axis ( ortogonal )
		Vector<double> nOef( 3 );

		nOef[ 0 ] = nAxi[ 1 ] * nRef[ 2 ] - nAxi[ 2 ] * nRef[ 1 ];
		nOef[ 1 ] = nAxi[ 2 ] * nRef[ 0 ] - nAxi[ 0 ] * nRef[ 2 ];
		nOef[ 2 ] = nAxi[ 0 ] * nRef[ 1 ] - nAxi[ 1 ] * nRef[ 0 ];

		nOef *= ( 1.0 / sqrt( nOef[ 0 ]*nOef[ 0 ] + nOef[ 1 ]*nOef[ 1 ] + nOef[ 2 ]*nOef[ 2 ] ) );

		// Transform to cartesian coordinates and sent to vertice nodes
		for( int i=0; i<mNumVerts; i++ )
		{
			// Position vector 
			Vector<double> OO_P( 3 );

			OO_P[ 0 ] = mNodes[ i ]->X() - OO_A[ 0 ];
			OO_P[ 1 ] = mNodes[ i ]->Y() - OO_A[ 1 ];
			OO_P[ 2 ] = mNodes[ i ]->Z() - OO_A[ 2 ];
               
			// Toroidal vector ( nTor = nAxi x OO_P )
			Vector<double> nTor( 3 );

			nTor[ 0 ] = nAxi[ 1 ] * OO_P[ 2 ] - nAxi[ 2 ] * OO_P[ 1 ];
			nTor[ 1 ] = nAxi[ 2 ] * OO_P[ 0 ] - nAxi[ 0 ] * OO_P[ 2 ];
			nTor[ 2 ] = nAxi[ 0 ] * OO_P[ 1 ] - nAxi[ 1 ] * OO_P[ 0 ];

			nTor *= ( 1.0 / sqrt( nTor[ 0 ]*nTor[ 0 ] + nTor[ 1 ]*nTor[ 1 ] + nTor[ 2 ]*nTor[ 2 ] ) );

			// Radial vector ( nRad = nTor x nAxi )
			Vector<double> nRad( 3 );

			nRad[ 0 ] = nTor[ 1 ] * nAxi[ 2 ] - nTor[ 2 ] * nAxi[ 1 ];
			nRad[ 1 ] = nTor[ 2 ] * nAxi[ 0 ] - nTor[ 0 ] * nAxi[ 2 ];
			nRad[ 2 ] = nTor[ 0 ] * nAxi[ 1 ] - nTor[ 1 ] * nAxi[ 0 ];
			 
			nRad *= ( 1.0 / sqrt( nRad[ 0 ]*nRad[ 0 ] + nRad[ 1 ]*nRad[ 1 ] + nRad[ 2 ]*nRad[ 2 ] ) );

			// Inner radial vector ( InnR = OO_P - Ro * nRad )
			Vector<double> InnR( 3 );

			InnR[ 0 ] = OO_P[ 0 ] - Ro * nRad[ 0 ];
			InnR[ 1 ] = OO_P[ 1 ] - Ro * nRad[ 1 ];
			InnR[ 2 ] = OO_P[ 2 ] - Ro * nRad[ 2 ];

			// Poloidal vector ( nPol = nInR x nTor )
			Vector<double> nPol( 3 );

			nPol[ 0 ] = InnR[ 1 ] * nTor[ 2 ] - InnR[ 2 ] * nTor[ 1 ];
			nPol[ 1 ] = InnR[ 2 ] * nTor[ 0 ] - InnR[ 0 ] * nTor[ 2 ];
			nPol[ 2 ] = InnR[ 0 ] * nTor[ 1 ] - InnR[ 1 ] * nTor[ 0 ];
			 
			nPol *= ( 1.0 / sqrt( nPol[ 0 ]*nPol[ 0 ] + nPol[ 1 ]*nPol[ 1 ] + nPol[ 2 ]*nPol[ 2 ] ) );

			// Poloidal ( Theta ) and toroidal ( Phi ) currents
			std::complex<double> J_Theta = cJpol;
			std::complex<double> J_Phi   = cJtor;

		    // Add plasma mode M 
			if( (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 5 ] == 1 )
			{
			    // Plasma mode index
		        int M = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 4 ];

				// Minor radius
			    double rt = sqrt( InnR[ 0 ]*InnR[ 0 ] + InnR[ 1 ]*InnR[ 1 ] + InnR[ 2 ]*InnR[ 2 ] );

			    // Theta angle
				double Y_Theta = InnR[ 0 ] * nAxi[ 0 ] + InnR[ 1 ] * nAxi[ 1 ] + InnR[ 2 ] * nAxi[ 2 ];
				double X_Theta = InnR[ 0 ] * nRad[ 0 ] + InnR[ 1 ] * nRad[ 1 ] + InnR[ 2 ] * nRad[ 2 ]; 
				
				double Theta = atan2( Y_Theta, X_Theta ); 

				// Phi angle
				double XY_Prj = OO_P[ 0 ] * nRad[ 0 ] + OO_P[ 1 ] * nRad[ 1 ] + OO_P[ 2 ] * nRad[ 2 ];

				double Y_Phi = XY_Prj * ( nRad[ 0 ] * nOef[ 0 ] + nRad[ 1 ] * nOef[ 1 ] + nRad[ 2 ] * nOef[ 2 ] );
				double X_Phi = XY_Prj * ( nRad[ 0 ] * nRef[ 0 ] + nRad[ 1 ] * nRef[ 1 ] + nRad[ 2 ] * nRef[ 2 ] ); 
				
				double Phi = atan2( Y_Phi, X_Phi ); 
			    
				// Plasma mode
			    J_Theta *= sin( M * Theta + Phi );
			    J_Phi   *= - M * ( ( Ro / rt ) + cos( Theta ) ) * sin( M * Theta + Phi ) - sin( Theta ) * cos( M * Theta + Phi );
			}

			cJx[ i ] += J_Theta * nPol[ 0 ] + J_Phi * nTor[ 0 ];
			cJy[ i ] += J_Theta * nPol[ 1 ] + J_Phi * nTor[ 1 ];
			cJz[ i ] += J_Theta * nPol[ 2 ] + J_Phi * nTor[ 2 ];
		}
	}

    //*****************************************************************************************************************
    //* - Imposed J currents on element nodes
    //*****************************************************************************************************************
    void JSource_1eg_FullWave::Get_J_Imposed_OnNodes( std::map< unsigned int, Vector<double> >& rJimp, 
		                                              std::map< unsigned int, Vector<double> >& iJimp )
	{
		// J sources on vertice nodes
		Vector< std::complex<double> > cJx( mNumVerts, std::complex<double>( 0.0, 0.0 ) );
		Vector< std::complex<double> > cJy( mNumVerts, std::complex<double>( 0.0, 0.0 ) );
		Vector< std::complex<double> > cJz( mNumVerts, std::complex<double>( 0.0, 0.0 ) );

		// Get cartesian and axisymmetric sources 
		if( (*mProperties)( COMPLEX_IBC_2o )[ 0 ] == 0 ) 
		{
		    Get_CartesianXYZ_Sources( cJx, cJy, cJz );
			Get_Axisymmetric_Sources( cJx, cJy, cJz );
		}
		// Get plasma mode currents sources
		else
		{
		    Get_PlasmaModeJs_Sources( cJx, cJy, cJz );
		}

		// Sent sources to map 
		for( int i=0; i<mNumVerts; i++ )
		{
			Vector<double> real_J_vector( 3 );
			
		    real_J_vector[ 0 ] = std::real( cJx[ i ] );
            real_J_vector[ 1 ] = std::real( cJy[ i ] );
            real_J_vector[ 2 ] = std::real( cJz[ i ] );

			Vector<double> imag_J_vector( 3 );

		    imag_J_vector[ 0 ] = std::imag( cJx[ i ] );
            imag_J_vector[ 1 ] = std::imag( cJy[ i ] );
            imag_J_vector[ 2 ] = std::imag( cJz[ i ] );

			int NodeId = mNodes[ i ]->Id();

			rJimp[ NodeId ] = real_J_vector;
			iJimp[ NodeId ] = imag_J_vector;
		}
	}

    //*****************************************************************************************************************
    //* - Residual vector
    //*****************************************************************************************************************
    void JSource_1eg_FullWave::GetResidualVector( Vector< std::complex<double> >& ResidualVector ) 
    {
		// Initializing residual vector
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		// J sources on vertice nodes
		Vector< std::complex<double> > cJx( mNumVerts, std::complex<double>( 0.0, 0.0 ) );
		Vector< std::complex<double> > cJy( mNumVerts, std::complex<double>( 0.0, 0.0 ) );
		Vector< std::complex<double> > cJz( mNumVerts, std::complex<double>( 0.0, 0.0 ) );

		// Get cartesian and axisymmetric sources 
		if( (*mProperties)( COMPLEX_IBC_2o )[ 0 ] == 0 )
		{
		    Get_CartesianXYZ_Sources( cJx, cJy, cJz );
			Get_Axisymmetric_Sources( cJx, cJy, cJz );
		}
		// Get plasma mode currents sources
		else
		{
		    Get_PlasmaModeJs_Sources( cJx, cJy, cJz );
		}

		// If potentials are off then the source term is equal to jwJ
		if( mPotentialsOn == false )
		{
		    std::complex<double> jw( 0.0, (*mProperties)( FREQUENCY ) );
		    
			cJx *= jw;
			cJy *= jw; 
			cJz *= jw; 
		}

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

		// 1st order Lagrangian base functions
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // 1st order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_1st( Nx, Ny, Nz, cX, cY, cZ );

		// Edge element source vector
		for( int i=0; i<mNumEdges; i++ )
		{
		    std::complex<double> vIntg_Ni_cJa( 0.0, 0.0 );
		    
		    for( int gp=0; gp<numGaussPoints; gp++ )
		    {
		    	vIntg_Ni_cJa += W[ gp ] * Nx[ i ][ gp ] * ( N[ 0 ][ gp ] * cJx[ 0 ] + N[ 1 ][ gp ] * cJx[ 1 ] + N[ 2 ][ gp ] * cJx[ 2 ] + N[ 3 ][ gp ] * cJx[ 3 ] );
				vIntg_Ni_cJa += W[ gp ] * Ny[ i ][ gp ] * ( N[ 0 ][ gp ] * cJy[ 0 ] + N[ 1 ][ gp ] * cJy[ 1 ] + N[ 2 ][ gp ] * cJy[ 2 ] + N[ 3 ][ gp ] * cJy[ 3 ] );
				vIntg_Ni_cJa += W[ gp ] * Nz[ i ][ gp ] * ( N[ 0 ][ gp ] * cJz[ 0 ] + N[ 1 ][ gp ] * cJz[ 1 ] + N[ 2 ][ gp ] * cJz[ 2 ] + N[ 3 ][ gp ] * cJz[ 3 ] );
		    }
			
            ResidualVector[ i + mNumNodes + 4 ] += vIntg_Ni_cJa;
		}
    }
} 