
#include "GenericRobin_1eg_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //**********************************************************************************************************************
    //*******************************************  Quasi-Static Mode  ******************************************************
    //**********************************************************************************************************************    
    // 
    // COMPLEX_IBC[ 0 ] = Real part of the generic Robin coefficient
    // COMPLEX_IBC[ 1 ] = Imag part of the generic Robin coefficient
    // 
    // COMPLEX_IBC_2o[ 0 ] = Flux module
    // COMPLEX_IBC_2o[ 1 ] = Flux phase
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 0
    //  
    //**********************************************************************************************************************
    //*******************************************  Plane Wave Mode  ********************************************************
    //**********************************************************************************************************************    
    //
    // COMPLEX_IBC[ 0 ] = Real part of the electric conductivity in S/m
    // COMPLEX_IBC[ 1 ] = Imag part of the electric conductivity in S/m
    // 
    // COMPLEX_IBC_2o[ 0 ] = Real part of the relative electric permittivity
    // COMPLEX_IBC_2o[ 1 ] = Imag part of the relative electric permittivity
    //  
    // COMPLEX_NEUMANN_FLOW[ 0 ] = Real part of the relative magnetic permeability
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Imag part of the relative magnetic permeability
    // 
    // COMPLEX_NEUMANN_FLOW[ 3 ] = Module of the E field parallel to the polarization vector
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Phase  of the E field parallel to the polarization vector
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 1 in single   plane wave  mode  
    //                           = 2 in multiple plane waves mode 
    //                            
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = X component of the polarization vector
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Y component of the polarization vector
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Z component of the polarization vector
	//
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = X component of the wave vector
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Y component of the wave vector
    // SINUSOIDAL_SURFACE_CURRENT[ 5 ] = Z component of the wave vector
    // 
    //**********************************************************************************************************************
    //*******************************************  Gaussian Beam Mode  *****************************************************
    //**********************************************************************************************************************    
    // 
    // COMPLEX_IBC[ 0 ] = Module of the E field perpendicular to the polarization and wave vectors at the waist center
    // COMPLEX_IBC[ 1 ] = Phase  of the E field perpendicular to the polarization and wave vectors at the waist center
    // 
    // COMPLEX_IBC_2o[ 0 ] = Gaussian beam waist radius on the direction parallel to the polarization vector
    // COMPLEX_IBC_2o[ 1 ] = Gaussian beam waist radius on the direction perpendicular to the polarization and wave vectors
    //   
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of the Gaussian beam waist
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of the Gaussian beam waist
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of the Gaussian beam waist 
    // 
    // COMPLEX_NEUMANN_FLOW[ 3 ] = Module of the E field parallel to the polarization vector at the waist center
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Phase  of the E field parallel to the polarization vector at the waist center
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 3
	// 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = X component of the polarization vector 
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Y component of the polarization vector 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Z component of the polarization vector
	//
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = X component of the wave vector 
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Y component of the wave vector 
    // SINUSOIDAL_SURFACE_CURRENT[ 5 ] = Z component of the wave vector 
    // 
    //**********************************************************************************************************************
    
    //**********************************************************************************************************************
    //* - Global IDs of the DOFs
    //**********************************************************************************************************************	
    void GenericRobin_1eg_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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

            for( int i=0; i<3; i++ ) 
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

    //**********************************************************************************************************************
    // - Area of the element
    //**********************************************************************************************************************
    double GenericRobin_1eg_FullWave::Calculate_Area()
    {
		double na[ 3 ], v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ 2 ]->X() - mNodes[ 0 ]->X();
        v2[ 1 ] = mNodes[ 2 ]->Y() - mNodes[ 0 ]->Y();
        v2[ 2 ] = mNodes[ 2 ]->Z() - mNodes[ 0 ]->Z();
		    
        v1[ 0 ] = mNodes[ 1 ]->X() - mNodes[ 0 ]->X();
        v1[ 1 ] = mNodes[ 1 ]->Y() - mNodes[ 0 ]->Y();
        v1[ 2 ] = mNodes[ 1 ]->Z() - mNodes[ 0 ]->Z();

        // Area = 0.5 * ||v2 x v1||
        na[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        na[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        na[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

		return ( 0.5 * sqrt( na[ 0 ]*na[ 0 ] + na[ 1 ]*na[ 1 ] + na[ 2 ]*na[ 2 ] ) );
    }

    //**********************************************************************************************************************
    //* - Exterior normal 
    //**********************************************************************************************************************	
    void GenericRobin_1eg_FullWave::Calculate_Exterior_Normal( Vector<double>& Next )
    {
		double v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ 2 ]->X() - mNodes[ 0 ]->X();
        v2[ 1 ] = mNodes[ 2 ]->Y() - mNodes[ 0 ]->Y();
        v2[ 2 ] = mNodes[ 2 ]->Z() - mNodes[ 0 ]->Z();
		    
        v1[ 0 ] = mNodes[ 1 ]->X() - mNodes[ 0 ]->X();
        v1[ 1 ] = mNodes[ 1 ]->Y() - mNodes[ 0 ]->Y();
        v1[ 2 ] = mNodes[ 1 ]->Z() - mNodes[ 0 ]->Z();

        Next.resize( 3 );

        Next[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        Next[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        Next[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

        Next *= ( 1.0 / sqrt( Next[ 0 ]*Next[ 0 ] + Next[ 1 ]*Next[ 1 ] + Next[ 2 ]*Next[ 2 ] ) );
    }

    //**********************************************************************************************************************
    //* - jk / mu = j * w * sqrt ( ep / mu )
    //**********************************************************************************************************************	
    std::complex<double> GenericRobin_1eg_FullWave::Calculate_jK_div_Mu()
    {
        // Gaussian beam 
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;
            
            return std::complex<double>( 0.0, mFreq * sqrt( eo / mo ) );
        }
        // Plane wave 
        else
        {
            if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] == 0.0 ) 
            {
                return std::complex<double>( 0.0, 0.0 );
            }            
            
            double wf = mFreq;
            double w2 = mFreq * mFreq;

            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;            
            
            double sg_real = (*mProperties)( COMPLEX_IBC )[ 0 ];
            double sg_imag = (*mProperties)( COMPLEX_IBC )[ 1 ];

            double ep_real = (*mProperties)( COMPLEX_IBC_2o )[ 0 ] * eo; 
            double ep_imag = (*mProperties)( COMPLEX_IBC_2o )[ 1 ] * eo; 

            double mu_real = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] * mo; 
            double mu_imag = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ] * mo; 

            std::complex<double> cUnit( 0.0                        , 1.0                         );
            std::complex<double> cMu  ( mu_real                    , mu_imag                     );
            std::complex<double> cEpw2( ep_real * w2 - sg_imag * wf, ep_imag * w2 + sg_real * wf );

            return ( cUnit * sqrt( cEpw2 / cMu ) );
        }
    }

    //**********************************************************************************************************************
    //* - jk = j * w * sqrt ( ep * mu )
    //**********************************************************************************************************************	
    std::complex<double> GenericRobin_1eg_FullWave::Calculate_jK()
    {
        // Gaussian beam 
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;
            
            return std::complex<double>( 0.0, mFreq * sqrt( eo * mo ) );
        }
        // Plane wave 
        else
        {
            double wf = mFreq;
            double w2 = mFreq * mFreq;

            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;            
            
            double sg_real = (*mProperties)( COMPLEX_IBC )[ 0 ];
            double sg_imag = (*mProperties)( COMPLEX_IBC )[ 1 ];

            double ep_real = (*mProperties)( COMPLEX_IBC_2o )[ 0 ] * eo; 
            double ep_imag = (*mProperties)( COMPLEX_IBC_2o )[ 1 ] * eo; 

            double mu_real = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] * mo; 
            double mu_imag = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ] * mo; 

            std::complex<double> cUnit( 0.0                        , 1.0                         );
            std::complex<double> cMu  ( mu_real                    , mu_imag                     );
            std::complex<double> cEpw2( ep_real * w2 - sg_imag * wf, ep_imag * w2 + sg_real * wf );
            
            return ( cUnit * sqrt( cEpw2 * cMu ) );
        }
    }

    //*******************************************************************************************************
    // - Get X, Y axis of the co-planar coordinate system 
    //*******************************************************************************************************
    void GenericRobin_1eg_FullWave::Get_CoPlanar_XY_Axis( Vector<double>& X_axis,  Vector<double>& Y_axis )
	{
        // Z axis of the co-planar coordinate system 
        Vector<double> Z_axis; 
        
        Calculate_Exterior_Normal( Z_axis );

        // X axis of the co-planar coordinate system 
        X_axis.resize( 3 );

        X_axis[ 0 ] = X( 2 ) - X( 1 ); 
		X_axis[ 1 ] = Y( 2 ) - Y( 1 ); 
		X_axis[ 2 ] = Z( 2 ) - Z( 1 );

        X_axis *= ( 1.0 / sqrt( X_axis[ 0 ]*X_axis[ 0 ] + X_axis[ 1 ]*X_axis[ 1 ] + X_axis[ 2 ]*X_axis[ 2 ] ) ); 

        // Y axis of the co-planar coordinate system ( Y_axis = X_axis x Z_axis )
        Y_axis.resize( 3 );

        Y_axis[ 0 ] = X_axis[ 1 ] * Z_axis[ 2 ] - X_axis[ 2 ] * Z_axis[ 1 ];
        Y_axis[ 1 ] = X_axis[ 2 ] * Z_axis[ 0 ] - X_axis[ 0 ] * Z_axis[ 2 ];
        Y_axis[ 2 ] = X_axis[ 0 ] * Z_axis[ 1 ] - X_axis[ 1 ] * Z_axis[ 0 ];

        Y_axis *= ( 1.0 / sqrt( Y_axis[ 0 ]*Y_axis[ 0 ] + Y_axis[ 1 ]*Y_axis[ 1 ] + Y_axis[ 2 ]*Y_axis[ 2 ] ) ); 
	}

    //*******************************************************************************************************
    // - Change coordinate system from ( X, Y, Z ) to ( pX, pY ) 
    //*******************************************************************************************************
    void GenericRobin_1eg_FullWave::Get_Vertices_Surface_Coordinates( Vector<double>& pX,  Vector<double>& pY )
    {
        // XY axis of the co-planar coordinate system 
        Vector<double> X_axis, Y_axis;

		Get_CoPlanar_XY_Axis( X_axis, Y_axis );

        // Position vectors of vertice nodes
        Vector<double> R_2( 3 ); 
        
        R_2[ 0 ] = X( 2 ) - X( 1 ); 
        R_2[ 1 ] = Y( 2 ) - Y( 1 ); 
        R_2[ 2 ] = Z( 2 ) - Z( 1 );
        
        Vector<double> R_3( 3 ); 

        R_3[ 0 ] = X( 3 ) - X( 1 ); 
        R_3[ 1 ] = Y( 3 ) - Y( 1 ); 
        R_3[ 2 ] = Z( 3 ) - Z( 1 );

        // Co-planar coordinates of vertice nodes
        pX.resize( mNumVerts ); 
        pY.resize( mNumVerts );

        pX[ 0 ] = 0.0; 
        pY[ 0 ] = 0.0;

        pX[ 1 ] = R_2[ 0 ] * X_axis[ 0 ] + R_2[ 1 ] * X_axis[ 1 ] + R_2[ 2 ] * X_axis[ 2 ]; 
        pY[ 1 ] = R_2[ 0 ] * Y_axis[ 0 ] + R_2[ 1 ] * Y_axis[ 1 ] + R_2[ 2 ] * Y_axis[ 2 ];

        pX[ 2 ] = R_3[ 0 ] * X_axis[ 0 ] + R_3[ 1 ] * X_axis[ 1 ] + R_3[ 2 ] * X_axis[ 2 ]; 
        pY[ 2 ] = R_3[ 0 ] * Y_axis[ 0 ] + R_3[ 1 ] * Y_axis[ 1 ] + R_3[ 2 ] * Y_axis[ 2 ];
    }

    //*******************************************************************************************************
    // - Lagrange 1st order base derivatives - DN[ X, Y ][ 0, 1, 2 ]
    //*******************************************************************************************************
    void GenericRobin_1eg_FullWave::Calculate_DN( double DN[ 2 ][ 3 ] )
    {
        double Cte = 1.0 / ( 2.0 * mArea ); 
        
        // Plane coordinates
        Vector<double> pX, pY; 
        
        Get_Vertices_Surface_Coordinates( pX, pY );

        // dNi/dx
        DN[ 0 ][ 0 ] = Cte * ( pY[ 1 ] - pY[ 2 ] );
        DN[ 0 ][ 1 ] = Cte * ( pY[ 2 ] - pY[ 0 ] );
        DN[ 0 ][ 2 ] = Cte * ( pY[ 0 ] - pY[ 1 ] );
                               
        // dNi/dy              
        DN[ 1 ][ 0 ] = Cte * ( pX[ 2 ] - pX[ 1 ] ); 
        DN[ 1 ][ 1 ] = Cte * ( pX[ 0 ] - pX[ 2 ] );
        DN[ 1 ][ 2 ] = Cte * ( pX[ 1 ] - pX[ 0 ] );
    }

    //*******************************************************************************************************
    // - Edge lengths with sign
    //*******************************************************************************************************
    void GenericRobin_1eg_FullWave::Calculate_Edge_Lengths( Vector<double>& eLength )
    {
		// Edge coordinates vector
		Vector< Vector<double> > R( mNumEdges ); 
		
		for( int i=0; i<mNumEdges; i++ ) 
		{
			R[ i ].resize( 3 );
		}

		// Edge: [0 1]
		R[ 0 ][ 0 ] = X( 2 ) - X( 1 ); 
		R[ 0 ][ 1 ] = Y( 2 ) - Y( 1 ); 
		R[ 0 ][ 2 ] = Z( 2 ) - Z( 1 );

		// Edge: [1 2]
		R[ 1 ][ 0 ] = X( 3 ) - X( 2 ); 
		R[ 1 ][ 1 ] = Y( 3 ) - Y( 2 ); 
		R[ 1 ][ 2 ] = Z( 3 ) - Z( 2 );

		// Edge: [2 0]
		R[ 2 ][ 0 ] = X( 1 ) - X( 3 ); 
		R[ 2 ][ 1 ] = Y( 1 ) - Y( 3 ); 
		R[ 2 ][ 2 ] = Z( 1 ) - Z( 3 );

		// Edge lengths
		for( int i=0; i<mNumEdges; i++ )
		{
			eLength[ i ] = sqrt( R[ i ][ 0 ]*R[ i ][ 0 ] + R[ i ][ 1 ]*R[ i ][ 1 ] + R[ i ][ 2 ]*R[ i ][ 2 ] );
		}
		
		// Check directions
		if( mNodes[ 0 ]->Id() > mNodes[ 1 ]->Id() ) eLength[ 0 ] *= -1.0;
		if( mNodes[ 1 ]->Id() > mNodes[ 2 ]->Id() ) eLength[ 1 ] *= -1.0;
		if( mNodes[ 2 ]->Id() > mNodes[ 0 ]->Id() ) eLength[ 2 ] *= -1.0;
	}

    //*******************************************************************************************************
    // - 1st order Whitney base for 2D triangles
    //*******************************************************************************************************
	void GenericRobin_1eg_FullWave::Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY )
	{
		double L1, L2, L3;
		
		int numPoints = cX.size();

		Nx.Resize( mNumEdges, numPoints );
		Ny.Resize( mNumEdges, numPoints );
		
		double DN[ 2 ][ 3 ]; Calculate_DN( DN );

		double dL1_dX = DN[ 0 ][ 0 ], dL1_dY = DN[ 1 ][ 0 ];
		double dL2_dX = DN[ 0 ][ 1 ], dL2_dY = DN[ 1 ][ 1 ];
		double dL3_dX = DN[ 0 ][ 2 ], dL3_dY = DN[ 1 ][ 2 ];

		Vector<double> eLength( mNumEdges ); Calculate_Edge_Lengths( eLength );

		for( int p=0; p<numPoints; p++ ) 
        {
            L2 = cX[ p ];
            L3 = cY[ p ];
            L1 = 1.0 - L2 - L3;

			// Edge: [0 1]
			Nx[ 0 ][ p ] = ( L1 * dL2_dX - L2 * dL1_dX ) * eLength[ 0 ];
            Ny[ 0 ][ p ] = ( L1 * dL2_dY - L2 * dL1_dY ) * eLength[ 0 ];
		    
			// Edge: [1 2]
			Nx[ 1 ][ p ] = ( L2 * dL3_dX - L3 * dL2_dX ) * eLength[ 1 ];
            Ny[ 1 ][ p ] = ( L2 * dL3_dY - L3 * dL2_dY ) * eLength[ 1 ];
		
			// Edge: [2 0]
			Nx[ 2 ][ p ] = ( L3 * dL1_dX - L1 * dL3_dX ) * eLength[ 2 ];
            Ny[ 2 ][ p ] = ( L3 * dL1_dY - L1 * dL3_dY ) * eLength[ 2 ];
        }
	}

    //**********************************************************************************************************************
    //* - Residual vector and quasi-static flux values ( needed for total J sources calculations )
    //**********************************************************************************************************************
    void GenericRobin_1eg_FullWave::GetResidualVector_Flux( std::map< unsigned int, std::complex<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

        if( ( mPotentialsOn == true ) && ( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 0.0 ) ) 
		{
            double Intg_Ni = mArea / 3.0;
            
            Vector<double> Coef = (*mProperties)( COMPLEX_IBC    );
            Vector<double> Flux = (*mProperties)( COMPLEX_IBC_2o );

            std::complex<double> cFlux( Flux[ 0 ] * cos( Flux[ 1 ] ), Flux[ 0 ] * sin( Flux[ 1 ] ) );

            for( int i=0; i<3; i++ ) 
            { 
                if( ( Coef[ 0 ] == 0.0 ) && ( Coef[ 1 ] == 0.0 ) ) 
                {
                    FluxValue[ mNodes[ i ]->Id() ] = cFlux;
                }
                
                ResidualVector[ i + 3*mNumNodes ] += cFlux * Intg_Ni;
            }
		}
    }

    //**********************************************************************************************************************
    //* - Residual vector for a single plane wave
    //**********************************************************************************************************************
    void GenericRobin_1eg_FullWave::GetResidualVector_Single_Wave( Vector< std::complex<double> >& ResidualVector )
    {
	    // Module and phase of the primary field ( E/A )
        double modul_F = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ]; 
        double phase_F = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];

		if( modul_F == 0.0 ) return;

		// E/A polarization vector ( normalized )
        Vector<double> Fpol( 3 );

        Fpol[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 0 ];
        Fpol[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 1 ];
        Fpol[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 2 ];

        Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );

        // Complex E/A vector
        Vector< std::complex<double> > cF( 3 );

        cF[ 0 ] = Fpol[ 0 ]; 
        cF[ 1 ] = Fpol[ 1 ]; 
        cF[ 2 ] = Fpol[ 2 ];

        cF *= std::complex<double>( modul_F * cos( phase_F ), modul_F * sin( phase_F ) ); 

		if( mPotentialsOn ) 
		{
			cF *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		}

        // Wave propagation vector ( normalized )
        Vector<double> Kdir( 3 );

        Kdir[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 3 ];
        Kdir[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 4 ];
        Kdir[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 5 ];

        Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );

		// Element exterior normal
		Vector<double> Next; 
		
		Calculate_Exterior_Normal( Next );

		// k*n
		double KdirNext = Kdir[ 0 ]*Next[ 0 ] + Kdir[ 1 ]*Next[ 1 ] + Kdir[ 2 ]*Next[ 2 ];

		// Element coefficient matrix ( n x n x E - n x k x E )
        Matrix<double> CMatrix( 3, 3 );

		CMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] - 1.0 ) - ( Kdir[ 0 ]*Next[ 0 ] - KdirNext );
		CMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ]       ) - ( Kdir[ 0 ]*Next[ 1 ]            );
		CMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ]       ) - ( Kdir[ 0 ]*Next[ 2 ]            );
										  
		CMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ]       ) - ( Kdir[ 1 ]*Next[ 0 ]            );
		CMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] - 1.0 ) - ( Kdir[ 1 ]*Next[ 1 ] - KdirNext );
		CMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ]       ) - ( Kdir[ 1 ]*Next[ 2 ]            );
										  
		CMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ]       ) - ( Kdir[ 2 ]*Next[ 0 ]            );
		CMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ]       ) - ( Kdir[ 2 ]*Next[ 1 ]            );
		CMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] - 1.0 ) - ( Kdir[ 2 ]*Next[ 2 ] - KdirNext );

		// Complex value jK
        std::complex<double> jK = Calculate_jK();

		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.0 * mArea;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // XY axis of the co-planar coordinate system 
        Vector<double> X_axis, Y_axis;

		Get_CoPlanar_XY_Axis( X_axis, Y_axis );

		// Input field on Gauss points 
		Vector< std::complex<double> > cFx( numGaussPoints ); 
        Vector< std::complex<double> > cFy( numGaussPoints ); 

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
            double alpha = cX[ gp ];
			double beta  = cY[ gp ];
			
			double rX = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
            double rY = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		    double rZ = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);
            
            double KdirR = Kdir[ 0 ] * rX + Kdir[ 1 ] * rY + Kdir[ 2 ] * rZ;

            Vector< std::complex<double> > gpcF = cF;  gpcF *= std::exp( jK * KdirR );

			std::complex<double> cFx3D = CMatrix[ 0 ][ 0 ] * gpcF[ 0 ] + CMatrix[ 0 ][ 1 ] * gpcF[ 1 ] + CMatrix[ 0 ][ 2 ] * gpcF[ 2 ]; 
            std::complex<double> cFy3D = CMatrix[ 1 ][ 0 ] * gpcF[ 0 ] + CMatrix[ 1 ][ 1 ] * gpcF[ 1 ] + CMatrix[ 1 ][ 2 ] * gpcF[ 2 ]; 
            std::complex<double> cFz3D = CMatrix[ 2 ][ 0 ] * gpcF[ 0 ] + CMatrix[ 2 ][ 1 ] * gpcF[ 1 ] + CMatrix[ 2 ][ 2 ] * gpcF[ 2 ]; 

			cFx[ gp ] = cFx3D * X_axis[ 0 ] + cFy3D * X_axis[ 1 ] + cFz3D * X_axis[ 2 ]; 
			cFy[ gp ] = cFx3D * Y_axis[ 0 ] + cFy3D * Y_axis[ 1 ] + cFz3D * Y_axis[ 2 ];  
		}

        // 1st order Whitney base functions
        Matrix<double> Nx, Ny; 
		
		Whitney2D_Ni_1st( Nx, Ny, cX, cY );

		// Complex value jK/Mu
        std::complex<double> jK_div_Mu = Calculate_jK_div_Mu();

	    // Residual vector
		for( int i=0; i<mNumEdges; i++ )
        {
            std::complex<double> sIntg_FNi( 0.0, 0.0 );

			for( int gp=0; gp<numGaussPoints; gp++ )
			{
				sIntg_FNi += W[ gp ] * ( Nx[ i ][ gp ] * cFx[ gp ] + Ny[ i ][ gp ] * cFy[ gp ] );
			}

            ResidualVector[ i + mNumNodes + mNumVerts ] += jK_div_Mu * sIntg_FNi;
        }  
    }

    //**********************************************************************************************************************
    //* - Residual vector for multiple plane waves
    //**********************************************************************************************************************
    void GenericRobin_1eg_FullWave::GetResidualVector_Multpl_Waves( Vector< Vector<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
	{
        // If plane waves set is empty then return
		if( FluxValue.size() == 0 ) return; 

	    // Complex value jK
        std::complex<double> jK = Calculate_jK();

		// Complex value jK/Mu
        std::complex<double> jK_div_Mu = Calculate_jK_div_Mu();

		// Element exterior normal
		Vector<double> Next; 
		
		Calculate_Exterior_Normal( Next );

		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.0 * mArea;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order Whitney base functions
        Matrix<double> Nx, Ny; 
		
		Whitney2D_Ni_1st( Nx, Ny, cX, cY );

        // XY axis of the co-planar coordinate system 
        Vector<double> X_axis, Y_axis;

		Get_CoPlanar_XY_Axis( X_axis, Y_axis );

		// Iterator over plane waves
		Vector< Vector<double> >::iterator pw_it;

        for( pw_it = FluxValue.begin(); pw_it != FluxValue.end(); ++pw_it )
		{
            // Module and phase of the primary field ( E/A )
            double modul_F = (*pw_it)[ 0 ]; 
            double phase_F = (*pw_it)[ 1 ];

			if( modul_F == 0.0 ) continue;
		    
		    // E/A polarization vector ( normalized )
            Vector<double> Fpol( 3 );
		    
            Fpol[ 0 ] = (*pw_it)[ 2 ];
            Fpol[ 1 ] = (*pw_it)[ 3 ];
            Fpol[ 2 ] = (*pw_it)[ 4 ];
		    
            Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );
		    
            // Complex E/A vector
            Vector< std::complex<double> > cF( 3 );
		    
            cF[ 0 ] = Fpol[ 0 ]; 
            cF[ 1 ] = Fpol[ 1 ]; 
            cF[ 2 ] = Fpol[ 2 ];
		    
            cF *= std::complex<double>( modul_F * cos( phase_F ), modul_F * sin( phase_F ) ); 
		    
		    if( mPotentialsOn ) 
		    {
		    	cF *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		    }
		    
            // Wave propagation vector ( normalized )
            Vector<double> Kdir( 3 );
		    
            Kdir[ 0 ] = (*pw_it)[ 5 ];
            Kdir[ 1 ] = (*pw_it)[ 6 ];
            Kdir[ 2 ] = (*pw_it)[ 7 ];
		    
            Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );

		    // k*n
		    double KdirNext = Kdir[ 0 ]*Next[ 0 ] + Kdir[ 1 ]*Next[ 1 ] + Kdir[ 2 ]*Next[ 2 ];
            
		    // Element coefficient matrix ( n x n x E - n x k x E )
            Matrix<double> CMatrix( 3, 3 );
            
		    CMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] - 1.0 ) - ( Kdir[ 0 ]*Next[ 0 ] - KdirNext );
		    CMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ]       ) - ( Kdir[ 0 ]*Next[ 1 ]            );
		    CMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ]       ) - ( Kdir[ 0 ]*Next[ 2 ]            );
		    								  
		    CMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ]       ) - ( Kdir[ 1 ]*Next[ 0 ]            );
		    CMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] - 1.0 ) - ( Kdir[ 1 ]*Next[ 1 ] - KdirNext );
		    CMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ]       ) - ( Kdir[ 1 ]*Next[ 2 ]            );
		    								  
		    CMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ]       ) - ( Kdir[ 2 ]*Next[ 0 ]            );
		    CMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ]       ) - ( Kdir[ 2 ]*Next[ 1 ]            );
		    CMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] - 1.0 ) - ( Kdir[ 2 ]*Next[ 2 ] - KdirNext );

		    // Input field on Gauss points 
		    Vector< std::complex<double> > cFx( numGaussPoints ); 
            Vector< std::complex<double> > cFy( numGaussPoints ); 
            
		    for( int gp=0; gp<numGaussPoints; gp++ )
		    {
                double alpha = cX[ gp ];
		    	double beta  = cY[ gp ];
		    	
		    	double rX = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
                double rY = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		        double rZ = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);
                
                double KdirR = Kdir[ 0 ] * rX + Kdir[ 1 ] * rY + Kdir[ 2 ] * rZ;
            
                Vector< std::complex<double> > gpcF = cF;  gpcF *= std::exp( jK * KdirR );
            
		    	std::complex<double> cFx3D = CMatrix[ 0 ][ 0 ] * gpcF[ 0 ] + CMatrix[ 0 ][ 1 ] * gpcF[ 1 ] + CMatrix[ 0 ][ 2 ] * gpcF[ 2 ]; 
                std::complex<double> cFy3D = CMatrix[ 1 ][ 0 ] * gpcF[ 0 ] + CMatrix[ 1 ][ 1 ] * gpcF[ 1 ] + CMatrix[ 1 ][ 2 ] * gpcF[ 2 ]; 
                std::complex<double> cFz3D = CMatrix[ 2 ][ 0 ] * gpcF[ 0 ] + CMatrix[ 2 ][ 1 ] * gpcF[ 1 ] + CMatrix[ 2 ][ 2 ] * gpcF[ 2 ]; 
            
		    	cFx[ gp ] = cFx3D * X_axis[ 0 ] + cFy3D * X_axis[ 1 ] + cFz3D * X_axis[ 2 ]; 
		    	cFy[ gp ] = cFx3D * Y_axis[ 0 ] + cFy3D * Y_axis[ 1 ] + cFz3D * Y_axis[ 2 ];  
		    }

	        // Residual vector
		    for( int i=0; i<mNumEdges; i++ )
            {
                std::complex<double> sIntg_FNi( 0.0, 0.0 );
            
		    	for( int gp=0; gp<numGaussPoints; gp++ )
		    	{
		    		sIntg_FNi += W[ gp ] * ( Nx[ i ][ gp ] * cFx[ gp ] + Ny[ i ][ gp ] * cFy[ gp ] );
		    	}
            
                ResidualVector[ i + mNumNodes + mNumVerts ] += jK_div_Mu * sIntg_FNi;
            }  
        }
	}

    //**********************************************************************************************************************
    //* - Residual vector for Gaussian beams
    //**********************************************************************************************************************
    void GenericRobin_1eg_FullWave::GetResidualVector_Gauss_Beams( Vector< std::complex<double> >& ResidualVector )
    {
        // Gaussian beam waist radius 
        double GBR_pol = (*mProperties)( COMPLEX_IBC_2o )[ 0 ]; 
        double GBR_per = (*mProperties)( COMPLEX_IBC_2o )[ 1 ];

		if( ( GBR_pol <= 0.0 ) || ( GBR_per <= 0.0 ) ) return;

		// Field polarization vector ( normalized )
        Vector<double> Fpol( 3 );

        Fpol[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 0 ];
        Fpol[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 1 ];
        Fpol[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 2 ];

        Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );

        // Wave propagation vector ( normalized )
        Vector<double> Kdir( 3 );

        Kdir[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 3 ];
        Kdir[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 4 ];
        Kdir[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 5 ];

        Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );

        // Remove any component of the polarization vector along the propagation direction (Gram-Schmidt),
        // so the waist-plane axes ( Fpol, Fper ) are orthogonal to Kdir even if the input vectors are not
        // exactly perpendicular. The negated comparison also catches NaN from zero-length input vectors.
        double FpolKdir = Fpol[ 0 ]*Kdir[ 0 ] + Fpol[ 1 ]*Kdir[ 1 ] + Fpol[ 2 ]*Kdir[ 2 ];

        Fpol[ 0 ] -= FpolKdir * Kdir[ 0 ];
        Fpol[ 1 ] -= FpolKdir * Kdir[ 1 ];
        Fpol[ 2 ] -= FpolKdir * Kdir[ 2 ];

        double normFpol = sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] );

        if( !( normFpol > 1.0e-12 ) ) return;

        Fpol *= ( 1.0 / normFpol );

        // Fper = Kdir X Fpol ( normalized )
        Vector<double> Fper( 3 );

        Fper[ 0 ] = Kdir[ 1 ]*Fpol[ 2 ] - Kdir[ 2 ]*Fpol[ 1 ];
        Fper[ 1 ] = Kdir[ 2 ]*Fpol[ 0 ] - Kdir[ 0 ]*Fpol[ 2 ];
        Fper[ 2 ] = Kdir[ 0 ]*Fpol[ 1 ] - Kdir[ 1 ]*Fpol[ 0 ];

        Fper *= ( 1.0 / sqrt( Fper[ 0 ]*Fper[ 0 ] + Fper[ 1 ]*Fper[ 1 ] + Fper[ 2 ]*Fper[ 2 ] ) );

        // Field parallel to the polarization vector
        Vector< std::complex<double> > Fpar( 3 );

        Fpar[ 0 ] = Fpol[ 0 ]; 
        Fpar[ 1 ] = Fpol[ 1 ]; 
        Fpar[ 2 ] = Fpol[ 2 ];

        double mod_Fpar = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ]; 
        double pha_Fpar = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];
        
        Fpar *= std::complex<double>( mod_Fpar * cos( pha_Fpar ), mod_Fpar * sin( pha_Fpar ) );

        // Field perpendicular to the polarization and wave vectors
        Vector< std::complex<double> > Fppw( 3 );

        Fppw[ 0 ] = Fper[ 0 ]; 
        Fppw[ 1 ] = Fper[ 1 ]; 
        Fppw[ 2 ] = Fper[ 2 ];

        double mod_Fppw = (*mProperties)( COMPLEX_IBC )[ 0 ]; 
        double pha_Fppw = (*mProperties)( COMPLEX_IBC )[ 1 ];
        
        Fppw *= std::complex<double>( mod_Fppw * cos( pha_Fppw ), mod_Fppw * sin( pha_Fppw ) );

        // Total field on waist central point
        Vector< std::complex<double> > Ftot( 3 );

        Ftot[ 0 ] = Fpar[ 0 ] + Fppw[ 0 ]; 
        Ftot[ 1 ] = Fpar[ 1 ] + Fppw[ 1 ]; 
        Ftot[ 2 ] = Fpar[ 2 ] + Fppw[ 2 ];

		if( mPotentialsOn ) 
		{
			Ftot *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		}

		// Complex wave vector jK 
        std::complex<double> jK = Calculate_jK();

        // Rayleight range
        double ZR_pol = 0.5 * std::imag( jK ) * GBR_pol * GBR_pol;
        double ZR_per = 0.5 * std::imag( jK ) * GBR_per * GBR_per;

        // Gaussian beam waist center coordinates
        double Xo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ];			
        double Yo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ];
        double Zo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 2 ]; 

		// Element exterior normal
		Vector<double> Next; 
		
		Calculate_Exterior_Normal( Next );

		// k*n
		double KdirNext = Kdir[ 0 ]*Next[ 0 ] + Kdir[ 1 ]*Next[ 1 ] + Kdir[ 2 ]*Next[ 2 ];

		// Element coefficient matrix ( n x n x E - n x k x E )
        Matrix<double> CMatrix( 3, 3 );

		CMatrix[ 0 ][ 0 ] = ( Next[ 0 ]*Next[ 0 ] - 1.0 ) - ( Kdir[ 0 ]*Next[ 0 ] - KdirNext );
		CMatrix[ 0 ][ 1 ] = ( Next[ 0 ]*Next[ 1 ]       ) - ( Kdir[ 0 ]*Next[ 1 ]            );
		CMatrix[ 0 ][ 2 ] = ( Next[ 0 ]*Next[ 2 ]       ) - ( Kdir[ 0 ]*Next[ 2 ]            );
										  
		CMatrix[ 1 ][ 0 ] = ( Next[ 1 ]*Next[ 0 ]       ) - ( Kdir[ 1 ]*Next[ 0 ]            );
		CMatrix[ 1 ][ 1 ] = ( Next[ 1 ]*Next[ 1 ] - 1.0 ) - ( Kdir[ 1 ]*Next[ 1 ] - KdirNext );
		CMatrix[ 1 ][ 2 ] = ( Next[ 1 ]*Next[ 2 ]       ) - ( Kdir[ 1 ]*Next[ 2 ]            );
										  
		CMatrix[ 2 ][ 0 ] = ( Next[ 2 ]*Next[ 0 ]       ) - ( Kdir[ 2 ]*Next[ 0 ]            );
		CMatrix[ 2 ][ 1 ] = ( Next[ 2 ]*Next[ 1 ]       ) - ( Kdir[ 2 ]*Next[ 1 ]            );
		CMatrix[ 2 ][ 2 ] = ( Next[ 2 ]*Next[ 2 ] - 1.0 ) - ( Kdir[ 2 ]*Next[ 2 ] - KdirNext );

		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.0 * mArea;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // XY axis of the co-planar coordinate system 
        Vector<double> X_axis, Y_axis;

		Get_CoPlanar_XY_Axis( X_axis, Y_axis );

		// Field vector on Gauss points 
		Vector< std::complex<double> > cFx( numGaussPoints ); 
        Vector< std::complex<double> > cFy( numGaussPoints ); 

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
            // Gauss point coordinates
            double alpha = cX[ gp ];
			double beta  = cY[ gp ];
			
			double Xg = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
            double Yg = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		    double Zg = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);

            // Gauss point position from waist center
            Vector<double> Rg( 3 );

            Rg[ 0 ] = Xg - Xo; 
            Rg[ 1 ] = Yg - Yo; 
            Rg[ 2 ] = Zg - Zo;

            // Gauss point distance from waist plane
            double Zd = Kdir[ 0 ]*Rg[ 0 ] + Kdir[ 1 ]*Rg[ 1 ] + Kdir[ 2 ]*Rg[ 2 ];

            // Accumulated phase 
            std::complex<double> AccPhase = jK * Zd;

            // Gouy phase shift
            std::complex<double> GouyPhase( 0.0, 0.5 * ( std::atan( Zd / ZR_pol ) + std::atan( Zd / ZR_per ) ) );

            // Q parameters
            std::complex<double> Qpol( Zd, -ZR_pol );
            std::complex<double> Qper( Zd, -ZR_per );

            // Node coordinates projection on waist plane
            double Xpol = Fpol[ 0 ]*Rg[ 0 ] + Fpol[ 1 ]*Rg[ 1 ] + Fpol[ 2 ]*Rg[ 2 ];
            double Xper = Fper[ 0 ]*Rg[ 0 ] + Fper[ 1 ]*Rg[ 1 ] + Fper[ 2 ]*Rg[ 2 ];

            // Q parameters phase
            std::complex<double> QpPhase = jK * 0.5 * ( ( Xpol * Xpol / Qpol ) + ( Xper * Xper / Qper ) );

            // Field normalization: sqrt( w0_pol*w0_per / ( w_pol*w_per ) ) = sqrt( ZR_pol*ZR_per / ( |Qpol|*|Qper| ) ) 
            double Beta = std::sqrt( ZR_pol * ZR_per / ( std::abs( Qpol ) * std::abs( Qper ) ) );

            // Field value on Gauss point
            Vector< std::complex<double> > cFg;

            cFg  = Ftot; 
            cFg *= Beta * std::exp( AccPhase - GouyPhase + QpPhase );

            // Edge field value on Gauss point
			std::complex<double> cFx3D = CMatrix[ 0 ][ 0 ] * cFg[ 0 ] + CMatrix[ 0 ][ 1 ] * cFg[ 1 ] + CMatrix[ 0 ][ 2 ] * cFg[ 2 ]; 
            std::complex<double> cFy3D = CMatrix[ 1 ][ 0 ] * cFg[ 0 ] + CMatrix[ 1 ][ 1 ] * cFg[ 1 ] + CMatrix[ 1 ][ 2 ] * cFg[ 2 ]; 
            std::complex<double> cFz3D = CMatrix[ 2 ][ 0 ] * cFg[ 0 ] + CMatrix[ 2 ][ 1 ] * cFg[ 1 ] + CMatrix[ 2 ][ 2 ] * cFg[ 2 ]; 

			cFx[ gp ] = cFx3D * X_axis[ 0 ] + cFy3D * X_axis[ 1 ] + cFz3D * X_axis[ 2 ]; 
			cFy[ gp ] = cFx3D * Y_axis[ 0 ] + cFy3D * Y_axis[ 1 ] + cFz3D * Y_axis[ 2 ];  
		}

        // 1st order Whitney base functions
        Matrix<double> Nx, Ny; 
		
		Whitney2D_Ni_1st( Nx, Ny, cX, cY );

		// Complex value jK/Mu
        std::complex<double> jK_div_Mu = Calculate_jK_div_Mu();

	    // Residual vector
		for( int i=0; i<mNumEdges; i++ )
        {
            std::complex<double> sIntg_FNi( 0.0, 0.0 );

			for( int gp=0; gp<numGaussPoints; gp++ )
			{
				sIntg_FNi += W[ gp ] * ( Nx[ i ][ gp ] * cFx[ gp ] + Ny[ i ][ gp ] * cFy[ gp ] );
			}

            ResidualVector[ i + mNumNodes + mNumVerts ] += jK_div_Mu * sIntg_FNi;
        } 
    }

    //**********************************************************************************************************************
    //* - Residual vector for plane waves and Gaussian beams
    //**********************************************************************************************************************
    void GenericRobin_1eg_FullWave::GetResidualVector_Flux( Vector< Vector<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
    {
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 1.0 )
        {
            GetResidualVector_Single_Wave( ResidualVector );
        }
		else if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 2.0 )
		{
            GetResidualVector_Multpl_Waves( FluxValue, ResidualVector );
		}
        else if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            GetResidualVector_Gauss_Beams( ResidualVector );
        }
    }

    //**********************************************************************************************************************
    //* - Stiffness matrix
    //**********************************************************************************************************************
    void GenericRobin_1eg_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Setting element stiffness matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );
        
        // Robin coefficient
        std::complex<double> cRobinCoeff( 0.0, 0.0 );

        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 0.0 )
        {
		    cRobinCoeff -= std::complex<double>( (*mProperties)( COMPLEX_IBC )[ 0 ], (*mProperties)( COMPLEX_IBC )[ 1 ] );
        }
        else
        {
            cRobinCoeff -= Calculate_jK_div_Mu();
        }

        // Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.0 * mArea;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order Whitney base functions
        Matrix<double> Nx, Ny;

		Whitney2D_Ni_1st( Nx, Ny, cX, cY );

		// Edge element matrix
		Matrix< std::complex<double> > EdgeEMatrix( mNumEdges, mNumEdges );

		// Upper diagonal 
		for( int i=0; i<mNumEdges; i++ )
		{
			for( int j=i; j<mNumEdges; j++ )
			{
				double sIntg_NiNj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					sIntg_NiNj += W[ gp ] * ( Nx[ i ][ gp ] * Nx[ j ][ gp ] + Ny[ i ][ gp ] * Ny[ j ][ gp ] );
				}

				EdgeEMatrix[ i ][ j ] = cRobinCoeff * sIntg_NiNj;
			}
		}

		// Lower diagonal 
        for( int i=0; i<mNumEdges; i++ )
        {
            for( int j=i+1; j<mNumEdges; j++ )
            {
                EdgeEMatrix[ j ][ i ] = EdgeEMatrix[ i ][ j ];
            }
        }        

		// Edge element matrix inclusion into stiffness matrix
		for( int i=0; i<mNumEdges; i++ )
        {
            for( int j=0; j<mNumEdges; j++ )
            {
                StiffMatrix[ i + mNumNodes + mNumVerts ][ j + mNumNodes + mNumVerts ] = EdgeEMatrix[ i ][ j ];
            }
        }     
    }
} 