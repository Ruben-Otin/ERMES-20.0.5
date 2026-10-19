
#include "CoaxialPortTEM_1eg_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //*******************************************************************************************************
    //* COMPLEX_IBC[ 0 ] : Coaxial port identification number. 
	//* 
	//* COMPLEX_IBC_2o[ 0 ] : Module of a complex number multiplying the input TEM field.
    //* COMPLEX_IBC_2o[ 1 ] : Phase  of a complex number multiplying the input TEM field.
    //*
    //* COMPLEX_NEUMANN_FLOW[ 0 ] : X coordinate of the coaxial center.
    //* COMPLEX_NEUMANN_FLOW[ 1 ] : Y coordinate of the coaxial center.
    //* COMPLEX_NEUMANN_FLOW[ 2 ] : Z coordinate of the coaxial center.
	//* 
    //* COMPLEX_NEUMANN_FLOW[ 3 ] : Diameter of the inner    cylinder of the coaxial waveguide.
    //* COMPLEX_NEUMANN_FLOW[ 4 ] : Diameter of the exterior cylinder of the coaxial waveguide.
	// 
    //* SINUSOIDAL_SURFACE_CURRENT[ 0 ] : Relative electric permittivity of the coaxial port.
    //* SINUSOIDAL_SURFACE_CURRENT[ 1 ] : Relative magnetic permeability of the coaxial port.
    //*******************************************************************************************************

    //*******************************************************************************************************
    //* - Global indexing of the DOFs
    //*******************************************************************************************************	
    void CoaxialPortTEM_1eg_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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

            for( int i=0; i<mNumVerts; i++ ) 
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

	//*******************************************************************************************************
    //* - Area of the element
    //*******************************************************************************************************
    double CoaxialPortTEM_1eg_FullWave::Calculate_Area() 
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

    //*******************************************************************************************************
    //* - TEM coaxial waveguide propagation constant = j * ( k / mu )
    //*******************************************************************************************************
    std::complex<double> CoaxialPortTEM_1eg_FullWave::Calculate_KMuCoeff() 
	{
		// Material properties of the coaxial port
		Vector<double> CC_MatProp = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT ); 

		double epr = CC_MatProp[ 0 ];
		double mur = CC_MatProp[ 1 ];

		double pi = 3.14159265358979;
		double eo = 8.8541878176e-12;
		double mo = pi * 4.0e-7;

		double ep = epr * eo;
		double mu = mur * mo;

        double w = (*mProperties)( FREQUENCY );

		std::complex<double> cKMu( 0.0, w * sqrt( ep / mu ) );
	
		return cKMu;
	}

    //*******************************************************************************************************
    //* - Normalization coefficient to make Integral_Surf_Port( E_TEM x H_TEM ) = 1.0
    //*******************************************************************************************************
    std::complex<double> CoaxialPortTEM_1eg_FullWave::Calculate_NrmCoeff() 
	{
		// Material properties of the coaxial port
		Vector<double> CC_MatProp = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT ); 

		double epr = CC_MatProp[ 0 ];
		double mur = CC_MatProp[ 1 ];	
		
		double pi = 3.14159265358979;
		double eo = 8.8541878176e-12;
		double mo = pi * 4.0e-7;

		double ep = epr * eo;
		double mu = mur * mo;

		double ctcSup = sqrt( mu / ep );

		// Dimensions of the coaxial port
		Vector<double> CC_Diameters = (*mProperties)( COMPLEX_NEUMANN_FLOW ); 

		double a = CC_Diameters[ 3 ];
		double b = CC_Diameters[ 4 ];

		double ctcInf = 2.0 * pi * log( b / a );

		std::complex<double> cNormCoeff( sqrt( ctcSup / ctcInf ), 0.0 );
	
		return cNormCoeff;
	}

    //*******************************************************************************************************
    //* - r0 vector ( distance to coordinates origin )
    //*******************************************************************************************************
	void CoaxialPortTEM_1eg_FullWave::Calculate_r0( double* r0 ) 
	{
		Vector<double> cNF = (*mProperties)( COMPLEX_NEUMANN_FLOW ); 
		
		r0[ 0 ] = cNF[ 0 ];
		r0[ 1 ] = cNF[ 1 ];
	    r0[ 2 ] = cNF[ 2 ];
	}

    //*******************************************************************************************************
    //* - Coaxial port identification number
    //*******************************************************************************************************
	double CoaxialPortTEM_1eg_FullWave::GetPortNumber()
	{
		Vector<double> cIBC = (*mProperties)( COMPLEX_IBC );

		return cIBC[ 0 ];
	}

    //*******************************************************************************************************
    //* - Exterior normal 
    //*******************************************************************************************************	
    void CoaxialPortTEM_1eg_FullWave::Calculate_Exterior_Normal( Vector<double>& Next )
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

    //*******************************************************************************************************
    // - Get X, Y axis of the co-planar coordinate system 
    //*******************************************************************************************************
    void CoaxialPortTEM_1eg_FullWave::Get_CoPlanar_XY_Axis( Vector<double>& X_axis,  Vector<double>& Y_axis )
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
    void CoaxialPortTEM_1eg_FullWave::Get_Vertices_Surface_Coordinates( Vector<double>& pX,  Vector<double>& pY )
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
    void CoaxialPortTEM_1eg_FullWave::Calculate_DN( double DN[ 2 ][ 3 ] )
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
    void CoaxialPortTEM_1eg_FullWave::Calculate_Edge_Lengths( Vector<double>& eLength )
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
	void CoaxialPortTEM_1eg_FullWave::Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY )
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

    //*******************************************************************************************************
    //* - Input field vector on Gauss points
    //*******************************************************************************************************
    void CoaxialPortTEM_1eg_FullWave::Input_Field_Vector( Vector<double>& Fx, Vector<double>& Fy, Vector<double>& cX, Vector<double>& cY )
	{
        // XY axis of the co-planar coordinate system 
        Vector<double> X_axis, Y_axis;

		Get_CoPlanar_XY_Axis( X_axis, Y_axis );

		// Port vector
		double r0[ 3 ]; Calculate_r0( r0 ); 

	    // Input field vector on Gauss points
		int numGaussPoints = cX.size();

		Fx.resize( numGaussPoints ); 
		Fy.resize( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            double alpha = cX[ gp ];
			double beta  = cY[ gp ];
			
			double rX = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
            double rY = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		    double rZ = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);
			 
			double rhoX = rX - r0[ 0 ];
			double rhoY = rY - r0[ 1 ];
			double rhoZ = rZ - r0[ 2 ];
			
			double InvRho2 = 1.0 / ( rhoX * rhoX + rhoY * rhoY + rhoZ * rhoZ );
			
			double Fx3D = InvRho2 * rhoX; 
			double Fy3D = InvRho2 * rhoY;  
			double Fz3D = InvRho2 * rhoZ; 

			Fx[ gp ] = Fx3D * X_axis[ 0 ] + Fy3D * X_axis[ 1 ] + Fz3D * X_axis[ 2 ]; 
			Fy[ gp ] = Fx3D * Y_axis[ 0 ] + Fy3D * Y_axis[ 1 ] + Fz3D * Y_axis[ 2 ];  
        }
	}

    //*******************************************************************************************************
    //* - Residual vector 
    //*******************************************************************************************************
    void CoaxialPortTEM_1eg_FullWave::GetResidualVector( Vector< std::complex<double> >& ResidualVector ) 
    {
        // Initializing residual vector
		ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );
		
		// If power is cero then return
		if( (*mProperties)( COMPLEX_IBC_2o )[ 0 ] == 0.0 )
		{
			return;
		}

		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

        // 1st order Whitney base functions
        Matrix<double> Nx, Ny; 
		
		Whitney2D_Ni_1st( Nx, Ny, cX, cY );

		// Input field on Gauss points ( F/r )
		Vector<double> Fx, Fy; 
		
		Input_Field_Vector( Fx, Fy, cX, cY );

        // Jacobian
        double Jacob = 2.0 * mArea;

		// TEM wave propagation cte: j*( k/mu ) 
		std::complex<double> cKMuCoeff = Calculate_KMuCoeff();

		// Normalization cte: sqrt( sqrt( mu/ep ) / 2 * pi * log( b/a ) )
		std::complex<double> cNrmCoeff = Calculate_NrmCoeff();

		// Integration cte
        std::complex<double> cIntgCte = -2.0 * Jacob * ( cKMuCoeff * cNrmCoeff );

        // Scale factor on E_input
		double Mod_Ein = (*mProperties)( COMPLEX_IBC_2o )[ 0 ];
	    double Pha_Ein = (*mProperties)( COMPLEX_IBC_2o )[ 1 ];

		std::complex<double> InputPortCoeff( Mod_Ein * cos( Pha_Ein ), Mod_Ein * sin( Pha_Ein ) );

		cIntgCte *= InputPortCoeff; 

		// Multiply integration cte by 1/jw for the AV potential formulation
		if( mPotentialsOn )
		{
		    std::complex<double> jw( 0.0, (*mProperties)( FREQUENCY ) );

		    cIntgCte *= ( 1.0 / jw ); 
		}

	    // Residual vector
		for( int i=0; i<mNumEdges; i++ )
        {
            double sIntg_FNi = 0.0;

			for( int gp=0; gp<numGaussPoints; gp++ )
			{
				sIntg_FNi += W[ gp ] * ( Nx[ i ][ gp ] * Fx[ gp ] + Ny[ i ][ gp ] * Fy[ gp ] );
			}

            ResidualVector[ i + mNumNodes + mNumVerts ] = cIntgCte * sIntg_FNi;
        }  
    }

    //*******************************************************************************************************
    //* - Stiffness matrix
    //*******************************************************************************************************	
    void CoaxialPortTEM_1eg_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Initializing matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );
        
		// - j * ( k / mu ) 
		std::complex<double> mcKMu = -Calculate_KMuCoeff(); 

        // Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.0 * mArea;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[gp] *= detJ; 
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

				EdgeEMatrix[ i ][ j ] = mcKMu * sIntg_NiNj;
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

    //*******************************************************************************************************
    //* - TEM projection on a plane 
    //*******************************************************************************************************
	std::complex<double> CoaxialPortTEM_1eg_FullWave::GetProjection()
	{
		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

        // 1st order Whitney base functions
		Matrix<double> Nx, Ny;

		Whitney2D_Ni_1st( Nx, Ny, cX, cY );

		// XY axis of the co-planar coordinate system 
        Vector<double> X_axis, Y_axis;

		Get_CoPlanar_XY_Axis( X_axis, Y_axis );

		// Port vector
		double r0[ 3 ]; Calculate_r0( r0 ); 

        // Ft on edges
		Vector< std::complex<double> > cFt( mNumEdges );

		if( mPotentialsOn )
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }
		}
		else
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }
		}

		// Inverse rho on Gauss points
		Vector<double> InvRho( numGaussPoints );

		// Fr on Gauss points
		Vector< std::complex<double> > cFr( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
			// Inverse rho on Gauss point
            double alpha = cX[ gp ];
			double beta  = cY[ gp ];

			double rX = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
            double rY = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		    double rZ = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);

			double rhoX = rX - r0[ 0 ];
			double rhoY = rY - r0[ 1 ];
			double rhoZ = rZ - r0[ 2 ];

			InvRho[ gp ] = 1.0 / sqrt( rhoX * rhoX + rhoY * rhoY + rhoZ * rhoZ );

			// Fr on Gauss point
		    std::complex<double> cFx( 0.0, 0.0 );
			std::complex<double> cFy( 0.0, 0.0 );
		
			for( int e=0; e<mNumEdges; e++ )
			{
		        cFx += cFt[ e ] * Nx[ e ][ gp ];
		        cFy += cFt[ e ] * Ny[ e ][ gp ];
			}

			std::complex<double> cFx3D = cFx * X_axis[ 0 ] + cFy * Y_axis[ 0 ];
			std::complex<double> cFy3D = cFx * X_axis[ 1 ] + cFy * Y_axis[ 1 ];
			std::complex<double> cFz3D = cFx * X_axis[ 2 ] + cFy * Y_axis[ 2 ];

			cFr[ gp ] = InvRho[ gp ] * ( cFx3D * rhoX + cFy3D * rhoY + cFz3D * rhoZ );
        }

		if( mPotentialsOn )
		{
            cFr *= std::complex<double>( 0.0, (*mProperties)( FREQUENCY ) );
		}

		// Integral ( Fr / rho ) over element surface
		std::complex<double> TEMProject( 0.0, 0.0 );

		for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
			TEMProject += W[ gp ] * InvRho[ gp ] * cFr[ gp ];  
        }

		return TEMProject;
	}

    //*******************************************************************************************************
    //* - Normalization term ( useful when using PMC symmetries )
    //*******************************************************************************************************
	std::complex<double> CoaxialPortTEM_1eg_FullWave::GetNormalization()
	{
		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		// Port vector
		double r0[ 3 ]; Calculate_r0( r0 ); 

        // Surface integral of inverse rho2
		double intg_InvRho2 = 0.0;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            double alpha = cX[ gp ];
			double beta  = cY[ gp ];

			double rX = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
            double rY = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		    double rZ = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);

			double rhoX = rX - r0[ 0 ];
			double rhoY = rY - r0[ 1 ];
			double rhoZ = rZ - r0[ 2 ];

			double InvRho2 = 1.0 / ( rhoX*rhoX + rhoY*rhoY + rhoZ*rhoZ );

			intg_InvRho2 += W[ gp ] * InvRho2;
        }
		
		// Normalization coefficient
		std::complex<double> Normalization( intg_InvRho2, 0.0 );

		std::complex<double> cNormCoeff = Calculate_NrmCoeff();

		Normalization *= cNormCoeff;

		return Normalization;
	}
} 