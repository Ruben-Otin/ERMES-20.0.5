
#include "RWPortTE10_1eg_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //*******************************************************************************************************
    //* COMPLEX_IBC[ 0 ] : Rectangular waveguide port identification number. 
	//* 
	//* COMPLEX_IBC_2o[ 0 ] : Module of a complex number multiplying the input TE10 field.
    //* COMPLEX_IBC_2o[ 1 ] : Phase  of a complex number multiplying the input TE10 field.
	//* 
    //* COMPLEX_NEUMANN_FLOW[ 0 ] : X coordinate of the origin corner.
    //* COMPLEX_NEUMANN_FLOW[ 1 ] : Y coordinate of the origin corner.
    //* COMPLEX_NEUMANN_FLOW[ 2 ] : Z coordinate of the origin corner.
	//* 
    //* SINUSOIDAL_SURFACE_CURRENT[ 0 ] : X coordinate of the upper corner.
    //* SINUSOIDAL_SURFACE_CURRENT[ 1 ] : Y coordinate of the upper corner.
	//* SINUSOIDAL_SURFACE_CURRENT[ 2 ] : Z coordinate of the upper corner.
	//* 
    //* SINUSOIDAL_SURFACE_CURRENT[ 3 ] : X coordinate of the right corner.
	//* SINUSOIDAL_SURFACE_CURRENT[ 4 ] : Y coordinate of the right corner.
    //* SINUSOIDAL_SURFACE_CURRENT[ 5 ] : Z coordinate of the right corner.
    //*******************************************************************************************************
	
    //*******************************************************************************************************
    //* - Global indexing of the DOFs
    //*******************************************************************************************************	
    void RWPortTE10_1eg_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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

    //*******************************************************************************************************
    //* - Area of the element
    //*******************************************************************************************************
    double RWPortTE10_1eg_FullWave::Calculate_Area() 
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
    //* - TE10 propagation coefficient
    //*******************************************************************************************************
    std::complex<double> RWPortTE10_1eg_FullWave::Calculate_BetaMuCoef() 
	{
		double w = (*mProperties)( FREQUENCY );
		double a = Calculate_Width();

		double pi = 3.14159265358979;
		double eo = 8.8541878176e-12;
		double mu = pi * 4.0e-7;

		double ko = w * sqrt( mu * eo );
		double kc = pi / a;

		std::complex<double> cBetaMu;
	
		if( ko < kc ) 
        {
            cBetaMu = std::complex<double>( -sqrt( kc*kc - ko*ko ) / mu, 0.0 );
        }
		else
        {
            cBetaMu = std::complex<double>( 0.0, sqrt( ko*ko - kc*kc ) / mu );
		}

		return cBetaMu;
	}

    //*******************************************************************************************************
    //* - Normalization coefficient to make Integral_Surf_Port( E_TE10 x H_TE10 ) = 1.0
    //*******************************************************************************************************
    std::complex<double> RWPortTE10_1eg_FullWave::Calculate_N10() 
	{
	    std::complex<double> cN10;

		std::complex<double> cBetaMu = Calculate_BetaMuCoef();

        double rBetaMu = -std::real( cBetaMu );
		double iBetaMu =  std::imag( cBetaMu );

        if( rBetaMu == 0.0 ) 
        {
            cN10 = std::complex<double>( 1.0 / sqrt( iBetaMu ), 0.0 );
        }
        else
        {
            double dummy = 0.5 * sqrt( 2.0 / rBetaMu );

            cN10 = std::complex<double>( dummy, -dummy );
        }

        double w = (*mProperties)( FREQUENCY );
		double a = Calculate_Width ();
		double h = Calculate_Height();

        cN10 *= -sqrt( 2.0 * w / ( a * h ) );

		return cN10;
	}

    //*******************************************************************************************************
    //* - Waveguide width
    //*******************************************************************************************************
	double RWPortTE10_1eg_FullWave::Calculate_Width()
	{
		// Corners coordinates
		Vector<double> cNF = (*mProperties)( COMPLEX_NEUMANN_FLOW );
		Vector<double> SSC = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT );

		// Horizontal vector
		double t[ 3 ];

		t[ 0 ] = SSC[ 3 ] - cNF[ 0 ];
		t[ 1 ] = SSC[ 4 ] - cNF[ 1 ];
		t[ 2 ] = SSC[ 5 ] - cNF[ 2 ];

		// Waveguide width
		return sqrt( t[ 0 ]*t[ 0 ] + t[ 1 ]*t[ 1 ] + t[ 2 ]*t[ 2 ] );
	}

    //*******************************************************************************************************
    //* - Waveguide height
    //*******************************************************************************************************
	double RWPortTE10_1eg_FullWave::Calculate_Height()
	{
		// Corners coordinates
		Vector<double> cNF = (*mProperties)( COMPLEX_NEUMANN_FLOW );
		Vector<double> SSC = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT );

		// Vertical vector
		double b[ 3 ];
		
		b[ 0 ] = SSC[ 0 ] - cNF[ 0 ];
		b[ 1 ] = SSC[ 1 ] - cNF[ 1 ];
		b[ 2 ] = SSC[ 2 ] - cNF[ 2 ];

		// Waveguide height
		return sqrt( b[ 0 ]*b[ 0 ] + b[ 1 ]*b[ 1 ] + b[ 2 ]*b[ 2 ] );
	}
		
    //*******************************************************************************************************
    //* - b vector ( vertical ) 
    //*******************************************************************************************************
	void RWPortTE10_1eg_FullWave::Calculate_b( double* b ) 
	{
		Vector<double> cNF = (*mProperties)( COMPLEX_NEUMANN_FLOW );
		Vector<double> SSC = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT );
	
		b[ 0 ] = SSC[ 0 ] - cNF[ 0 ];
		b[ 1 ] = SSC[ 1 ] - cNF[ 1 ];
		b[ 2 ] = SSC[ 2 ] - cNF[ 2 ];

		double height = Calculate_Height();
		double b_norm = 1.0 / height;

        b[ 0 ] *= b_norm; 
        b[ 1 ] *= b_norm; 
        b[ 2 ] *= b_norm;
	}

    //*******************************************************************************************************
    //* - t vector ( horizontal ) 
    //*******************************************************************************************************
	void RWPortTE10_1eg_FullWave::Calculate_t( double* t ) 
	{
		Vector<double> cNF = (*mProperties)( COMPLEX_NEUMANN_FLOW );
		Vector<double> SSC = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT );

		t[ 0 ] = SSC[ 3 ] - cNF[ 0 ];
		t[ 1 ] = SSC[ 4 ] - cNF[ 1 ];
		t[ 2 ] = SSC[ 5 ] - cNF[ 2 ];

		double width = Calculate_Width();
		double tnorm = 1.0 / width;

        t[ 0 ] *= tnorm; 
        t[ 1 ] *= tnorm; 
        t[ 2 ] *= tnorm;
	}

    //*******************************************************************************************************
    //* - r0 vector ( distance to origin )
    //*******************************************************************************************************
	void RWPortTE10_1eg_FullWave::Calculate_r0( double* r0 ) 
	{
		Vector<double> cNF = (*mProperties)( COMPLEX_NEUMANN_FLOW );
		
		r0[ 0 ] = cNF[ 0 ];
		r0[ 1 ] = cNF[ 1 ];
	    r0[ 2 ] = cNF[ 2 ];
	}

	//*******************************************************************************************************
    //* - Port identification number
    //*******************************************************************************************************
	double RWPortTE10_1eg_FullWave::GetPortNumber()
	{
		Vector<double> cIBC = (*mProperties)( COMPLEX_IBC );

		return cIBC[ 0 ];
	}

    //*******************************************************************************************************
    //* - Exterior normal 
    //*******************************************************************************************************	
    void RWPortTE10_1eg_FullWave::Calculate_Exterior_Normal( Vector<double>& Next )
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
    void RWPortTE10_1eg_FullWave::Get_CoPlanar_XY_Axis( Vector<double>& X_axis,  Vector<double>& Y_axis )
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
    void RWPortTE10_1eg_FullWave::Get_Vertices_Surface_Coordinates( Vector<double>& pX,  Vector<double>& pY )
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
    void RWPortTE10_1eg_FullWave::Calculate_DN( double DN[ 2 ][ 3 ] )
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
    void RWPortTE10_1eg_FullWave::Calculate_Edge_Lengths( Vector<double>& eLength )
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
	void RWPortTE10_1eg_FullWave::Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY )
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
    void RWPortTE10_1eg_FullWave::Input_Field_Vector( Vector<double>& Fx, Vector<double>& Fy, Vector<double>& cX, Vector<double>& cY )
	{
        // XY axis of the co-planar coordinate system 
        Vector<double> X_axis, Y_axis;

		Get_CoPlanar_XY_Axis( X_axis, Y_axis );

		// kc = ( pi / a )
		double pi = 3.141592653589793;
		double a  = Calculate_Width();
		double kc = pi / a;

		// Port vectors
		double r0[ 3 ]; Calculate_r0( r0 );
		double t [ 3 ]; Calculate_t ( t  ); 
		double b [ 3 ]; Calculate_b ( b  ); 

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

			double Xt = ( rX - r0[ 0 ] ) * t[ 0 ] + ( rY - r0[ 1 ] ) * t[ 1 ] + ( rZ - r0[ 2 ] ) * t[ 2 ]; 
			
			double Fb3D = sin( kc * Xt );  

			double Fx3D = b[ 0 ] * Fb3D; 
			double Fy3D = b[ 1 ] * Fb3D;   
			double Fz3D = b[ 2 ] * Fb3D;  

			Fx[ gp ] = Fx3D * X_axis[ 0 ] + Fy3D * X_axis[ 1 ] + Fz3D * X_axis[ 2 ]; 
			Fy[ gp ] = Fx3D * Y_axis[ 0 ] + Fy3D * Y_axis[ 1 ] + Fz3D * Y_axis[ 2 ];  
        }
	}

    //*******************************************************************************************************
    //* - Residual vector 
    //*******************************************************************************************************
    void RWPortTE10_1eg_FullWave::GetResidualVector( Vector< std::complex<double> >& ResidualVector ) 
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

		// TE10 propagation coefficient
		std::complex<double> cBetaMu = Calculate_BetaMuCoef();

		// Normalization constant
		std::complex<double> cN10 = Calculate_N10();

		// Integration constant
		std::complex<double> cIntgCte = -2.0 * Jacob * ( cBetaMu * cN10 );

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
    void RWPortTE10_1eg_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Initializing stiffness matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		// Negative TE10 propagation coefficient
		std::complex<double> mcBetaMu = -Calculate_BetaMuCoef(); 

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

				EdgeEMatrix[ i ][ j ] = mcBetaMu * sIntg_NiNj;
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
    //* - TE10 projection on a plane  
    //*******************************************************************************************************
	std::complex<double> RWPortTE10_1eg_FullWave::GetProjection()
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

		// Port vectors
		double r0[ 3 ]; Calculate_r0( r0 );
		double t [ 3 ]; Calculate_t ( t  ); 
		double b [ 3 ]; Calculate_b ( b  ); 

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

		// Fb in Gauss points
		Vector< std::complex<double> > cFb( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
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
		
			cFb[ gp ] = cFx3D * b[ 0 ] + cFy3D * b[ 1 ] + cFz3D * b[ 2 ];
		}
		
		if( mPotentialsOn )
		{
            cFb *= std::complex<double>( 0.0, (*mProperties)( FREQUENCY ) );
		}

		// kc = ( pi / a )
		double pi = 3.141592653589793;
		double a  = Calculate_Width();
		double kc = pi / a;

		// Integral over port surface
		std::complex<double> TE10Project( 0.0, 0.0 );

		for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            double alpha = cX[ gp ];
			double beta  = cY[ gp ];

			double rX = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
            double rY = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		    double rZ = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);

			double Xt = ( rX - r0[ 0 ] ) * t[ 0 ] + ( rY - r0[ 1 ] ) * t[ 1 ] + ( rZ - r0[ 2 ] ) * t[ 2 ]; 

			TE10Project += W[ gp ] * sin( kc * Xt ) * cFb[ gp ];  
        }

		return TE10Project;
	}

    //*******************************************************************************************************
    //* - Normalization term ( useful when using PMC symmetries )
    //*******************************************************************************************************
	std::complex<double> RWPortTE10_1eg_FullWave::GetNormalization()
	{
		// Gauss points
		Vector<double> cX, cY, W;

		int nGaussPoints = GaussPoints2D_Order12( cX, cY, W );
		
		// kc = ( pi / a )
		double pi = 3.141592653589793;
		double a  = Calculate_Width();
		double kc = pi / a;

		// Port vectors
		double r0[ 3 ]; Calculate_r0( r0 );
		double t [ 3 ]; Calculate_t ( t  ); 

        // Surface integral to normalize in case of symmetry
		double intg_Sin2 = 0.0;

		for( int gp=0; gp<nGaussPoints; gp++ ) 
        {
            double alpha = cX[ gp ];
			double beta  = cY[ gp ];

			double rX = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
            double rY = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		    double rZ = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);

			double Xt = ( rX - r0[ 0 ] ) * t[ 0 ] + ( rY - r0[ 1 ] ) * t[ 1 ] + ( rZ - r0[ 2 ] ) * t[ 2 ]; 

			intg_Sin2 += W[ gp ] * sin( kc * Xt ) * sin( kc * Xt );  
        }
		
		// Normalization coefficient
		std::complex<double> Normalization( intg_Sin2, 0.0 );

		std::complex<double> cN10 = Calculate_N10();

		Normalization *= cN10;
		
		return Normalization;
	}
} 