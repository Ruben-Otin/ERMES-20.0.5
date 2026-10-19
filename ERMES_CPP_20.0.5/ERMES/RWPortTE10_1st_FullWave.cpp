
#include "RWPortTE10_1st_FullWave.h"
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
    void RWPortTE10_1st_FullWave::GetEquationIdVector( Vector<int>& EquationId )   
    {
        EquationId.resize( mNumDofs );

		if( mPotentialsOn )
		{
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                EquationId[ i               ] = mNodes[ i ]->pDofcAx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcAy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcAz()->EquationId();
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
    double RWPortTE10_1st_FullWave::Calculate_Area() 
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
    std::complex<double> RWPortTE10_1st_FullWave::Calculate_BetaMuCoef() 
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
    std::complex<double> RWPortTE10_1st_FullWave::Calculate_N10() 
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
	double RWPortTE10_1st_FullWave::Calculate_Width()
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
	double RWPortTE10_1st_FullWave::Calculate_Height()
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
	void RWPortTE10_1st_FullWave::Calculate_b( double* b ) 
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
	void RWPortTE10_1st_FullWave::Calculate_t( double* t ) 
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
	void RWPortTE10_1st_FullWave::Calculate_r0( double* r0 ) 
	{
		Vector<double> cNF = (*mProperties)( COMPLEX_NEUMANN_FLOW );
		
		r0[ 0 ] = cNF[ 0 ];
		r0[ 1 ] = cNF[ 1 ];
	    r0[ 2 ] = cNF[ 2 ];
	}

	//*******************************************************************************************************
    //* - Port identification number
    //*******************************************************************************************************
	double RWPortTE10_1st_FullWave::GetPortNumber()
	{
		Vector<double> cIBC = (*mProperties)( COMPLEX_IBC );

		return cIBC[ 0 ];
	}

    //*******************************************************************************************************
    //* - Residual vector 
    //*******************************************************************************************************
    void RWPortTE10_1st_FullWave::GetResidualVector( Vector< std::complex<double> >& ResidualVector ) 
    {
        // Initializing residual vector
		ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		// If power is cero then return
		if( (*mProperties)( COMPLEX_IBC_2o )[ 0 ] == 0.0 )
		{
			return;
		}

		// kc = ( pi / a )
		double pi = 3.141592653589793;
		double a  = Calculate_Width();
		double kc = pi / a;

		// Port vectors
		double r0[ 3 ]; Calculate_r0( r0 );
		double t [ 3 ]; Calculate_t ( t  ); 
		double b [ 3 ]; Calculate_b ( b  ); 
		 
		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		// Surface 1st order base
		Matrix<double> N;

		Lagrange2D_Ni_1st( N, cX, cY );

		// Input fields on Gauss points
		Vector<double> Fb( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            double alpha = cX[ gp ];
			double beta  = cY[ gp ];

			double rX = alpha * ( X(2) - X(1) ) + beta * ( X(3) - X(1) ) + X(1);
            double rY = alpha * ( Y(2) - Y(1) ) + beta * ( Y(3) - Y(1) ) + Y(1);
		    double rZ = alpha * ( Z(2) - Z(1) ) + beta * ( Z(3) - Z(1) ) + Z(1);

			double Xt = ( rX - r0[ 0 ] ) * t[ 0 ] + ( rY - r0[ 1 ] ) * t[ 1 ] + ( rZ - r0[ 2 ] ) * t[ 2 ]; 
			
			Fb[ gp ] = sin( kc * Xt );  
        }

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
		for( int i=0; i<mNumNodes; i++ )
        {
            double intg_FbNi = 0.0;

            for( int gp=0; gp<numGaussPoints; gp++ ) 
			{
				intg_FbNi += W[ gp ] * Fb[ gp ] * N[ i ][ gp ];
			}

            ResidualVector[ i               ] = cIntgCte * ( b[ 0 ] * intg_FbNi );
            ResidualVector[ i +   mNumNodes ] = cIntgCte * ( b[ 1 ] * intg_FbNi );
			ResidualVector[ i + 2*mNumNodes ] = cIntgCte * ( b[ 2 ] * intg_FbNi );
        }  
    }

    //*******************************************************************************************************
    //* - Stiffness matrix
    //*******************************************************************************************************	
    void RWPortTE10_1st_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Initializing matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		// Negative TE10 propagation coefficient
		std::complex<double> mcBetaMu = -Calculate_BetaMuCoef(); 

		// NiNj surface integral
        Matrix<double> intg_NiNj; 
		
		intg_NiNj.Resize( mNumNodes, mNumNodes, 0.0 );

        for( int i=0; i<mNumNodes; i++ )
        {
        	for( int j=0; j<mNumNodes; j++ )
        	{
        	    if( i != j )
                {
					intg_NiNj[ i ][ j ] = mArea / 12.0;
                }
        	    else
                {
					intg_NiNj[ i ][ j ] = mArea / 6.0;		         
                }
        	}
        }

        // Stiffness matrix
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                StiffMatrix[ i               ][ j               ] = mcBetaMu * intg_NiNj[ i ][ j ];
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = mcBetaMu * intg_NiNj[ i ][ j ];
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = mcBetaMu * intg_NiNj[ i ][ j ];
            }
        }
    }

    //*******************************************************************************************************
    //* - TE10 projection on a plane 
    //*******************************************************************************************************
	std::complex<double> RWPortTE10_1st_FullWave::GetProjection()
	{
		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		// Lagrange 1st order base
		Matrix<double> N;

		Lagrange2D_Ni_1st( N, cX, cY );

		// Port vectors
		double r0[ 3 ]; Calculate_r0( r0 );
		double t [ 3 ]; Calculate_t ( t  ); 
		double b [ 3 ]; Calculate_b ( b  ); 

        // F on nodes
		Vector< std::complex<double> > cFx( mNumNodes );
		Vector< std::complex<double> > cFy( mNumNodes );
		Vector< std::complex<double> > cFz( mNumNodes );

		if( mPotentialsOn )
		{
		    for( int i=0; i<mNumNodes; i++ )
		    {
				cFx[ i ] = (*mProperties)( cAx, *mNodes[ i ] );
				cFy[ i ] = (*mProperties)( cAy, *mNodes[ i ] );
				cFz[ i ] = (*mProperties)( cAz, *mNodes[ i ] );
		    }
		}
		else
		{
		    for( int i=0; i<mNumNodes; i++ )
		    {
				cFx[ i ] = (*mProperties)( cEx, *mNodes[ i ] );
				cFy[ i ] = (*mProperties)( cEy, *mNodes[ i ] );
				cFz[ i ] = (*mProperties)( cEz, *mNodes[ i ] );
		    }
		}

		// Fb in Gauss points
		Vector< std::complex<double> > cFb( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
			std::complex<double> gpcFx( 0.0, 0.0 );
			std::complex<double> gpcFy( 0.0, 0.0 );
			std::complex<double> gpcFz( 0.0, 0.0 );
			
			for( int i=0; i<mNumNodes; i++ ) 
			{
				gpcFx += N[ i ][ gp ] * cFx[ i ];
				gpcFy += N[ i ][ gp ] * cFy[ i ];
				gpcFz += N[ i ][ gp ] * cFz[ i ];
			}
		
			cFb[ gp ] = b[ 0 ] * gpcFx + b[ 1 ] * gpcFy + b[ 2 ] * gpcFz;
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
	std::complex<double> RWPortTE10_1st_FullWave::GetNormalization()
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