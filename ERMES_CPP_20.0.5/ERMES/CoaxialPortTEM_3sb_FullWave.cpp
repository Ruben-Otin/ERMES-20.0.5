
#include "CoaxialPortTEM_3sb_FullWave.h"
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
    void CoaxialPortTEM_3sb_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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
    double CoaxialPortTEM_3sb_FullWave::Calculate_Area() 
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
    std::complex<double> CoaxialPortTEM_3sb_FullWave::Calculate_KMuCoeff() 
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
    std::complex<double> CoaxialPortTEM_3sb_FullWave::Calculate_NrmCoeff() 
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
	void CoaxialPortTEM_3sb_FullWave::Calculate_r0( double* r0 ) 
	{
		Vector<double> cNF = (*mProperties)( COMPLEX_NEUMANN_FLOW ); 
		
		r0[ 0 ] = cNF[ 0 ];
		r0[ 1 ] = cNF[ 1 ];
	    r0[ 2 ] = cNF[ 2 ];
	}

    //*******************************************************************************************************
    //* - Coaxial port identification number
    //*******************************************************************************************************
	double CoaxialPortTEM_3sb_FullWave::GetPortNumber()
	{
		Vector<double> cIBC = (*mProperties)( COMPLEX_IBC );

		return cIBC[ 0 ];
	}

    //*******************************************************************************************************
    //* - Residual vector 
    //*******************************************************************************************************
    void CoaxialPortTEM_3sb_FullWave::GetResidualVector( Vector< std::complex<double> >& ResidualVector ) 
    {
        // Initializing residual vector
		ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );
		
		// If power is cero then return
		if( (*mProperties)( COMPLEX_IBC_2o )[ 0 ] == 0.0 )
		{
			return;
		}

		// Port vector
		double r0[ 3 ]; Calculate_r0( r0 ); 

		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		// Surface 3sb base
		Matrix<double> N;

		Lagrange2D_Ni_3sb( N, cX, cY );

		// Input fields on Gauss points
		Vector<double> Fx( numGaussPoints );
		Vector<double> Fy( numGaussPoints );
		Vector<double> Fz( numGaussPoints );

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

			Fx[ gp ] = InvRho2 * rhoX; 
			Fy[ gp ] = InvRho2 * rhoY;  
			Fz[ gp ] = InvRho2 * rhoZ; 
        }

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
		for( int i=0; i<mNumNodes; i++ )
        {
            double intg_FxNi = 0.0;
			double intg_FyNi = 0.0;
			double intg_FzNi = 0.0;

            for( int gp=0; gp<numGaussPoints; gp++ ) 
			{
				intg_FxNi += W[ gp ] * Fx[ gp ] * N[ i ][ gp ];
				intg_FyNi += W[ gp ] * Fy[ gp ] * N[ i ][ gp ];
				intg_FzNi += W[ gp ] * Fz[ gp ] * N[ i ][ gp ];
			}

            ResidualVector[ i               ] = cIntgCte * intg_FxNi;
            ResidualVector[ i +   mNumNodes ] = cIntgCte * intg_FyNi;
			ResidualVector[ i + 2*mNumNodes ] = cIntgCte * intg_FzNi;
        }  
    }

    //*******************************************************************************************************
    //* - Stiffness matrix
    //*******************************************************************************************************	
    void CoaxialPortTEM_3sb_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Initializing stiffness matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		// - j * ( k / mu ) 
		std::complex<double> mcKMu = -Calculate_KMuCoeff(); 
        
        // Surface integrals
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

        Matrix<double> N; 
		
		Lagrange2D_Ni_3sb( N, cX, cY );

        Matrix<double> intg_NiNj; 
		
		intg_NiNj.Resize( mNumNodes, mNumNodes, 0.0 );

        double Jacob = 2.0 * mArea;

        for( int i=0; i<mNumNodes; i++ )
        {
        	for( int j=0; j<mNumNodes; j++ )
        	{
                for( int gp=0; gp<numGaussPoints; gp++ ) 
                {
                    intg_NiNj[ i ][ j ] += Jacob * W[ gp ] * N[ i ][ gp ] * N[ j ][ gp ];
                }
        	}
        }

		// Stiffness matrix
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                StiffMatrix[ i               ][ j               ] = mcKMu * intg_NiNj[ i ][ j ];
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = mcKMu * intg_NiNj[ i ][ j ];
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = mcKMu * intg_NiNj[ i ][ j ];
            }
        }     
    }

    //*******************************************************************************************************
    //* - TEM projection on a plane 
    //*******************************************************************************************************
	std::complex<double> CoaxialPortTEM_3sb_FullWave::GetProjection()
	{
		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		// Surface 3sb base
		Matrix<double> N;

		Lagrange2D_Ni_3sb( N, cX, cY );

		// Port vector
		double r0[ 3 ]; Calculate_r0( r0 ); 

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
			std::complex<double> gpcFx( 0.0, 0.0 );
			std::complex<double> gpcFy( 0.0, 0.0 );
			std::complex<double> gpcFz( 0.0, 0.0 );
			
			for( int i=0; i<mNumNodes; i++ ) 
			{
				gpcFx += N[ i ][ gp ] * cFx[ i ];
				gpcFy += N[ i ][ gp ] * cFy[ i ];
				gpcFz += N[ i ][ gp ] * cFz[ i ];
			}

			cFr[ gp ] = InvRho[ gp ] * ( gpcFx * rhoX + gpcFy * rhoY + gpcFz * rhoZ );
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
	std::complex<double> CoaxialPortTEM_3sb_FullWave::GetNormalization()
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

			double InvRho2 = 1.0 / ( rhoX * rhoX + rhoY * rhoY + rhoZ * rhoZ );

			intg_InvRho2 += W[ gp ] * InvRho2;
        }
		
		// Normalization coefficient
		std::complex<double> Normalization( intg_InvRho2, 0.0 );

		std::complex<double> cNormCoeff = Calculate_NrmCoeff();

		Normalization *= cNormCoeff;

		return Normalization;
	}
} 