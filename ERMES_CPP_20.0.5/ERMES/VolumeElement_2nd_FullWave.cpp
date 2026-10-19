
#include "VolumeElement_2nd_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
    //******************************************************************************************************************************************
    //* - Global index of the DOFs
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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
    double VolumeElement_2nd_FullWave::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //******************************************************************************************************************************************
    // - 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::Calculate_DN( double DN[ 3 ][ 4 ] )
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
    // - dN/dx, dN/dy, dN/dz on cX, cY, cZ points
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::Derivatives_On_Points( Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz,
         						                            Vector<double>&   cX, Vector<double>&   cY, Vector<double>&   cZ )
    {
		int nPoints = cX.size();

		dNdx.Resize( mNumNodes, nPoints );
		dNdy.Resize( mNumNodes, nPoints );
		dNdz.Resize( mNumNodes, nPoints );
		
		double DN[ 3 ][ 4 ];

        Calculate_DN( DN );

		double dL1dx, dL2dx, dL3dx, dL4dx;
        double dL1dy, dL2dy, dL3dy, dL4dy;
		double dL1dz, dL2dz, dL3dz, dL4dz;

        dL1dx = DN[ 0 ][ 0 ], dL2dx = DN[ 0 ][ 1 ], dL3dx = DN[ 0 ][ 2 ], dL4dx = DN[ 0 ][ 3 ];
        dL1dy = DN[ 1 ][ 0 ], dL2dy = DN[ 1 ][ 1 ], dL3dy = DN[ 1 ][ 2 ], dL4dy = DN[ 1 ][ 3 ];
		dL1dz = DN[ 2 ][ 0 ], dL2dz = DN[ 2 ][ 1 ], dL3dz = DN[ 2 ][ 2 ], dL4dz = DN[ 2 ][ 3 ];

		// dNi/dx, dNi/dy, dNi/dz
        for( int p=0; p<nPoints; p++ ) 
        {
            double L2 = cX[ p ];
            double L3 = cY[ p ];
			double L4 = cZ[ p ];
            double L1 = 1.0 - L2 - L3 - L4;

			// Ni = Li * ( 2 * Li - 1.00 )
			dNdx[ 0 ][ p ] = dL1dx * ( 4.0 * L1 - 1.0 );
			dNdy[ 0 ][ p ] = dL1dy * ( 4.0 * L1 - 1.0 );
			dNdz[ 0 ][ p ] = dL1dz * ( 4.0 * L1 - 1.0 );
				  
			dNdx[ 1 ][ p ] = dL2dx * ( 4.0 * L2 - 1.0 );
			dNdy[ 1 ][ p ] = dL2dy * ( 4.0 * L2 - 1.0 );
			dNdz[ 1 ][ p ] = dL2dz * ( 4.0 * L2 - 1.0 );
				  
			dNdx[ 2 ][ p ] = dL3dx * ( 4.0 * L3 - 1.0 );
			dNdy[ 2 ][ p ] = dL3dy * ( 4.0 * L3 - 1.0 );
			dNdz[ 2 ][ p ] = dL3dz * ( 4.0 * L3 - 1.0 );
				  
			dNdx[ 3 ][ p ] = dL4dx * ( 4.0 * L4 - 1.0 );
			dNdy[ 3 ][ p ] = dL4dy * ( 4.0 * L4 - 1.0 );
			dNdz[ 3 ][ p ] = dL4dz * ( 4.0 * L4 - 1.0 );

			// Ni = 4 * Lj * Lk
			dNdx[ 4 ][ p ] = 4.0 * ( L2 * dL1dx + L1 * dL2dx );
			dNdy[ 4 ][ p ] = 4.0 * ( L2 * dL1dy + L1 * dL2dy );
			dNdz[ 4 ][ p ] = 4.0 * ( L2 * dL1dz + L1 * dL2dz );

			dNdx[ 5 ][ p ] = 4.0 * ( L2 * dL3dx + L3 * dL2dx );
			dNdy[ 5 ][ p ] = 4.0 * ( L2 * dL3dy + L3 * dL2dy );
			dNdz[ 5 ][ p ] = 4.0 * ( L2 * dL3dz + L3 * dL2dz );

			dNdx[ 6 ][ p ] = 4.0 * ( L1 * dL3dx + L3 * dL1dx );
			dNdy[ 6 ][ p ] = 4.0 * ( L1 * dL3dy + L3 * dL1dy );
			dNdz[ 6 ][ p ] = 4.0 * ( L1 * dL3dz + L3 * dL1dz );

			dNdx[ 7 ][ p ] = 4.0 * ( L1 * dL4dx + L4 * dL1dx );
			dNdy[ 7 ][ p ] = 4.0 * ( L1 * dL4dy + L4 * dL1dy );
			dNdz[ 7 ][ p ] = 4.0 * ( L1 * dL4dz + L4 * dL1dz );

			dNdx[ 8 ][ p ] = 4.0 * ( L2 * dL4dx + L4 * dL2dx );
			dNdy[ 8 ][ p ] = 4.0 * ( L2 * dL4dy + L4 * dL2dy );
			dNdz[ 8 ][ p ] = 4.0 * ( L2 * dL4dz + L4 * dL2dz );

			dNdx[ 9 ][ p ] = 4.0 * ( L3 * dL4dx + L4 * dL3dx );
			dNdy[ 9 ][ p ] = 4.0 * ( L3 * dL4dy + L4 * dL3dy );
			dNdz[ 9 ][ p ] = 4.0 * ( L3 * dL4dz + L4 * dL3dz );
        }		
	}

    //******************************************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
                                                                  Matrix< std::complex<double> >& EleStiffMatrix, 
                                                                  Vector< std::complex<double> >& ResidualVector )
    {
	    ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		if( mPotentialsOn ) 
		{
            for( int j=0; j<4; j++ ) 
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
    // - Element stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
        // Resizing stiffness matrix
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );	

		// Get non-regularized matrix
		if( mPeso == 0.0 ) 
		{
		    Get_Curl_Curl_Matrix( StiffMatrix ); 
		}
		// Get regularized matrix
		else
		{
		    Get_Curl_Div_Matrix( StiffMatrix ); 
		}
		
		// If AV-potentials are activated then add V contribution
		if( mPotentialsOn ) 
		{
		    Get_AV_VA_VV_Matrix( StiffMatrix ); 
		}
	}

    //******************************************************************************************************************************************
    // - Element regularized stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::Get_Curl_Div_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
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

		std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> cMuInv = 1.0 / cMu;

		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 2nd order basis functions 
		Matrix<double> N;

		Lagrange3D_Ni_2nd( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		Derivatives_On_Points( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Upper diagonal [Kxx], [Kyy], [Kzz]
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=i; j<mNumNodes; j++ )
			{
				double intg_CurlDiv = 0.0; 
				double intg_NiNj    = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_CurlDiv += W[ gp ] * ( dNdx[ i ][ gp ] * dNdx[ j ][ gp ] + dNdy[ i ][ gp ] * dNdy[ j ][ gp ] + dNdz[ i ][ gp ] * dNdz[ j ][ gp ] );
					intg_NiNj    += W[ gp ] * (    N[ i ][ gp ] *    N[ j ][ gp ] );
				}

				StiffMatrix[ i               ][ j               ] = ( cMuInv * intg_CurlDiv ) - ( w2cEps * intg_NiNj );
				StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( cMuInv * intg_CurlDiv ) - ( w2cEps * intg_NiNj );
				StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( cMuInv * intg_CurlDiv ) - ( w2cEps * intg_NiNj );
			}
		}

		// [Kxy], [Kxz], [Kyz] 
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=0; j<mNumNodes; j++ )
			{
				double intg_CurlDiv_XY = 0.0;
				double intg_CurlDiv_XZ = 0.0;
				double intg_CurlDiv_YZ = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_CurlDiv_XY += W[ gp ] * ( dNdx[ i ][ gp ] * dNdy[ j ][ gp ] - dNdy[ i ][ gp ] * dNdx[ j ][ gp ] );
					intg_CurlDiv_XZ += W[ gp ] * ( dNdx[ i ][ gp ] * dNdz[ j ][ gp ] - dNdz[ i ][ gp ] * dNdx[ j ][ gp ] );
					intg_CurlDiv_YZ += W[ gp ] * ( dNdy[ i ][ gp ] * dNdz[ j ][ gp ] - dNdz[ i ][ gp ] * dNdy[ j ][ gp ] );
				}
	            
				StiffMatrix[ i             ][ j +   mNumNodes ] = ( cMuInv * intg_CurlDiv_XY );
				StiffMatrix[ i             ][ j + 2*mNumNodes ] = ( cMuInv * intg_CurlDiv_XZ );
				StiffMatrix[ i + mNumNodes ][ j + 2*mNumNodes ] = ( cMuInv * intg_CurlDiv_YZ );
			}
		}

		// Lower diagonal ( [EE] | [AA] )
        for( int i=0; i<30; i++ )
        {
            for( int j=i+1; j<30; j++ )
            {
                StiffMatrix[ j ][ i ] = StiffMatrix[ i ][ j ];
            }
        }        
	}

    //******************************************************************************************************************************************
    // - Element non-regularized stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::Get_Curl_Curl_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
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

		std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> cMuInv = 1.0 / cMu;

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // Basis functions 2nd order 
		Matrix<double> N;

		Lagrange3D_Ni_2nd( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		Derivatives_On_Points( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Upper diagonal [Kxx], [Kyy], [Kzz]
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=i; j<mNumNodes; j++ )
			{
				double intg_CurlCurl_XX = 0.0; 
				double intg_CurlCurl_YY = 0.0; 
				double intg_CurlCurl_ZZ = 0.0; 
				double intg_NiNj        = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_CurlCurl_XX += W[ gp ] * ( dNdy[ i ][ gp ] * dNdy[ j ][ gp ] + dNdz[ i ][ gp ] * dNdz[ j ][ gp ] );
					intg_CurlCurl_YY += W[ gp ] * ( dNdx[ i ][ gp ] * dNdx[ j ][ gp ] + dNdz[ i ][ gp ] * dNdz[ j ][ gp ] );
					intg_CurlCurl_ZZ += W[ gp ] * ( dNdx[ i ][ gp ] * dNdx[ j ][ gp ] + dNdy[ i ][ gp ] * dNdy[ j ][ gp ] );
					intg_NiNj        += W[ gp ] * (    N[ i ][ gp ] *    N[ j ][ gp ] );
				}

				StiffMatrix[ i               ][ j               ] = ( cMuInv * intg_CurlCurl_XX ) - ( w2cEps * intg_NiNj );
				StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( cMuInv * intg_CurlCurl_YY ) - ( w2cEps * intg_NiNj );
				StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( cMuInv * intg_CurlCurl_ZZ ) - ( w2cEps * intg_NiNj );
			}
		}

		// [Kxy], [Kxz], [Kyz] 
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=0; j<mNumNodes; j++ )
			{
				double intg_CurlCurl_XY = 0.0;
				double intg_CurlCurl_XZ = 0.0;
				double intg_CurlCurl_YZ = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_CurlCurl_XY += W[ gp ] * ( dNdy[ i ][ gp ] * dNdx[ j ][ gp ] );
					intg_CurlCurl_XZ += W[ gp ] * ( dNdz[ i ][ gp ] * dNdx[ j ][ gp ] );
					intg_CurlCurl_YZ += W[ gp ] * ( dNdz[ i ][ gp ] * dNdy[ j ][ gp ] );
				}
	            
				StiffMatrix[ i             ][ j +   mNumNodes ] = ( -cMuInv * intg_CurlCurl_XY );
				StiffMatrix[ i             ][ j + 2*mNumNodes ] = ( -cMuInv * intg_CurlCurl_XZ );
				StiffMatrix[ i + mNumNodes ][ j + 2*mNumNodes ] = ( -cMuInv * intg_CurlCurl_YZ );
			}
		}

		// Lower diagonal ( [EE] | [AA] )
        for( int i=0; i<30; i++ )
        {
            for( int j=i+1; j<30; j++ )
            {
                StiffMatrix[ j ][ i ] = StiffMatrix[ i ][ j ];
            }
        }   
	}

    //******************************************************************************************************************************************
    // - AV-potentials contribution [AV], [VA], [VV] to stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::Get_AV_VA_VV_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
        // Element material properties
		double wf = (*mProperties)( FREQUENCY );

        double w2 = wf * wf;
        double eo = 8.8541878176e-12;
	    
		double sgm_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
		double sgm_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        double eps_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double eps_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		std::complex<double> w2cEps;

		if( mComplxFreqOn )
		{
		    std::complex<double> jUnit( 0.0, 1.0 );
		    
		    std::complex<double> sLaplace  = mComplxFreq;
		    std::complex<double> sLaplace2 = sLaplace * sLaplace;
		    
            w2cEps = sLaplace * ( sgm_real + jUnit * sgm_imag ) - sLaplace2 * ( eps_real + jUnit * eps_imag );		 
		}
		else
		{
            w2cEps = std::complex<double>( w2 * eps_real - wf * sgm_imag, w2 * eps_imag + wf * sgm_real );
		}

		std::complex<double> nw2cEpsVol = -mVolume * w2cEps;

		Vector< std::complex<double> > nw2cEps_intgVolNi( mNumNodes );

		for( int i=0; i< 4; i++ ) nw2cEps_intgVolNi[ i ] = -nw2cEpsVol / 20.0;
		for( int i=4; i<10; i++ ) nw2cEps_intgVolNi[ i ] =  nw2cEpsVol /  5.0;

		// Derivatives 1st order basis
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		// [Kvv]
		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
				StiffMatrix[ i + 30 ][ j + 30 ] = nw2cEpsVol * ( DN[ 0 ][ i ]*DN[ 0 ][ j ] + DN[ 1 ][ i ]*DN[ 1 ][ j ] + DN[ 2 ][ i ]*DN[ 2 ][ j ] );
        	}
        }

		// [Kxv], [Kyv], [Kzv]
		for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
        		StiffMatrix[ i               ][ j + 3*mNumNodes ] = nw2cEps_intgVolNi[ i ] * DN[ 0 ][ j ];
        		StiffMatrix[ i +   mNumNodes ][ j + 3*mNumNodes ] = nw2cEps_intgVolNi[ i ] * DN[ 1 ][ j ];
        		StiffMatrix[ i + 2*mNumNodes ][ j + 3*mNumNodes ] = nw2cEps_intgVolNi[ i ] * DN[ 2 ][ j ];
        	}
        }

        // [Kvx], [Kvy], [Kvz]
		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<mNumNodes; j++ ) 
        	{
        		StiffMatrix[ i + 3*mNumNodes ][ j               ] = nw2cEps_intgVolNi[ j ] * DN[ 0 ][ i ];
        		StiffMatrix[ i + 3*mNumNodes ][ j +   mNumNodes ] = nw2cEps_intgVolNi[ j ] * DN[ 1 ][ i ];
        		StiffMatrix[ i + 3*mNumNodes ][ j + 2*mNumNodes ] = nw2cEps_intgVolNi[ j ] * DN[ 2 ][ i ];
        	}
        }
	}

    //******************************************************************************************************************************************
    // - Rotational on nodes
    //******************************************************************************************************************************************
	void VolumeElement_2nd_FullWave::Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes )
	{
		// Nodes natural coordinates
		Vector<double> cX( mNumNodes ); 
		Vector<double> cY( mNumNodes );
		Vector<double> cZ( mNumNodes );

		// Corners 
		cX[ 0 ] = 0.0; cY[ 0 ] = 0.0; cZ[ 0 ] = 0.0;
		cX[ 1 ] = 1.0; cY[ 1 ] = 0.0; cZ[ 1 ] = 0.0;
		cX[ 2 ] = 0.0; cY[ 2 ] = 1.0; cZ[ 2 ] = 0.0;
		cX[ 3 ] = 0.0; cY[ 3 ] = 0.0; cZ[ 3 ] = 1.0;
		
		// Edges 
		cX[ 4 ] = 0.5; cY[ 4 ] = 0.0; cZ[ 4 ] = 0.0; 
		cX[ 5 ] = 0.5; cY[ 5 ] = 0.5; cZ[ 5 ] = 0.0; 
		cX[ 6 ] = 0.0; cY[ 6 ] = 0.5; cZ[ 6 ] = 0.0; 
		cX[ 7 ] = 0.0; cY[ 7 ] = 0.0; cZ[ 7 ] = 0.5; 
		cX[ 8 ] = 0.5; cY[ 8 ] = 0.0; cZ[ 8 ] = 0.5;
		cX[ 9 ] = 0.0; cY[ 9 ] = 0.5; cZ[ 9 ] = 0.5; 

		// Derivatives on nodes
        Matrix<double> dNdx, dNdy, dNdz;

		Derivatives_On_Points( dNdx, dNdy, dNdz, cX, cY, cZ );
        
		// Get E / A in nodes
		Vector< std::complex<double> > cFx( mNumNodes );
		Vector< std::complex<double> > cFy( mNumNodes );
		Vector< std::complex<double> > cFz( mNumNodes );
		    
		if( mPotentialsOn )
		{
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        cFx[ n ] = (*mProperties)( cAx, *mNodes[ n ] );
		        cFy[ n ] = (*mProperties)( cAy, *mNodes[ n ] );
		        cFz[ n ] = (*mProperties)( cAz, *mNodes[ n ] );
		    }
		}
		else
		{
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        cFx[ n ] = (*mProperties)( cEx, *mNodes[ n ] );
		        cFy[ n ] = (*mProperties)( cEy, *mNodes[ n ] );
		        cFz[ n ] = (*mProperties)( cEz, *mNodes[ n ] );
		    }		
		}

		// Rotational on nodes
		rot_OnNodes.resize( mNumNodes );

		for( int n=0; n<mNumNodes; n++ )
		{
		    Vector< std::complex<double> > cRot( 3, std::complex<double>( 0.0, 0.0 ) );
		 
		    for( int i=0; i<mNumNodes; i++ )
		    {
		        cRot[ 0 ] += cFz[ i ]*dNdy[ i ][ n ] - cFy[ i ]*dNdz[ i ][ n ];
			    cRot[ 1 ] += cFx[ i ]*dNdz[ i ][ n ] - cFz[ i ]*dNdx[ i ][ n ];
			    cRot[ 2 ] += cFy[ i ]*dNdx[ i ][ n ] - cFx[ i ]*dNdy[ i ][ n ];
		    }					   

		    rot_OnNodes[ n ] = cRot;
		}
	}
	  
    //******************************************************************************************************************************************
    // - Rotational on Gauss points
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs )
	{
		// Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );

		// Derivatives on Gauss points
        Matrix<double> dNdx, dNdy, dNdz;

		Derivatives_On_Points( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Get E / A in nodes
		Vector< std::complex<double> > cFx( mNumNodes );
		Vector< std::complex<double> > cFy( mNumNodes );
		Vector< std::complex<double> > cFz( mNumNodes );
		    
		if( mPotentialsOn )
		{
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        cFx[ n ] = (*mProperties)( cAx, *mNodes[ n ] );
		        cFy[ n ] = (*mProperties)( cAy, *mNodes[ n ] );
		        cFz[ n ] = (*mProperties)( cAz, *mNodes[ n ] );
		    }
		}
		else
		{
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        cFx[ n ] = (*mProperties)( cEx, *mNodes[ n ] );
		        cFy[ n ] = (*mProperties)( cEy, *mNodes[ n ] );
		        cFz[ n ] = (*mProperties)( cEz, *mNodes[ n ] );
		    }		
		}

		// Rotational on Gauss points
		rot_OnGPs.resize( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
            Vector< std::complex<double> > cRot( 3, std::complex<double>( 0.0, 0.0 ) );
		 
		    for( int i=0; i<mNumNodes; i++ )
		    {
		        cRot[ 0 ] += cFz[ i ]*dNdy[ i ][ gp ] - cFy[ i ]*dNdz[ i ][ gp ];
			    cRot[ 1 ] += cFx[ i ]*dNdz[ i ][ gp ] - cFz[ i ]*dNdx[ i ][ gp ];
			    cRot[ 2 ] += cFy[ i ]*dNdx[ i ][ gp ] - cFx[ i ]*dNdy[ i ][ gp ];
		    }

		    rot_OnGPs[ gp ] = cRot;
		}
	}

    //******************************************************************************************************************************************
    // - Electric field on nodes
    //******************************************************************************************************************************************
	void VolumeElement_2nd_FullWave::Calculate_E_field_OnNodes( cVector2Type& Ef_OnNodes )
	{
		Ef_OnNodes.resize( mNumNodes );

		for( int n=0; n<mNumNodes; n++ )
		{
		    Ef_OnNodes[ n ].resize( 3 );
		}

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

            for( int i=0; i<4; i++ )
            {
                cGradVs[ 0 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 0 ][ i ];
                cGradVs[ 1 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 1 ][ i ];
                cGradVs[ 2 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 2 ][ i ];
            }

		    for( int n=0; n<mNumNodes; n++ )
		    {
		        Ef_OnNodes[ n ][ 0 ] = jw * ( (*mProperties)( cAx, *mNodes[ n ] ) + cGradVs[ 0 ] );
		    	Ef_OnNodes[ n ][ 1 ] = jw * ( (*mProperties)( cAy, *mNodes[ n ] ) + cGradVs[ 1 ] );
		    	Ef_OnNodes[ n ][ 2 ] = jw * ( (*mProperties)( cAz, *mNodes[ n ] ) + cGradVs[ 2 ] );
		    }	
		}
		else
		{
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        Ef_OnNodes[ n ][ 0 ] = (*mProperties)( cEx, *mNodes[ n ] );
		    	Ef_OnNodes[ n ][ 1 ] = (*mProperties)( cEy, *mNodes[ n ] );
		    	Ef_OnNodes[ n ][ 2 ] = (*mProperties)( cEz, *mNodes[ n ] );
		    }			
		}
	}

    //******************************************************************************************************************************************
    // - Electric field on Gauss points
    //******************************************************************************************************************************************
    void VolumeElement_2nd_FullWave::Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs )
	{
		// Get E field on nodes
		cVector2Type cE_OnNodes;
		
		Calculate_E_field_OnNodes( cE_OnNodes );
		
		// Get GiD internal Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );

        // Get basis functions values on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2nd( N, cX, cY, cZ );

	    // Electric field on Gauss points
		Ef_OnGPs.resize( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
			Vector< std::complex<double> > cE_OnGP( 3, std::complex<double>( 0.0, 0.0 ) );

			for( int n=0; n<mNumNodes; n++ )
		    {       
				cE_OnGP[ 0 ] += cE_OnNodes[ n ][ 0 ] * N[ n ][ gp ];
			    cE_OnGP[ 1 ] += cE_OnNodes[ n ][ 1 ] * N[ n ][ gp ];
			    cE_OnGP[ 2 ] += cE_OnNodes[ n ][ 2 ] * N[ n ][ gp ];
			}

			Ef_OnGPs[ gp ] = cE_OnGP;
		}	
	}

    //******************************************************************************************************************************************
    // - Inner GiD Gauss points 
    //******************************************************************************************************************************************
    int VolumeElement_2nd_FullWave::GetInnerGiDGaussPoints( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, int numResultsOnGPs )
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