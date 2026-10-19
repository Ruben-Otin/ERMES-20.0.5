
#include "VolumeElement_1st_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{ 
    //***************************************************************************************************************************
    // - Global index of the DOFs
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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

    //***************************************************************************************************************************
    // - Volume of the element
    //***************************************************************************************************************************
    double VolumeElement_1st_FullWave::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //***************************************************************************************************************************
    // - 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::Calculate_DN( double DN[ 3 ][ 4 ] )
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

    //***************************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
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

    //***************************************************************************************************************************
    // - Stiffness matrix
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
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
		    Get_Curl_Div_Matrix ( StiffMatrix ); 
		}
		
		// If AV-potentials are activated then add V contribution
		if( mPotentialsOn ) 
		{
		    Get_AV_VA_VV_Matrix ( StiffMatrix ); 
		}
	}

    //***************************************************************************************************************************
    // - Regularized stiffness matrix
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::Get_Curl_Div_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
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

		std::complex<double> w2cEpsVol = mVolume * w2cEps;
	    std::complex<double> cMuInvVol = mVolume / cMu;
		
        std::complex<double> NiNj_Matrix[ 4 ][ 4 ];

        for( int i=0; i<4; i++ )
        {
        	for( int j=0; j<4; j++ )
        	{
        	    if( i!=j )
				{
					NiNj_Matrix[ i ][ j ] = w2cEpsVol / 20.0;
				}
        	    else   
				{
					NiNj_Matrix[ i ][ j ] = w2cEpsVol / 10.0;	
				}
        	}
        }
        
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        std::complex<double> eK_ij;
        
        for( int i=0; i<4; i++ )
        {
        	for( int j=0; j<4; j++ ) 
        	{
        		// [Kxx], [Kyy], [Kzz]  
        		eK_ij = cMuInvVol * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] ) - NiNj_Matrix[ i ][ j ];
        
        		StiffMatrix[ i     ][ j     ] = eK_ij;
        		StiffMatrix[ i + 4 ][ j + 4 ] = eK_ij;
        		StiffMatrix[ i + 8 ][ j + 8 ] = eK_ij;
        
        		// [Kxy], [Kyx]  
        		eK_ij = cMuInvVol * ( DN[ 0 ][ i ] * DN[ 1 ][ j ] - DN[ 1 ][ i ] * DN[ 0 ][ j ] );
        
        		StiffMatrix[ i     ][ j + 4 ] =  eK_ij;
        		StiffMatrix[ i + 4 ][ j     ] = -eK_ij;
        
        		// [Kxz], [Kzx] 
        		eK_ij = cMuInvVol * ( DN[ 0 ][ i ] * DN[ 2 ][ j ] - DN[ 2 ][ i ] * DN[ 0 ][ j ] );             
        
        		StiffMatrix[ i     ][ j + 8 ] =  eK_ij;
        		StiffMatrix[ i + 8 ][ j     ] = -eK_ij;
        
        		// [Kyz], [Kzy]  
        		eK_ij = cMuInvVol * ( DN[ 1 ][ i ] * DN[ 2 ][ j ] - DN[ 2 ][ i ] * DN[ 1 ][ j ] );
            
        		StiffMatrix[ i + 4 ][ j + 8 ] =  eK_ij;
        		StiffMatrix[ i + 8 ][ j + 4 ] = -eK_ij;
        	}
        }
	}

    //***************************************************************************************************************************
    // - Non-regularized stiffness matrix
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::Get_Curl_Curl_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{
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

		std::complex<double> w2cEpsVol = mVolume * w2cEps;
	    std::complex<double> cMuInvVol = mVolume / cMu;
		
        std::complex<double> NiNj_Matrix[ 4 ][ 4 ];

        for( int i=0; i<4; i++ )
        {
        	for( int j=0; j<4; j++ )
        	{
        	    if( i!=j ) 
				{
					NiNj_Matrix[ i ][ j ] = w2cEpsVol / 20.0;
				}
        	    else 
				{
					NiNj_Matrix[ i ][ j ] = w2cEpsVol / 10.0;	
				}
        	}
        }
        
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		// [Kxx], [Kyy], [Kzz]  
		for( int i=0; i<4; i++ )
        {
            for( int j=i; j<4; j++ ) 
        	{
        		StiffMatrix[ i     ][ j     ] = cMuInvVol * ( DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] ) - NiNj_Matrix[ i ][ j ];
        		StiffMatrix[ i + 4 ][ j + 4 ] = cMuInvVol * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] ) - NiNj_Matrix[ i ][ j ];
        		StiffMatrix[ i + 8 ][ j + 8 ] = cMuInvVol * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] ) - NiNj_Matrix[ i ][ j ];
        	}
        }

        // [Kxy], [Kxz], [Kyz]
        for( int i=0; i<4; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
        		StiffMatrix[ i     ][ j + 4 ] = -cMuInvVol * ( DN[ 1 ][ i ] * DN[ 0 ][ j ] );
        		StiffMatrix[ i     ][ j + 8 ] = -cMuInvVol * ( DN[ 2 ][ i ] * DN[ 0 ][ j ] );
        		StiffMatrix[ i + 4 ][ j + 8 ] = -cMuInvVol * ( DN[ 2 ][ i ] * DN[ 1 ][ j ] );
        	}
        }

		// Lower diagonal
		for( int i=0; i<12; i++ )
		{
			for( int j=i+1; j<12; j++ )
			{
				StiffMatrix[ j ][ i ] = StiffMatrix[ i ][ j ];
			}
		} 
	}

    //***************************************************************************************************************************
    // - AV-potentials contribution [AV], [VA], [VV] to stiffness matrix
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::Get_AV_VA_VV_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
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

		std::complex<double> nw2cEpsVol     = -mVolume * w2cEps;
		std::complex<double> nw2cEpsVolDiv4 =  nw2cEpsVol / 4.0;

        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
				// [Kxv], [Kyv], [Kzv]
        		StiffMatrix[ i      ][ j + 12 ] = nw2cEpsVolDiv4 * DN[ 0 ][ j ];
        		StiffMatrix[ i + 4  ][ j + 12 ] = nw2cEpsVolDiv4 * DN[ 1 ][ j ];
        		StiffMatrix[ i + 8  ][ j + 12 ] = nw2cEpsVolDiv4 * DN[ 2 ][ j ];
																	   
 		 	 	// [Kvx], [Kvy], [Kvz]								   
           		StiffMatrix[ i + 12 ][ j      ] = nw2cEpsVolDiv4 * DN[ 0 ][ i ];
        		StiffMatrix[ i + 12 ][ j + 4  ] = nw2cEpsVolDiv4 * DN[ 1 ][ i ];
        		StiffMatrix[ i + 12 ][ j + 8  ] = nw2cEpsVolDiv4 * DN[ 2 ][ i ];

				// [Kvv]
				StiffMatrix[ i + 12 ][ j + 12 ] = nw2cEpsVol * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );
        	}
        }
	}

    //***************************************************************************************************************************
    // - Rotational on nodes
    //***************************************************************************************************************************
	void VolumeElement_1st_FullWave::Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes )
	{
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

		// Derivatives on nodes
		double DN[3][4]; Calculate_DN( DN );

		// Calculate rotational 
		Vector< std::complex<double> > cRot( 3, std::complex<double>( 0.0, 0.0 ) );

		for( int i=0; i<mNumNodes; i++ )
		{
		    cRot[ 0 ] += cFz[ i ] * DN[ 1 ][ i ] - cFy[ i ] * DN[ 2 ][ i ];
		    cRot[ 1 ] += cFx[ i ] * DN[ 2 ][ i ] - cFz[ i ] * DN[ 0 ][ i ];
		    cRot[ 2 ] += cFy[ i ] * DN[ 0 ][ i ] - cFx[ i ] * DN[ 1 ][ i ];
		}	

		// Rotational on nodes
		rot_OnNodes.resize( mNumNodes );

		for( int n=0; n<mNumNodes; n++ )
		{
		    rot_OnNodes[ n ] = cRot;
		}
	}
	  
    //***************************************************************************************************************************
    // - Rotational on Gauss points
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs )
	{
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

		// Calculate derivatives on Gauss points
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		// Calculate rotational 
		Vector< std::complex<double> > cRot( 3, std::complex<double>( 0.0, 0.0 ) );

		for( int i=0; i<mNumNodes; i++ )
		{
		    cRot[ 0 ] += cFz[ i ] * DN[ 1 ][ i ] - cFy[ i ] * DN[ 2 ][ i ];
		    cRot[ 1 ] += cFx[ i ] * DN[ 2 ][ i ] - cFz[ i ] * DN[ 0 ][ i ];
		    cRot[ 2 ] += cFy[ i ] * DN[ 0 ][ i ] - cFx[ i ] * DN[ 1 ][ i ];
		}	

		// Number of Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );

		// Rotational on Gauss points
		rot_OnGPs.resize( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
		    rot_OnGPs[ gp ] = cRot;
		}
	}

    //***************************************************************************************************************************
    // - Electric field on nodes
    //***************************************************************************************************************************
	void VolumeElement_1st_FullWave::Calculate_E_field_OnNodes( cVector2Type& Ef_OnNodes )
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

    //***************************************************************************************************************************
    // - Electric field on Gauss points
    //***************************************************************************************************************************
    void VolumeElement_1st_FullWave::Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs )
	{
		// Get E field on nodes
		cVector2Type cE_OnNodes;

		Calculate_E_field_OnNodes( cE_OnNodes );
		
		// Get GiD internal Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );

        // Get basis functions values on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

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

    //***************************************************************************************************************************
    // - Inner GiD Gauss points 
    //***************************************************************************************************************************
    int VolumeElement_1st_FullWave::GetInnerGiDGaussPoints( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, int numResultsOnGPs )
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