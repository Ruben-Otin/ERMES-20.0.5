
#include "VolumeElement_1st_ColdPlasma.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
    //********************************************************************************************************************************
    // - Global index of the DOFs
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::GetEquationIdVector( Vector<int>& EquationId )
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

    //********************************************************************************************************************************
    // - Volume of the element
    //********************************************************************************************************************************
    double VolumeElement_1st_ColdPlasma::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //********************************************************************************************************************************
    //* - Set element as IHL or plasma 
    //********************************************************************************************************************************	
    void VolumeElement_1st_ColdPlasma::SetPlasmaIHL() 
    { 
        // Material properties
        double sg_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
        double sg_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        double ep_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY );
        double ep_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY );

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY );
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY );

        // If the material is vaccum then it is considered as plasma. 
        // We only consider the first one written in the plasma file *-11.dat
        // This allows the use of different plasmas volumes (for visualization). 
        if( ( sg_real == 0.0 ) && ( sg_imag == 0.0 ) &&
            ( ep_real == 1.0 ) && ( ep_imag == 0.0 ) && 
            ( mu_real == 1.0 ) && ( mu_imag == 0.0 )  )
        {
            mIs_IHL = false;        
        }
        else
        {
            mIs_IHL = true;        
        }
    }

    //********************************************************************************************************************************
    // - 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::Calculate_DN( double DN[ 3 ][ 4 ] )
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

    //********************************************************************************************************************************
    //* - Permittivity tensor on nodes 
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::PermittivityTensor_OnNodes( Vector< Matrix< std::complex<double> > >& TEn )
    {
        // Permittivity tensor on nodes
        TEn.resize( mNumNodes );
        
        // Constants
        double wf = (*mProperties)( FREQUENCY );
        double eo = 8.8541878176e-12;
        
        std::complex<double> cZero( 0.0, 0.0 );
        std::complex<double> cUnit( 0.0, 1.0 );

        if( mIs_IHL )
        {
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                TEn[ i ].Resize( 3, 3 );

                double sgm_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
                double sgm_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

                double eps_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
                double eps_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

                std::complex<double> cEps( eps_real - ( sgm_imag / wf ), eps_imag + ( sgm_real / wf ) );

                // IHL permittivity tensor 
                Matrix< std::complex<double> > K( 3, 3 );

                K[ 0 ][ 0 ] = cEps ; K[ 0 ][ 1 ] = cZero; K[ 0 ][ 2 ] = cZero; 
                K[ 1 ][ 0 ] = cZero; K[ 1 ][ 1 ] = cEps ; K[ 1 ][ 2 ] = cZero; 
                K[ 2 ][ 0 ] = cZero; K[ 2 ][ 1 ] = cZero; K[ 2 ][ 2 ] = cEps ;

                TEn[ i ] = K;
            }
        }
        else
        {
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                mpColdPlasma->Get_PermittivityTensor_OnNode( TEn[ i ], mNodes[ i ], wf );
            }
        }
    }

    //********************************************************************************************************************************
    //* - Permittivity tensor on Gauss points ( extrapolating TEn to GPs with N )
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp,  
                                                                         Vector< Matrix< std::complex<double> > >& TEn , 
                                                                         Matrix< double >& N )
    {
        // Number of Gauss points
        int numGaussPoints = N[ 0 ].size();

        // Initializing vector of 3x3 permittivity tensors 
        TEgp.resize( numGaussPoints );

        // Permittivity tensor on Gauss points  
        for( int gp=0; gp<numGaussPoints; gp++ )
        {
            TEgp[ gp ].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );

            for( int ki=0; ki<3; ki++ )
            {
                for( int kj=0; kj<3; kj++ )
                {
                    for( int i=0; i<mNumNodes; i++ ) 
                    {
                        TEgp[ gp ][ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * N[ i ][ gp ];
                    }
                }
            }
        }
    }

    //********************************************************************************************************************************
    //* - Permittivity tensor on Gauss points ( extrapolating TEn to GPs with N )
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, Matrix< double >& N )
    {
        // Permitivitty tensor on nodes
        Vector< Matrix< std::complex<double> > > TEn;

        PermittivityTensor_OnNodes( TEn );

        // Number of Gauss points
        int numGaussPoints = N[ 0 ].size();

        // Initializing vector of 3x3 permittivity tensors 
        TEgp.resize( numGaussPoints );

        // Permittivity tensor on Gauss points  
        for( int gp=0; gp<numGaussPoints; gp++ )
        {
            TEgp[ gp ].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );

            for( int ki=0; ki<3; ki++ )
            {
                for( int kj=0; kj<3; kj++ )
                {
                    for( int i=0; i<mNumNodes; i++ ) 
                    {
                        TEgp[ gp ][ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * N[ i ][ gp ];
                    }
                }
            }
        }
    }

    //********************************************************************************************************************************
    //* - Derivatives of the permittivity tensor on Gauss points
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::PermittivityTensor_Derivatives( Matrix< std::complex<double> >& dTEdx, 
                                                                       Matrix< std::complex<double> >& dTEdy,
                                                                       Matrix< std::complex<double> >& dTEdz,
                                                                       Vector< Matrix< std::complex<double> > >& TEn )
    {
        // 3x3 permittivity tensor derivatives
        dTEdx.Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );
        dTEdy.Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );
        dTEdz.Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );

        // 1st order base derivatives
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        // TE derivatives: dTE/dx = Sum_i{ TEni * dNi/dx } 
        for( int ki=0; ki<3; ki++ )
        {
            for( int kj=0; kj<3; kj++ )
            {
                for( int i=0; i<mNumNodes; i++ )
                {
                    dTEdx[ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * DN[ 0 ][ i ];
                    dTEdy[ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * DN[ 1 ][ i ];
                    dTEdz[ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * DN[ 2 ][ i ];
                }
            }
        }
    }

    //******************************************************************************************************************************** 
    //* - Inverse of the regularization weight tau
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::Get_invTau_OnGaussPoints( Vector< std::complex<double> >& invTau, Vector< Matrix< std::complex<double> > >& TEgp )
    {
        int numGaussPoints = TEgp.size();

        invTau.resize( numGaussPoints );
        
        for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            std::complex<double> Tr_EcE;
            
            Tr_EcE = ( TEgp[ gp ][ 0 ][ 0 ] * std::conj( TEgp[ gp ][ 0 ][ 0 ] ) + 
                       TEgp[ gp ][ 0 ][ 1 ] * std::conj( TEgp[ gp ][ 0 ][ 1 ] ) + 
                       TEgp[ gp ][ 0 ][ 2 ] * std::conj( TEgp[ gp ][ 0 ][ 2 ] ) + 
                             
                       TEgp[ gp ][ 1 ][ 0 ] * std::conj( TEgp[ gp ][ 1 ][ 0 ] ) + 
                       TEgp[ gp ][ 1 ][ 1 ] * std::conj( TEgp[ gp ][ 1 ][ 1 ] ) + 
                       TEgp[ gp ][ 1 ][ 2 ] * std::conj( TEgp[ gp ][ 1 ][ 2 ] ) + 
                             
                       TEgp[ gp ][ 2 ][ 0 ] * std::conj( TEgp[ gp ][ 2 ][ 0 ] ) + 
                       TEgp[ gp ][ 2 ][ 1 ] * std::conj( TEgp[ gp ][ 2 ][ 1 ] ) + 
                       TEgp[ gp ][ 2 ][ 2 ] * std::conj( TEgp[ gp ][ 2 ][ 2 ] ) );             
            
            invTau[ gp ] = 3.0 / Tr_EcE;
        }    
    }

    //********************************************************************************************************************************
    //* - Regularization terms Div_eNi and aDiv_eE
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::DivD_a_DivD_OnGaussPoints( Vector< Matrix<std::complex<double> > >& Div_ce_N, 
                                                                  Vector< Matrix<std::complex<double> > >& aDiv_e_E, 
                                                                  Vector< Matrix<std::complex<double> > >& TEn ,
                                                                  Vector< Matrix<std::complex<double> > >& TEgp,  
                                                                  Matrix<double>& N )
    {
        // Number of Gauss points
        int numGaussPoints = N[ 0 ].size();

        // Permittivity tensor derivatives
        Matrix< std::complex<double> > dTEdx, dTEdy, dTEdz;

        PermittivityTensor_Derivatives( dTEdx, dTEdy, dTEdz, TEn );

        // Regularization constant
        Vector< std::complex<double> > invTau;

        Get_invTau_OnGaussPoints( invTau, TEgp );

        // 1st order base derivatives
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        // Div( conj([e])*Ni )
        Div_ce_N.resize( numGaussPoints );

        // invTau * Div( [e]*E )
        aDiv_e_E.resize( numGaussPoints );

        // Div_eNi and aDiv_eE
        for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            Div_ce_N[ gp ].Resize( 3, mNumNodes, std::complex<double>( 0.0, 0.0 ) );
            aDiv_e_E[ gp ].Resize( 3, mNumNodes, std::complex<double>( 0.0, 0.0 ) );

            for( int c=0; c<3; c++ )
            {
                for( int i=0; i<mNumNodes; i++ )
                {
                    std::complex<double> gpSum( 0.0, 0.0 );

                    gpSum += TEgp[ gp ][ 0 ][ c ] * DN[ 0 ][ i ] + N[ i ][ gp ] * dTEdx[ 0 ][ c ];
                    gpSum += TEgp[ gp ][ 1 ][ c ] * DN[ 1 ][ i ] + N[ i ][ gp ] * dTEdy[ 1 ][ c ];
                    gpSum += TEgp[ gp ][ 2 ][ c ] * DN[ 2 ][ i ] + N[ i ][ gp ] * dTEdz[ 2 ][ c ];

                    Div_ce_N[ gp ][ c ][ i ] = std::conj( gpSum );
                    aDiv_e_E[ gp ][ c ][ i ] = invTau[ gp ] * gpSum;
                }
            }
        }
    }

    //***************************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //***************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
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

    //********************************************************************************************************************************
	//* - Stiffness matrix
	//********************************************************************************************************************************
	void VolumeElement_1st_ColdPlasma::GetStiffnessMatrix( Matrix<std::complex<double> >& StiffMatrix )
	{
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );	

		if( mPeso == 0.0 ) 
		{
		    Get_Non_Regularized_Matrix( StiffMatrix ); 
		}
		else if( mPotentialsOn ) 
	    {
	        Get_Regularized_Matrix_AV( StiffMatrix ); 
	    }
        else
        {
            Get_Regularized_Matrix_Ef( StiffMatrix ); 
        }	
	}

    //********************************************************************************************************************************
    //* - Non-regularized stiffness matrix ( Ef and AV formulation )
    //********************************************************************************************************************************	
    void VolumeElement_1st_ColdPlasma::Get_Non_Regularized_Matrix( Matrix< std::complex<double> >& StiffMatrix )
    {  
        // Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        double pi = 3.14159265358979;
        double mo = pi * 4.0e-7;

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> invMu = 1.0 / cMu;
		
		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // 1st order base derivatives
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        // Permittivity tensor on Gauss points
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, N );

        // [Kxx], [Kyy], [Kzz]
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=0; j<mNumNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl_XX = mVolume * ( DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );
                std::complex<double> intg_Curl_Curl_YY = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );
                std::complex<double> intg_Curl_Curl_ZZ = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] );

		        std::complex<double> intg_Ni_TEp_Nj_XX( 0.0, 0.0 );
                std::complex<double> intg_Ni_TEp_Nj_YY( 0.0, 0.0 );
                std::complex<double> intg_Ni_TEp_Nj_ZZ( 0.0, 0.0 );

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Ni_TEp_Nj_XX += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_YY += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 1 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 2 ] * N[ j ][ gp ] );
				}
	           
				StiffMatrix[ i               ][ j               ] = ( invMu * intg_Curl_Curl_XX ) - ( w2 * intg_Ni_TEp_Nj_XX );
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( invMu * intg_Curl_Curl_YY ) - ( w2 * intg_Ni_TEp_Nj_YY );
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Curl_ZZ ) - ( w2 * intg_Ni_TEp_Nj_ZZ );
			}
		}

        // [Kxy], [Kxz], [Kyz] 
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=0; j<mNumNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl_XY = mVolume * ( - DN[ 1 ][ i ] * DN[ 0 ][ j ] ); 
                std::complex<double> intg_Curl_Curl_XZ = mVolume * ( - DN[ 2 ][ i ] * DN[ 0 ][ j ] ); 
                std::complex<double> intg_Curl_Curl_YZ = mVolume * ( - DN[ 2 ][ i ] * DN[ 1 ][ j ] ); 

		        std::complex<double> intg_Ni_TEp_Nj_XY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_XZ( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YZ( 0.0, 0.0 ); 

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Ni_TEp_Nj_XY += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 1 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_XZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 2 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_YZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 2 ] * N[ j ][ gp ] );
				}
	                
				StiffMatrix[ i             ][ j +   mNumNodes ] = ( invMu * intg_Curl_Curl_XY ) - ( w2 * intg_Ni_TEp_Nj_XY );
                StiffMatrix[ i             ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Curl_XZ ) - ( w2 * intg_Ni_TEp_Nj_XZ );
                StiffMatrix[ i + mNumNodes ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Curl_YZ ) - ( w2 * intg_Ni_TEp_Nj_YZ );
			}
		}

        // [Kyx], [Kzx], [Kzy] 
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=0; j<mNumNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl_YX = mVolume * ( - DN[ 0 ][ i ] * DN[ 1 ][ j ] );
                std::complex<double> intg_Curl_Curl_ZX = mVolume * ( - DN[ 0 ][ i ] * DN[ 2 ][ j ] );
                std::complex<double> intg_Curl_Curl_ZY = mVolume * ( - DN[ 1 ][ i ] * DN[ 2 ][ j ] );

		        std::complex<double> intg_Ni_TEp_Nj_YX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZY( 0.0, 0.0 ); 

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Ni_TEp_Nj_YX += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZX += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZY += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 1 ] * N[ j ][ gp ] );
				}
	                
				StiffMatrix[ i +   mNumNodes ][ j             ] = ( invMu * intg_Curl_Curl_YX ) - ( w2 * intg_Ni_TEp_Nj_YX );
                StiffMatrix[ i + 2*mNumNodes ][ j             ] = ( invMu * intg_Curl_Curl_ZX ) - ( w2 * intg_Ni_TEp_Nj_ZX );
                StiffMatrix[ i + 2*mNumNodes ][ j + mNumNodes ] = ( invMu * intg_Curl_Curl_ZY ) - ( w2 * intg_Ni_TEp_Nj_ZY );
			}
		}

        // [AV], [VA], [VV]
        if( mPotentialsOn )
        {
            Add_AV_VA_VV_Matrix( StiffMatrix, TEgp ); 
        }
    }

    //********************************************************************************************************************************
    //* - Regularized stiffness matrix ( Ef formulation )
    //********************************************************************************************************************************	
    void VolumeElement_1st_ColdPlasma::Get_Regularized_Matrix_Ef( Matrix<std::complex<double> >& StiffMatrix )
    {  
        // Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        double pi = 3.14159265358979;
        double mo = pi * 4.0e-7;

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> invMu = 1.0 / cMu;

		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // 1st order base derivatives
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        // Permitivitty tensor on nodes
        Vector< Matrix< std::complex<double> > > TEn;

        PermittivityTensor_OnNodes( TEn );

        // Permittivity tensor on Gauss points
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, TEn, N );

        // Div( conj(e)*N ) on Gauss points
        Vector< Matrix< std::complex<double> > > Div_ce_N;

        // ( 1.0 / ( 1/3 * tr( conj(e)*e ) ) * Div( e*E ) on Gauss points
        Vector< Matrix< std::complex<double> > > aDiv_e_E;

        DivD_a_DivD_OnGaussPoints( Div_ce_N, aDiv_e_E, TEn, TEgp, N );

        // [Kxx], [Kyy], [Kzz]
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl_XX = mVolume * ( DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] ); 
                std::complex<double> intg_Curl_Curl_YY = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] ); 
                std::complex<double> intg_Curl_Curl_ZZ = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] ); 

                std::complex<double> intg_DivD_DivD_XX( 0.0, 0.0 ); 
                std::complex<double> intg_DivD_DivD_YY( 0.0, 0.0 );
                std::complex<double> intg_DivD_DivD_ZZ( 0.0, 0.0 );

                std::complex<double> intg_Ni_TEp_Nj_XX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZZ( 0.0, 0.0 ); 

                for( int gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_DivD_DivD_XX += W[ gp ] * ( Div_ce_N[ gp ][ 0 ][ i ] * aDiv_e_E[ gp ][ 0 ][ j ] );
                    intg_DivD_DivD_YY += W[ gp ] * ( Div_ce_N[ gp ][ 1 ][ i ] * aDiv_e_E[ gp ][ 1 ][ j ] );
                    intg_DivD_DivD_ZZ += W[ gp ] * ( Div_ce_N[ gp ][ 2 ][ i ] * aDiv_e_E[ gp ][ 2 ][ j ] );

                    intg_Ni_TEp_Nj_XX += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_YY += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 1 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 2 ] * N[ j ][ gp ] );
                }

                StiffMatrix[ i               ][ j               ] = ( invMu * ( intg_Curl_Curl_XX + intg_DivD_DivD_XX ) ) - ( w2 * intg_Ni_TEp_Nj_XX );
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( invMu * ( intg_Curl_Curl_YY + intg_DivD_DivD_YY ) ) - ( w2 * intg_Ni_TEp_Nj_YY );
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( invMu * ( intg_Curl_Curl_ZZ + intg_DivD_DivD_ZZ ) ) - ( w2 * intg_Ni_TEp_Nj_ZZ );
            }
        }

        // [Kxy], [Kxz], [Kyz] 
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl_XY = mVolume * ( - DN[ 1 ][ i ] * DN[ 0 ][ j ] ); 
                std::complex<double> intg_Curl_Curl_XZ = mVolume * ( - DN[ 2 ][ i ] * DN[ 0 ][ j ] ); 
                std::complex<double> intg_Curl_Curl_YZ = mVolume * ( - DN[ 2 ][ i ] * DN[ 1 ][ j ] ); 

                std::complex<double> intg_DivD_DivD_XY( 0.0, 0.0 ); 
                std::complex<double> intg_DivD_DivD_XZ( 0.0, 0.0 ); 
                std::complex<double> intg_DivD_DivD_YZ( 0.0, 0.0 ); 

                std::complex<double> intg_Ni_TEp_Nj_XY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_XZ( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YZ( 0.0, 0.0 ); 

                for ( int gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_DivD_DivD_XY += W[ gp ] * ( Div_ce_N[ gp ][ 0 ][ i ] * aDiv_e_E[ gp ][ 1 ][ j ] );
                    intg_DivD_DivD_XZ += W[ gp ] * ( Div_ce_N[ gp ][ 0 ][ i ] * aDiv_e_E[ gp ][ 2 ][ j ] );
                    intg_DivD_DivD_YZ += W[ gp ] * ( Div_ce_N[ gp ][ 1 ][ i ] * aDiv_e_E[ gp ][ 2 ][ j ] );

                    intg_Ni_TEp_Nj_XY += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 1 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_XZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 2 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_YZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 2 ] * N[ j ][ gp ] );
                }

                StiffMatrix[ i             ][ j +   mNumNodes ] = ( invMu * ( intg_Curl_Curl_XY + intg_DivD_DivD_XY ) ) - ( w2 * intg_Ni_TEp_Nj_XY );
                StiffMatrix[ i             ][ j + 2*mNumNodes ] = ( invMu * ( intg_Curl_Curl_XZ + intg_DivD_DivD_XZ ) ) - ( w2 * intg_Ni_TEp_Nj_XZ );
                StiffMatrix[ i + mNumNodes ][ j + 2*mNumNodes ] = ( invMu * ( intg_Curl_Curl_YZ + intg_DivD_DivD_YZ ) ) - ( w2 * intg_Ni_TEp_Nj_YZ );
            }
        }

        // [Kyx], [Kzx], [Kzy] 
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl_YX = mVolume * ( - DN[ 0 ][ i ] * DN[ 1 ][ j ] );
                std::complex<double> intg_Curl_Curl_ZX = mVolume * ( - DN[ 0 ][ i ] * DN[ 2 ][ j ] );
                std::complex<double> intg_Curl_Curl_ZY = mVolume * ( - DN[ 1 ][ i ] * DN[ 2 ][ j ] );

                std::complex<double> intg_DivD_DivD_YX( 0.0, 0.0 ); 
                std::complex<double> intg_DivD_DivD_ZX( 0.0, 0.0 ); 
                std::complex<double> intg_DivD_DivD_ZY( 0.0, 0.0 ); 

                std::complex<double> intg_Ni_TEp_Nj_YX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZY( 0.0, 0.0 ); 

                for( int gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_DivD_DivD_YX += W[ gp ] * ( Div_ce_N[ gp ][ 1 ][ i ] * aDiv_e_E[ gp ][ 0 ][ j ] );
                    intg_DivD_DivD_ZX += W[ gp ] * ( Div_ce_N[ gp ][ 2 ][ i ] * aDiv_e_E[ gp ][ 0 ][ j ] );
                    intg_DivD_DivD_ZY += W[ gp ] * ( Div_ce_N[ gp ][ 2 ][ i ] * aDiv_e_E[ gp ][ 1 ][ j ] );

                    intg_Ni_TEp_Nj_YX += W[gp] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZX += W[gp] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZY += W[gp] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 1 ] * N[ j ][ gp ] );
                }

                StiffMatrix[ i +   mNumNodes ][ j             ] = ( invMu * ( intg_Curl_Curl_YX + intg_DivD_DivD_YX ) ) - ( w2 * intg_Ni_TEp_Nj_YX );
                StiffMatrix[ i + 2*mNumNodes ][ j             ] = ( invMu * ( intg_Curl_Curl_ZX + intg_DivD_DivD_ZX ) ) - ( w2 * intg_Ni_TEp_Nj_ZX );
                StiffMatrix[ i + 2*mNumNodes ][ j + mNumNodes ] = ( invMu * ( intg_Curl_Curl_ZY + intg_DivD_DivD_ZY ) ) - ( w2 * intg_Ni_TEp_Nj_ZY );
            }
        }
    }

    //********************************************************************************************************************************
    //* - Regularized stiffness matrix ( AV formulation )
    //********************************************************************************************************************************	
    void VolumeElement_1st_ColdPlasma::Get_Regularized_Matrix_AV( Matrix<std::complex<double> >& StiffMatrix )
    {  
        // Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        double pi = 3.14159265358979;
        double mo = pi * 4.0e-7;

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> invMu = 1.0 / cMu;

		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // 1st order base derivatives
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        // Permittivity tensor on Gauss points
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, N );

        // [Kxx], [Kyy], [Kzz]
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Curl_Div_XYZ = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );

                std::complex<double> intg_Ni_TEp_Nj_XX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZZ( 0.0, 0.0 ); 

                for( int gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Ni_TEp_Nj_XX += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_YY += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 1 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 2 ] * N[ j ][ gp ] );
                }

                StiffMatrix[ i               ][ j               ] = ( invMu * intg_Curl_Div_XYZ ) - ( w2 * intg_Ni_TEp_Nj_XX );
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( invMu * intg_Curl_Div_XYZ ) - ( w2 * intg_Ni_TEp_Nj_YY );
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Div_XYZ ) - ( w2 * intg_Ni_TEp_Nj_ZZ );
            }
        }

        // [Kxy], [Kxz], [Kyz] 
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Curl_Div_XY = mVolume * ( DN[ 0 ][ i ] * DN[ 1 ][ j ] - DN[ 1 ][ i ] * DN[ 0 ][ j ] );
				std::complex<double> intg_Curl_Div_XZ = mVolume * ( DN[ 0 ][ i ] * DN[ 2 ][ j ] - DN[ 2 ][ i ] * DN[ 0 ][ j ] );
				std::complex<double> intg_Curl_Div_YZ = mVolume * ( DN[ 1 ][ i ] * DN[ 2 ][ j ] - DN[ 2 ][ i ] * DN[ 1 ][ j ] );

                std::complex<double> intg_Ni_TEp_Nj_XY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_XZ( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YZ( 0.0, 0.0 ); 

                for( int gp=0; gp<numGaussPoints; gp++ )
                {
					intg_Ni_TEp_Nj_XY += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 1 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_XZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 2 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_YZ += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 2 ] * N[ j ][ gp ] );
                }

                StiffMatrix[ i             ][ j +   mNumNodes ] = ( invMu * intg_Curl_Div_XY ) - ( w2 * intg_Ni_TEp_Nj_XY );
                StiffMatrix[ i             ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Div_XZ ) - ( w2 * intg_Ni_TEp_Nj_XZ );
                StiffMatrix[ i + mNumNodes ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Div_YZ ) - ( w2 * intg_Ni_TEp_Nj_YZ );
            }
        }

        // [Kyx], [Kzx], [Kzy] 
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Curl_Div_YX = mVolume * ( DN[ 1 ][ i ] * DN[ 0 ][ j ] - DN[ 0 ][ i ] * DN[ 1 ][ j ] );
				std::complex<double> intg_Curl_Div_ZX = mVolume * ( DN[ 2 ][ i ] * DN[ 0 ][ j ] - DN[ 0 ][ i ] * DN[ 2 ][ j ] );
				std::complex<double> intg_Curl_Div_ZY = mVolume * ( DN[ 2 ][ i ] * DN[ 1 ][ j ] - DN[ 1 ][ i ] * DN[ 2 ][ j ] );

                std::complex<double> intg_Ni_TEp_Nj_YX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZY( 0.0, 0.0 ); 

                for( int gp=0; gp<numGaussPoints; gp++ )
                {
					intg_Ni_TEp_Nj_YX += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZX += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 0 ] * N[ j ][ gp ] );
                    intg_Ni_TEp_Nj_ZY += W[ gp ] * ( N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 1 ] * N[ j ][ gp ] );
                }

                StiffMatrix[ i +   mNumNodes ][ j             ] = ( invMu * intg_Curl_Div_YX ) - ( w2 * intg_Ni_TEp_Nj_YX );
                StiffMatrix[ i + 2*mNumNodes ][ j             ] = ( invMu * intg_Curl_Div_ZX ) - ( w2 * intg_Ni_TEp_Nj_ZX );
                StiffMatrix[ i + 2*mNumNodes ][ j + mNumNodes ] = ( invMu * intg_Curl_Div_ZY ) - ( w2 * intg_Ni_TEp_Nj_ZY );
            }
        }

        // [AV], [VA], [VV]
        Add_AV_VA_VV_Matrix( StiffMatrix, TEgp ); 
    }

    //********************************************************************************************************************************
    // - AV-potentials contribution [AV], [VA], [VV] to stiffness matrix
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::Add_AV_VA_VV_Matrix( Matrix< std::complex<double> >& StiffMatrix, Vector< Matrix< std::complex<double> > >& TEgp ) 
	{ 
        // Problem frequency
		double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;
        
		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;
        
		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

		// 1st order basis derivatives 
        double DN[ 3 ][ 4 ]; 
        
        Calculate_DN( DN );

        // [Kvv]
		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
                std::complex<double> intg_dNi_TEp_dNj( 0.0, 0.0 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_dNi_TEp_dNj += W[ gp ] * DN[ 0 ][ i ] * ( TEgp[ gp ][ 0 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 0 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 0 ][ 2 ] * DN[ 2 ][ j ] );
                    intg_dNi_TEp_dNj += W[ gp ] * DN[ 1 ][ i ] * ( TEgp[ gp ][ 1 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 1 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 1 ][ 2 ] * DN[ 2 ][ j ] );
                    intg_dNi_TEp_dNj += W[ gp ] * DN[ 2 ][ i ] * ( TEgp[ gp ][ 2 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 2 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 2 ][ 2 ] * DN[ 2 ][ j ] );
				}
	                
				StiffMatrix[ i + 3*mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_dNi_TEp_dNj );
        	}
        }

		// [Kxv], [Kyv], [Kzv]
		for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
                std::complex<double> intg_NiX_TEp_dNj( 0.0, 0.0 ); 
                std::complex<double> intg_NiY_TEp_dNj( 0.0, 0.0 ); 
                std::complex<double> intg_NiZ_TEp_dNj( 0.0, 0.0 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_NiX_TEp_dNj += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 0 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 0 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 0 ][ 2 ] * DN[ 2 ][ j ] );
                    intg_NiY_TEp_dNj += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 1 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 1 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 1 ][ 2 ] * DN[ 2 ][ j ] );
                    intg_NiZ_TEp_dNj += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 2 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 2 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 2 ][ 2 ] * DN[ 2 ][ j ] );
				}
	                
				StiffMatrix[ i               ][ j + 3*mNumNodes ] = ( - w2 * intg_NiX_TEp_dNj );
                StiffMatrix[ i +   mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_NiY_TEp_dNj );
                StiffMatrix[ i + 2*mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_NiZ_TEp_dNj );
        	}
        }

		// [Kvx], [Kvy], [Kvz]
		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<mNumNodes; j++ ) 
        	{
                std::complex<double> intg_dNi_TEp_NjX( 0.0, 0.0 ); 
                std::complex<double> intg_dNi_TEp_NjY( 0.0, 0.0 ); 
                std::complex<double> intg_dNi_TEp_NjZ( 0.0, 0.0 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_dNi_TEp_NjX += W[ gp ] * N[ j ][ gp ] * ( TEgp[ gp ][ 0 ][ 0 ] * DN[ 0 ][ i ] + TEgp[ gp ][ 1 ][ 0 ] * DN[ 1 ][ i ] + TEgp[ gp ][ 2 ][ 0 ] * DN[ 2 ][ i ] );
                    intg_dNi_TEp_NjY += W[ gp ] * N[ j ][ gp ] * ( TEgp[ gp ][ 0 ][ 1 ] * DN[ 0 ][ i ] + TEgp[ gp ][ 1 ][ 1 ] * DN[ 1 ][ i ] + TEgp[ gp ][ 2 ][ 1 ] * DN[ 2 ][ i ] );
                    intg_dNi_TEp_NjZ += W[ gp ] * N[ j ][ gp ] * ( TEgp[ gp ][ 0 ][ 2 ] * DN[ 0 ][ i ] + TEgp[ gp ][ 1 ][ 2 ] * DN[ 1 ][ i ] + TEgp[ gp ][ 2 ][ 2 ] * DN[ 2 ][ i ] );
				}
	                
				StiffMatrix[ i + 3*mNumNodes ][ j               ] = ( - w2 * intg_dNi_TEp_NjX );
                StiffMatrix[ i + 3*mNumNodes ][ j +   mNumNodes ] = ( - w2 * intg_dNi_TEp_NjY );
                StiffMatrix[ i + 3*mNumNodes ][ j + 2*mNumNodes ] = ( - w2 * intg_dNi_TEp_NjZ );
        	}
        }
	}
} 