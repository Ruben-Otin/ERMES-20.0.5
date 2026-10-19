
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
                EquationId[ i               ] = mNodes[i]->pDofcAx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[i]->pDofcAy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[i]->pDofcAz()->EquationId();
				EquationId[ i + 3*mNumNodes ] = mNodes[i]->pDofcVs()->EquationId();
            }				
		}
		else
		{
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                EquationId[ i               ] = mNodes[i]->pDofcEx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[i]->pDofcEy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[i]->pDofcEz()->EquationId();
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
         
        return fabs( Det / 6.00 );
    }

    //********************************************************************************************************************************
    //* - Set element as IHL material or plasma 
    //********************************************************************************************************************************	
    void VolumeElement_1st_ColdPlasma::SetPlasmaIHL() 
    { 
        // Material properties
        double sigma   = (*mProperties)(  IHL_ELECTRIC_CONDUCTIVITY );

        double ep_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY );
        double ep_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY );

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY );
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY );

        // If the material is vaccum then it is considered as plasma. 
        // We only consider the first one written in the plasma file *-8.dat
        // This allows the use of different plasmas volumes (for visualization). 
        if ( ( sigma == 0.0 ) && ( ep_real == 1.0 ) && ( ep_imag == 0.0 ) && ( mu_real == 1.0 ) && ( mu_imag == 0.0 ) )
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
    void VolumeElement_1st_ColdPlasma::Calculate_DN( double DN[3][4] )
    {
        double Cte = 1.00 / ( 6.00 * mVolume ); 

        // dNi/dx
        DN[0][0] = Cte * ( Y(3)*Z(2) + Y(4)*Z(3) + Y(2)*Z(4) - Y(3)*Z(4) - Y(2)*Z(3) - Y(4)*Z(2) );
        DN[0][1] = Cte * ( Y(3)*Z(4) + Y(4)*Z(1) + Y(1)*Z(3) - Y(3)*Z(1) - Y(1)*Z(4) - Y(4)*Z(3) );
        DN[0][2] = Cte * ( Y(2)*Z(1) + Y(1)*Z(4) + Y(4)*Z(2) - Y(2)*Z(4) - Y(4)*Z(1) - Y(1)*Z(2) );
        DN[0][3] = Cte * ( Y(2)*Z(3) + Y(3)*Z(1) + Y(1)*Z(2) - Y(2)*Z(1) - Y(3)*Z(2) - Y(1)*Z(3) );
        
        // dNi/dy
        DN[1][0] = Cte * ( X(3)*Z(4) + X(4)*Z(2) + X(2)*Z(3) - X(3)*Z(2) - X(2)*Z(4) - X(4)*Z(3) ); 
        DN[1][1] = Cte * ( X(3)*Z(1) + X(1)*Z(4) + X(4)*Z(3) - X(1)*Z(3) - X(3)*Z(4) - X(4)*Z(1) );
        DN[1][2] = Cte * ( X(2)*Z(4) + X(4)*Z(1) + X(1)*Z(2) - X(2)*Z(1) - X(4)*Z(2) - X(1)*Z(4) );
        DN[1][3] = Cte * ( X(1)*Z(3) + X(3)*Z(2) + X(2)*Z(1) - X(3)*Z(1) - X(1)*Z(2) - X(2)*Z(3) );
        
        // dNi/dz
        DN[2][0] = Cte * ( X(4)*Y(3) + X(2)*Y(4) + X(3)*Y(2) - X(2)*Y(3) - X(4)*Y(2) - X(3)*Y(4) );
        DN[2][1] = Cte * ( X(3)*Y(4) + X(4)*Y(1) + X(1)*Y(3) - X(3)*Y(1) - X(4)*Y(3) - X(1)*Y(4) ); 
        DN[2][2] = Cte * ( X(2)*Y(1) + X(4)*Y(2) + X(1)*Y(4) - X(4)*Y(1) - X(2)*Y(4) - X(1)*Y(2) ); 
        DN[2][3] = Cte * ( X(2)*Y(3) + X(3)*Y(1) + X(1)*Y(2) - X(2)*Y(1) - X(3)*Y(2) - X(1)*Y(3) );    
    }

    //********************************************************************************************************************************
    // - Calculation of dN/dx, dN/dy, dN/dz on cX, cY, cZ points
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::DerivativesOnPoints( Matrix< double >& dNdx,
		                                                    Matrix< double >& dNdy,
													        Matrix< double >& dNdz,
													        Vector< double >&   cX, 
													        Vector< double >&   cY, 
													        Vector< double >&   cZ )
    {
		int numPoints = cX.size();

		dNdx.Resize( mNumNodes, numPoints );
		dNdy.Resize( mNumNodes, numPoints );
		dNdz.Resize( mNumNodes, numPoints );
		
		double DN[3][4]; 
		
		Calculate_DN( DN );

		double dL1dx, dL2dx, dL3dx, dL4dx;
        double dL1dy, dL2dy, dL3dy, dL4dy;
		double dL1dz, dL2dz, dL3dz, dL4dz;

        dL1dx = DN[0][0], dL2dx = DN[0][1], dL3dx = DN[0][2], dL4dx = DN[0][3];
        dL1dy = DN[1][0], dL2dy = DN[1][1], dL3dy = DN[1][2], dL4dy = DN[1][3];
		dL1dz = DN[2][0], dL2dz = DN[2][1], dL3dz = DN[2][2], dL4dz = DN[2][3];

		// dNi/dx, dNi/dy, dNi/dz
        for( int p=0; p<numPoints; p++ ) 
        {
			dNdx[0][p] = dL1dx;
			dNdy[0][p] = dL1dy;
			dNdz[0][p] = dL1dz;

			dNdx[1][p] = dL2dx;
			dNdy[1][p] = dL2dy;
			dNdz[1][p] = dL2dz;

			dNdx[2][p] = dL3dx;
			dNdy[2][p] = dL3dy;
			dNdz[2][p] = dL3dz;

			dNdx[3][p] = dL4dx;
			dNdy[3][p] = dL4dy;
			dNdz[3][p] = dL4dz;
        }
	}

    //******************************************************************************************************************************** 
    //* - Calculates rotation matrix 
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::Rotation_Matrix( Matrix< std::complex<double> >& T, double Bx, double By, double Bz )
    {    
         // Unit vectors
         Vector< double > n( 3 ), t( 3 ), b( 3 );

         // b = B / |B|
         double Bnorm = std::sqrt( Bx*Bx + By*By + Bz*Bz );

         b[0] = Bx / Bnorm; 
         b[1] = By / Bnorm; 
         b[2] = Bz / Bnorm;

         // n _|_ b
         double nnorm = std::sqrt( b[0]*b[0] + b[1]*b[1] );

         if ( nnorm > 0.00 )
         {
             n[0] = b[1] / nnorm; 
             n[1] =-b[0] / nnorm; 
             n[2] = 0.00;
         }
         else
         {
             n[0] = b[2];
             n[1] = 0.00; 
             n[2] = 0.00;
         }

         // t = b x n
         t[0] = b[1]*n[2] - b[2]*n[1];
         t[1] = b[2]*n[0] - b[0]*n[2];
         t[2] = b[0]*n[1] - b[1]*n[0];

         // Rotation matrix 
         // T = [ nx tx bx
         //       ny ty by 
         //       nz tz bz ]
         T.Resize( 3, 3 );

         for ( int row=0; row<3; row++ )
         {
             T[row][0] = std::complex<double>( n[row] ); 
             T[row][1] = std::complex<double>( t[row] ); 
             T[row][2] = std::complex<double>( b[row] );       
         }
    }

    //********************************************************************************************************************************
    //* - Calculates permittivity tensor on nodes 
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

        if ( mIs_IHL )
        {
            for ( int i=0; i<mNumNodes; i++ ) 
            { 
                TEn[i].Resize( 3, 3 );

                double sigma    = (*mProperties)(  IHL_ELECTRIC_CONDUCTIVITY )     ;
                double eps_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
                double eps_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

                std::complex<double> cEps( eps_real, eps_imag + ( sigma / wf ) );

                // IHL permittivity tensor 
                Matrix< std::complex<double> > K( 3, 3 );

                K[0][0] = cEps ; K[0][1] = cZero; K[0][2] = cZero; 
                K[1][0] = cZero; K[1][1] = cEps ; K[1][2] = cZero; 
                K[2][0] = cZero; K[2][1] = cZero; K[2][2] = cEps ;

                TEn[i] = K;
            }
        }
        else
        {
            for ( int i=0; i<mNumNodes; i++ ) 
            { 
                TEn[i].Resize( 3, 3 );

                // Plasma parameters on nodes
                Vector< double > Bext; 

                Vector< std::complex<double> > SDPRL;

                mpColdPlasma->Get_Bext_SDPRL_InNode( Bext, SDPRL, mNodes[i], wf );

                double Bx = Bext[0];
                double By = Bext[1];
                double Bz = Bext[2];

                std::complex<double> S = SDPRL[0];
                std::complex<double> D = SDPRL[1];
                std::complex<double> P = SDPRL[2];
               
                // Permittivity tensor in b = B/|B| coordinates
                Matrix<std::complex<double> > K( 3, 3 );

                K[0][0] = S      ; K[0][1] = -cUnit*D; K[0][2] = cZero; 
                K[1][0] = cUnit*D; K[1][1] =  S      ; K[1][2] = cZero; 
                K[2][0] = cZero  ; K[2][1] =  cZero  ; K[2][2] = P    ;

                K *= eo;

                // Rotation matrix 
                Matrix< std::complex<double> > T; 

                Rotation_Matrix( T, Bx, By, Bz );

                // TEn = T * K * Tt 
                // Be carefull!!, T.Transpose() changes T permanently
                TEn[i] = T * K;
                TEn[i] = TEn[i] * T.Transpose();
            }
        }
    }

    //********************************************************************************************************************************
    //* - Calculates permittivity tensor on Gauss points ( extrapolating TEn to GPs with N )
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp,  
                                                                         Vector< Matrix< std::complex<double> > >& TEn , 
                                                                         Matrix< double >& N )
    {
        // Number of Gauss points
        int numGaussPoints = N[0].size();

        // Initializing vector of 3x3 permittivity tensors 
        TEgp.resize( numGaussPoints );

        // Permittivity tensor on Gauss points  
        for ( int gp=0; gp<numGaussPoints; gp++ )
        {
            TEgp[gp].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );

            for ( int ki=0; ki<3; ki++ )
            {
                for ( int kj=0; kj<3; kj++ )
                {
                    for ( int i=0; i<mNumNodes; i++ ) 
                    {
                        TEgp[gp][ki][kj] += TEn[i][ki][kj] * N[i][gp];
                    }
                }
            }
        }
    }

    //********************************************************************************************************************************
    //* - Calculate derivatives of the permittivity tensor on Gauss points
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::PermittivityTensor_Derivatives( Vector< Matrix< std::complex<double> > >& dTEdx, 
                                                                       Vector< Matrix< std::complex<double> > >& dTEdy,
                                                                       Vector< Matrix< std::complex<double> > >& dTEdz,
                                                                       Vector< Matrix< std::complex<double> > >& TEn  ,
                                                                       Matrix< double >& dNdx, 
                                                                       Matrix< double >& dNdy, 
                                                                       Matrix< double >& dNdz )
    {
        // Number of Gauss points
        int numGaussPoints = dNdx[0].size();

        // Initializing vector of 3x3 permittivity tensors derivatives
        dTEdx.resize( numGaussPoints );
        dTEdy.resize( numGaussPoints );
        dTEdz.resize( numGaussPoints );

        // TE derivatives: (dTE/dx)|gp = Sum_i{ TEni * (dNi/dx)|gp } 
        for ( int gp=0; gp<numGaussPoints; gp++ )
        { 
            dTEdx[gp].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );
            dTEdy[gp].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );
            dTEdz[gp].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );

            for ( int ki=0; ki<3; ki++ )
            {
                for ( int kj=0; kj<3; kj++ )
                {
                    for ( int i=0; i<mNumNodes; i++ )
                    {
                        dTEdx[gp][ki][kj] += TEn[i][ki][kj] * dNdx[i][gp];
                        dTEdy[gp][ki][kj] += TEn[i][ki][kj] * dNdy[i][gp];
                        dTEdz[gp][ki][kj] += TEn[i][ki][kj] * dNdz[i][gp];
                    }
                }
            }
        }
    }

    //******************************************************************************************************************************** 
    //* - Inverse of the regularization weight tau
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::Get_invTau_InGaussPoints( Vector<         std::complex<double>   >& invTau, 
                                                                 Vector< Matrix< std::complex<double> > >& TEgp   )
    {
        int numGaussPoints = TEgp.size();

        invTau.resize( numGaussPoints );
        
        for ( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            std::complex<double> Tr_EcE;
            
            Tr_EcE = ( TEgp[gp][0][0] * std::conj( TEgp[gp][0][0] ) + 
                       TEgp[gp][0][1] * std::conj( TEgp[gp][0][1] ) + 
                       TEgp[gp][0][2] * std::conj( TEgp[gp][0][2] ) + 
                       
                       TEgp[gp][1][0] * std::conj( TEgp[gp][1][0] ) + 
                       TEgp[gp][1][1] * std::conj( TEgp[gp][1][1] ) + 
                       TEgp[gp][1][2] * std::conj( TEgp[gp][1][2] ) + 
                       
                       TEgp[gp][2][0] * std::conj( TEgp[gp][2][0] ) + 
                       TEgp[gp][2][1] * std::conj( TEgp[gp][2][1] ) + 
                       TEgp[gp][2][2] * std::conj( TEgp[gp][2][2] ) );             
            
            invTau[gp] = 3.0 / Tr_EcE;
        }    
    }

    //********************************************************************************************************************************
    //* - Calculation of the regularization terms Div_eNi and aDiv_eE
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::DivD_a_DivD_OnGaussPoints( Vector< Matrix<std::complex<double> > >& Div_ce_N, 
                                                                  Vector< Matrix<std::complex<double> > >& aDiv_e_E, 
                                                                  Vector< Matrix<std::complex<double> > >& TEn ,
                                                                  Vector< Matrix<std::complex<double> > >& TEgp,  
                                                                  Matrix< double >& N, 
                                                                  Matrix< double >& dNdx, 
                                                                  Matrix< double >& dNdy, 
                                                                  Matrix< double >& dNdz )
    {
        // Number of Gauss points
        int numGaussPoints = N[0].size();

        // Permittivity tensor derivatives
        Vector< Matrix< std::complex<double> > > dTEdx;
        Vector< Matrix< std::complex<double> > > dTEdy;
        Vector< Matrix< std::complex<double> > > dTEdz;

        PermittivityTensor_Derivatives( dTEdx, dTEdy, dTEdz, TEn, dNdx, dNdy, dNdz );

        // Computing: Div( conj([e])*Ni )
        Div_ce_N.resize( numGaussPoints );

        // Computing: invTau * Div( [e]*E )
        aDiv_e_E.resize( numGaussPoints );

        // Regularization constant
        Vector< std::complex<double> > invTau;

        Get_invTau_InGaussPoints( invTau, TEgp );

        for ( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            Div_ce_N[gp].Resize( 3, mNumNodes, std::complex<double>( 0.0, 0.0 ) );
            aDiv_e_E[gp].Resize( 3, mNumNodes, std::complex<double>( 0.0, 0.0 ) );

            for ( int c=0; c<3; c++ )
            {
                for ( int i=0; i<mNumNodes; i++ )
                {
                    std::complex<double> gpSum ( 0.00, 0.00 );

                    gpSum += TEgp[gp][0][c] * dNdx[i][gp] + N[i][gp] * dTEdx[gp][0][c];
                    gpSum += TEgp[gp][1][c] * dNdy[i][gp] + N[i][gp] * dTEdy[gp][1][c];
                    gpSum += TEgp[gp][2][c] * dNdz[i][gp] + N[i][gp] * dTEdz[gp][2][c];

                    Div_ce_N[gp][c][i] = std::conj( gpSum );

                    aDiv_e_E[gp][c][i] = invTau[gp] * gpSum;
                }
            }
        }
    }

    //********************************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue    , 
                                                                    Matrix  <               std::complex<double> >& EleStiffMatrix,  
                                                                    Vector  <               std::complex<double> >& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		if ( mPotentialsOn ) 
		{
            for ( int j=0; j<4; j++ ) 
            { 
                if ( FixedValue.find( mNodes[j]->Id() ) != FixedValue.end() )
                {    
                    std::complex<double> cVoltage = FixedValue[ mNodes[j]->Id() ];
		    
                    for( int i=0; i<mNumDofs; i++ ) 
                    {
                        ResidualVector[ i ] -= EleStiffMatrix[ i ][ j + 3*mNumNodes ] * cVoltage;
                    }
                }
            }
		}
    }

    //********************************************************************************************************************************
	//* - Element stiffness matrix
	//********************************************************************************************************************************
	void VolumeElement_1st_ColdPlasma::GetStiffnessMatrix( Matrix<std::complex<double> >& StiffMatrix )
	{
        // Resizing stiffness matrix
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );	

		// Get non-regularized matrix
		if ( mPeso == 0.00 ) 
		{
		    Get_NonRegularized_Matrix( StiffMatrix ); 
		}
		// Get regularized matrix
		else
		{
		    if ( mPotentialsOn ) 
            {
                Get_Regularized_AV_Matrix( StiffMatrix ); 
            }
            else                 
            {
                Get_Regularized_Ef_Matrix( StiffMatrix ); 
            }
		}
	}

    //********************************************************************************************************************************
    //* - Get non-regularized stiffness matrix
    //********************************************************************************************************************************	
    void VolumeElement_1st_ColdPlasma::Get_NonRegularized_Matrix( Matrix< std::complex<double> >& StiffMatrix )
    {  
        // Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        double pi = 3.1415926536;
	    double mo = pi * 4.00e-7; 

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu   ( mu_real, mu_imag );
		std::complex<double> invMu = 1.00 / cMu;
		
		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[gp] *= detJ; 
		}

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		DerivativesOnPoints( dNdx, dNdy, dNdz, cX, cY, cZ );

        // Permitivitty tensor on nodes
        Vector< Matrix< std::complex<double> > > TEn;

        PermittivityTensor_OnNodes( TEn );

        // Permittivity tensor on Gauss points
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, TEn, N );

        // Upper diagonal [Kxx], [Kyy], [Kzz]
		for ( int i=0; i<mNumNodes; i++ )
		{
			for ( int j=i; j<mNumNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl_XX( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_Curl_YY( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_Curl_ZZ( 0.00, 0.00 ); 

		        std::complex<double> intg_Ni_TEp_Nj_XX( 0.00, 0.00 );
                std::complex<double> intg_Ni_TEp_Nj_YY( 0.00, 0.00 );
                std::complex<double> intg_Ni_TEp_Nj_ZZ( 0.00, 0.00 );

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_Curl_XX += W[gp] * ( dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
                    intg_Curl_Curl_YY += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
                    intg_Curl_Curl_ZZ += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] );

					intg_Ni_TEp_Nj_XX += W[gp] * ( N[i][gp] * TEgp[gp][0][0] * N[j][gp] );
                    intg_Ni_TEp_Nj_YY += W[gp] * ( N[i][gp] * TEgp[gp][1][1] * N[j][gp] );
                    intg_Ni_TEp_Nj_ZZ += W[gp] * ( N[i][gp] * TEgp[gp][2][2] * N[j][gp] );
				}
	           
				StiffMatrix[ i               ][ j               ] = ( invMu * intg_Curl_Curl_XX ) - ( w2 * intg_Ni_TEp_Nj_XX );
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( invMu * intg_Curl_Curl_YY ) - ( w2 * intg_Ni_TEp_Nj_YY );
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Curl_ZZ ) - ( w2 * intg_Ni_TEp_Nj_ZZ );
			}
		}

        // [Kxy], [Kxz], [Kyz] 
		for ( int i=0; i<mNumNodes; i++ )
		{
			for ( int j=0; j<mNumNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl_XY( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_Curl_XZ( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_Curl_YZ( 0.00, 0.00 ); 

		        std::complex<double> intg_Ni_TEp_Nj_XY( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj_XZ( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj_YZ( 0.00, 0.00 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_Curl_XY += W[gp] * ( - dNdy[i][gp] * dNdx[j][gp] );
                    intg_Curl_Curl_XZ += W[gp] * ( - dNdz[i][gp] * dNdx[j][gp] );
                    intg_Curl_Curl_YZ += W[gp] * ( - dNdz[i][gp] * dNdy[j][gp] );

					intg_Ni_TEp_Nj_XY += W[gp] * ( N[i][gp] * TEgp[gp][0][1] * N[j][gp] );
                    intg_Ni_TEp_Nj_XZ += W[gp] * ( N[i][gp] * TEgp[gp][0][2] * N[j][gp] );
                    intg_Ni_TEp_Nj_YZ += W[gp] * ( N[i][gp] * TEgp[gp][1][2] * N[j][gp] );
				}
	                
				StiffMatrix[ i             ][ j +   mNumNodes ] = ( invMu * intg_Curl_Curl_XY ) - ( w2 * intg_Ni_TEp_Nj_XY );
                StiffMatrix[ i             ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Curl_XZ ) - ( w2 * intg_Ni_TEp_Nj_XZ );
                StiffMatrix[ i + mNumNodes ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Curl_YZ ) - ( w2 * intg_Ni_TEp_Nj_YZ );
			}
		}

        if ( mPotentialsOn )
        {
            // Upper diagonal [Kvv]
		    for( int i=0; i<mNumNodes; i++ )
            {
                for( int j=i; j<mNumNodes; j++ ) 
            	{
                    std::complex<double> intg_dNi_TEp_dNj( 0.00, 0.00 ); 
            
		    		for ( int gp=0; gp<numGaussPoints; gp++ )
		    		{
		    			intg_dNi_TEp_dNj += W[gp] * dNdx[i][gp] * ( TEgp[gp][0][0] * dNdx[j][gp] + TEgp[gp][0][1] * dNdy[j][gp] + TEgp[gp][0][2] * dNdz[j][gp] );
                        intg_dNi_TEp_dNj += W[gp] * dNdy[i][gp] * ( TEgp[gp][1][0] * dNdx[j][gp] + TEgp[gp][1][1] * dNdy[j][gp] + TEgp[gp][1][2] * dNdz[j][gp] );
                        intg_dNi_TEp_dNj += W[gp] * dNdz[i][gp] * ( TEgp[gp][2][0] * dNdx[j][gp] + TEgp[gp][2][1] * dNdy[j][gp] + TEgp[gp][2][2] * dNdz[j][gp] );
		    		}
	                    
		    		StiffMatrix[ i + 3*mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_dNi_TEp_dNj );
            	}
            }
            
		    // [Kxv], [Kyv], [Kzv]
		    for( int i=0; i<mNumNodes; i++ )
            {
                for( int j=0; j<mNumNodes; j++ ) 
            	{
                    std::complex<double> intg_NiX_TEp_dNj( 0.00, 0.00 ); 
                    std::complex<double> intg_NiY_TEp_dNj( 0.00, 0.00 ); 
                    std::complex<double> intg_NiZ_TEp_dNj( 0.00, 0.00 ); 
            
		    		for ( int gp=0; gp<numGaussPoints; gp++ )
		    		{
		    			intg_NiX_TEp_dNj += W[gp] * N[i][gp] * ( TEgp[gp][0][0] * dNdx[j][gp] + TEgp[gp][0][1] * dNdy[j][gp] + TEgp[gp][0][2] * dNdz[j][gp] );
                        intg_NiY_TEp_dNj += W[gp] * N[i][gp] * ( TEgp[gp][1][0] * dNdx[j][gp] + TEgp[gp][1][1] * dNdy[j][gp] + TEgp[gp][1][2] * dNdz[j][gp] );
                        intg_NiZ_TEp_dNj += W[gp] * N[i][gp] * ( TEgp[gp][2][0] * dNdx[j][gp] + TEgp[gp][2][1] * dNdy[j][gp] + TEgp[gp][2][2] * dNdz[j][gp] );
		    		}
	                    
		    		StiffMatrix[ i               ][ j + 3*mNumNodes ] = ( - w2 * intg_NiX_TEp_dNj );
                    StiffMatrix[ i +   mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_NiY_TEp_dNj );
                    StiffMatrix[ i + 2*mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_NiZ_TEp_dNj );
            	}
            }
        }

        // Lower diagonal 
		Generate_Lower_Diagonal( StiffMatrix );
    }

    //********************************************************************************************************************************
    //* - Get regularized stiffness matrix for the electric field formulation
    //********************************************************************************************************************************	
    void VolumeElement_1st_ColdPlasma::Get_Regularized_Ef_Matrix( Matrix<std::complex<double> >& StiffMatrix )
    {  
        // Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        double pi = 3.1415926536;
	    double mo = pi * 4.00e-7; 

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu   ( mu_real, mu_imag );
		std::complex<double> invMu = 1.00 / cMu;

		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[gp] *= detJ; 
		}

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		DerivativesOnPoints( dNdx, dNdy, dNdz, cX, cY, cZ );

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

        DivD_a_DivD_OnGaussPoints( Div_ce_N, aDiv_e_E, TEn, TEgp, N, dNdx, dNdy, dNdz );

        // Upper diagonal [Kxx], [Kyy], [Kzz]
        for ( int i=0; i<mNumNodes; i++ )
        {
            for ( int j=i; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl_XX( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_Curl_YY( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_Curl_ZZ( 0.00, 0.00 ); 

                std::complex<double> intg_DivD_DivD_XX( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD_YY( 0.00, 0.00 );
                std::complex<double> intg_DivD_DivD_ZZ( 0.00, 0.00 );

                std::complex<double> intg_Ni_TEp_Nj_XX( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj_YY( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZZ( 0.00, 0.00 ); 

                for ( int gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Curl_Curl_XX += W[gp] * ( dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
                    intg_Curl_Curl_YY += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
                    intg_Curl_Curl_ZZ += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] );

                    intg_DivD_DivD_XX += W[gp] * ( Div_ce_N[gp][0][i] * aDiv_e_E[gp][0][j] );
                    intg_DivD_DivD_YY += W[gp] * ( Div_ce_N[gp][1][i] * aDiv_e_E[gp][1][j] );
                    intg_DivD_DivD_ZZ += W[gp] * ( Div_ce_N[gp][2][i] * aDiv_e_E[gp][2][j] );

                    intg_Ni_TEp_Nj_XX += W[gp] * ( N[i][gp] * TEgp[gp][0][0] * N[j][gp] );
                    intg_Ni_TEp_Nj_YY += W[gp] * ( N[i][gp] * TEgp[gp][1][1] * N[j][gp] );
                    intg_Ni_TEp_Nj_ZZ += W[gp] * ( N[i][gp] * TEgp[gp][2][2] * N[j][gp] );
                }

                StiffMatrix[ i               ][ j               ] = ( invMu * ( intg_Curl_Curl_XX + intg_DivD_DivD_XX ) ) - ( w2 * intg_Ni_TEp_Nj_XX );
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( invMu * ( intg_Curl_Curl_YY + intg_DivD_DivD_YY ) ) - ( w2 * intg_Ni_TEp_Nj_YY );
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( invMu * ( intg_Curl_Curl_ZZ + intg_DivD_DivD_ZZ ) ) - ( w2 * intg_Ni_TEp_Nj_ZZ );
            }
        }

        // [Kxy], [Kxz], [Kyz] 
        for ( int i=0; i<mNumNodes; i++ )
        {
            for ( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl_XY( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_Curl_XZ( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_Curl_YZ( 0.00, 0.00 ); 

                std::complex<double> intg_DivD_DivD_XY( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD_XZ( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD_YZ( 0.00, 0.00 ); 

                std::complex<double> intg_Ni_TEp_Nj_XY( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj_XZ( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj_YZ( 0.00, 0.00 ); 

                for ( int gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Curl_Curl_XY += W[gp] * ( - dNdy[i][gp] * dNdx[j][gp] );
                    intg_Curl_Curl_XZ += W[gp] * ( - dNdz[i][gp] * dNdx[j][gp] );
                    intg_Curl_Curl_YZ += W[gp] * ( - dNdz[i][gp] * dNdy[j][gp] );

                    intg_DivD_DivD_XY += W[gp] * ( Div_ce_N[gp][0][i] * aDiv_e_E[gp][1][j] );
                    intg_DivD_DivD_XZ += W[gp] * ( Div_ce_N[gp][0][i] * aDiv_e_E[gp][2][j] );
                    intg_DivD_DivD_YZ += W[gp] * ( Div_ce_N[gp][1][i] * aDiv_e_E[gp][2][j] );

                    intg_Ni_TEp_Nj_XY += W[gp] * ( N[i][gp] * TEgp[gp][0][1] * N[j][gp] );
                    intg_Ni_TEp_Nj_XZ += W[gp] * ( N[i][gp] * TEgp[gp][0][2] * N[j][gp] );
                    intg_Ni_TEp_Nj_YZ += W[gp] * ( N[i][gp] * TEgp[gp][1][2] * N[j][gp] );
                }

                StiffMatrix[ i             ][ j +   mNumNodes ] = ( invMu * ( intg_Curl_Curl_XY + intg_DivD_DivD_XY ) ) - ( w2 * intg_Ni_TEp_Nj_XY );
                StiffMatrix[ i             ][ j + 2*mNumNodes ] = ( invMu * ( intg_Curl_Curl_XZ + intg_DivD_DivD_XZ ) ) - ( w2 * intg_Ni_TEp_Nj_XZ );
                StiffMatrix[ i + mNumNodes ][ j + 2*mNumNodes ] = ( invMu * ( intg_Curl_Curl_YZ + intg_DivD_DivD_YZ ) ) - ( w2 * intg_Ni_TEp_Nj_YZ );
            }
        }

        // Lower diagonal 
		Generate_Lower_Diagonal( StiffMatrix );
    }

    //********************************************************************************************************************************
    // - Get regularized stiffness matrix for the AV-potential formulation
    //********************************************************************************************************************************
    void VolumeElement_1st_ColdPlasma::Get_Regularized_AV_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
        // Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        std::complex<double> jw( 0.0, -wf );

        double pi = 3.1415926536;
	    double mo = pi * 4.00e-7; 

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu   ( mu_real, mu_imag );
		std::complex<double> invMu = 1.00 / cMu;
		
		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[gp] *= detJ; 
		}

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		DerivativesOnPoints( dNdx, dNdy, dNdz, cX, cY, cZ );

        // Permitivitty tensor on nodes
        Vector< Matrix< std::complex<double> > > TEn;

        PermittivityTensor_OnNodes( TEn );

        // Permittivity tensor on Gauss points
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, TEn, N );

        // Upper diagonal [Kxx], [Kyy], [Kzz]
		for ( int i=0; i<mNumNodes; i++ )
		{
			for ( int j=i; j<mNumNodes; j++ )
			{
                std::complex<double> intg_Curl_Div_XYZ( 0.00, 0.00 ); 

		        std::complex<double> intg_Ni_TEp_Nj_XX( 0.00, 0.00 );
                std::complex<double> intg_Ni_TEp_Nj_YY( 0.00, 0.00 );
                std::complex<double> intg_Ni_TEp_Nj_ZZ( 0.00, 0.00 );

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
                    intg_Curl_Div_XYZ += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );

					intg_Ni_TEp_Nj_XX += W[gp] * ( N[i][gp] * TEgp[gp][0][0] * N[j][gp] );
                    intg_Ni_TEp_Nj_YY += W[gp] * ( N[i][gp] * TEgp[gp][1][1] * N[j][gp] );
                    intg_Ni_TEp_Nj_ZZ += W[gp] * ( N[i][gp] * TEgp[gp][2][2] * N[j][gp] );
				}
	           
				StiffMatrix[ i               ][ j               ] = ( invMu * intg_Curl_Div_XYZ ) - ( w2 * intg_Ni_TEp_Nj_XX );
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( invMu * intg_Curl_Div_XYZ ) - ( w2 * intg_Ni_TEp_Nj_YY );
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_Div_XYZ ) - ( w2 * intg_Ni_TEp_Nj_ZZ );
			}
		}

        // [Kxy], [Kxz], [Kyz] 
		for ( int i=0; i<mNumNodes; i++ )
		{
			for ( int j=0; j<mNumNodes; j++ )
			{
                std::complex<double> intg_Curl_DivA_XY( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_DivA_XZ( 0.00, 0.00 ); 
                std::complex<double> intg_Curl_DivA_YZ( 0.00, 0.00 ); 

		        std::complex<double> intg_Ni_TEp_Nj_XY( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj_XZ( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj_YZ( 0.00, 0.00 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_DivA_XY += W[gp] * ( dNdx[i][gp] * dNdy[j][gp] - dNdy[i][gp] * dNdx[j][gp] );
                    intg_Curl_DivA_XZ += W[gp] * ( dNdx[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdx[j][gp] );
                    intg_Curl_DivA_YZ += W[gp] * ( dNdy[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdy[j][gp] );

					intg_Ni_TEp_Nj_XY += W[gp] * ( N[i][gp] * TEgp[gp][0][1] * N[j][gp] );
                    intg_Ni_TEp_Nj_XZ += W[gp] * ( N[i][gp] * TEgp[gp][0][2] * N[j][gp] );
                    intg_Ni_TEp_Nj_YZ += W[gp] * ( N[i][gp] * TEgp[gp][1][2] * N[j][gp] );
				}
	                
				StiffMatrix[ i             ][ j +   mNumNodes ] = ( invMu * intg_Curl_DivA_XY ) - ( w2 * intg_Ni_TEp_Nj_XY );
                StiffMatrix[ i             ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_DivA_XZ ) - ( w2 * intg_Ni_TEp_Nj_XZ );
                StiffMatrix[ i + mNumNodes ][ j + 2*mNumNodes ] = ( invMu * intg_Curl_DivA_YZ ) - ( w2 * intg_Ni_TEp_Nj_YZ );
			}
		}

        // Upper diagonal [Kvv]
		for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=i; j<mNumNodes; j++ ) 
        	{
                std::complex<double> intg_dNi_TEp_dNj( 0.00, 0.00 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_dNi_TEp_dNj += W[gp] * dNdx[i][gp] * ( TEgp[gp][0][0] * dNdx[j][gp] + TEgp[gp][0][1] * dNdy[j][gp] + TEgp[gp][0][2] * dNdz[j][gp] );
                    intg_dNi_TEp_dNj += W[gp] * dNdy[i][gp] * ( TEgp[gp][1][0] * dNdx[j][gp] + TEgp[gp][1][1] * dNdy[j][gp] + TEgp[gp][1][2] * dNdz[j][gp] );
                    intg_dNi_TEp_dNj += W[gp] * dNdz[i][gp] * ( TEgp[gp][2][0] * dNdx[j][gp] + TEgp[gp][2][1] * dNdy[j][gp] + TEgp[gp][2][2] * dNdz[j][gp] );
				}
	                
				StiffMatrix[ i + 3*mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_dNi_TEp_dNj );
        	}
        }

		// [Kxv], [Kyv], [Kzv]
		for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ ) 
        	{
                std::complex<double> intg_NiX_TEp_dNj( 0.00, 0.00 ); 
                std::complex<double> intg_NiY_TEp_dNj( 0.00, 0.00 ); 
                std::complex<double> intg_NiZ_TEp_dNj( 0.00, 0.00 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_NiX_TEp_dNj += W[gp] * N[i][gp] * ( TEgp[gp][0][0] * dNdx[j][gp] + TEgp[gp][0][1] * dNdy[j][gp] + TEgp[gp][0][2] * dNdz[j][gp] );
                    intg_NiY_TEp_dNj += W[gp] * N[i][gp] * ( TEgp[gp][1][0] * dNdx[j][gp] + TEgp[gp][1][1] * dNdy[j][gp] + TEgp[gp][1][2] * dNdz[j][gp] );
                    intg_NiZ_TEp_dNj += W[gp] * N[i][gp] * ( TEgp[gp][2][0] * dNdx[j][gp] + TEgp[gp][2][1] * dNdy[j][gp] + TEgp[gp][2][2] * dNdz[j][gp] );
				}
	                
				StiffMatrix[ i               ][ j + 3*mNumNodes ] = ( - w2 * intg_NiX_TEp_dNj );
                StiffMatrix[ i +   mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_NiY_TEp_dNj );
                StiffMatrix[ i + 2*mNumNodes ][ j + 3*mNumNodes ] = ( - w2 * intg_NiZ_TEp_dNj );
        	}
        }

        // Lower diagonal 
		Generate_Lower_Diagonal( StiffMatrix );
	}

    //********************************************************************************************************************************
	//* - Generates lower diagonal from the upper diagonal
	//********************************************************************************************************************************
	void VolumeElement_1st_ColdPlasma::Generate_Lower_Diagonal( Matrix<std::complex<double> >& StiffMatrix )
	{
		if ( mIs_IHL )
        {
            for ( int i=0; i<mNumDofs; i++ )
            {
                for ( int j=i+1; j<mNumDofs; j++ )
                {
                    StiffMatrix[j][i] = StiffMatrix[i][j];
                }
            }
        }
        else
        {
		    for ( int i=0; i<mNumDofs; i++ )
		    {
			    for ( int j=i+1; j<mNumDofs; j++ )
			    {
				    StiffMatrix[j][i] = std::conj( StiffMatrix[i][j] );
			    }
		    }
        }
	}
} 