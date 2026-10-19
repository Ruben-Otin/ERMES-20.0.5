
#include "FarField_1st_ColdPlasma.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //********************************************************************************************************************************
    //* - Global index of the DOFs
    //********************************************************************************************************************************	
    void FarField_1st_ColdPlasma::GetEquationIdVector( Vector<int>& EquationId )
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
    //* - Area of the element
    //********************************************************************************************************************************	
    double FarField_1st_ColdPlasma::Calculate_Area() 
    {
		double na[3], v1[3], v2[3];

        v2[0] = mNodes[2]->X() - mNodes[0]->X();
        v2[1] = mNodes[2]->Y() - mNodes[0]->Y();
        v2[2] = mNodes[2]->Z() - mNodes[0]->Z();

        v1[0] = mNodes[1]->X() - mNodes[0]->X();
        v1[1] = mNodes[1]->Y() - mNodes[0]->Y();
        v1[2] = mNodes[1]->Z() - mNodes[0]->Z();

        // Area = 0.5 * ( v2 x v1 )
        na[0] = v2[1]*v1[2] - v2[2]*v1[1];
        na[1] = v2[2]*v1[0] - v2[0]*v1[2];
        na[2] = v2[0]*v1[1] - v2[1]*v1[0];

        return ( 0.5 * sqrt( na[0]*na[0] + na[1]*na[1] + na[2]*na[2] ) );
    }

    //********************************************************************************************************************************
    //* - Exterior normal 
    //********************************************************************************************************************************	
    void FarField_1st_ColdPlasma::Calculate_ExtNormal( Vector<double>& extN )
    {
        double v1[3], v2[3];

        v2[0] = mNodes[2]->X() - mNodes[0]->X();
        v2[1] = mNodes[2]->Y() - mNodes[0]->Y();
        v2[2] = mNodes[2]->Z() - mNodes[0]->Z();

        v1[0] = mNodes[1]->X() - mNodes[0]->X();
        v1[1] = mNodes[1]->Y() - mNodes[0]->Y();
        v1[2] = mNodes[1]->Z() - mNodes[0]->Z();

        extN.resize(3);

        extN[0] = v2[1]*v1[2] - v2[2]*v1[1];
        extN[1] = v2[2]*v1[0] - v2[0]*v1[2];
        extN[2] = v2[0]*v1[1] - v2[1]*v1[0];

        double inv_extN_norm = 1.00 / sqrt( extN[0]*extN[0] + extN[1]*extN[1] + extN[2]*extN[2] );

		extN *= inv_extN_norm;
    }

    //********************************************************************************************************************************
    //* - Returns the wave type applied to the cold plasma far field surface
    //********************************************************************************************************************************
    String FarField_1st_ColdPlasma::Get_Wave_Type()
    {
        String sWaveType;

        Vector<double> vGRC = (*mProperties)( COMPLEX_IBC );

        if      ( vGRC[0] == 1.0 ) sWaveType = "FW";
        else if ( vGRC[0] == 2.0 ) sWaveType = "RW";
        else if ( vGRC[0] == 3.0 ) sWaveType = "LW";
        else if ( vGRC[0] == 4.0 ) sWaveType = "PW";
        else                       sWaveType = "FW";

        return sWaveType;
    }

    //******************************************************************************************************************************** 
    //* -  Calculates rotation matrix 
    //********************************************************************************************************************************
    void FarField_1st_ColdPlasma::Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz )
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

        if ( nnorm > 0.0 )
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
    //* - Calculates permittivity tensor on Gauss points ( extrapolating TEn to GPs with N )
    //********************************************************************************************************************************
    void FarField_1st_ColdPlasma::PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp,
                                                                    Matrix< double >& N )
    {
        // Permittivity tensor on nodes
        Vector< Matrix< std::complex<double> > > TEn( mNumNodes );
        
        // Constants
        double eo = 8.8541878176e-12;

        std::complex<double> cZero( 0.0, 0.0 );
        std::complex<double> cUnit( 0.0, 1.0 );

        for ( int i=0; i<mNumNodes; i++ ) 
        { 
            TEn[i].Resize( 3, 3 );

            // Plasma parameters on nodes            
            Vector< double > Bext;

            Vector< std::complex<double> > SDPRL;

            mpColdPlasma->Get_Bext_SDPRL_InNode( Bext, SDPRL, mNodes[i], mFreq );

            double Bx = Bext[0];
            double By = Bext[1];
            double Bz = Bext[2];

            std::complex<double> S = SDPRL[0];
            std::complex<double> D = SDPRL[1];
            std::complex<double> P = SDPRL[2];

            // Permittivity tensor in b = B/|B| coordinates
            Matrix< std::complex<double> > K( 3, 3 );

            K[0][0] = S      ; K[0][1] = -cUnit*D; K[0][2] = cZero; 
            K[1][0] = cUnit*D; K[1][1] =  S      ; K[1][2] = cZero; 
            K[2][0] = cZero  ; K[2][1] =  cZero  ; K[2][2] = P    ;

            K *= eo;

            // Rotation matrix 
            Matrix< std::complex<double> > T; 

            Rotation_Matrix( T, Bx, By, Bz );

            // TEn = T * K * Tt 
            // Be carefull!!, T.Transpose() change T permanently
            TEn[i] = T * K;
            TEn[i] = TEn[i] * T.Transpose();
        }

        // Number of Gauss points
        int numGaussPoints = N[0].size();

        // Initializing vector of 3x3 permittivity tensors 
        TEgp.resize( numGaussPoints );

        // Permittivity tensor on Gauss points  
        for ( int gp=0; gp<numGaussPoints; gp++ )
        {
            TEgp[gp].Resize( 3, 3 );

            for ( int ki=0; ki<3; ki++ )
            {
                for ( int kj=0; kj<3; kj++ )
                {
                    TEgp[gp][ki][kj] = cZero;

                    for ( int i=0; i<mNumNodes; i++ ) 
                    {
                        TEgp[gp][ki][kj] += TEn[i][ki][kj] * N[i][gp];
                    }
                }
            }
        }
    }

    //********************************************************************************************************************************
    //* - Calculates wave vector tensor on Gauss points ( extrapolating TK to GPs with N )
    //********************************************************************************************************************************
    void FarField_1st_ColdPlasma::jwK_OnGaussPoints( Vector< std::complex<double> >& jwKgp, 
                                                     Matrix< double >& N )
    {
        // Wave vectors on nodes
        Vector< std::complex<double> > jwKn( mNumNodes );

        // Constants
        double eo = 8.8541878176e-12;
        double pi = 3.14159265358979;
        double mu = pi * 4.0e-7;

        std::complex<double> cZero( 0.0,   0.0 );
        std::complex<double> rUnit( 1.0,   0.0 );
        std::complex<double>   j_w( 0.0, mFreq );

        for ( int i=0; i<mNumNodes; i++ ) 
        { 
            // Plasma parameters on nodes
            Vector< double > Bext; 

            Vector< std::complex<double> > SDPRL;

            mpColdPlasma->Get_Bext_SDPRL_InNode( Bext, SDPRL, mNodes[i], mFreq );

            double Bx = Bext[0];
            double By = Bext[1];
            double Bz = Bext[2];

            std::complex<double> S = SDPRL[0];
            std::complex<double> D = SDPRL[1];
            std::complex<double> P = SDPRL[2];
            std::complex<double> R = SDPRL[3];
            std::complex<double> L = SDPRL[4];

            // Wave type
            String sWaveType = Get_Wave_Type();

            // Wave vector jw * K_pla
            std::complex<double> K_pla;

            if ( sWaveType == "FW" )
            {
                if ( abs( S ) > 0.0 ) 
                {
                    K_pla = sqrt( mu * eo * ( R * L / S ) );
                }
                else   
                {
                    K_pla = sqrt( mu * eo * ( rUnit ) );
                }
            }
            else if ( sWaveType == "RW" )
            {
                K_pla = sqrt( mu * eo * R );
            }
            else if ( sWaveType == "LW" )
            {
                K_pla = sqrt( mu * eo * L );
            }
            else if ( sWaveType == "PW" )
            {
                K_pla = sqrt( mu * eo * P );
            }

            // Nodal wave vector
            jwKn[i] = j_w * K_pla;
        }	

        // Number of Gauss points
        int numGaussPoints = N[0].size();

        // Wave vector on gauss points
        jwKgp.resize( numGaussPoints );

        for ( int gp=0; gp<numGaussPoints; gp++ )
        {
            jwKgp[gp] = cZero;

            for ( int i=0; i<mNumNodes; i++ ) 
            {
                jwKgp[gp] += jwKn[i] * N[i][gp];
            }
        }
    }

    //******************************************************************************************************************************** 
    //* - Inverse of the regularization weight tau
    //********************************************************************************************************************************
    void FarField_1st_ColdPlasma::Get_invTau_InGaussPoints( Vector<         std::complex<double>   >& invTau, 
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
    //* - Element stiffness matrix
    //********************************************************************************************************************************	
    void FarField_1st_ColdPlasma::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );	

		if ( mPotentialsOn ) 
        {
            Get_Stiffness_Matrix_AV( StiffMatrix ); 
        }
        else                 
        {
            Get_Stiffness_Matrix_Ef( StiffMatrix ); 
        }
    }
		
    //********************************************************************************************************************************
    //* - Element stiffness matrix for the electric field formulation
    //********************************************************************************************************************************	
    void FarField_1st_ColdPlasma::Get_Stiffness_Matrix_Ef( Matrix< std::complex<double> >& StiffMatrix ) 
    {
        // Constants
        double pi = 3.14159265358979;
        double mu = pi * 4.0e-7;

        std::complex<double> invMu( 1.0 / mu, 0.0 );
 
        // Gauss points and weights
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order10( cX, cY, W );

		double detJ = 2.0 * Calculate_Area();

		for ( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            W[gp] *= detJ; 
        }

        // Lagrangian base on Gauss points
		Matrix<double> N;

		Lagrange2D_Ni_1st( N, cX, cY );

        // Permittivity tensor on Gauss points
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, N );

        // Wave vector tensor on Gauss points 
        Vector< std::complex<double> > jwK;

        jwK_OnGaussPoints( jwK, N );

        // Regularization constant on Gauss points
        Vector< std::complex<double> > invTau;

        Get_invTau_InGaussPoints( invTau, TEgp );

        // Exterior normal 
        Vector<double> extN; 
        
        Calculate_ExtNormal( extN );

        // Stiffness matrix
        for ( int i=0; i<mNumNodes; i++ )
        {
            for ( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> eNTE_X, eNTE_Y, eNTE_Z;

                std::complex<double> intg_nxNi_jwK_nxNj_XX( 0.0, 0.0 );
                std::complex<double> intg_nxNi_jwK_nxNj_YY( 0.0, 0.0 );
                std::complex<double> intg_nxNi_jwK_nxNj_ZZ( 0.0, 0.0 );

                std::complex<double> intg_nxNi_jwK_nxNj_XY( 0.0, 0.0 );
                std::complex<double> intg_nxNi_jwK_nxNj_XZ( 0.0, 0.0 );
                std::complex<double> intg_nxNi_jwK_nxNj_YZ( 0.0, 0.0 );

                std::complex<double> intg_neNi_rct_neNj_XX( 0.0, 0.0 );
                std::complex<double> intg_neNi_rct_neNj_YY( 0.0, 0.0 );
                std::complex<double> intg_neNi_rct_neNj_ZZ( 0.0, 0.0 );

                std::complex<double> intg_neNi_rct_neNj_XY( 0.0, 0.0 );
                std::complex<double> intg_neNi_rct_neNj_XZ( 0.0, 0.0 );
                std::complex<double> intg_neNi_rct_neNj_YZ( 0.0, 0.0 );

                std::complex<double> intg_neNi_rct_neNj_YX( 0.0, 0.0 );
                std::complex<double> intg_neNi_rct_neNj_ZX( 0.0, 0.0 );
                std::complex<double> intg_neNi_rct_neNj_ZY( 0.0, 0.0 );
               
                for ( int gp=0; gp<numGaussPoints; gp++ ) 
                {
                    eNTE_X = extN[0] * TEgp[gp][0][0] + extN[1] * TEgp[gp][1][0] + extN[2] * TEgp[gp][2][0];
                    eNTE_Y = extN[0] * TEgp[gp][0][1] + extN[1] * TEgp[gp][1][1] + extN[2] * TEgp[gp][2][1];
                    eNTE_Z = extN[0] * TEgp[gp][0][2] + extN[1] * TEgp[gp][1][2] + extN[2] * TEgp[gp][2][2];

                    intg_nxNi_jwK_nxNj_XX += W[gp] * jwK[gp] * ( extN[0] * extN[0] - 1 ) * N[i][gp] * N[j][gp];
                    intg_nxNi_jwK_nxNj_YY += W[gp] * jwK[gp] * ( extN[1] * extN[1] - 1 ) * N[i][gp] * N[j][gp];
                    intg_nxNi_jwK_nxNj_ZZ += W[gp] * jwK[gp] * ( extN[2] * extN[2] - 1 ) * N[i][gp] * N[j][gp];
                                                                          
                    intg_nxNi_jwK_nxNj_XY += W[gp] * jwK[gp] * ( extN[0] * extN[1]     ) * N[i][gp] * N[j][gp];
                    intg_nxNi_jwK_nxNj_XZ += W[gp] * jwK[gp] * ( extN[0] * extN[2]     ) * N[i][gp] * N[j][gp];
                    intg_nxNi_jwK_nxNj_YZ += W[gp] * jwK[gp] * ( extN[1] * extN[2]     ) * N[i][gp] * N[j][gp];

                    intg_neNi_rct_neNj_XX += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_X ) * eNTE_X * N[i][gp] * N[j][gp];
                    intg_neNi_rct_neNj_YY += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_Y ) * eNTE_Y * N[i][gp] * N[j][gp];
                    intg_neNi_rct_neNj_ZZ += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_Z ) * eNTE_Z * N[i][gp] * N[j][gp];
                                                                                                
                    intg_neNi_rct_neNj_XY += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_X ) * eNTE_Y * N[i][gp] * N[j][gp];
                    intg_neNi_rct_neNj_XZ += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_X ) * eNTE_Z * N[i][gp] * N[j][gp];
                    intg_neNi_rct_neNj_YZ += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_Y ) * eNTE_Z * N[i][gp] * N[j][gp];
                                                                                                                              
                    intg_neNi_rct_neNj_YX += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_Y ) * eNTE_X * N[i][gp] * N[j][gp];
                    intg_neNi_rct_neNj_ZX += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_Z ) * eNTE_X * N[i][gp] * N[j][gp];
                    intg_neNi_rct_neNj_ZY += W[gp] * jwK[gp] * invTau[gp] * std::conj( eNTE_Z ) * eNTE_Y * N[i][gp] * N[j][gp];
                }

                StiffMatrix[ i               ][ j               ] = invMu * ( intg_nxNi_jwK_nxNj_XX - intg_neNi_rct_neNj_XX ); 
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = invMu * ( intg_nxNi_jwK_nxNj_YY - intg_neNi_rct_neNj_YY ); 
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = invMu * ( intg_nxNi_jwK_nxNj_ZZ - intg_neNi_rct_neNj_ZZ ); 

                StiffMatrix[ i               ][ j +   mNumNodes ] = invMu * ( intg_nxNi_jwK_nxNj_XY - intg_neNi_rct_neNj_XY ); 
                StiffMatrix[ i               ][ j + 2*mNumNodes ] = invMu * ( intg_nxNi_jwK_nxNj_XZ - intg_neNi_rct_neNj_XZ ); 
                StiffMatrix[ i +   mNumNodes ][ j + 2*mNumNodes ] = invMu * ( intg_nxNi_jwK_nxNj_YZ - intg_neNi_rct_neNj_YZ ); 

                StiffMatrix[ i +   mNumNodes ][ j               ] = invMu * ( intg_nxNi_jwK_nxNj_XY - intg_neNi_rct_neNj_YX ); 
                StiffMatrix[ i + 2*mNumNodes ][ j               ] = invMu * ( intg_nxNi_jwK_nxNj_XZ - intg_neNi_rct_neNj_ZX ); 
                StiffMatrix[ i + 2*mNumNodes ][ j +   mNumNodes ] = invMu * ( intg_nxNi_jwK_nxNj_YZ - intg_neNi_rct_neNj_ZY ); 
            }
        }
    }

    //********************************************************************************************************************************
    //* - Element stiffness matrix for the AV-potential formulation
    //********************************************************************************************************************************	
    void FarField_1st_ColdPlasma::Get_Stiffness_Matrix_AV( Matrix< std::complex<double> >& StiffMatrix ) 
    {
        // Constants
        double pi = 3.14159265358979;
        double mu = pi * 4.0e-7;

        std::complex<double> invMu( 1.0 / mu, 0.0 );
 
        // Gauss points and weights
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order10( cX, cY, W );

		double detJ = 2.0 * Calculate_Area();

		for ( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            W[gp] *= detJ; 
        }

        // Lagrangian base on Gauss points
		Matrix<double> N;

		Lagrange2D_Ni_1st( N, cX, cY );

        // Wave vector tensor on Gauss points 
        Vector< std::complex<double> > jwK;

        jwK_OnGaussPoints( jwK, N );

        // Stiffness matrix
        for ( int i=0; i<mNumNodes; i++ )
        {
            for ( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Ni_jwK_Nj( 0.0, 0.0 );
                
                for ( int gp=0; gp<numGaussPoints; gp++ ) 
                {
                    intg_Ni_jwK_Nj += W[gp] * jwK[gp] * N[i][gp] * N[j][gp];
                }

                StiffMatrix[ i               ][ j               ] = ( - invMu * intg_Ni_jwK_Nj ); 
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = ( - invMu * intg_Ni_jwK_Nj ); 
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = ( - invMu * intg_Ni_jwK_Nj ); 
            }
        }
    }
} 