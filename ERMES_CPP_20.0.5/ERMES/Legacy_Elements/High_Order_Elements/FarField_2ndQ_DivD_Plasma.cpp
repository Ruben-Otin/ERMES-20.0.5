
#include "FarField_2ndQ_DivD_Plasma.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //**********************************************************************************
    //* - Global index of the nodes
    //**********************************************************************************	
    void FarField_2ndQ_DivD_Plasma::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize(numDofs);

        for (int i=0; i<numNodes; i++) 
        { 
            EquationId[i             ] = mNodes[i]->pDofcEx()->EquationId();
            EquationId[i+ numNodes   ] = mNodes[i]->pDofcEy()->EquationId();
			EquationId[i+(numNodes*2)] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

    //***********************************************************************************
    //* - Calculation of the Jacobian determinant
    //***********************************************************************************
    void FarField_2ndQ_DivD_Plasma::Calculate_detJ( Vector<double>& detJ, 
                                                    Vector<double>&   cX, 
                                                    Vector<double>&   cY )
    {
		Vector<double> dNda( numNodes );
		Vector<double> dNdb( numNodes );
        
        int numPoints = cX.size();

		detJ.resize( numPoints );

        for (int p=0; p<numPoints; p++) 
        {
			// dNi/da ,dNi/db
			dNda[0] =   4.00*cX[p] + 4.00*cY[p] - 3.00;
			dNdb[0] =   4.00*cX[p] + 4.00*cY[p] - 3.00;

			dNda[1] =   4.00*cX[p] - 1.00;
			dNdb[1] =   0.00;

			dNda[2] =   0.00;
			dNdb[2] =   4.00*cY[p] - 1.00;

			dNda[3] = - 8.00*cX[p] - 4.00*cY[p] + 4.00;
			dNdb[3] = - 4.00*cX[p];

			dNda[4] =   4.00*cY[p];
			dNdb[4] =   4.00*cX[p];

			dNda[5] = - 4.00*cY[p];
			dNdb[5] = - 8.00*cY[p] - 4.00*cX[p] + 4.00;

			double dXda = 0.00;	double dXdb = 0.00;
			double dYda = 0.00;	double dYdb = 0.00;
			double dZda = 0.00;	double dZdb = 0.00;

			for (int i=0; i<numNodes; i++) 
			{
				dXda += dNda[i] * mNodes[i]->X();  dXdb += dNdb[i] * mNodes[i]->X();
				dYda += dNda[i] * mNodes[i]->Y();  dYdb += dNdb[i] * mNodes[i]->Y();         
				dZda += dNda[i] * mNodes[i]->Z();  dZdb += dNdb[i] * mNodes[i]->Z();
			}

			Vector<double> normal(3);

			normal[0] = dZda*dYdb - dYda*dZdb;
			normal[1] = dXda*dZdb - dZda*dXdb;
			normal[2] = dYda*dXdb - dXda*dYdb;

			detJ[p] = sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00));
		}
	}

    //************************************************************************************************* 
    //* -  Calculates the rotation matrix 
    //*************************************************************************************************
    void FarField_2ndQ_DivD_Plasma::Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz  )
    {    
        // Unit vectors
        Vector<double> n(3), t(3), b(3);

        // b = B / |B|
        double Bnorm = std::sqrt( Bx*Bx + By*By + Bz*Bz );

        b[0] = Bx / Bnorm; 
        b[1] = By / Bnorm; 
        b[2] = Bz / Bnorm;

        // n _|_ b
        double nnorm = std::sqrt( b[0]*b[0] + b[1]*b[1] );

        if (nnorm > 0.00)
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
        T.Resize(3,3);

        for ( int row=0; row<3; row++ )
        {
            T[row][0] = std::complex<double>( n[row] ); 
            T[row][1] = std::complex<double>( t[row] ); 
            T[row][2] = std::complex<double>( b[row] );       
        }
    }

    //********************************************************************************************************************
    //* - Calculates permittivity tensor on Gauss points ( extrapolating TEn with N )
    //********************************************************************************************************************
    void FarField_2ndQ_DivD_Plasma::PermittivityTensor_OnGaussPoints( Vector<Matrix<std::complex<double> > >& TEgp, 
                                                                      Matrix<double>& N )
    {
        // Permittivity tensor on nodes
        Vector<Matrix<std::complex<double> > > TEn;
        
        TEn.resize( numNodes );

        // Complex constants
        std::complex<double> cZero( 0.0, 0.0 );
        std::complex<double> cUnit( 0.0, 1.0 );

        for (int i=0; i<numNodes; i++) 
        { 
            TEn[i].Resize(3,3);

            // Plasma parameters on nodes            
            Vector<double> Bext;
            Vector<std::complex<double> > SDPRL;

            mpColdPlasma->Get_Bext_SDPRL_InNode( Bext, SDPRL, mNodes[i], freq );

            double Bx = Bext[0];
            double By = Bext[1];
            double Bz = Bext[2];

            std::complex<double> S  = SDPRL[0];
            std::complex<double> D  = SDPRL[1];
            std::complex<double> P  = SDPRL[2];

            // Permittivity tensor in b = B/|B| coordinates
            Matrix<std::complex<double> > K(3,3);

            K[0][0] = S      ; K[0][1] = -cUnit*D; K[0][2] = cZero; 
            K[1][0] = cUnit*D; K[1][1] =  S      ; K[1][2] = cZero; 
            K[2][0] = cZero  ; K[2][1] =  cZero  ; K[2][2] = P    ;

            K *= eo;

            // Rotation matrix 
            Matrix<std::complex<double> > T; 

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

                    for ( int i=0; i<numNodes; i++ ) 
                    {
                        TEgp[gp][ki][kj] += TEn[i][ki][kj] * N[i][gp];
                    }// end ( int i=0; i<numNodes; i++ ) 
                }// end ( int kj=0; kj<3; kj++ )
            }// end ( int ki=0; ki<3; ki++ )
        }// end ( int gp=0; gp<numGaussPoints; gp++ )
    }

    //********************************************************************************************************************
    //* - Returns a string with the plasma wave type applied to the surface
    //********************************************************************************************************************
    String FarField_2ndQ_DivD_Plasma::Get_Wave_Type()
    {
        String sWaveType;

        Vector<double> vGRC = (*mProperties)(COMPLEX_IBC);

        if      (vGRC[0] == 1.0) sWaveType = "FW";
        else if (vGRC[0] == 2.0) sWaveType = "RW";
        else if (vGRC[0] == 3.0) sWaveType = "LW";
        else if (vGRC[0] == 4.0) sWaveType = "PW";
        else                     sWaveType = "FW";

        return sWaveType;
    }

    //********************************************************************************************************************
    //* - Calculates wave vector tensor on Gauss points ( extrapolating TK with N )
    //********************************************************************************************************************
    void FarField_2ndQ_DivD_Plasma::jwK_OnGaussPoints( Vector<std::complex<double> >& jwKgp, Matrix<double>& N )
    {
        // Wave vectors on nodes
        Vector<std::complex<double> > jwKn;

        jwKn.resize( numNodes );

        // Complex constants
        std::complex<double> cZero( 0.0,  0.0 );
        std::complex<double> rUnit( 1.0,  0.0 );
        std::complex<double> j_w  ( 0.0, freq );

        for (int i=0; i<numNodes; i++) 
        { 
            // Plasma parameters on nodes
            Vector<double> Bext; 
            Vector<std::complex<double> > SDPRL;

            mpColdPlasma->Get_Bext_SDPRL_InNode( Bext, SDPRL, mNodes[i], freq );

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

            // Wave vector jw*K_pla
            std::complex<double> K_pla;

            if ( sWaveType == "FW" )
            {
                if (abs(S) > 0.0) K_pla = sqrt( mu * eo * ( R * L / S ) );
                else              K_pla = sqrt( mu * eo * ( rUnit     ) );
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

            for ( int i=0; i<numNodes; i++ ) 
            {
                jwKgp[gp] += jwKn[i] * N[i][gp];
            }
        }
    }

    //****************************************************************************************
    //* - Calculates exterior normal on Gauss points
    //****************************************************************************************	
    void FarField_2ndQ_DivD_Plasma::ExtNormal_OnGaussPoints( Vector<Vector<double> >& eNgp, 
                                                             Vector<double>& cX, 
                                                             Vector<double>& cY )
    {
        // Surface derivatives
        Vector<double> dNda( numNodes );
        Vector<double> dNdb( numNodes );

        int numGaussPoints = cX.size();

        eNgp.resize( numGaussPoints );

        for (int gp=0; gp<numGaussPoints; gp++) 
        {
            // dNi/da ,dNi/db
            dNda[0] =   4.00*cX[gp] + 4.00*cY[gp] - 3.00;
            dNdb[0] =   4.00*cX[gp] + 4.00*cY[gp] - 3.00;

            dNda[1] =   4.00*cX[gp] - 1.00;
            dNdb[1] =   0.00;

            dNda[2] =   0.00;
            dNdb[2] =   4.00*cY[gp] - 1.00;

            dNda[3] = - 8.00*cX[gp] - 4.00*cY[gp] + 4.00;
            dNdb[3] = - 4.00*cX[gp];

            dNda[4] =   4.00*cY[gp];
            dNdb[4] =   4.00*cX[gp];

            dNda[5] = - 4.00*cY[gp];
            dNdb[5] = - 8.00*cY[gp] - 4.00*cX[gp] + 4.00;

            double dXda = 0.00;	double dXdb = 0.00;
            double dYda = 0.00;	double dYdb = 0.00;
            double dZda = 0.00;	double dZdb = 0.00;

            for (int i=0; i<numNodes; i++) 
            {
                dXda += dNda[i] * mNodes[i]->X();  dXdb += dNdb[i] * mNodes[i]->X();
                dYda += dNda[i] * mNodes[i]->Y();  dYdb += dNdb[i] * mNodes[i]->Y();         
                dZda += dNda[i] * mNodes[i]->Z();  dZdb += dNdb[i] * mNodes[i]->Z();
            }

            Vector<double> normal(3);

            normal[0] = dZda*dYdb - dYda*dZdb;
            normal[1] = dXda*dZdb - dZda*dXdb;
            normal[2] = dYda*dXdb - dXda*dYdb;

            normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));

            eNgp[gp].resize(3);

            eNgp[gp] = normal;
        }
    }

    //****************************************************************************************
    //* - Calculation of the stiffness matrix
    //****************************************************************************************	
    void FarField_2ndQ_DivD_Plasma::GetStiffMatrix( Matrix<std::complex<double> >& StiffMatrix ) 
    {
        // Useful constants
        std::complex<double> invMu( 1.00/mu, 0.00 );
        std::complex<double> cZero( 0.00   , 0.00 );
        
        // Loop indexes
		int i, j, gp;

        // Resizing Stiffness matrix
        StiffMatrix.Resize( numDofs, numDofs, cZero );

        // Gauss points and weights
		Vector<double> cX; 
		Vector<double> cY;
		Vector<double>  W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

        // Jacobian
		Vector<double> detJ;

        Calculate_detJ( detJ, cX, cY );

		for ( gp=0; gp<numGaussPoints; gp++ ) W[gp] *= detJ[gp]; 

        // Lagrangian base on Gauss points
		Matrix<double> N;

		Lagrange2D_Ni_2nd( N, cX, cY );

        // Permittivity tensor on Gauss points
        Vector<Matrix<std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, N );

         // Wave vector tensor on Gauss points (interpolated with N)
        Vector<std::complex<double> > jwK;

        jwK_OnGaussPoints( jwK, N );

        // Regularization constant on Gauss points
        Vector<std::complex<double> > invNormTr2( numGaussPoints );

        for ( gp=0; gp<numGaussPoints; gp++ ) 
        {
            std::complex<double> NormTraza2 = ( TEgp[gp][0][0] * std::conj( TEgp[gp][0][0] ) + TEgp[gp][0][1] * std::conj( TEgp[gp][0][1] ) + TEgp[gp][0][2] * std::conj( TEgp[gp][0][2] ) + 
                                                TEgp[gp][1][0] * std::conj( TEgp[gp][1][0] ) + TEgp[gp][1][1] * std::conj( TEgp[gp][1][1] ) + TEgp[gp][1][2] * std::conj( TEgp[gp][1][2] ) + 
                                                TEgp[gp][2][0] * std::conj( TEgp[gp][2][0] ) + TEgp[gp][2][1] * std::conj( TEgp[gp][2][1] ) + TEgp[gp][2][2] * std::conj( TEgp[gp][2][2] ) );             

            invNormTr2[gp] = 3.0 / ( NormTraza2 );
        }

        // Exterior normal on Gauss points
        Vector<Vector<double> > eNgp;
        
        ExtNormal_OnGaussPoints( eNgp, cX, cY );
  
        //////////////////////////////////////     [Kxx]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;
                
                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );
               
                int ci = 0; int cj = 0;
                
                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj] - 1;

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]) ;
                }

                StiffMatrix[i][j] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        //////////////////////////////////////     [Kyy]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;

                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );

                int ci = 1; int cj = 1;

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj] - 1;

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]) ;
                }

                StiffMatrix[i+numNodes][j+numNodes] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        //////////////////////////////////////     [Kzz]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;

                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );

                int ci = 2; int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj] - 1;

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]) ;
                }

                StiffMatrix[i+(numNodes*2)][j+(numNodes*2)] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        //////////////////////////////////////     [Kxy]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;

                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );

                int ci = 0; int cj = 1;

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj];

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]) ;
                }

                StiffMatrix[i][j+numNodes] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        //////////////////////////////////////     [Kxz]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;

                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );

                int ci = 0; int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj];

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]) ;
                }

                StiffMatrix[i][j+(numNodes*2)] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        //////////////////////////////////////     [Kyz]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;

                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );

                int ci = 1; int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj];

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]) ;
                }

                StiffMatrix[i+numNodes][j+(numNodes*2)] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        ////////////////////////////////////////     [Kyx]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;

                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );

                int ci = 1; int cj = 0;

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj];

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]);
                }

                StiffMatrix[i+numNodes][j] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        //////////////////////////////////////     [Kzx]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;

                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );

                int ci = 2; int cj = 0;

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj];

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]);
                }

                StiffMatrix[i+(2*numNodes)][j] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        //////////////////////////////////////     [Kzy]     //////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                double nxn; 

                std::complex<double> nTEi;
                std::complex<double> nTEj;

                std::complex<double> int_nxNi_jwK_nxNj( 0.00, 0.00 );
                std::complex<double> int_neNi_rct_neNj( 0.00, 0.00 );

                int ci = 2; int cj = 1;

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    nxn = eNgp[gp][ci]*eNgp[gp][cj];

                    int_nxNi_jwK_nxNj += W[gp] * jwK[gp] * nxn * ( N[i][gp] * N[j][gp] );

                    nTEi = eNgp[gp][0] * TEgp[gp][0][ci] + eNgp[gp][1] * TEgp[gp][1][ci] + eNgp[gp][2] * TEgp[gp][2][ci];
                    nTEj = eNgp[gp][0] * TEgp[gp][0][cj] + eNgp[gp][1] * TEgp[gp][1][cj] + eNgp[gp][2] * TEgp[gp][2][cj];

                    int_neNi_rct_neNj += W[gp] * ( jwK[gp] * invNormTr2[gp] ) * ( std::conj(nTEi) *  N[i][gp] ) * ( nTEj * N[j][gp]) ;
                }

                StiffMatrix[i+(2*numNodes)][j+numNodes] = invMu * ( int_nxNi_jwK_nxNj - int_neNi_rct_neNj ); 
            }
        }

        ////////////////////////////// Lower diagonal block ///////////////////////////////
        //for ( i=0; i<numDofs; i++ )
        //{
        //    for ( j=i+1; j<numDofs; j++ )
        //    {
        //        StiffMatrix[j][i] = StiffMatrix[i][j];
        //    }
        //}
    }
 
} /* end namespace Kratos */