
#include "FarField_2ndQ_Plasma.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //**********************************************************************************
    //* - Global index of the nodes
    //**********************************************************************************	
    void FarField_2ndQ_Plasma::GetEquationIdVector( Vector<int>& EquationId )
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
    void FarField_2ndQ_Plasma::Calculate_detJ( Vector<double>& detJ, 
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
    void FarField_2ndQ_Plasma::Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz  )
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

    //**********************************************************************************************************************
    //* - Calculates wave vector tensor on nodes 
    //**********************************************************************************************************************
    void FarField_2ndQ_Plasma::Tensor_K_OnNodes( Vector<Matrix<std::complex<double> > >& TKn )
    {
        // Permittivity tensor on nodes
        TKn.resize( numNodes );

        // Complex constants
        std::complex<double> cZero( 0.0,  0.0 );
        std::complex<double> cUnit( 0.0,  1.0 );
        std::complex<double> rUnit( 1.0,  0.0 );
        std::complex<double> j_w  ( 0.0, freq );

        for (int i=0; i<numNodes; i++) 
        { 
            TKn[i].Resize( 3, 3, cZero );

            // Plasma parameters on nodes
            Vector<double> PlasmaPn = mpColdPlasma->Get_AllPlasmaParameters_InNode( mNodes[i], freq );

            double Bx = PlasmaPn[1];
            double By = PlasmaPn[2];
            double Bz = PlasmaPn[3];

            double S = PlasmaPn[4];
            double D = PlasmaPn[5];
            double P = PlasmaPn[6];

            double R = PlasmaPn[7];
            double L = PlasmaPn[8];

            // Wave vector components
            double K_per2;
            double K_par2;

            std::complex<double> K_per;
            std::complex<double> K_par; 

            // Perpendicular wave vector K 
            if      ( mpColdPlasma->Far_Field_Wave_Type() == 'R' ) K_per2 = mu * eo * R;
            else if ( mpColdPlasma->Far_Field_Wave_Type() == 'L' ) K_per2 = mu * eo * L;
            else                                                   K_per2 = mu * eo * R;

            if ( K_per2 < 0.0 ) K_per = cUnit * sqrt(-K_per2);
            else                K_per = rUnit * sqrt(+K_per2);

            // Parallel wave vector K 
            K_par2 = K_per2; //mu * eo * P;

            if ( K_par2 < 0.0 ) K_par = cUnit * sqrt(-K_par2);
            else                K_par = rUnit * sqrt(+K_par2);

            // Wave vector tensor in b = B/|B| coordinates
            Matrix<std::complex<double> > Kb; 
            Kb.Resize( 3, 3, cZero );

            Kb[0][0] = j_w * K_per;   
            Kb[1][1] = j_w * K_per;  
            Kb[2][2] = j_w * K_par;

            // Rotation matrix 
            Matrix<std::complex<double> > T; 

            Rotation_Matrix( T, Bx, By, Bz );

            // TKn = T * Kb * Tt 
            // Be carefull!!, T.Transpose() change T permanently
            TKn[i] = T * Kb;
            TKn[i] = TKn[i] * T.Transpose();
        }	
    }

    //********************************************************************************************************************
    //* - Calculates wave vector tensor on Gauss points ( extrapolating TK with N )
    //********************************************************************************************************************
    void FarField_2ndQ_Plasma::Tensor_K_OnGaussPoints( Vector<Matrix<std::complex<double> > >& TKgp, Matrix<double>& N )
    {
        // Extracting permitivitty tensor on nodes
        Vector<Matrix<std::complex<double> > > TKn;

        Tensor_K_OnNodes( TKn );

        // Number of Gauss points
        int numGaussPoints = N[0].size();

        // Initializing vector of 3x3 permittivity tensors 
        TKgp.resize( numGaussPoints );

        // Complex constants
        std::complex<double> cZero( 0.0, 0.0 );

        // Permittivity tensor on Gauss points  
        for ( int gp=0; gp<numGaussPoints; gp++ )
        {
            TKgp[gp].Resize( 3, 3, cZero );

            for ( int ki=0; ki<3; ki++ )
            {
                for ( int kj=0; kj<3; kj++ )
                {
                    for ( int i=0; i<numNodes; i++ ) 
                    {
                        TKgp[gp][ki][kj] += TKn[i][ki][kj] * N[i][gp];
                    }// end ( int i=0; i<numNodes; i++ ) 
                }// end ( int kj=0; kj<3; kj++ )
            }// end ( int ki=0; ki<3; ki++ )
        }// end ( int gp=0; gp<numGaussPoints; gp++ )
    }

    //****************************************************************************************
    //* - Calculation of the stiffness matrix
    //****************************************************************************************	
    void FarField_2ndQ_Plasma::GetStiffMatrix( Matrix<std::complex<double> >& StiffMatrix ) 
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

		int numGaussPoints = GaussPoints2D_Order8( cX, cY, W );

        // Jacobian
		Vector<double> detJ;

        Calculate_detJ( detJ, cX, cY );

		for ( gp=0; gp<numGaussPoints; gp++ ) W[gp] *= detJ[gp]; 

        // Lagrangian base on Gauss points
		Matrix<double> N;

		Lagrange2D_Ni_2nd( N, cX, cY );

        // Resizing Stiffness matrix
        StiffMatrix.Resize( numDofs, numDofs );

        // Wave vector tensor on Gauss points (interpolated with N)
        Vector<Matrix<std::complex<double> > > TK;

        Tensor_K_OnGaussPoints( TK, N );

        /////////////////////////// [Kxx], [Kyy], [Kzz] ///////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=i; j<numNodes; j++ )
            {
                std::complex<double> int_NiKNj_XX( 0.00, 0.00 );
                std::complex<double> int_NiKNj_YY( 0.00, 0.00 );
                std::complex<double> int_NiKNj_ZZ( 0.00, 0.00 );

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    int_NiKNj_XX += W[gp] * ( N[i][gp] * TK[gp][0][0] * N[j][gp] );
                    int_NiKNj_YY += W[gp] * ( N[i][gp] * TK[gp][1][1] * N[j][gp] );
                    int_NiKNj_ZZ += W[gp] * ( N[i][gp] * TK[gp][2][2] * N[j][gp] );  
                }

                StiffMatrix[i             ][j             ] = -invMu * int_NiKNj_XX;
                StiffMatrix[i+ numNodes   ][j+ numNodes   ] = -invMu * int_NiKNj_YY;
                StiffMatrix[i+(numNodes*2)][j+(numNodes*2)] = -invMu * int_NiKNj_ZZ;       
            }
        }
        //////////////////////////////// [Kxy] //////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> int_NiKNj_XY( 0.00, 0.00 );

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    int_NiKNj_XY += W[gp] * ( N[i][gp] * TK[gp][0][1] * N[j][gp] );
                }

                StiffMatrix[i][j+numNodes] = -invMu * int_NiKNj_XY;
            }
        }
        //////////////////////////////// [Kxz] //////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> int_NiKNj_XZ( 0.00, 0.00 );

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    int_NiKNj_XZ += W[gp] * ( N[i][gp] * TK[gp][0][2] * N[j][gp] );
                }

                StiffMatrix[i][j+(numNodes*2)] = -invMu * int_NiKNj_XZ;
            }
        }
        //////////////////////////////// [Kyz] //////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> int_NiKNj_YZ( 0.00, 0.00 );

                for ( gp=0; gp<numGaussPoints; gp++ ) 
                {
                    int_NiKNj_YZ += W[gp] * ( N[i][gp] * TK[gp][1][2] * N[j][gp] );    
                }

                StiffMatrix[i+numNodes][j+(numNodes*2)] = -invMu * int_NiKNj_YZ;
            }
        }
        //////////////////////////// Lower diagonal block ///////////////////////////////
        for ( i=0; i<numDofs; i++ )
        {
            for ( j=i+1; j<numDofs; j++ )
            {
                StiffMatrix[j][i] = StiffMatrix[i][j];
            }
        }
    }
 
} /* end namespace Kratos */