
#include "GenericRobin_3sb_ColdPlasma.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //*************************************************************************************************************************************
    //* - Global index of the nodes
    //*************************************************************************************************************************************	
    void GenericRobin_3sb_ColdPlasma::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( numDofs );

        for ( int i=0; i<numNodes; i++ ) 
        { 
            EquationId[i             ] = mNodes[i]->pDofcEx()->EquationId();
            EquationId[i + numNodes  ] = mNodes[i]->pDofcEy()->EquationId();
			EquationId[i + numNodes*2] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

    //************************************************************************************************************************************* 
    //* -  Calculates the rotation matrix 
    //*************************************************************************************************************************************
    void GenericRobin_3sb_ColdPlasma::Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz )
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

    //*************************************************************************************************************************************
    //* - Calculates permittivity tensor on Gauss points ( extrapolating TEn with N )
    //*************************************************************************************************************************************
    void GenericRobin_3sb_ColdPlasma::PermittivityTensor_OnGaussPoints( Vector<Matrix<std::complex<double> > >& TEgp, Matrix<double>& N )
    {
        // Permittivity tensor on nodes
        Vector<Matrix<std::complex<double> > > TEn;
        
        TEn.resize( numNodes );

        // Complex constants
        std::complex<double> cZero( 0.0, 0.0 );
        std::complex<double> cUnit( 0.0, 1.0 );

        double eo = 8.8541878176e-12;

        for (int i=0; i<numNodes; i++) 
        { 
            TEn[i].Resize(3,3);

            // Plasma parameters on nodes            
            Vector<double> Bext;
            Vector< std::complex<double> > SDPRL;

            mpColdPlasma->Get_Bext_SDPRL_InNode( Bext, SDPRL, mNodes[i], freq );

            double Bx = Bext[0];
            double By = Bext[1];
            double Bz = Bext[2];

            std::complex<double> S = SDPRL[0];
            std::complex<double> D = SDPRL[1];
            std::complex<double> P = SDPRL[2];

            // Permittivity tensor in b = B/|B| coordinates
            Matrix< std::complex<double> > K(3,3);

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

                    for ( int i=0; i<numNodes; i++ ) 
                    {
                        TEgp[gp][ki][kj] += TEn[i][ki][kj] * N[i][gp];
                    }// end ( int i=0; i<numNodes; i++ ) 
                }// end ( int kj=0; kj<3; kj++ )
            }// end ( int ki=0; ki<3; ki++ )
        }// end ( int gp=0; gp<numGaussPoints; gp++ )
    }

    //*************************************************************************************************************************************
    //* - Calculates the area of the element
    //*************************************************************************************************************************************	
    double GenericRobin_3sb_ColdPlasma::Calculate_Area() 
    {
		double n[3], v1[3], v2[3];

        v2[0] = mNodes[2]->X() - mNodes[0]->X();
        v2[1] = mNodes[2]->Y() - mNodes[0]->Y();
        v2[2] = mNodes[2]->Z() - mNodes[0]->Z();

        v1[0] = mNodes[1]->X() - mNodes[0]->X();
        v1[1] = mNodes[1]->Y() - mNodes[0]->Y();
        v1[2] = mNodes[1]->Z() - mNodes[0]->Z();

        //area = 0.5*(v2 x v1)
        n[0] = v2[1]*v1[2] - v2[2]*v1[1];
        n[1] = v2[2]*v1[0] - v2[0]*v1[2];
        n[2] = v2[0]*v1[1] - v2[1]*v1[0];

        double area = 0.5 * sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

		return area;
    }

    //*************************************************************************************************************************************
    //* - Calculation of the stiffness matrix
    //*************************************************************************************************************************************	
    void GenericRobin_3sb_ColdPlasma::GetStiffnessMatrix( Matrix<std::complex<double> >& StiffMatrix ) 
    {
        // Loop indexes
        int i, j, gp;

        // Initializing elemental stiff matrix
        StiffMatrix.Resize( numDofs, numDofs );

        // Gaussian weights and points
		Vector<double> cX; 
		Vector<double> cY;
		Vector<double>  W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

        // Element base in Gauss points
		Matrix<double> N;

		Lagrange2D_Ni_3sb( N, cX, cY );

        // Permittivity tensor on Gauss points (interpolated with Ni)
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, N );

        // Jacobian 
		double detJ = 2.00 * Calculate_Area();
        
        // Plasma sheath coefficients and integral constants
        Vector<double> vPSC = (*mProperties)(COMPLEX_IBC);

        std::complex<double> cCoef( -detJ*vPSC[0], -detJ*vPSC[1] );

        //////////////////////////////////////// Block [Kxx] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=i; j<numNodes; j++ )
            {
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 0; 
                int cj = 0;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Ni_TEp_Nj += W[gp] * ( N[i][gp] * TEgp[gp][ci][cj] * N[j][gp] );
                }

                StiffMatrix[i][j] = cCoef * intg_Ni_TEp_Nj ;
            }
        }
        //////////////////////////////////////// Block [Kyy] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=i; j<numNodes; j++ )
            {
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 1; 
                int cj = 1;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Ni_TEp_Nj += W[gp] * ( N[i][gp] * TEgp[gp][ci][cj] * N[j][gp] );
                }

                StiffMatrix[i+numNodes][j+numNodes] = cCoef * intg_Ni_TEp_Nj ;
            }
        }
        //////////////////////////////////////// Block [Kzz] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=i; j<numNodes; j++ )
            {
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 2; 
                int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Ni_TEp_Nj += W[gp] * ( N[i][gp] * TEgp[gp][ci][cj] * N[j][gp] );
                }

                StiffMatrix[i+(2*numNodes)][j+(2*numNodes)] = cCoef * intg_Ni_TEp_Nj ;
            }
        }
        //////////////////////////////////////// Block [Kxy] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 0; 
                int cj = 1;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Ni_TEp_Nj += W[gp] * ( N[i][gp] * TEgp[gp][ci][cj] * N[j][gp] );
                }

                StiffMatrix[i][j+numNodes] = cCoef * intg_Ni_TEp_Nj ;
            }
        }
        /////////////////////////////////////// Block [Kxz] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 0; 
                int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Ni_TEp_Nj += W[gp] * ( N[i][gp] * TEgp[gp][ci][cj] * N[j][gp] );
                }

                StiffMatrix[i][j+(2*numNodes)] = cCoef * intg_Ni_TEp_Nj ;
            }
        }
        ////////////////////////////////////// Block [Kyz] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 );   
                
                int ci = 1; 
                int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Ni_TEp_Nj += W[gp] * ( N[i][gp] * TEgp[gp][ci][cj] * N[j][gp] );
                }

                StiffMatrix[i+numNodes][j+(2*numNodes)] = cCoef * intg_Ni_TEp_Nj ;
            }
        }

        // Lower diagonal
        for ( i=0; i<numDofs; i++ )
		{
            for ( j=i+1; j<numDofs; j++ )
            {
                StiffMatrix[j][i] = std::conj( StiffMatrix[i][j] );
            }
		}
    }
} 