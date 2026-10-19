
#include "JSource_3th.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
   /**********************************************************************************
    * -Global index of the nodes
    **********************************************************************************/	
    void JSource_3th::GetEquationIdVector(Vector<int>& EquationId)
    {
        const int numNodes = 20;
		const int numDofs  = numNodes*3;
		
		EquationId.resize(numDofs);

        for (int i=0; i<numNodes; i++) 
        { 
            EquationId[i              ] = mNodes[i]->pDofcEx()->EquationId();
            EquationId[i + numNodes   ] = mNodes[i]->pDofcEy()->EquationId();
			EquationId[i +(numNodes*2)] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

   /**********************************************************
    * - Calculate element volume
    **********************************************************/
    double JSource_3th::Calculate_Volume()
    {
         double det;

         det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
               - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
               + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
               - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
         return fabs(det/6.00);
    }

   /**********************************************************************************
    * - Calculation of the residual vector
    **********************************************************************************/
    void JSource_3th::GetResidualVector(Vector<std::complex<double> >& ResidualVector) 
    {
        const int numNodes = 20;
		const int numDofs  = numNodes*3;
		
		int i,j;
		
		ResidualVector.resize(numDofs);

		for (i=0; i<numDofs; i++) ResidualVector[i] = std::complex<double>(0.00,0.00);

        double freq = (*mProperties)(FREQUENCY);

        std::vector<double> vJ  = (*mProperties)(SINUSOIDAL_SURFACE_CURRENT);
		std::vector<double> vJa = (*mProperties)(COMPLEX_IBC);

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int gp;

		int nGaussPoints = GaussPoints3D_Order5(cX, cY, cZ, W);

		Matrix<double> N;

		Lagrange3D_Ni_3rd(N, cX, cY, cZ);

        double jacob = 6.00*Calculate_Volume();

		if ((vJ[0] != 0.00) || (vJ[2] != 0.00) || (vJ[4] != 0.00))
		{
			// j*w*J
			std::complex<double> cJx(-freq * vJ[0] * sin(vJ[1]), freq * vJ[0] * cos(vJ[1]));
			std::complex<double> cJy(-freq * vJ[2] * sin(vJ[3]), freq * vJ[2] * cos(vJ[3]));
			std::complex<double> cJz(-freq * vJ[4] * sin(vJ[5]), freq * vJ[4] * cos(vJ[5]));

			double intNi;

			for (i = 0; i < numNodes; i++)
			{
				intNi = 0.00;
				for (gp = 0; gp < nGaussPoints; gp++) intNi += W[gp] * N[i][gp];

				ResidualVector[i                 ] += cJx * (jacob * intNi);
				ResidualVector[i +  numNodes     ] += cJy * (jacob * intNi);
				ResidualVector[i + (numNodes * 2)] += cJz * (jacob * intNi);
			}
		}

		if (vJa[0] != 0.00)
		{
			// j*w*Ja
			std::complex<double> cJa(-freq * vJa[0] * sin(vJa[1]), freq * vJa[0] * cos(vJa[1]));

			std::vector<std::complex<double> > cJxc(numNodes);
			std::vector<std::complex<double> > cJyc(numNodes);
			std::vector<std::complex<double> > cJzc(numNodes);

			std::vector<double> Jaxis = (*mProperties)(COMPLEX_NEUMANN_FLOW);

			for (i = 0; i < numNodes; i++)
			{
				// Position vector 
				double px = mNodes[i]->X() - Jaxis[0];
				double py = mNodes[i]->Y() - Jaxis[1];
				double pz = mNodes[i]->Z() - Jaxis[2];

				// Rotation axis
				double ax = Jaxis[3];
				double ay = Jaxis[4];
				double az = Jaxis[5];

				// r = a x p 
				double rx = ay * pz - az * py;
				double ry = az * px - ax * pz;
				double rz = ax * py - ay * px;

				double norm = sqrt(rx * rx + ry * ry + rz * rz);

				if (norm > 0.00)
				{
					rx = rx / norm;
					ry = ry / norm;
					rz = rz / norm;
				}
				else
				{
					rx = 0.00;
					ry = 0.00;
					rz = 0.00;
				}

				cJxc[i] = cJa * rx;
				cJyc[i] = cJa * ry;
				cJzc[i] = cJa * rz;
			}

			double intNiNj;

			for (i = 0; i < numNodes; i++)
			{
				for (j = 0; j < numNodes; j++)
				{
					intNiNj = 0.00;

					for (gp = 0; gp < nGaussPoints; gp++) intNiNj += W[gp] * N[i][gp] * N[j][gp];

					ResidualVector[i                 ] += cJxc[j] * (jacob * intNiNj);
					ResidualVector[i +  numNodes     ] += cJyc[j] * (jacob * intNiNj);
					ResidualVector[i + (numNodes * 2)] += cJzc[j] * (jacob * intNiNj);
				}

			}//end for (i=0; i<numNodes; i++)

		}// end if ( vJa[0] != 0.00 )

    }// end function
   
} 