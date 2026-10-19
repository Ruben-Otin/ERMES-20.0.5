
#include "JSource_2ndQ.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
   /**********************************************************************************
    * -Global index of the nodes
    **********************************************************************************/	
    void JSource_2ndQ::GetEquationIdVector(Vector<int>& EquationId)
    {
        const int numNodes = 10;
		const int numDofs  = numNodes*3;
		
		EquationId.resize(numDofs);

        for (int i=0; i<numNodes; i++) 
        { 
            EquationId[i              ] = mNodes[i]->pDofcEx()->EquationId();
            EquationId[i + numNodes   ] = mNodes[i]->pDofcEy()->EquationId();
			EquationId[i +(numNodes*2)] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

   /***********************************************************************************
    * - Calculation of dN/da, dN/db, dN/du on cX,cY,cZ points
    ***********************************************************************************/
    void JSource_2ndQ::NaturalDerivatives3D_2ndQ(Matrix<double>& dNda,
		                                         Matrix<double>& dNdb,
												 Matrix<double>& dNdu,
												 std::vector<double>& cX, 
												 std::vector<double>& cY, 
												 std::vector<double>& cZ)
    {
		int p;
		int numNodes  = 10;
		int numPoints = cX.size();

		dNda.Resize(numNodes, numPoints);
		dNdb.Resize(numNodes, numPoints);
		dNdu.Resize(numNodes, numPoints);

		double L1,L2,L3,L4;

		//dNi/da ,dNi/db, dNi/du
        for ( p=0; p<numPoints; p++ ) 
        {
            L2 = cX[p];
            L3 = cY[p];
			L4 = cZ[p];
            L1 = 1.00 - L2 - L3 - L4;

			dNda[0][p] =  1.00 - 4.00 * L1;
			dNdb[0][p] =  1.00 - 4.00 * L1;
			dNdu[0][p] =  1.00 - 4.00 * L1;

			dNda[1][p] =  4.00 * L2 - 1.00;
			dNdb[1][p] =  0.00;
			dNdu[1][p] =  0.00;

			dNda[2][p] =  0.00;
			dNdb[2][p] =  4.00 * L3 - 1.00;
			dNdu[2][p] =  0.00;

			dNda[3][p] =  0.00;
			dNdb[3][p] =  0.00;
			dNdu[3][p] =  4.00 * L4 - 1.00;

			dNda[4][p] = -4.00 * L2 + 4.00 * L1;
			dNdb[4][p] = -4.00 * L2;
			dNdu[4][p] = -4.00 * L2;

			dNda[5][p] =  4.00 * L3;
			dNdb[5][p] =  4.00 * L2;
			dNdu[5][p] =  0.00;

			dNda[6][p] = -4.00 * L3 ;
			dNdb[6][p] = -4.00 * L3 + 4.00 * L1;
			dNdu[6][p] = -4.00 * L3 ;

			dNda[7][p] = -4.00 * L4;
			dNdb[7][p] = -4.00 * L4;
			dNdu[7][p] = -4.00 * L4 + 4.00 * L1;

			dNda[8][p] =  4.00 * L4;
			dNdb[8][p] =  0.00;
			dNdu[8][p] =  4.00 * L2;

			dNda[9][p] =  0.00;
			dNdb[9][p] =  4.00 * L4;
			dNdu[9][p] =  4.00 * L3;
		}
	}

   /***********************************************************************************
    * - Calculation of the Jacobian determinant
    ***********************************************************************************/
    void JSource_2ndQ::Calculate_detJ(std::vector<double>& detJ,
		                              std::vector<double>& cX, 
									  std::vector<double>& cY, 
									  std::vector<double>& cZ)
    {
		int i, p;
		int numNodes  = 10;
		int numPoints = cX.size();

		detJ.resize(numPoints, 0.00);

		Matrix<double> dNda, dNdb, dNdu;

		NaturalDerivatives3D_2ndQ(dNda, dNdb, dNdu, cX, cY, cZ);

		for ( p=0; p<numPoints; p++ ) 
        {
			double dXda = 0.00;	double dXdb = 0.00; double dXdu = 0.00;
			double dYda = 0.00;	double dYdb = 0.00; double dYdu = 0.00;
			double dZda = 0.00;	double dZdb = 0.00; double dZdu = 0.00;
			
			for ( i=0; i<numNodes; i++ ) 
			{
				dXda += dNda[i][p] * mNodes[i]->X();  dXdb += dNdb[i][p] * mNodes[i]->X();  dXdu += dNdu[i][p] * mNodes[i]->X();
				dYda += dNda[i][p] * mNodes[i]->Y();  dYdb += dNdb[i][p] * mNodes[i]->Y();  dYdu += dNdu[i][p] * mNodes[i]->Y();
				dZda += dNda[i][p] * mNodes[i]->Z();  dZdb += dNdb[i][p] * mNodes[i]->Z();  dZdu += dNdu[i][p] * mNodes[i]->Z();
			}
	
			detJ[p] = dXda*dYdb*dZdu + dXdb*dZda*dYdu + dXdu*dYda*dZdb - dXda*dZdb*dYdu - dXdb*dYda*dZdu - dXdu*dZda*dYdb;
		}
	}

   /**********************************************************************************
    * - Calculation of the residual vector
    **********************************************************************************/
    void JSource_2ndQ::GetResidualVector(Vector<std::complex<double> >& ResidualVector) 
    {
        int i, j, gp;
		
		const int numNodes = 10;
		const int numDofs  = numNodes*3;
		
		ResidualVector.resize(numDofs);

		for ( i=0; i<numDofs; i++ ) ResidualVector[i] = std::complex<double>(0.00,0.00);

        double freq = (*mProperties)(FREQUENCY);

        std::vector<double> vJ  = (*mProperties)(SINUSOIDAL_SURFACE_CURRENT);
		std::vector<double> vJa = (*mProperties)(COMPLEX_IBC);

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order4(cX, cY, cZ, W);

		std::vector<double> detJ;

        Calculate_detJ(detJ, cX, cY, cZ);

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ[gp]; 

		Matrix<double> N;

		Lagrange3D_Ni_2nd(N, cX, cY, cZ);

		if ((vJ[0] != 0.00) || (vJ[2] != 0.00) || (vJ[4] != 0.00))
		{
			// j*w*J
			std::complex<double> cJx(-freq * vJ[0] * sin(vJ[1]), freq * vJ[0] * cos(vJ[1]));
			std::complex<double> cJy(-freq * vJ[2] * sin(vJ[3]), freq * vJ[2] * cos(vJ[3]));
			std::complex<double> cJz(-freq * vJ[4] * sin(vJ[5]), freq * vJ[4] * cos(vJ[5]));

			double intNi;

			for (i=0; i<numNodes; i++)
			{
				intNi = 0.00;

				for (gp=0; gp<nGaussPoints; gp++) intNi += W[gp] * N[i][gp];

				ResidualVector[i             ] += cJx * intNi;
				ResidualVector[i+ numNodes   ] += cJy * intNi;
				ResidualVector[i+(numNodes*2)] += cJz * intNi;
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

					ResidualVector[i                 ] += cJxc[j] * intNiNj;
					ResidualVector[i +  numNodes     ] += cJyc[j] * intNiNj;
					ResidualVector[i + (numNodes * 2)] += cJzc[j] * intNiNj;
				}

			}//end for (i=0; i<numNodes; i++)

		}// end if ( vJa[0] != 0.00 )

    }// end function
    
} 