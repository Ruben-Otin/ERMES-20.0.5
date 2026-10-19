
#include "GenericRobin_2ndQ.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
   /**********************************************************************************
    * - Global index of the nodes
    **********************************************************************************/	
    void GenericRobin_2ndQ::GetEquationIdVector(Vector<int>& EquationId)
    {
        const int numNodes = 6;
		const int numDofs  = numNodes*3;
		
		EquationId.resize(numDofs);

        for (int i=0; i<numNodes; i++) 
        { 
            EquationId[i             ] = mNodes[i]->pDofcEx()->EquationId();
            EquationId[i+ numNodes   ] = mNodes[i]->pDofcEy()->EquationId();
			EquationId[i+(numNodes*2)] = mNodes[i]->pDofcEz()->EquationId();
        }
    }
	
   /***********************************************************************************
    * - Calculation of the Jacobian determinant
    ***********************************************************************************/
    void GenericRobin_2ndQ::Calculate_detJ(std::vector<double>& detJ,
		                                   std::vector<double>& cX, 
									       std::vector<double>& cY)
    {
		int i, p;
		int numNodes  = 6;
		int numPoints = cX.size();

		std::vector<double> dNda(numNodes);
		std::vector<double> dNdb(numNodes);

		detJ.resize(numPoints);

        for (p=0; p<numPoints; p++) 
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

			for (i=0; i<numNodes; i++) 
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

   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void GenericRobin_2ndQ::GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix) 
    {
		const int numNodes = 6;
		const int numDofs  = numNodes*3;

		StiffMatrix.Resize(numDofs,numDofs);

        int i, j, gp;

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		std::vector<double> detJ;

        Calculate_detJ(detJ, cX, cY);

		for (gp=0; gp<nGaussPoints; gp++) W[gp] *= detJ[gp]; 

		Matrix<double> N;

		Lagrange2D_Ni_2nd(N, cX, cY);

		double intNiNj;

		Vector<double> vGRC = (*mProperties)(COMPLEX_IBC);

		std::complex<double> ccteNN(-vGRC[0], -vGRC[1]);
		std::complex<double> ccero ( 0.00   ,  0.00   );

		/////////////////////////// Bloque Superior ///////////////////////////////
		for (i=0; i<numNodes; i++)
        {
            for (j=i; j<numNodes; j++)
            {
                intNiNj = 0.00;

                for (gp=0;gp<nGaussPoints;gp++) intNiNj += W[gp]*N[i][gp]*N[j][gp];
                
				//[Kxx],[Kxy],[Kxz]
                StiffMatrix[i][j             ] = ccteNN * intNiNj;
				StiffMatrix[i][j+ numNodes   ] = ccero;
				StiffMatrix[i][j+(numNodes*2)] = ccero;

				//[Kyy],[Kyz]
				StiffMatrix[i+ numNodes][j+ numNodes   ] = ccteNN * intNiNj;
				StiffMatrix[i+ numNodes][j+(numNodes*2)] = ccero;

				//[Kzz]
				StiffMatrix[i+(numNodes*2)][j+(numNodes*2)] = ccteNN * intNiNj;
            }
        }

		////////////////////////// Bloque Simetrico ///////////////////////////////
		for (i=0; i<numDofs; i++)
        {
            for (j=i+1; j<numDofs; j++)
            {
                StiffMatrix[j][i] =  StiffMatrix[i][j];
            }
        }
    }
 
} /* end namespace Kratos */ 