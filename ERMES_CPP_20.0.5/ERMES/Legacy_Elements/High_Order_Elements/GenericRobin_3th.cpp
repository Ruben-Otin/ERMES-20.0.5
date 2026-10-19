
#include "GenericRobin_3th.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
   /**********************************************************************************
    * - Global index of the nodes
    **********************************************************************************/	
    void GenericRobin_3th::GetEquationIdVector(Vector<int>& EquationId)
    {
        const int numNodes = 10;
		const int numDofs  = numNodes*3;
		
		EquationId.resize(numDofs);

        for (int i=0; i<numNodes; i++) 
        { 
            EquationId[i             ] = mNodes[i]->pDofcEx()->EquationId();
            EquationId[i+ numNodes   ] = mNodes[i]->pDofcEy()->EquationId();
			EquationId[i+(numNodes*2)] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void GenericRobin_3th::GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix) 
    {
		const int numNodes = 10;
		const int numDofs  = numNodes*3;

		StiffMatrix.Resize(numDofs,numDofs);

        int i,j;

		////////////////////// Calculating area /////////////////////
		double n[3],v1[3],v2[3];

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

		double area = 0.5*sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double>  W;

		int gp;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		Matrix<double> N;

		Lagrange2D_Ni_3rd(N, cX, cY);

		double intNiNj;

		double jacob = 2.00*area;

		Vector<double> vGRC = (*mProperties)(COMPLEX_IBC);

		std::complex<double> ccteNN(-jacob*vGRC[0],-jacob*vGRC[1]);

		std::complex<double> ccero(0.00,0.00);

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