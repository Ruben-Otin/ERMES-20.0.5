
#include "FarField_ERMES.h"
#include "../includes/variables.h"


namespace Kratos
{  	


   /**********************************************************************************
    * - Set the problem´s frequency
    **********************************************************************************/	
    void FarField_ERMES::SetFrequency(double NewFrequency)
    {
        freq = NewFrequency;
    }

   /**********************************************************************************
    * - Global index of the nodes
    **********************************************************************************/	
    void FarField_ERMES::GetEquationIdVector(int* EqIdVector_9)
    {
        for (int i=0; i<3; i++) 
        {  
            EqIdVector_9[i  ] = mNodes[i]->pDofcEx()->EquationId();
            EqIdVector_9[i+3] = mNodes[i]->pDofcEy()->EquationId();
            EqIdVector_9[i+6] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

   /**********************************************************************************
    * - Calculation of the residual vector (Linear interpolation)
    **********************************************************************************/
    void FarField_ERMES::GetResidualVector(std::complex<double>* ResidualVector_9) 
    {
        std::complex<double> zero(0.00,0.00);

        for (int i=0; i<9; i++) 
        {
            (*ResidualVector_9) = zero;
             ResidualVector_9++;
        }
    }

   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void FarField_ERMES::GetStiffMatrix9x9(std::complex<double> StiffMatrix_9x9[9][9])
    {
		int i,j;
        
        ////////////////////// Calculating normal vector /////////////////////

		double n[3],v1[3],v2[3];

        v2[0] = mNodes[2]->X() - mNodes[0]->X();
        v2[1] = mNodes[2]->Y() - mNodes[0]->Y();
        v2[2] = mNodes[2]->Z() - mNodes[0]->Z();

        v1[0] = mNodes[1]->X() - mNodes[0]->X();
        v1[1] = mNodes[1]->Y() - mNodes[0]->Y();
        v1[2] = mNodes[1]->Z() - mNodes[0]->Z();

        //normal = (v2 x v1)/mod_normal
        n[0] = v2[1]*v1[2] - v2[2]*v1[1];
        n[1] = v2[2]*v1[0] - v2[0]*v1[2];
        n[2] = v2[0]*v1[1] - v2[1]*v1[0];

		double area = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

        double  norma = 1.00/area;
        n[0] *= norma; 
        n[1] *= norma; 
        n[2] *= norma;

		area = 0.5*area;

        ////////////////// Calculating NiNj matrix ////////////////////

        double pi = 3.1415926535897932384626433832795,
			   mu = pi*4.00e-7,
               eo = 8.854e-12;

		double koXareaDivmu = -area*freq*sqrt(eo/mu);

        std::complex<double> cCteDiv12(0.00, koXareaDivmu/12.00),
                             cCteDiv6 (0.00, koXareaDivmu/6.00);

        std::complex<double> NiNj_matrix[3][3];

        for (i=0; i<3; i++)
        {
            for (j=0; j<3; j++)
            {
                  if (i!=j) 
                      NiNj_matrix[i][j] = cCteDiv12;
                  else      
                      NiNj_matrix[i][j] = cCteDiv6;
            }
        }

		////////////// Calculating local stiffnes matrix //////////////

        std::complex<double> integral_NiNj;
        std::complex<double> zero(0.00,0.00);

        //Filling stiffness matrix
        for(i=0; i<3; i++)
        {
           for(j=0; j<3; j++)
           {
              integral_NiNj = NiNj_matrix[i][j];

			  // [Kxx]
              StiffMatrix_9x9[i  ][j  ] = integral_NiNj;
              // [Kxy]
              StiffMatrix_9x9[i  ][j+3] = zero; 
              // [Kxz]
              StiffMatrix_9x9[i  ][j+6] = zero;

              // [Kyx]
              StiffMatrix_9x9[i+3][j  ] = zero; 
              // [Kyy]
              StiffMatrix_9x9[i+3][j+3] = integral_NiNj;
              // [Kyz]
              StiffMatrix_9x9[i+3][j+6] = zero;  

              // [Kzx]
              StiffMatrix_9x9[i+6][j  ] = zero; 
              // [Kzy]
              StiffMatrix_9x9[i+6][j+3] = zero;  
			  // [Kzz]
              StiffMatrix_9x9[i+6][j+6] = integral_NiNj;
	  
           }
        }

    }
 
} /* end namespace Kratos */ 