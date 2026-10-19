
#include "J_Source_ERMES.h"
#include "../includes/variables.h"


namespace Kratos
{   
   /**********************************************************************************
    * -Global index of the nodes
    **********************************************************************************/	
    void J_Source_ERMES::GetEquationIdVector(int* EqIdVector_12)
    {
        for (int i=0; i<4; i++) 
        { 
            EqIdVector_12[i  ] = mNodes[i]->pDofcEx()->EquationId();
            EqIdVector_12[i+4] = mNodes[i]->pDofcEy()->EquationId();
            EqIdVector_12[i+8] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

   /**********************************************************************************
	* - Get magnitude vector
	**********************************************************************************/
    void J_Source_ERMES::GetMagnitudeVector(std::complex<double>* MagnitudVector) 
    {  
        for (int i=0; i<4; i++) 
        {
            MagnitudVector[i  ] = (*mProperties)(cEx, *mNodes[i]);
            MagnitudVector[i+4] = (*mProperties)(cEy, *mNodes[i]);
            MagnitudVector[i+8] = (*mProperties)(cEz, *mNodes[i]);
        }
    }

   /**********************************************************
    * - Calculate element volume
    **********************************************************/
    double J_Source_ERMES::Calculate_Volume()
    {
        double det;

        det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs(det/6.00);
    }

  
   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void J_Source_ERMES::GetStiffMatrix12x12(std::complex<double> StiffMatrix_12x12[12][12])
    {
        std::complex<double> zero(0.00,0.00);

        for(int i=0; i<12; i++)
        {
            for(int j=0; j<12; j++) 
            {
                StiffMatrix_12x12[i][j] = zero;
            }
        }
    }


   /**********************************************************************************
    * - Calculation of the residual vector
    **********************************************************************************/
    void J_Source_ERMES::GetResidualVector(std::complex<double>* ResidualVector_12) 
    {
        double freq = (*mProperties)(FREQUENCY);

        std::vector<double> vJ = (*mProperties)(SINUSOIDAL_SURFACE_CURRENT);
          
        // j*w*J
        std::complex<double> cJx( -freq*vJ[0]*sin(vJ[1]), freq*vJ[0]*cos(vJ[1]) );
        std::complex<double> cJy( -freq*vJ[2]*sin(vJ[3]), freq*vJ[2]*cos(vJ[3]) );
        std::complex<double> cJz( -freq*vJ[4]*sin(vJ[5]), freq*vJ[4]*cos(vJ[5]) );

        int i;

        double volumeDiv4 = Calculate_Volume()/4.00;

        std::complex<double> jwJIntNi;

        jwJIntNi = cJx*volumeDiv4;
        for (i=0; i<4; i++)  ResidualVector_12[i] = jwJIntNi;

        jwJIntNi = cJy*volumeDiv4;
        for (i=4; i<8; i++)  ResidualVector_12[i] = jwJIntNi;

        jwJIntNi = cJz*volumeDiv4;
        for (i=8; i<12; i++) ResidualVector_12[i] = jwJIntNi;
    }

   
    
} /* end namespace Kratos */ 