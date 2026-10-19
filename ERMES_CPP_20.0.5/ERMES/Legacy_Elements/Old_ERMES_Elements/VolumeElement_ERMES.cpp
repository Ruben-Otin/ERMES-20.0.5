

#include "VolumeElement_ERMES.h"
#include "../includes/variables.h"


namespace Kratos
{   
   /**********************************************************************************
    * -Global index of the nodes
    **********************************************************************************/	
    void VolumeElement_ERMES::GetEquationIdVector(int* EqIdVector_12)
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
    void VolumeElement_ERMES::GetMagnitudeVector(std::complex<double>* MagnitudVector) 
    {  
        for (int i=0; i<4; i++) 
        {
            MagnitudVector[i  ] = (*mProperties)(cEx, *mNodes[i]);
            MagnitudVector[i+4] = (*mProperties)(cEy, *mNodes[i]);
            MagnitudVector[i+8] = (*mProperties)(cEz, *mNodes[i]);
        }
    }

   /**********************************************************************************
    * - Calculation of the residual vector
    **********************************************************************************/
    void VolumeElement_ERMES::GetResidualVector(std::complex<double>* ResidualVector_12) 
    {
        std::complex<double> zero(0.00,0.00);

        for (int i=0; i<12; i++) 
        {
            (*ResidualVector_12) = zero;
             ResidualVector_12++;
        }
    }

   /**********************************************************
    * - Calculate element volume
    **********************************************************/
    double VolumeElement_ERMES::Calculate_Volume()
    {
        double det;

        det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs(det/6.00);
    }

   /**********************************************************
    * - Calculation of dN/dz column vector
    **********************************************************/
    void VolumeElement_ERMES::Calculate_DN(double DN[3][4])
    {
        double elementVolume = Calculate_Volume();

        double IntgCte = 1.00/(6.00*elementVolume); 

        //dN/dx
        DN[0][0] = IntgCte*( Y(3)*Z(2) + Y(4)*Z(3) + Y(2)*Z(4) - Y(3)*Z(4) - Y(2)*Z(3) - Y(4)*Z(2) );
        DN[0][1] = IntgCte*( Y(3)*Z(4) + Y(4)*Z(1) + Y(1)*Z(3) - Y(3)*Z(1) - Y(1)*Z(4) - Y(4)*Z(3) );
        DN[0][2] = IntgCte*( Y(2)*Z(1) + Y(1)*Z(4) + Y(4)*Z(2) - Y(2)*Z(4) - Y(4)*Z(1) - Y(1)*Z(2) );
        DN[0][3] = IntgCte*( Y(2)*Z(3) + Y(3)*Z(1) + Y(1)*Z(2) - Y(2)*Z(1) - Y(3)*Z(2) - Y(1)*Z(3) );
        
        //dN/dy
        DN[1][0] = IntgCte*( X(3)*Z(4) + X(4)*Z(2) + X(2)*Z(3) - X(3)*Z(2) - X(2)*Z(4) - X(4)*Z(3) ); 
        DN[1][1] = IntgCte*( X(3)*Z(1) + X(1)*Z(4) + X(4)*Z(3) - X(1)*Z(3) - X(3)*Z(4) - X(4)*Z(1) );
        DN[1][2] = IntgCte*( X(2)*Z(4) + X(4)*Z(1) + X(1)*Z(2) - X(2)*Z(1) - X(4)*Z(2) - X(1)*Z(4) );
        DN[1][3] = IntgCte*( X(1)*Z(3) + X(3)*Z(2) + X(2)*Z(1) - X(3)*Z(1) - X(1)*Z(2) - X(2)*Z(3) );
        
        //dN/dz
        DN[2][0] = IntgCte*( X(4)*Y(3) + X(2)*Y(4) + X(3)*Y(2) - X(2)*Y(3) - X(4)*Y(2) - X(3)*Y(4) );
        DN[2][1] = IntgCte*( X(3)*Y(4) + X(4)*Y(1) + X(1)*Y(3) - X(3)*Y(1) - X(4)*Y(3) - X(1)*Y(4) ); 
        DN[2][2] = IntgCte*( X(2)*Y(1) + X(4)*Y(2) + X(1)*Y(4) - X(4)*Y(1) - X(2)*Y(4) - X(1)*Y(2) ); 
        DN[2][3] = IntgCte*( X(2)*Y(3) + X(3)*Y(1) + X(1)*Y(2) - X(2)*Y(1) - X(3)*Y(2) - X(1)*Y(3) );    
    }

   
    
   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void VolumeElement_ERMES::GetStiffMatrix12x12(std::complex<double> StiffMatrix_12x12[12][12])
    {
        int i,j;
        
        double freq    = (*mProperties)(FREQUENCY),
               sigma   = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY),
               epsylon = (*mProperties)(IHL_ELECTRICAL_PERMITTIVITY),
               mu      = (*mProperties)(IHL_MAGNETIC_PERMEABILITY);

        double volume  = Calculate_Volume();

        double rcte_aux = volume*freq*freq*epsylon;
        double icte_aux = volume*freq*sigma;

        std::complex<double> cCteDiv20(rcte_aux/20.00, icte_aux/20.00);
        std::complex<double> cCteDiv10(rcte_aux/10.00, icte_aux/10.00);

        std::complex<double> NiNj_matrix[4][4];

        for (i=0; i<4; i++)
        {
            for (j=0; j<4; j++)
            {
                 if (i!=j) 
                     NiNj_matrix[i][j] = cCteDiv20;
                 else      
                     NiNj_matrix[i][j] = cCteDiv10;		         
            }
        }

        double DN[3][4];
        Calculate_DN(DN);

		double matrix_ij;

        std::complex<double> cNiNj;

		double cte_aux = volume/mu;

        for(i=0; i<4; i++)
        {
            for(j=0; j<4; j++) 
            {
                /////////////  Bloques: [Kxx], [Kyy], [Kzz]  //////////////////////////
                cNiNj     = -NiNj_matrix[i][j];
                matrix_ij = cte_aux * ( DN[0][i]*DN[0][j] + DN[1][i]*DN[1][j] + DN[2][i]*DN[2][j] ) + std::real(cNiNj);
 
                StiffMatrix_12x12[i  ][j  ] = std::complex<double>(matrix_ij, std::imag(cNiNj));
                StiffMatrix_12x12[i+4][j+4] = std::complex<double>(matrix_ij, std::imag(cNiNj));
                StiffMatrix_12x12[i+8][j+8] = std::complex<double>(matrix_ij, std::imag(cNiNj));

                /////////////  Bloques: [Kxy], [Kyx]  //////////////////////////
                matrix_ij = cte_aux * ( DN[0][i]*DN[1][j] - DN[1][i]*DN[0][j] );
          
                StiffMatrix_12x12[i  ][j+4] = std::complex<double>( matrix_ij, 0.00);
                StiffMatrix_12x12[i+4][j  ] = std::complex<double>(-matrix_ij, 0.00);

                /////////////  Bloques: [Kxz] , [Kzx]  //////////////////////////
                matrix_ij = cte_aux * ( DN[0][i]*DN[2][j] - DN[2][i]*DN[0][j] );             
          
                StiffMatrix_12x12[i  ][j+8] = std::complex<double>( matrix_ij, 0.00);
                StiffMatrix_12x12[i+8][j  ] = std::complex<double>(-matrix_ij, 0.00);

                /////////////  Bloques: [Kyz], [Kzy]  //////////////////////////
                matrix_ij = cte_aux * ( DN[1][i]*DN[2][j] - DN[2][i]*DN[1][j] );
                
                StiffMatrix_12x12[i+4][j+8] = std::complex<double>( matrix_ij, 0.00);
                StiffMatrix_12x12[i+8][j+4] = std::complex<double>(-matrix_ij, 0.00);
            }
        }
    }

  

} /* end namespace Kratos */ 