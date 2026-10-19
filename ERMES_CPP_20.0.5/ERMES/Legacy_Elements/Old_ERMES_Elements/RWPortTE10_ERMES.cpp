
#include "RWPortTE10_ERMES.h"
#include "../includes/variables.h"


namespace Kratos
{  	
	
   /**********************************************************************************
	* - Get magnitude vector
	**********************************************************************************/
    void RWPortTE10_ERMES::GetMagnitudeVector(std::complex<double>* MagnitudVector) 
    {       
        for (int i=0; i<3; i++) 
        {
            MagnitudVector[i  ] = (*mProperties)(cEx, *mNodes[i]);
            MagnitudVector[i+3] = (*mProperties)(cEy, *mNodes[i]);
            MagnitudVector[i+6] = (*mProperties)(cEz, *mNodes[i]);
        }
    }

   /**********************************************************************************
    * -Global index of the nodes
    **********************************************************************************/	
    void RWPortTE10_ERMES::GetEquationIdVector(int* EqIdVector_9)
    {
        for (int i=0; i<3; i++) 
        {  
            EqIdVector_9[i  ] = mNodes[i]->pDofcEx()->EquationId();
            EqIdVector_9[i+3] = mNodes[i]->pDofcEy()->EquationId();
            EqIdVector_9[i+6] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

   /**********************************************************************************
    **********************************************************************************
	**
    ** FREQUENCY      : w
    **
	** COMPLEX_IBC[0] : tipo de elemento:
	**                   - 0 : puerto de entrada TE10
    **                   - 1 : puerto de salida 1
	**				     - 2 : puerto de salida 2
    **
	** COMPLEX_IBC[1]     : coord. x de la esquina origen
	** COMPLEX_IBC_2o [0] : coord. y
	** COMPLEX_IBC_2o [1] : coord. z
    **
    ** COMPLEX_NEUMANN_FLOW[0] : coord. x de la esquina superior
	** COMPLEX_NEUMANN_FLOW[1] : coord. y
	** COMPLEX_NEUMANN_FLOW[2] : coord. z
    **
	** COMPLEX_NEUMANN_FLOW[3] : coord. x de la esquina derecha
	** COMPLEX_NEUMANN_FLOW[4] : coord. y 
	** COMPLEX_NEUMANN_FLOW[5] : coord. z
    **
    **********************************************************************************
	**********************************************************************************/

   /**********************************************************************************
    * - Set up variables for the new framework
    **********************************************************************************/
    void RWPortTE10_ERMES::SetUpVariables() 
	{
		////////////////////// Calculating normal vector /////////////////////

		double v1[3],v2[3];

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

        double areaDiv12 = area/12.00;
        double areaDiv6  = area/6.00;

        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                  if (i!=j) 
                      NiNj_matrix[i][j] = areaDiv12;
                  else      
                      NiNj_matrix[i][j] = areaDiv6;
            }
        }

		////////////////// Calculating beta/mu coefficient ////////////////////

		double pi   = 3.1415926535897932384626433832795,
               freq = (*mProperties)(FREQUENCY),
			   mu   = pi*4.00e-7,
		       ko   = freq*sqrt(8.854e-12*mu);

		Vector<double> cIBC   = (*mProperties)(COMPLEX_IBC),
                       cIBC2o = (*mProperties)(COMPLEX_IBC_2o),
	                   cNF    = (*mProperties)(COMPLEX_NEUMANN_FLOW);

		double t[3];

		t[0] = cNF[3] - cIBC[1];
		t[1] = cNF[4] - cIBC2o[0];
		t[2] = cNF[5] - cIBC2o[1];

		double	kc = pi / sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
	
		if (ko < kc) 
        {
            cbetaMu = std::complex<double> (-sqrt(kc*kc - ko*ko)/mu, 0.00);
        }
		else
        {
            cbetaMu = std::complex<double> (0.00, sqrt(ko*ko - kc*kc)/mu);
		}
	}

   /**********************************************************************************
    * - Calculation of the residual vector (Linear interpolation)
    **********************************************************************************/
    void RWPortTE10_ERMES::GetResidualVector(std::complex<double>* ResidualVector_9) 
    {
	   Vector<double> puerto = (*mProperties)(COMPLEX_IBC);

       if(puerto[0] == 0.00)
	   {
			////////////////////// Calculating vector b and t /////////////////////

			int i, j;

			Vector<double> cIBC   = (*mProperties)(COMPLEX_IBC),
						   cIBC2o = (*mProperties)(COMPLEX_IBC_2o),
						   cNF    = (*mProperties)(COMPLEX_NEUMANN_FLOW);

			// New coordinate origin
			double r0[3];

			r0[0] = cIBC[1];
			r0[1] = cIBC2o[0];
			r0[2] = cIBC2o[1];

			// Vertical vector
			double b[3];
			// Horizontal vector
			double t[3];

			for (i=0; i<3; i++)
			{
				b[i] = cNF[i  ] - r0[i];
				t[i] = cNF[i+3] - r0[i];
			}

			// Waveguide height
			double h = sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
			// Waveguide width
			double a = sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );

            double  norma = 1.00/h;
            b[0] *= norma; 
            b[1] *= norma; 
            b[2] *= norma;

            norma = 1.00/a;
            t[0] *= norma; 
            t[1] *= norma; 
            t[2] *= norma;

			// b component of the electric field.
			double Eb[3];

			// Relative X coordinate of the node.
			double xt;
            // kc = (pi/a)
			double kc = 3.1415926535897932384626433832795/a; 

            for (i=0; i<3; i++)
			{ 
				xt    = (mNodes[i]->X() - r0[0])*t[0] + 
					    (mNodes[i]->Y() - r0[1])*t[1] + 
						(mNodes[i]->Z() - r0[2])*t[2]; 

                Eb[i] = sin(kc*xt);
			}

			// Suma en j de Ej*int(NiNj)
            double Ex_NiNj[3] = {0.00},
				   Ey_NiNj[3] = {0.00},
				   Ez_NiNj[3] = {0.00};

            for (i=0; i<3; i++)
			{ 
			    for (j=0; j<3; j++)
				{ 
                    double aux_ij = NiNj_matrix[i][j];

					Ex_NiNj[i] += b[0]*Eb[j]*aux_ij;
				    Ey_NiNj[i] += b[1]*Eb[j]*aux_ij;
				    Ez_NiNj[i] += b[2]*Eb[j]*aux_ij;
				}
			}


            std::complex<double> cN10(0.00,0.00);

            double rbetaMu = -std::real(cbetaMu),
                   ibetaMu =  std::imag(cbetaMu);
 
            if ( (rbetaMu == 0.00) ) 
            {
                cN10 = std::complex<double> (1.00/sqrt(ibetaMu), 0.00);
            }
            else
            {
                double dummy = 0.5*sqrt(2.00/rbetaMu);
                cN10 = std::complex<double> (dummy,-dummy);
            }

            double freq = (*mProperties)(FREQUENCY);

            cN10 *= -sqrt(2.00*freq/(a*h));

            std::complex<double> cCte = -2.00*(cbetaMu*cN10);


			for(i=0; i<3; i++) 
			{
				// {bx}
				ResidualVector_9[i  ] = cCte*Ex_NiNj[i];
				// {by}
				ResidualVector_9[i+3] = cCte*Ey_NiNj[i];
      	        // {bz}
				ResidualVector_9[i+6] = cCte*Ez_NiNj[i];
			}

        }
        else
        {
            std::complex<double> zero(0.00,0.00);

            for (int i=0; i<9; i++) 
            {
                (*ResidualVector_9) = zero;
                 ResidualVector_9++;
            }
        }

    }


   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void RWPortTE10_ERMES::GetStiffMatrix9x9(std::complex<double> StiffMatrix_9x9[9][9])
    {
		SetUpVariables();

        std::complex<double> zero(0.00,0.00);
        std::complex<double> CteXintegral_NiNj;

        //Filling stiffness matrix
        for(int i=0; i<3; i++)
        {
           for(int j=0; j<3; j++)
           {
              CteXintegral_NiNj = -cbetaMu*NiNj_matrix[i][j];

			  // [Kxx]
              StiffMatrix_9x9[i  ][j  ] = CteXintegral_NiNj;
              // [Kxy]
              StiffMatrix_9x9[i  ][j+3] = zero; 
              // [Kxz]
              StiffMatrix_9x9[i  ][j+6] = zero;

              // [Kyx]
              StiffMatrix_9x9[i+3][j  ] = zero; 
              // [Kyy]
              StiffMatrix_9x9[i+3][j+3] = CteXintegral_NiNj;
              // [Kyz]
              StiffMatrix_9x9[i+3][j+6] = zero;  

              // [Kzx]
              StiffMatrix_9x9[i+6][j  ] = zero; 
              // [Kzy]
              StiffMatrix_9x9[i+6][j+3] = zero;  
			  // [Kzz]
              StiffMatrix_9x9[i+6][j+6] = CteXintegral_NiNj;
	  
           }
        }

    }
 
} /* end namespace Kratos */ 