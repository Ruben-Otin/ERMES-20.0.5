
#include "CoaxialPortTEM_2ndQ.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
   /**********************************************************************************
    **********************************************************************************
	**
    ** FREQUENCY      : w
    **
	** COMPLEX_IBC[0] : tipo de elemento:
	**                   - 0 : puerto de entrada TEM
    **                   - 1 : puerto de salida 1
	**				     - 2 : puerto de salida 2
    **
	** COMPLEX_IBC[1]     : Coord. x centro del coaxial.
	** COMPLEX_IBC_2o [0] : Coord. y centro del coaxial.
	** COMPLEX_IBC_2o [1] : Coord. z centro del coaxial.
    **
    ** COMPLEX_NEUMANN_FLOW[0] : Inner radius a.
	** COMPLEX_NEUMANN_FLOW[1] : Exterior radius b.
	** COMPLEX_NEUMANN_FLOW[2] : Electrical permittivity.
    **
	** COMPLEX_NEUMANN_FLOW[3] : Magnetic permeability.
	** COMPLEX_NEUMANN_FLOW[4] : 0.00
	** COMPLEX_NEUMANN_FLOW[5] : 0.00
    **
    **********************************************************************************
	**********************************************************************************/

   /**********************************************************************************
    * - Global index of the nodes
    **********************************************************************************/	
    void CoaxialPortTEM_2ndQ::GetEquationIdVector(Vector<int>& EquationId)
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

   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void CoaxialPortTEM_2ndQ::GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix) 
    {
		const int numNodes = 6;
		const int numDofs  = numNodes*3;

		StiffMatrix.Resize(numDofs,numDofs);

        int i,j;
		
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double>  W;

		int gp;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		std::vector<double> detJ;

        Calculate_detJ(detJ, cX, cY);

		for (gp=0; gp<nGaussPoints; gp++) W[gp] *= detJ[gp]; 

		Matrix<double> N;

		Lagrange2D_Ni_2nd(N, cX, cY);

		double intNiNj;

		std::complex<double> ccteNN = - CalculateKMuCoef();

		std::complex<double> ccero(0.00, 0.00);

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

   /***********************************************************************************
    * - Calculation of the Jacobian determinant
    ***********************************************************************************/
    void CoaxialPortTEM_2ndQ::Calculate_detJ(std::vector<double>& detJ,
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
    * - Returns current port number
    **********************************************************************************/
	double CoaxialPortTEM_2ndQ::getPortNumber()
	{
		Vector<double> cIBC = (*mProperties)(COMPLEX_IBC);

		return cIBC[0];
	}

   /**********************************************************************************
    * - Set up j*(k/mu) coef
    **********************************************************************************/
    std::complex<double> CoaxialPortTEM_2ndQ::CalculateKMuCoef() 
	{
		//Freq = w = 2*pi*v
		double freq = (*mProperties)(FREQUENCY);
		
		//Relative electrical permittivity and relative magnetic permeability
		Vector<double> cNF = (*mProperties)(COMPLEX_NEUMANN_FLOW); 

		double rep = cNF[2];
		double rmu = cNF[3];

		double pi = 3.1415926535897932384626433832795;

		double eo  = 8.8541878176e-12;
		double muo = pi*4.00e-7;

		double ep = rep * eo;
		double mu = rmu * muo;

		std::complex<double> cKMu( 0.00, freq*sqrt(ep/mu) );
	
		return cKMu;
	}

   /**********************************************************************************
    * - Set up normalization coef
    **********************************************************************************/
    std::complex<double> CoaxialPortTEM_2ndQ::CalculateNormalizationCoef() 
	{
		double freq = (*mProperties)(FREQUENCY);
		
		//a,b,relative electrical permittivity and relative magnetic permeability
		Vector<double> cNF = (*mProperties)(COMPLEX_NEUMANN_FLOW);
		
		double a = cNF[0];
		double b = cNF[1];

		double rep = cNF[2];
		double rmu = cNF[3];

		double pi = 3.1415926535897932384626433832795;

		double eo  = 8.8541878176e-12;
		double muo = pi*4.00e-7;

		double ep = rep * eo;
		double mu = rmu * muo;

		double ctcSup = sqrt(mu/ep);
		double ctcInf = 2*pi*log(b/a);

		double realNorm = sqrt(ctcSup/ctcInf);

		std::complex<double> cNormCoef(realNorm, 0.00);
	
		return cNormCoef;
	}

   /**********************************************************************************
    * - Calculation of the residual vector (Linear interpolation)
    **********************************************************************************/
	void CoaxialPortTEM_2ndQ::Calculate_r0(double* r0) 
	{
		Vector<double> cIBC   = (*mProperties)(COMPLEX_IBC),
					   cIBC2o = (*mProperties)(COMPLEX_IBC_2o); 
		
		r0[0] = cIBC  [1];
		r0[1] = cIBC2o[0];
	    r0[2] = cIBC2o[1];
	}

   /**********************************************************************************
    * - Calculation of the residual vector 
    **********************************************************************************/
    void CoaxialPortTEM_2ndQ::GetResidualVector(Vector<std::complex<double> >& ResidualVector) 
    {
        int i, gp;
		const int numNodes = 6;
		const int numDofs  = numNodes*3;
		
		ResidualVector.resize(numDofs);
		
		Vector<double> puerto = (*mProperties)(COMPLEX_IBC);
		
		//Output Port
		if (puerto[0] != 0.00)
		{
			std::complex<double> czero(0.00,0.00);

			for (int d=0; d<numDofs; d++) ResidualVector[d] = czero;

			return;
		}

		//Excitation Port
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		std::vector<double> detJ;

        Calculate_detJ(detJ, cX, cY);

		for (gp=0; gp<nGaussPoints; gp++) W[gp] *= detJ[gp]; 

		Matrix<double> N;

		Lagrange2D_Ni_2nd(N, cX, cY);

		std::vector<double> Ex(nGaussPoints),
			                Ey(nGaussPoints),
							Ez(nGaussPoints);

		double r0[3];
		Calculate_r0(r0);

		double rho[3];
		double invRho2;

		for ( gp=0; gp<nGaussPoints; gp++ ) 
        {
			double rX = 0.00;	
			double rY = 0.00;	
			double rZ = 0.00;	

			for ( i=0; i<numNodes; i++ ) 
			{
				rX += N[i][gp] * mNodes[i]->X();  
				rY += N[i][gp] * mNodes[i]->Y();  
				rZ += N[i][gp] * mNodes[i]->Z();  
			}

			rho[0] = rX - r0[0];
			rho[1] = rY - r0[1];
			rho[2] = rZ - r0[2];

			invRho2 = 1.00/(rho[0]*rho[0] + rho[1]*rho[1] + rho[2]*rho[2]);

			Ex[gp] = invRho2 * rho[0]; 
			Ey[gp] = invRho2 * rho[1];  
			Ez[gp] = invRho2 * rho[2]; 
        }

		std::complex<double> cKMu      = CalculateKMuCoef();
		std::complex<double> cNormCoef = CalculateNormalizationCoef();

        std::complex<double> cCte = -2.00 * ( cKMu * cNormCoef );
	
		for ( i=0; i<numNodes; i++ )
        {
            double intExNi = 0.00;
			double intEyNi = 0.00;
			double intEzNi = 0.00;

            for ( gp=0; gp<nGaussPoints; gp++ )
			{
				intExNi += W[gp]*Ex[gp]*N[i][gp];
				intEyNi += W[gp]*Ey[gp]*N[i][gp];
				intEzNi += W[gp]*Ez[gp]*N[i][gp];
			}

            ResidualVector[i             ] = cCte*intExNi;
            ResidualVector[i+ numNodes   ] = cCte*intEyNi;
			ResidualVector[i+(numNodes*2)] = cCte*intEzNi;
        } 
    }

   /**********************************************************************************
    * - Returns TEM projection on a plane surface 
    **********************************************************************************/
	std::complex<double> CoaxialPortTEM_2ndQ::getProjection()
	{
		int i, gp;
		const int numNodes = 6;

        std::vector<double> cX; 
		std::vector<double>	cY; 
		std::vector<double>  W;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		std::vector<double> detJ;

        Calculate_detJ(detJ, cX, cY);

		for (gp=0; gp<nGaussPoints; gp++) W[gp] *= detJ[gp]; 

		Matrix<double> N;

		Lagrange2D_Ni_2nd(N, cX, cY);

		std::vector<std::complex<double> > cEr(nGaussPoints);

		std::vector<double> invRho(nGaussPoints);

		std::complex<double> cE_x, cE_y, cE_z;

		double r0[3];   
		Calculate_r0(r0); 

		double r[3];

		double invModule;

		//Er and invRho in gauss points
		for ( gp=0; gp<nGaussPoints; gp++ ) 
        {
			std::complex<double> cExInGp(0.00,0.00);
			std::complex<double> cEyInGp(0.00,0.00);
			std::complex<double> cEzInGp(0.00,0.00);

			for ( i=0; i<numNodes; i++ ) 
			{
				cE_x = (*mProperties)(cEx, *mNodes[i]);
				cE_y = (*mProperties)(cEy, *mNodes[i]);
				cE_z = (*mProperties)(cEz, *mNodes[i]);

			    cExInGp += N[i][gp] * cE_x;
				cEyInGp += N[i][gp] * cE_y;
				cEzInGp += N[i][gp] * cE_z;
			}

			double rX = 0.00;	
			double rY = 0.00;	
			double rZ = 0.00;	

			for ( i=0; i<numNodes; i++ ) 
			{
				rX += N[i][gp] * mNodes[i]->X();  
				rY += N[i][gp] * mNodes[i]->Y();  
				rZ += N[i][gp] * mNodes[i]->Z();  
			}

			r[0] = rX - r0[0];
			r[1] = rY - r0[1];
			r[2] = rZ - r0[2];

			invModule  = 1.00 / sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

			invRho[gp] = invModule;

			r[0] *= invModule;
			r[1] *= invModule;
			r[2] *= invModule;

			cEr[gp] = cExInGp * r[0] + cEyInGp * r[1] + cEzInGp * r[2];
        }

		std::complex<double> TEMProject(0.00,0.00);

		//Integral ( Er * (1.00 / rho) ) over element surface
		for ( gp=0; gp<nGaussPoints; gp++ ) 
        {
			TEMProject += ( W[gp] * invRho[gp] ) * cEr[gp];  
        }

		double portNumber = getPortNumber();

        //Cuando los puertos estan dentro del volumen los elementos se repiten
		if ( int(portNumber) != portNumber ) TEMProject *= 0.5;

		return TEMProject;
	}

   /**********************************************************************************
    * - Returns a normalization term (useful when using PMC symmetries)
    **********************************************************************************/
	std::complex<double> CoaxialPortTEM_2ndQ::getNormalization()
	{
		int i, gp;
		const int numNodes = 6;

		std::vector<double> cX; 
		std::vector<double>	cY; 
		std::vector<double>  W;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		std::vector<double> detJ;
        Calculate_detJ(detJ, cX, cY);

		for (gp=0; gp<nGaussPoints; gp++) W[gp] *= detJ[gp]; 

		Matrix<double> N;
		Lagrange2D_Ni_2nd(N, cX, cY);

		double r0[3];
		Calculate_r0(r0);

		double rho[3];

		double intInvRho2 = 0.00;

		//Integral ( 1.00 / rho^2 ), para normalizar en caso de simetrías
		for ( gp=0; gp<nGaussPoints; gp++ ) 
        {
			double rX = 0.00;	
			double rY = 0.00;	
			double rZ = 0.00;	

			for ( i=0; i<numNodes; i++ ) 
			{
				rX += N[i][gp] * mNodes[i]->X();  
				rY += N[i][gp] * mNodes[i]->Y();  
				rZ += N[i][gp] * mNodes[i]->Z();  
			}

			rho[0] = rX - r0[0];
			rho[1] = rY - r0[1];
			rho[2] = rZ - r0[2];

			intInvRho2 += W[gp] / (rho[0]*rho[0] + rho[1]*rho[1] + rho[2]*rho[2]);
        }

		double portNumber = getPortNumber();

        //Cuando los puertos estan dentro del volumen los elementos se repiten
		if ( int(portNumber) != portNumber ) intInvRho2 *= 0.5;
		
		std::complex<double> Normalization(intInvRho2, 0.00);

		std::complex<double> cNormCoef = CalculateNormalizationCoef();

		Normalization *= cNormCoef;

		return Normalization;
	}
} 