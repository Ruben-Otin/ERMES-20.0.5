
#include "RWPortTE10_3th.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
   /**********************************************************************************
    * - Global index of the nodes
    **********************************************************************************/	
    void RWPortTE10_3th::GetEquationIdVector(Vector<int>& EquationId)
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
    * - Set up beta/mu coef
    **********************************************************************************/
    std::complex<double> RWPortTE10_3th::CalculateBetaMuCoef() 
	{
		double freq = (*mProperties)(FREQUENCY);
		double pi   = 3.1415926535897932384626433832795;
		double mu   = pi * (4.00e-7);
		double a    = CalculateWidth();

		double ko = freq * sqrt( mu * (8.8541878176e-12) );
		double kc = pi/a;

		std::complex<double> cbetaMu;
	
		if (ko < kc) 
        {
            cbetaMu = std::complex<double> (-sqrt(kc*kc - ko*ko)/mu, 0.00);
        }
		else
        {
            cbetaMu = std::complex<double> (0.00, sqrt(ko*ko - kc*kc)/mu);
		}

		return cbetaMu;
	}

   /**********************************************************************************
    * - Calculates waveguide width
    **********************************************************************************/
	double RWPortTE10_3th::CalculateWidth()
	{
		Vector<double> cIBC   = (*mProperties)(COMPLEX_IBC),
					   cIBC2o = (*mProperties)(COMPLEX_IBC_2o),
					   cNF    = (*mProperties)(COMPLEX_NEUMANN_FLOW);

		// Horizontal vector
		double t[3];

		t[0] = cNF[3] - cIBC  [1];
		t[1] = cNF[4] - cIBC2o[0];
		t[2] = cNF[5] - cIBC2o[1];

		// Waveguide width
		return sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
	}

   /**********************************************************************************
    * - Calculates waveguide height
    **********************************************************************************/
	double RWPortTE10_3th::CalculateHeight()
	{
		Vector<double> cIBC   = (*mProperties)(COMPLEX_IBC),
					   cIBC2o = (*mProperties)(COMPLEX_IBC_2o),
					   cNF    = (*mProperties)(COMPLEX_NEUMANN_FLOW);

		//Vertical vector
		double b[3];
		
		b[0] = cNF[0] - cIBC  [1];
		b[1] = cNF[1] - cIBC2o[0];
		b[2] = cNF[2] - cIBC2o[1];

		// Waveguide height
		return sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
	}
		
   /**********************************************************************************
    * - Calculates N10 constant
    **********************************************************************************/
    std::complex<double> RWPortTE10_3th::CalculateN10() 
	{
		std::complex<double> cN10;

		std::complex<double> cbetaMu = CalculateBetaMuCoef();

        double rbetaMu = -std::real(cbetaMu);
		double ibetaMu =  std::imag(cbetaMu);

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
		double a    = CalculateWidth();
		double h    = CalculateHeight();

        cN10 *= -sqrt(2.00*freq/(a*h));

		return cN10;
	}

   /**********************************************************************************
    * - Calculate area
    **********************************************************************************/
    double RWPortTE10_3th::CalculateArea() 
	{
		double n[3],v1[3],v2[3];

        v2[0] = mNodes[2]->X() - mNodes[0]->X();
        v2[1] = mNodes[2]->Y() - mNodes[0]->Y();
        v2[2] = mNodes[2]->Z() - mNodes[0]->Z();

        v1[0] = mNodes[1]->X() - mNodes[0]->X();
        v1[1] = mNodes[1]->Y() - mNodes[0]->Y();
        v1[2] = mNodes[1]->Z() - mNodes[0]->Z();

        //area = 0.5 * ||v2 x v1||
        n[0] = v2[1]*v1[2] - v2[2]*v1[1];
        n[1] = v2[2]*v1[0] - v2[0]*v1[2];
        n[2] = v2[0]*v1[1] - v2[1]*v1[0];

		return ( 0.5 * sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) );
	}

   /**********************************************************************************
    * - Calculates b (vertical) vector
    **********************************************************************************/
	void RWPortTE10_3th::Calculate_b(double* b) 
	{
		Vector<double> cIBC   = (*mProperties)(COMPLEX_IBC),
					   cIBC2o = (*mProperties)(COMPLEX_IBC_2o),
					   cNF    = (*mProperties)(COMPLEX_NEUMANN_FLOW); 
		
		b[0] = cNF[0] - cIBC  [1];
		b[1] = cNF[1] - cIBC2o[0];
		b[2] = cNF[2] - cIBC2o[1];

		double h     = CalculateHeight();
		double norma = 1.00/h;

        b[0] *= norma; 
        b[1] *= norma; 
        b[2] *= norma;
	}

   /**********************************************************************************
    * - Calculates t (horizontal) vector
    **********************************************************************************/
	void RWPortTE10_3th::Calculate_t(double* t) 
	{
		Vector<double> cIBC   = (*mProperties)(COMPLEX_IBC),
					   cIBC2o = (*mProperties)(COMPLEX_IBC_2o),
					   cNF    = (*mProperties)(COMPLEX_NEUMANN_FLOW); 
		
		t[0] = cNF[3] - cIBC  [1];
		t[1] = cNF[4] - cIBC2o[0];
		t[2] = cNF[5] - cIBC2o[1];

		double a     = CalculateWidth();
		double norma = 1.00/a;

        t[0] *= norma; 
        t[1] *= norma; 
        t[2] *= norma;
	}

   /**********************************************************************************
    * - Calculates r0 vector (distance to origin)
    **********************************************************************************/
	void RWPortTE10_3th::Calculate_r0(double* r0) 
	{
		Vector<double> cIBC   = (*mProperties)(COMPLEX_IBC),
					   cIBC2o = (*mProperties)(COMPLEX_IBC_2o); 
		
		r0[0] = cIBC  [1];
		r0[1] = cIBC2o[0];
	    r0[2] = cIBC2o[1];
	}

   /**********************************************************************************
    * - Calculation of the residual vector (Linear interpolation)
    **********************************************************************************/
    void RWPortTE10_3th::GetResidualVector(Vector<std::complex<double> >& ResidualVector) 
    {
        const int numNodes = 10;
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

		int gp;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		Matrix<double> N;

		Lagrange2D_Ni_3rd(N, cX, cY);

		std::vector<double> Eb(nGaussPoints);

		// kc = (pi/a)
		double pi = 3.1415926535897932384626433832795;
		double a  = CalculateWidth();
		double kc = pi/a;

		double t[3];   
		Calculate_t(t); 

		double r0[3];   
		Calculate_r0(r0); 

		double alpha, beta;
		double rX,rY,rZ;
		double Xt;

		double X1 = mNodes[0]->X(), Y1 = mNodes[0]->Y(), Z1 = mNodes[0]->Z();
		double X2 = mNodes[1]->X(), Y2 = mNodes[1]->Y(), Z2 = mNodes[1]->Z();
		double X3 = mNodes[2]->X(), Y3 = mNodes[2]->Y(), Z3 = mNodes[2]->Z();

		for ( gp=0; gp<nGaussPoints; gp++ ) 
        {
            alpha = cX[gp];
			beta  = cY[gp];

			rX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
            rY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		    rZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

			Xt = ( rX - r0[0] ) * t[0] + ( rY - r0[1] ) * t[1] + ( rZ - r0[2] ) * t[2]; 

			Eb[gp] = sin(kc*Xt);  
        }

		double b[3];   
		Calculate_b(b); 

		double area  = CalculateArea(); 
		double jacob = 2.00*area;

		std::complex<double> cbetaMu =  CalculateBetaMuCoef();
		std::complex<double> cN10    =  CalculateN10();

        std::complex<double> cCte = -2.00*jacob*(cbetaMu*cN10);
	
		for (int i=0; i<numNodes; i++)
        {
            double intEbNi = 0.00;
            for ( gp=0; gp<nGaussPoints; gp++ ) intEbNi += W[gp]*Eb[gp]*N[i][gp];

            ResidualVector[i             ] = cCte*(b[0]*intEbNi);
            ResidualVector[i+ numNodes   ] = cCte*(b[1]*intEbNi);
			ResidualVector[i+(numNodes*2)] = cCte*(b[2]*intEbNi);
        }  
    }

   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void RWPortTE10_3th::GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix) 
    {
		const int numNodes = 10;
		const int numDofs  = numNodes*3;

		StiffMatrix.Resize(numDofs,numDofs);

        int i,j;

		double area = CalculateArea(); 

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double>  W;

		int gp;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		Matrix<double> N;

		Lagrange2D_Ni_3rd(N, cX, cY);

		double intNiNj;

		double jacob = 2.00*area;

		std::complex<double> cbetaMu =  CalculateBetaMuCoef(); 
		std::complex<double> ccteNN  = -jacob*cbetaMu;
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

   /**********************************************************************************
    * - Returns current port number
    **********************************************************************************/
	double RWPortTE10_3th::getPortNumber()
	{
		Vector<double> cIBC = (*mProperties)(COMPLEX_IBC);

		return cIBC[0];
	}

   /**********************************************************************************
    * - Returns TE10 projection on a plane surface 
    **********************************************************************************/
	std::complex<double> RWPortTE10_3th::getProjection()
	{
		const int numNodes = 10;

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double>  W;

		int gp;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);

		Matrix<double> N;

		Lagrange2D_Ni_3rd(N, cX, cY);

		std::vector<std::complex<double> > cEb(nGaussPoints);

		std::complex<double> cE_x, cE_y, cE_z;

		double b[3];   
		Calculate_b(b);

		//Eb in gauss points
		for ( gp=0; gp<nGaussPoints; gp++ ) 
        {
			std::complex<double> cEbInGp(0.00,0.00);

			for (int i=0; i<numNodes; i++ ) 
			{
				cE_x = (*mProperties)(cEx, *mNodes[i]);
				cE_y = (*mProperties)(cEy, *mNodes[i]);
				cE_z = (*mProperties)(cEz, *mNodes[i]);

			    cEbInGp += N[i][gp] * ( b[0]*cE_x + b[1]*cE_y + b[2]*cE_z );
			}
           
			cEb[gp] = cEbInGp;
        }
		
		// kc = (pi/a)
		double pi = 3.1415926535897932384626433832795;
		double a  = CalculateWidth();
		double kc = pi/a;

		double t[3];   
		Calculate_t(t); 

		double r0[3];   
		Calculate_r0(r0); 

		double alpha, beta;
		double rX,rY,rZ;
		double Xt;

		double X1 = mNodes[0]->X(), Y1 = mNodes[0]->Y(), Z1 = mNodes[0]->Z();
		double X2 = mNodes[1]->X(), Y2 = mNodes[1]->Y(), Z2 = mNodes[1]->Z();
		double X3 = mNodes[2]->X(), Y3 = mNodes[2]->Y(), Z3 = mNodes[2]->Z();

		std::complex<double> TE10Project(0.00,0.00);

		//Integral ( Eb*sin(kc*xt) ) over element surface
		for ( gp=0; gp<nGaussPoints; gp++ ) 
        {
            alpha = cX[gp];
			beta  = cY[gp];

			rX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
            rY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		    rZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

			Xt = ( rX - r0[0] ) * t[0] + ( rY - r0[1] ) * t[1] + ( rZ - r0[2] ) * t[2]; 

			TE10Project += ( W[gp] * sin(kc*Xt) ) * cEb[gp];  
        }

		double portNumber = getPortNumber();

        //Cuando los puertos estan dentro del volumen los elementos se repiten
		if ( int(portNumber) != portNumber ) TE10Project *= 0.5;

		return TE10Project;
	}

   /**********************************************************************************
    * - Returns a normalization term (useful when using PMC symmetries)
    **********************************************************************************/
	std::complex<double> RWPortTE10_3th::getNormalization()
	{
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints2D_Order12(cX, cY, W);
		
		// kc = (pi/a)
		double pi = 3.1415926535897932384626433832795;
		double a  = CalculateWidth();
		double kc = pi/a;

		double t[3];   
		Calculate_t(t); 

		double r0[3];   
		Calculate_r0(r0); 

		double alpha, beta;
		double rX,rY,rZ;
		double Xt;

		double X1 = mNodes[0]->X(), Y1 = mNodes[0]->Y(), Z1 = mNodes[0]->Z();
		double X2 = mNodes[1]->X(), Y2 = mNodes[1]->Y(), Z2 = mNodes[1]->Z();
		double X3 = mNodes[2]->X(), Y3 = mNodes[2]->Y(), Z3 = mNodes[2]->Z();

		double intSin2 = 0.00;

		//Integral ( sin(kc*xt)*sin(kc*xt) ), para normalizar en caso de simetrías
		for (int gp=0; gp<nGaussPoints; gp++ ) 
        {
            alpha = cX[gp];
			beta  = cY[gp];

			rX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
            rY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		    rZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

			Xt = ( rX - r0[0] ) * t[0] + ( rY - r0[1] ) * t[1] + ( rZ - r0[2] ) * t[2]; 

			intSin2 += W[gp] * sin(kc*Xt) * sin(kc*Xt);  
        }

		double portNumber = getPortNumber();

        //Cuando los puertos estan dentro del volumen los elementos se repiten
		if ( int(portNumber) != portNumber ) intSin2 *= 0.5;
		
		std::complex<double> Normalization(intSin2,0.00);

		std::complex<double> cN10 = CalculateN10();

		Normalization *= cN10;

		return Normalization;
	}

} /* end namespace Kratos */ 