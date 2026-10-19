
#include "VolumeElement_2nd.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
   /**********************************************************************************
    * -Global index of the nodes
    **********************************************************************************/	
    void VolumeElement_2nd::GetEquationIdVector(Vector<int>& EquationId)
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

   /**********************************************************
    * - Calculate element volume
    **********************************************************/
    double VolumeElement_2nd::Calculate_Volume()
    {
        double det;

        det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs(det/6.00);
    }

   /**********************************************************
    * - Calculation of dN/dx, dN/dy, dN/dz column vector
    **********************************************************/
    void VolumeElement_2nd::Calculate_DN(double DN[3][4])
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

   /***********************************************************************************
    * - Calculation of dN/dx, dN/dy, dN/dz on cX,cY,cZ points
    ***********************************************************************************/
    void VolumeElement_2nd::LagrangeDerivatives3D_2nd(Matrix<double>& dNdx,
		                                              Matrix<double>& dNdy,
													  Matrix<double>& dNdz,
													  std::vector<double>& cX, 
													  std::vector<double>& cY, 
													  std::vector<double>& cZ)
    {
		int nNodes  = 10;
		int nPoints = cX.size();

		dNdx.Resize(nNodes, nPoints);
		dNdy.Resize(nNodes, nPoints);
		dNdz.Resize(nNodes, nPoints);
		
		double DN[3][4];
        Calculate_DN(DN);

		double dL1dx,dL2dx,dL3dx,dL4dx;
        double dL1dy,dL2dy,dL3dy,dL4dy;
		double dL1dz,dL2dz,dL3dz,dL4dz;

        dL1dx = DN[0][0], dL2dx = DN[0][1], dL3dx = DN[0][2], dL4dx = DN[0][3];
        dL1dy = DN[1][0], dL2dy = DN[1][1], dL3dy = DN[1][2], dL4dy = DN[1][3];
		dL1dz = DN[2][0], dL2dz = DN[2][1], dL3dz = DN[2][2], dL4dz = DN[2][3];

		double L1,L2,L3,L4;

		//dNi/dx ,dNi/dy, dNi/dz
        for (int p=0; p<nPoints; p++ ) 
        {
            L2 = cX[p];
            L3 = cY[p];
			L4 = cZ[p];
            L1 = 1.00 - L2 - L3 - L4;

			int nn;
			
			////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////   POINTS   ////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////
			double L;
			double dNdL;
			
			//////////////////////////////////////// Points: 0  /////////////////////////////////////////////////
			nn   = 0;
			L    = L1;
			dNdL = 4.00 * L - 1.00;
			
			dNdx[nn][p] = dNdL * dL1dx;
			dNdy[nn][p] = dNdL * dL1dy;
			dNdz[nn][p] = dNdL * dL1dz;

			//////////////////////////////////////// Points: 1  /////////////////////////////////////////////////
			nn   = 1;
			L    = L2;
			dNdL = 4.00 * L - 1.00;
			
			dNdx[nn][p] = dNdL * dL2dx;
			dNdy[nn][p] = dNdL * dL2dy;
			dNdz[nn][p] = dNdL * dL2dz;

			//////////////////////////////////////// Points: 2  /////////////////////////////////////////////////
			nn   = 2;
			L    = L3;
			dNdL = 4.00 * L - 1.00;
			
			dNdx[nn][p] = dNdL * dL3dx;
			dNdy[nn][p] = dNdL * dL3dy;
			dNdz[nn][p] = dNdL * dL3dz;

			//////////////////////////////////////// Points: 3  /////////////////////////////////////////////////
			nn   = 3;
			L    = L4;
			dNdL = 4.00 * L - 1.00;
			
			dNdx[nn][p] = dNdL * dL4dx;
			dNdy[nn][p] = dNdL * dL4dy;
			dNdz[nn][p] = dNdL * dL4dz;

			////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////   EDGES   /////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
			double Li,Lj;
            double dNdLi, dNdLj;
			double dLidx, dLidy, dLidz;
			double dLjdx, dLjdy, dLjdz;
			
            ////////////////////////////////////// Edge: [0 1] /////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L2; dLjdx = dL2dx; dLjdy = dL2dy; dLjdz = dL2dz;

			//nodo: 4
            nn    = 4;
			dNdLi = 4.00 * Lj;
			dNdLj = 4.00 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;
			
			////////////////////////////////////// Edge: [1 2] /////////////////////////////////////////////////
			Li = L2; dLidx = dL2dx; dLidy = dL2dy; dLidz = dL2dz;
			Lj = L3; dLjdx = dL3dx; dLjdy = dL3dy; dLjdz = dL3dz;

            //nodo: 5
            nn    = 5;
			dNdLi = 4.00 * Lj;
			dNdLj = 4.00 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [2 0] /////////////////////////////////////////////////
			Li = L3; dLidx = dL3dx; dLidy = dL3dy; dLidz = dL3dz;
			Lj = L1; dLjdx = dL1dx; dLjdy = dL1dy; dLjdz = dL1dz;

			//nodo: 6
            nn    = 6;
			dNdLi = 4.00 * Lj;
			dNdLj = 4.00 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [0 3] /////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 7
            nn    = 7;
			dNdLi = 4.00 * Lj;
			dNdLj = 4.00 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [1 3] /////////////////////////////////////////////////
			Li = L2; dLidx = dL2dx; dLidy = dL2dy; dLidz = dL2dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 8
            nn    = 8;
			dNdLi = 4.00 * Lj;
			dNdLj = 4.00 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [2 3] /////////////////////////////////////////////////
			Li = L3; dLidx = dL3dx; dLidy = dL3dy; dLidz = dL3dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 9
            nn    = 9;
			dNdLi = 4.00 * Lj;
			dNdLj = 4.00 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;	
        }
	}

   /**********************************************************************************
    * - Set singularities
    **********************************************************************************/	
    void VolumeElement_2nd::SetPeso(double Peso)
    {
        mPeso = Peso;
    }

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Distance to an edge
	////////////////////////////////////////////////////////////////////////////////////////////////////
	double VolumeElement_2nd::RtoS(double alpha, double beta, double gamma)
    {
		double X1 = mNodes[0]->X(), Y1 = mNodes[0]->Y();
		double X2 = mNodes[1]->X(), Y2 = mNodes[1]->Y();
		double X3 = mNodes[2]->X(), Y3 = mNodes[2]->Y();
		double X4 = mNodes[3]->X(), Y4 = mNodes[3]->Y();
        
		double X = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + gamma * ( X4 - X1 ) + X1;
        double Y = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + gamma * ( Y4 - Y1 ) + Y1;

        double X0 = 20.00;
        double Y0 = 5.00;

        return sqrt( (X-X0)*(X-X0) + (Y-Y0)*(Y-Y0) );
    }

   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void VolumeElement_2nd::GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix) 
    {  
		int i, j, gp;
		const int numNodes = 10;
		const int numDofs  = numNodes*3;

		StiffMatrix.Resize(numDofs,numDofs);

		// Material properties
		double eo = 8.8541878176e-12;
		double pi = 3.1415926535897932384626433832795;
		double mo = (4.00e-7)*pi;
        
        double freq     = (*mProperties)(FREQUENCY),
               sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY),
               eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * eo,
               eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * eo,
               mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mo,
               mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mo;

        std::complex<double> cEps( eps_real, eps_imag + (sigma/freq) );
        std::complex<double> cMu ( mu_real , mu_imag );
      
		// Gaussian derivation weights
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order4(cX, cY, cZ, W);

		double detJ = 6.00*Calculate_Volume();

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions and derivatives
		Matrix<double> N;

		Lagrange3D_Ni_2nd(N, cX, cY, cZ);

		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2nd(dNdx, dNdy, dNdz, cX, cY, cZ);

		if ( mPeso != 0.00 ) 
		{
			/////////////////////////////////////////////////////////////////////////////////////////////////////
			double intCurlDiv, intNN;

			std::complex<double> eKij;

			std::complex<double> cteCurlDiv = 1.00 / cMu;
			std::complex<double> cteNN      = (freq*freq) * cEps;  

			////////////////////////////////////// Bloque [Kxx],[Kyy],[Kzz] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=i; j<numNodes; j++ )
				{
					intCurlDiv = 0.00; 
					intNN      = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
						intNN      += W[gp] * (    N[i][gp] *    N[j][gp] );
					}

					eKij = (cteCurlDiv * intCurlDiv) - (cteNN * intNN);

					StiffMatrix[i             ][j             ] = eKij;
					StiffMatrix[i+   numNodes ][j+   numNodes ] = eKij;
					StiffMatrix[i+(2*numNodes)][j+(2*numNodes)] = eKij;
				}
			}

			////////////////////////////////////// Bloque [Kxy] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=0; j<numNodes; j++ )
				{
					intCurlDiv = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdy[j][gp] - dNdy[i][gp] * dNdx[j][gp] );
					}
	                
					StiffMatrix[i][j+numNodes] = (cteCurlDiv * intCurlDiv);
				}
			}

			////////////////////////////////////// Bloque [Kxz] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=0; j<numNodes; j++ )
				{
					intCurlDiv = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdx[j][gp] );
					}
	                
					StiffMatrix[i][j+(2*numNodes)] = (cteCurlDiv * intCurlDiv);
				}
			}

			////////////////////////////////////// Bloque [Kyz] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=0; j<numNodes; j++ )
				{
					intCurlDiv = 0.00; 

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurlDiv += W[gp] * ( dNdy[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdy[j][gp] );
					}
	                
					StiffMatrix[i+numNodes][j+(2*numNodes)] = (cteCurlDiv * intCurlDiv);
				}
			}
		}
		else
		{
			/////////////////////////////////////////////////////////////////////////////////////////////////////
			double intCurl, intNN;

			std::complex<double> cteCurl = 1.00 / cMu;
			std::complex<double> cteNN   = (freq*freq) * cEps;    

			////////////////////////////////////// Bloque [Kxx] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=i; j<numNodes; j++ )
				{
					intCurl = 0.00;
					intNN   = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
						intNN   += W[gp] * (    N[i][gp] *    N[j][gp] );
					}
	                
					StiffMatrix[i][j] = (cteCurl * intCurl) - (cteNN * intNN);
				}
			}

			////////////////////////////////////// Bloque [Kyy] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=i; j<numNodes; j++ )
				{
					intCurl = 0.00;
					intNN   = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
						intNN   += W[gp] * (    N[i][gp] *    N[j][gp] );
					}
	                
					StiffMatrix[i+numNodes][j+numNodes] = (cteCurl * intCurl) - (cteNN * intNN);
				}
			}

			////////////////////////////////////// Bloque [Kzz] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=i; j<numNodes; j++ )
				{
					intCurl = 0.00; 
					intNN   = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] );
						intNN   += W[gp] * (    N[i][gp] *    N[j][gp] );
					}
	                
					StiffMatrix[i+(2*numNodes)][j+(2*numNodes)] = (cteCurl * intCurl) - (cteNN * intNN);
				}
			}

			////////////////////////////////////// Bloque [Kxy] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=0; j<numNodes; j++ )
				{
					intCurl = 0.00; 

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdy[i][gp] * dNdx[j][gp] );
					}
	                
					StiffMatrix[i][j+numNodes] = - (cteCurl * intCurl);
				}
			}

			////////////////////////////////////// Bloque [Kxz] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=0; j<numNodes; j++ )
				{
					intCurl = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdz[i][gp] * dNdx[j][gp] );
					}
	                
					StiffMatrix[i][j+(2*numNodes)] = - (cteCurl * intCurl);
				}
			}

			////////////////////////////////////// Bloque [Kyz] ////////////////////////////////////////////////
			for ( i=0; i<numNodes; i++ )
			{
				for ( j=0; j<numNodes; j++ )
				{
					intCurl = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdz[i][gp] * dNdy[j][gp] );
					}
	                
					StiffMatrix[i+numNodes][j+(2*numNodes)] = - (cteCurl * intCurl);
				}
			}
		}

        //////////////////////////////////// Bloque Simetrico /////////////////////////////////////////////
        for ( i=0; i<numDofs; i++ )
        {
            for ( j=i+1; j<numDofs; j++ )
            {
                StiffMatrix[j][i] =  StiffMatrix[i][j];
            }
        }      
    }

   /*******************************************************************
    * - Calculate derivatives
    ********************************************************************/
    void VolumeElement_2nd::DerivativesInNodes(Matrix<double>& ndNdx, 
		                                       Matrix<double>& ndNdy, 
											   Matrix<double>& ndNdz)
    {
		const int numNodes = 10;
		
		std::vector<double> cX(numNodes); 
		std::vector<double>	cY(numNodes);
		std::vector<double> cZ(numNodes);

		//Points 0,1,2,3
		cX[0] = 0.00; cY[0] = 0.00; cZ[0] = 0.00;
		cX[1] = 1.00; cY[1] = 0.00; cZ[1] = 0.00;
		cX[2] = 0.00; cY[2] = 1.00; cZ[2] = 0.00;
		cX[3] = 0.00; cY[3] = 0.00; cZ[3] = 1.00;
		
		//Edge [0 1]
		cX[4] = 0.50; cY[4] = 0.00; cZ[4] = 0.00; 
		//Edge [1 2]
		cX[5] = 0.50; cY[5] = 0.50; cZ[5] = 0.00; 
		//Edge [2 0]
		cX[6] = 0.00; cY[6] = 0.50; cZ[6] = 0.00; 
		//Edge [0 3]
		cX[7] = 0.00; cY[7] = 0.00; cZ[7] = 0.50; 
		//Edge [1 3]
		cX[8] = 0.50; cY[8] = 0.00; cZ[8] = 0.50;
		//Edge [2 3]
		cX[9] = 0.00; cY[9] = 0.50; cZ[9] = 0.50; 

		LagrangeDerivatives3D_2nd(ndNdx, ndNdy, ndNdz, cX, cY, cZ);
	}

   /*******************************************************************
    * - Calculate derivatives
    ********************************************************************/
    void VolumeElement_2nd::CalculateDerivatives(Vector<std::complex<double> >& ncHx,
		                                         Vector<std::complex<double> >& ncHy,
												 Vector<std::complex<double> >& ncHz)
    {
        int i,n;

		const int numNodes = 10;
        
		double pi = 3.1415926535897932384626433832795;
		double mo = (4.00e-7)*pi;
        
        double freq     = (*mProperties)(FREQUENCY),
               mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mo,
               mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mo;

        std::complex<double> cMu(mu_real, mu_imag);
        std::complex<double> cUnit(0.00, 1.00);
        std::complex<double> inv_jwmu = 1.00/(freq*cUnit*cMu);

		Matrix<double> ndNdx; 
		Matrix<double> ndNdy; 
		Matrix<double> ndNdz;

		DerivativesInNodes(ndNdx, ndNdy, ndNdz);

		Vector<std::complex<double> > ncEx(numNodes);
		Vector<std::complex<double> > ncEy(numNodes);
		Vector<std::complex<double> > ncEz(numNodes);
		     
		for(n=0; n<numNodes; n++)
		{
			ncEx[n] = (*mProperties)(cEx, *mNodes[n]);
			ncEy[n] = (*mProperties)(cEy, *mNodes[n]);
			ncEz[n] = (*mProperties)(cEz, *mNodes[n]);
		}

		std::complex<double> czero(0.00,0.00);

        ncHx.resize(numNodes);
		ncHy.resize(numNodes);
		ncHz.resize(numNodes);

		for(n=0; n<numNodes; n++)
		{
			ncHx[n] = czero;
			ncHy[n] = czero;
			ncHz[n] = czero;

			for(i=0; i<numNodes; i++)
			{
				ncHx[n] += ndNdy[i][n]*ncEz[i] - ndNdz[i][n]*ncEy[i];
				ncHy[n] += ndNdz[i][n]*ncEx[i] - ndNdx[i][n]*ncEz[i];
				ncHz[n] += ndNdx[i][n]*ncEy[i] - ndNdy[i][n]*ncEx[i];
			}
		}
			
		ncHx *= inv_jwmu;
		ncHy *= inv_jwmu;
		ncHz *= inv_jwmu;       
    }

} /* end namespace Kratos */ 