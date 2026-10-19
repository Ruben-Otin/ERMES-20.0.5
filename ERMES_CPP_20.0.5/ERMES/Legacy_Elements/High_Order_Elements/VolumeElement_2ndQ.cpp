
#include "VolumeElement_2ndQ.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"


namespace Kratos
{   
   /**********************************************************************************
    * -Global index of the nodes
    **********************************************************************************/	
    void VolumeElement_2ndQ::GetEquationIdVector(Vector<int>& EquationId)
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
    * - Calculation of dN/dx, dN/dy, dN/dz on cX, cY, cZ points
    ***********************************************************************************/
    void VolumeElement_2ndQ::LagrangeDerivatives3D_2ndQ(Matrix<double>& dNdx,
		                                                Matrix<double>& dNdy,
													    Matrix<double>& dNdz,
													    std::vector<double>& cX, 
													    std::vector<double>& cY, 
													    std::vector<double>& cZ)
    {
		int i, p;
		int numNodes  = 10;
		int numPoints = cX.size();

		Matrix<double> dNda, dNdb, dNdu;

		NaturalDerivatives3D_2ndQ(dNda, dNdb, dNdu, cX, cY, cZ);

		Matrix<double> invJ(3,3);

		dNdx.Resize(numNodes, numPoints);
		dNdy.Resize(numNodes, numPoints);
		dNdz.Resize(numNodes, numPoints);

		//dNi/dx ,dNi/dy, dNi/dz
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
	
			double  detJ = dXda*dYdb*dZdu + dXdb*dZda*dYdu + dXdu*dYda*dZdb - dXda*dZdb*dYdu - dXdb*dYda*dZdu - dXdu*dZda*dYdb;
			double idetJ = 1.00/detJ;

			invJ[0][0] = idetJ*( dYdb*dZdu - dZdb*dYdu ); invJ[0][1] = idetJ*( dZda*dYdu - dYda*dZdu ); invJ[0][2] = idetJ*( dYda*dZdb - dZda*dYdb );
			invJ[1][0] = idetJ*( dZdb*dXdu - dXdb*dZdu ); invJ[1][1] = idetJ*( dXda*dZdu - dZda*dXdu ); invJ[1][2] = idetJ*( dZda*dXdb - dXda*dZdb );
			invJ[2][0] = idetJ*( dXdb*dYdu - dYdb*dXdu ); invJ[2][1] = idetJ*( dYda*dXdu - dXda*dYdu ); invJ[2][2] = idetJ*( dXda*dYdb - dYda*dXdb );

			for ( i=0; i<numNodes; i++ ) 
			{
				dNdx[i][p] = invJ[0][0] * dNda[i][p] + invJ[0][1] * dNdb[i][p] + invJ[0][2] * dNdu[i][p];
				dNdy[i][p] = invJ[1][0] * dNda[i][p] + invJ[1][1] * dNdb[i][p] + invJ[1][2] * dNdu[i][p];
				dNdz[i][p] = invJ[2][0] * dNda[i][p] + invJ[2][1] * dNdb[i][p] + invJ[2][2] * dNdu[i][p];
			}
		}
	}

   /***********************************************************************************
    * - Calculation of dN/da, dN/db, dN/du on cX,cY,cZ points
    ***********************************************************************************/
    void VolumeElement_2ndQ::NaturalDerivatives3D_2ndQ(Matrix<double>& dNda,
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
    void VolumeElement_2ndQ::Calculate_detJ(std::vector<double>& detJ,
		                                    std::vector<double>& cX, 
											std::vector<double>& cY, 
											std::vector<double>& cZ)
    {
		int i, p;
		int numNodes  = 10;
		int numPoints = cX.size();

		detJ.resize(numPoints,0.00);

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
    * - Set singularities
    **********************************************************************************/	
    void VolumeElement_2ndQ::SetPeso(double Peso)
    {
		mPeso = Peso;
    }

   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void VolumeElement_2ndQ::GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix) 
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

		std::vector<double> detJ;

        Calculate_detJ(detJ, cX, cY, cZ);

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ[gp]; 

        // Basis functions and derivatives
		Matrix<double> N;

		Lagrange3D_Ni_2nd(N, cX, cY, cZ);

		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2ndQ(dNdx, dNdy, dNdz, cX, cY, cZ);

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
    void VolumeElement_2ndQ::DerivativesInNodes(Matrix<double>& ndNdx, 
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

		LagrangeDerivatives3D_2ndQ(ndNdx, ndNdy, ndNdz, cX, cY, cZ);
	}

   /*******************************************************************
    * - Calculate derivatives
    ********************************************************************/
    void VolumeElement_2ndQ::CalculateDerivatives(Vector<std::complex<double> >& ncHx,
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
		     
		for( n=0; n<numNodes; n++ )
		{
			ncEx[n] = (*mProperties)(cEx, *mNodes[n]);
			ncEy[n] = (*mProperties)(cEy, *mNodes[n]);
			ncEz[n] = (*mProperties)(cEz, *mNodes[n]);
		}

		std::complex<double> czero(0.00,0.00);

        ncHx.resize(numNodes);
		ncHy.resize(numNodes);
		ncHz.resize(numNodes);

		for( n=0; n<numNodes; n++ )
		{
			ncHx[n] = czero;
			ncHy[n] = czero;
			ncHz[n] = czero;

			for( i=0; i<numNodes; i++ )
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