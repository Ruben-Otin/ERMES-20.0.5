
#include "VolumeElement_2ndQ_DivD_Plasma.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
    //************************************************************************************
    //* - Global index of the nodes
    //************************************************************************************	
    void VolumeElement_2ndQ_DivD_Plasma::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize(numDofs);

        for (int i=0; i<numNodes; i++) 
        { 
            EquationId[i              ] = mNodes[i]->pDofcEx()->EquationId();
            EquationId[i + numNodes   ] = mNodes[i]->pDofcEy()->EquationId();
			EquationId[i +(numNodes*2)] = mNodes[i]->pDofcEz()->EquationId();
        }
    }

    //************************************************************************************
    //* - Set volume element as IHL material or plasma
    //************************************************************************************	
    bool VolumeElement_2ndQ_DivD_Plasma::Set_Plasma_or_IHL() 
    { 
        // Material properties
        double sigma    = (*mProperties)( IHL_ELECTRIC_CONDUCTIVITY  );

        double epr_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY );
        double epr_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY );

        double mu_real  = (*mProperties)( REAL_MAGNETIC_PERMEABILITY );
        double mu_imag  = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY );

        // If the material is vaccum then it is considered as plasma. 
        // This allows the use of different plasmas volumes (good for PBC and visualization).
        // We only consider the first one written in the plasma file *-8.dat
        if ( (sigma == 0.0) && (epr_real == 1.0) && (epr_imag == 0.0) && (mu_real == 1.0) && (mu_imag == 0.0) )
        {
            Is_IHL = false;        
        }
        else
        {
            Is_IHL = true;        
        }

        return Is_IHL;
    }

    //*************************************************************************************
    //* - Calculation of dN/dx, dN/dy, dN/dz on cX, cY, cZ points
    //*************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::LagrangeDerivatives3D_2ndQ( Matrix<double>& dNdx,
																     Matrix<double>& dNdy,
																     Matrix<double>& dNdz,
                                                                     Vector<double>&   cX, 
                                                                     Vector<double>&   cY, 
                                                                     Vector<double>&   cZ )
    {
		int i, p;
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

    //*******************************************************************************************
    //* - Calculation of dN/da, dN/db, dN/du on cX,cY,cZ points
    //*******************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::NaturalDerivatives3D_2ndQ( Matrix<double>& dNda,
															        Matrix<double>& dNdb,
															        Matrix<double>& dNdu,
															        Vector<double>&   cX, 
															        Vector<double>&   cY, 
															        Vector<double>&   cZ )
    {
		int p;
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

    //************************************************************************************
    //* - Calculation of the Jacobian determinant
    //************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::Calculate_detJ( Vector<double>& detJ,
		                                                 Vector<double>&   cX, 
											             Vector<double>&   cY, 
											             Vector<double>&   cZ )
    {
		int i, p;
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

    //*********************************************************************************************
	//* - Generates the elemental stiffness matrix
	//*********************************************************************************************
	void VolumeElement_2ndQ_DivD_Plasma::GetStiffMatrix( Matrix<std::complex<double> >& StiffMatrix )
	{
		if (mPeso != 0.0)
        {
            Generate_UpperDiagonal_CurlDiv ( StiffMatrix );
        }
        else
        {
            Generate_UpperDiagonal_CurlCurl( StiffMatrix );
        }
        
		Generate_LowerDiagonal( StiffMatrix );
	}

    //************************************************************************************************************************
    //* - Generates the diagonal and the upper diagonal part of the elemental stiffness matrix (Un-gauged elements)
    //************************************************************************************************************************	
    void VolumeElement_2ndQ_DivD_Plasma::Generate_UpperDiagonal_CurlCurl( Matrix<std::complex<double> >& StiffMatrix )
    {  
        // Loop indexes
		int i, j, gp;

		// Initializing elemental stiff matrix
		StiffMatrix.Resize( numDofs, numDofs );

		// Magnetic material properties
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mo;

        std::complex<double> cMu ( mu_real , mu_imag );

		// Gaussian derivation weights
		Vector<double> cX; 
		Vector<double> cY;
		Vector<double> cZ;
		Vector<double>  W;

        int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		Vector<double> detJ;

        Calculate_detJ( detJ, cX, cY, cZ );

		for ( gp=0; gp<numGaussPoints; gp++ ) W[gp] *= detJ[gp];
		
        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2nd( N, cX, cY, cZ );

		// Cartesian derivatives of the basis functions on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2ndQ( dNdx, dNdy, dNdz, cX, cY, cZ );

        // Extracting permitivitty tensor on nodes
        Vector<Matrix<std::complex<double> > > TEn;

        PermittivityTensor_OnNodes( TEn );

        // Permittivity tensor on Gauss points
        Vector<Matrix<std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, TEn, N );

		// Filling stiffness matrix with cotributions from Curl(Ni)*(1/mu)*Curl(Nj) and freq2*Ni*[Ep]*Nj
		std::complex<double> invMu( 1.00 / cMu  );
		std::complex<double> freq2( freq * freq );    

		//////////////////////////////////////// Block [Kxx] ////////////////////////////////////////////////
		for ( i=0; i<numNodes; i++ )
		{
			for ( j=i; j<numNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
		        std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

				for ( gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_Curl += W[gp] * ( dNdy[i][gp] *    dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
					intg_Ni_TEp_Nj += W[gp] * (    N[i][gp] * TEgp[gp][0][0] *    N[j][gp]               );
				}
	           
				StiffMatrix[i][j] = ( invMu * intg_Curl_Curl ) - ( freq2 * intg_Ni_TEp_Nj );
			}
		}
		//////////////////////////////////////// Block [Kyy] ////////////////////////////////////////////////
		for ( i=0; i<numNodes; i++ )
		{
			for ( j=i; j<numNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
		        std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

				for ( gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_Curl += W[gp] * ( dNdx[i][gp] *    dNdx[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
					intg_Ni_TEp_Nj += W[gp] * (    N[i][gp] * TEgp[gp][1][1] *    N[j][gp]               );
				}
	                
				StiffMatrix[i+numNodes][j+numNodes] = ( invMu * intg_Curl_Curl ) - ( freq2 * intg_Ni_TEp_Nj );
			}
		}
		//////////////////////////////////////// Block [Kzz] ////////////////////////////////////////////////
		for ( i=0; i<numNodes; i++ )
		{
			for ( j=i; j<numNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
		        std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

				for ( gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_Curl += W[gp] * ( dNdx[i][gp] *    dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] );
					intg_Ni_TEp_Nj += W[gp] * (    N[i][gp] * TEgp[gp][2][2] *    N[j][gp]               );
				}
	                
				StiffMatrix[i+(2*numNodes)][j+(2*numNodes)] = ( invMu * intg_Curl_Curl ) - ( freq2 * intg_Ni_TEp_Nj );
			}
		}
		//////////////////////////////////////// Block [Kxy] ////////////////////////////////////////////////
		for ( i=0; i<numNodes; i++ )
		{
			for ( j=0; j<numNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
		        std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

				for ( gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_Curl += W[gp] * ( - dNdy[i][gp] *    dNdx[j][gp]            );
					intg_Ni_TEp_Nj += W[gp] * (      N[i][gp] * TEgp[gp][0][1] * N[j][gp] );
				}
	                
				StiffMatrix[i][j+numNodes] = ( invMu * intg_Curl_Curl ) - ( freq2 * intg_Ni_TEp_Nj );;
			}
		}
		/////////////////////////////////////// Block [Kxz] ////////////////////////////////////////////////
		for ( i=0; i<numNodes; i++ )
		{
			for ( j=0; j<numNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
		        std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

				for ( gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_Curl += W[gp] * ( - dNdz[i][gp] *    dNdx[j][gp]            );
					intg_Ni_TEp_Nj += W[gp] * (      N[i][gp] * TEgp[gp][0][2] * N[j][gp] );
				}
	                
				StiffMatrix[i][j+(2*numNodes)] = ( invMu * intg_Curl_Curl ) - ( freq2 * intg_Ni_TEp_Nj );
			}
		}
		////////////////////////////////////// Block [Kyz] ////////////////////////////////////////////////
		for ( i=0; i<numNodes; i++ )
		{
			for ( j=0; j<numNodes; j++ )
			{
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
		        std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

				for ( gp=0; gp<numGaussPoints; gp++ )
				{
					intg_Curl_Curl += W[gp] * ( - dNdz[i][gp] *    dNdy[j][gp]            );
					intg_Ni_TEp_Nj += W[gp] * (      N[i][gp] * TEgp[gp][1][2] * N[j][gp] );
				}
	                
				StiffMatrix[i+numNodes][j+(2*numNodes)] = ( invMu * intg_Curl_Curl ) - ( freq2 * intg_Ni_TEp_Nj );
			}
		}
    }

    //*****************************************************************************************************************************
    //* - Generates the diagonal and the upper diagonal part of the elemental stiffness matrix (Gauged elements - Div(E) = 0.0)
    //*****************************************************************************************************************************	
    void VolumeElement_2ndQ_DivD_Plasma::Generate_UpperDiagonal_CurlDiv( Matrix<std::complex<double> >& StiffMatrix )
    {  
        // Loop indexes
        int i, j, gp;

        // Initializing elemental stiff matrix
        StiffMatrix.Resize( numDofs, numDofs );

        // Magnetic material properties
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mo;

        std::complex<double> cMu ( mu_real , mu_imag );

        // Gaussian weights and points
        Vector<double> cX; 
        Vector<double> cY;
        Vector<double> cZ;
        Vector<double>  W;

        int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

        // Jacobian
        Vector<double> detJ;

        Calculate_detJ( detJ, cX, cY, cZ );

        for ( gp=0; gp<numGaussPoints; gp++ ) W[gp] *= detJ[gp];

        // Basis functions on Gauss points 
        Matrix<double> N;

        Lagrange3D_Ni_2nd( N, cX, cY, cZ );

        // Cartesian derivatives of the basis functions on Gauss points
        Matrix<double> dNdx, dNdy, dNdz;

        LagrangeDerivatives3D_2ndQ( dNdx, dNdy, dNdz, cX, cY, cZ );

        // Extracting permitivitty tensor on nodes
        Vector<Matrix<std::complex<double> > > TEn;

        PermittivityTensor_OnNodes( TEn );

        // Permittivity tensor on Gauss points (interpolated with Ni)
        Vector<Matrix<std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, TEn, N );

        // Div( conj(e)*N ) on Gauss points
        Vector<Matrix<std::complex<double> > > Div_ce_N;

        // ( 1.0 / ( conj(e)*e ) ) * Div( e*E ) on Gauss points
        Vector<Matrix<std::complex<double> > > aDiv_e_E;

        DivD_a_DivD_OnGaussPoints( Div_ce_N, aDiv_e_E, TEn, TEgp, N, dNdx, dNdy, dNdz );

        // Temporal variables and constants
        std::complex<double> invMu( 1.00 / cMu  );
        std::complex<double> freq2( freq * freq );  

        //////////////////////////////////////// Block [Kxx] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=i; j<numNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 0; 
                int cj = 0;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Curl_Curl += W[gp] * ( dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
                    intg_DivD_DivD += W[gp] * (       Div_ce_N[gp][ci][i] * aDiv_e_E[gp][cj][j]       );
                    intg_Ni_TEp_Nj += W[gp] * (        N[i][gp] * TEgp[gp][ci][cj] * N[j][gp]         );
                }

                StiffMatrix[i][j] = ( invMu * intg_Curl_Curl ) + ( invMu * intg_DivD_DivD ) - ( freq2 * intg_Ni_TEp_Nj );
            }
        }
        //////////////////////////////////////// Block [Kyy] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=i; j<numNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 1; 
                int cj = 1;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Curl_Curl += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
                    intg_DivD_DivD += W[gp] * (       Div_ce_N[gp][ci][i] * aDiv_e_E[gp][cj][j]       );
                    intg_Ni_TEp_Nj += W[gp] * (         N[i][gp] * TEgp[gp][ci][cj] * N[j][gp]        );
                }

                StiffMatrix[i+numNodes][j+numNodes] = ( invMu * intg_Curl_Curl ) + ( invMu * intg_DivD_DivD ) - ( freq2 * intg_Ni_TEp_Nj );
            }
        }
        //////////////////////////////////////// Block [Kzz] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=i; j<numNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 2; 
                int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Curl_Curl += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] );
                    intg_DivD_DivD += W[gp] * (       Div_ce_N[gp][ci][i] * aDiv_e_E[gp][cj][j]       );
                    intg_Ni_TEp_Nj += W[gp] * (         N[i][gp] * TEgp[gp][ci][cj] * N[j][gp]        );
                }

                StiffMatrix[i+(2*numNodes)][j+(2*numNodes)] = ( invMu * intg_Curl_Curl ) + ( invMu * intg_DivD_DivD ) - ( freq2 * intg_Ni_TEp_Nj );
            }
        }
        //////////////////////////////////////// Block [Kxy] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 0; 
                int cj = 1;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Curl_Curl += W[gp] * (       - dNdy[i][gp] * dNdx[j][gp]         );
                    intg_DivD_DivD += W[gp] * ( Div_ce_N[gp][ci][i] * aDiv_e_E[gp][cj][j] );
                    intg_Ni_TEp_Nj += W[gp] * (   N[i][gp] * TEgp[gp][ci][cj] * N[j][gp]  );
                }

                StiffMatrix[i][j+numNodes] = ( invMu * intg_Curl_Curl ) + ( invMu * intg_DivD_DivD ) - ( freq2 * intg_Ni_TEp_Nj );
            }
        }
        /////////////////////////////////////// Block [Kxz] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 ); 

                int ci = 0; 
                int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Curl_Curl += W[gp] * (       - dNdz[i][gp] * dNdx[j][gp]         );
                    intg_DivD_DivD += W[gp] * ( Div_ce_N[gp][ci][i] * aDiv_e_E[gp][cj][j] );
                    intg_Ni_TEp_Nj += W[gp] * (   N[i][gp] * TEgp[gp][ci][cj] * N[j][gp]  );
                }

                StiffMatrix[i][j+(2*numNodes)] = ( invMu * intg_Curl_Curl ) + ( invMu * intg_DivD_DivD ) - ( freq2 * intg_Ni_TEp_Nj );
            }
        }
        ////////////////////////////////////// Block [Kyz] ////////////////////////////////////////////////
        for ( i=0; i<numNodes; i++ )
        {
            for ( j=0; j<numNodes; j++ )
            {
                std::complex<double> intg_Curl_Curl( 0.00, 0.00 ); 
                std::complex<double> intg_DivD_DivD( 0.00, 0.00 ); 
                std::complex<double> intg_Ni_TEp_Nj( 0.00, 0.00 );   
                
                int ci = 1; 
                int cj = 2;

                for ( gp=0; gp<numGaussPoints; gp++ )
                {
                    intg_Curl_Curl += W[gp] * (       - dNdz[i][gp] * dNdy[j][gp]         );
                    intg_DivD_DivD += W[gp] * ( Div_ce_N[gp][ci][i] * aDiv_e_E[gp][cj][j] );
                    intg_Ni_TEp_Nj += W[gp] * (  N[i][gp] * TEgp[gp][ci][cj] * N[j][gp]   );
                }

                StiffMatrix[i+numNodes][j+(2*numNodes)] = ( invMu * intg_Curl_Curl ) + ( invMu * intg_DivD_DivD ) - ( freq2 * intg_Ni_TEp_Nj );
            }
        }
    }

    //***************************************************************************************************************
	//* - Generates the lower diagonal part of the stiff matrix from the upper diagonal
	//***************************************************************************************************************
	void VolumeElement_2ndQ_DivD_Plasma::Generate_LowerDiagonal( Matrix<std::complex<double> >& StiffMatrix )
	{
		int i, j;

        if ( Is_IHL )
        {
            for ( i=0; i<numDofs; i++)
            {
                for ( j=i+1; j<numDofs; j++ )
                {
                    StiffMatrix[j][i] = StiffMatrix[i][j];
                }
            }
        }
        else
        {
		    for ( i=0; i<numDofs; i++)
		    {
			    for ( j=i+1; j<numDofs; j++ )
			    {
				    StiffMatrix[j][i] = std::conj( StiffMatrix[i][j] );
			    }
		    }
        }
	}

    //********************************************************************************************************************
    //* - Calculate derivatives of the permittivity tensor on Gauss points
    //********************************************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::PermittivityTensor_Derivatives( Vector<Matrix<std::complex<double> > >& dTEdx, 
                                                                         Vector<Matrix<std::complex<double> > >& dTEdy,
                                                                         Vector<Matrix<std::complex<double> > >& dTEdz,
                                                                         Vector<Matrix<std::complex<double> > >& TEn  ,
                                                                         Matrix<double>& dNdx, 
                                                                         Matrix<double>& dNdy, 
                                                                         Matrix<double>& dNdz )
    {
        // Number of Gauss points
        int numGaussPoints = dNdx[0].size();

        // Complex cero
        std::complex<double> cZero(0.0, 0.0);

        // Initializing vector of 3x3 permittivity tensors derivatives
        dTEdx.resize( numGaussPoints );
        dTEdy.resize( numGaussPoints );
        dTEdz.resize( numGaussPoints );

        // Computing derivatives: (dTE/dx)|gp = Sum_i{ TEni * (dNi/dx)|gp } ...
        for ( int gp=0; gp<numGaussPoints; gp++ )
        { 
            dTEdx[gp].Resize( 3, 3 );
            dTEdy[gp].Resize( 3, 3 );
            dTEdz[gp].Resize( 3, 3 );

            for ( int ki=0; ki<3; ki++ )
            {
                for ( int kj=0; kj<3; kj++ )
                {
                    dTEdx[gp][ki][kj] = cZero;
                    dTEdy[gp][ki][kj] = cZero;
                    dTEdz[gp][ki][kj] = cZero;

                    for ( int i=0; i<numNodes; i++ )
                    {
                        dTEdx[gp][ki][kj] += TEn[i][ki][kj] * dNdx[i][gp];
                        dTEdy[gp][ki][kj] += TEn[i][ki][kj] * dNdy[i][gp];
                        dTEdz[gp][ki][kj] += TEn[i][ki][kj] * dNdz[i][gp];
                    }
                }
            }
        }
    }

    //*****************************************************************************************************************
    //* - Calculation of the regularization terms Div_eNi and aDiv_eE
    //*****************************************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::DivD_a_DivD_OnGaussPoints( Vector<Matrix<std::complex<double> > >& Div_ce_N, 
                                                                    Vector<Matrix<std::complex<double> > >& aDiv_e_E, 
                                                                    Vector<Matrix<std::complex<double> > >& TEn ,
                                                                    Vector<Matrix<std::complex<double> > >& TEgp,  
                                                                    Matrix<double>& N   , 
                                                                    Matrix<double>& dNdx, 
                                                                    Matrix<double>& dNdy, 
                                                                    Matrix<double>& dNdz )
    {
        // Number of Gauss points
        int numGaussPoints = N[0].size();

        // Complex cero
        std::complex<double> cZero(0.0, 0.0);

        // Permittivity tensor derivatives
        Vector<Matrix<std::complex<double> > > dTEdx;
        Vector<Matrix<std::complex<double> > > dTEdy;
        Vector<Matrix<std::complex<double> > > dTEdz;

        PermittivityTensor_Derivatives( dTEdx, dTEdy, dTEdz, TEn, dNdx, dNdy, dNdz );

        // Computing: Div( conj([e])*Ni )
        Div_ce_N.resize( numGaussPoints );

        // Computing: invTau * Div( [e]*E )
        aDiv_e_E.resize( numGaussPoints );

        // Regularization constants
        std::complex<double> NormTraza2;
        std::complex<double> invNormTr2;

        for ( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            Div_ce_N[gp].Resize( 3, numNodes );
            aDiv_e_E[gp].Resize( 3, numNodes );

            NormTraza2 = ( TEgp[gp][0][0] * std::conj( TEgp[gp][0][0] ) + TEgp[gp][0][1] * std::conj( TEgp[gp][0][1] ) + TEgp[gp][0][2] * std::conj( TEgp[gp][0][2] ) + 
                           TEgp[gp][1][0] * std::conj( TEgp[gp][1][0] ) + TEgp[gp][1][1] * std::conj( TEgp[gp][1][1] ) + TEgp[gp][1][2] * std::conj( TEgp[gp][1][2] ) + 
                           TEgp[gp][2][0] * std::conj( TEgp[gp][2][0] ) + TEgp[gp][2][1] * std::conj( TEgp[gp][2][1] ) + TEgp[gp][2][2] * std::conj( TEgp[gp][2][2] ) );             
            
            invNormTr2 = 3.0 / ( NormTraza2 );
          
            for ( int c=0; c<3; c++ )
            {
                for ( int i=0; i<numNodes; i++ )
                {
                    std::complex<double> gpSum ( 0.00, 0.00 );

                    gpSum += TEgp[gp][0][c] * dNdx[i][gp] + N[i][gp] * dTEdx[gp][0][c];
                    gpSum += TEgp[gp][1][c] * dNdy[i][gp] + N[i][gp] * dTEdy[gp][1][c];
                    gpSum += TEgp[gp][2][c] * dNdz[i][gp] + N[i][gp] * dTEdz[gp][2][c];

                    Div_ce_N[gp][c][i] = std::conj( gpSum );

                    aDiv_e_E[gp][c][i] = invNormTr2 * gpSum;
                }
            }
        }
    }

    //************************************************************************************************* 
    //* -  Calculates the rotation matrix 
    //*************************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz  )
    {    
         // Unit vectors
         Vector<double> n(3), t(3), b(3);

         // b = B / |B|
         double Bnorm = std::sqrt( Bx*Bx + By*By + Bz*Bz );

         b[0] = Bx / Bnorm; 
         b[1] = By / Bnorm; 
         b[2] = Bz / Bnorm;

         // n _|_ b
         double nnorm = std::sqrt( b[0]*b[0] + b[1]*b[1] );

         if (nnorm > 0.00)
         {
             n[0] = b[1] / nnorm; 
             n[1] =-b[0] / nnorm; 
             n[2] = 0.00;
         }
         else
         {
             n[0] = b[2];
             n[1] = 0.00; 
             n[2] = 0.00;
         }

         // t = b x n
         t[0] = b[1]*n[2] - b[2]*n[1];
         t[1] = b[2]*n[0] - b[0]*n[2];
         t[2] = b[0]*n[1] - b[1]*n[0];

         // Rotation matrix 
         // T = [ nx tx bx
         //       ny ty by 
         //       nz tz bz ]
         T.Resize(3,3);

         for ( int row=0; row<3; row++ )
         {
             T[row][0] = std::complex<double>( n[row] ); 
             T[row][1] = std::complex<double>( t[row] ); 
             T[row][2] = std::complex<double>( b[row] );       
         }
    }

    //**********************************************************************************************************************
    //* - Calculates permittivity tensor on nodes 
    //**********************************************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::PermittivityTensor_OnNodes( Vector<Matrix<std::complex<double> > >& TEn )
    {
        // Permittivity tensor on nodes
        TEn.resize( numNodes );
        
        // Complex constants
        std::complex<double> cZero( 0.0, 0.0 );
        std::complex<double> cUnit( 0.0, 1.0 );

        if ( Is_IHL )
        {
            for (int i=0; i<numNodes; i++) 
            { 
                TEn[i].Resize(3,3);

                double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
                double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
                double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

                std::complex<double> cEps( eps_real, eps_imag + (sigma/freq) );

                // IHL permittivity tensor 
                Matrix<std::complex<double> > K(3,3);

                K[0][0] = cEps ; K[0][1] = cZero; K[0][2] = cZero; 
                K[1][0] = cZero; K[1][1] = cEps ; K[1][2] = cZero; 
                K[2][0] = cZero; K[2][1] = cZero; K[2][2] = cEps ;

                TEn[i] = K;
            }
        }
        else
        {
            for (int i=0; i<numNodes; i++) 
            { 
                TEn[i].Resize(3,3);

                // Plasma parameters on nodes
                Vector<double> Bext; 
                Vector<std::complex<double> > SDPRL;

                mpColdPlasma->Get_Bext_SDPRL_InNode( Bext, SDPRL, mNodes[i], freq );

                double Bx = Bext[0];
                double By = Bext[1];
                double Bz = Bext[2];

                std::complex<double> S = SDPRL[0];
                std::complex<double> D = SDPRL[1];
                std::complex<double> P = SDPRL[2];
               
                // Permittivity tensor in b = B/|B| coordinates
                Matrix<std::complex<double> > K(3,3);

                K[0][0] = S      ; K[0][1] = -cUnit*D; K[0][2] = cZero; 
                K[1][0] = cUnit*D; K[1][1] =  S      ; K[1][2] = cZero; 
                K[2][0] = cZero  ; K[2][1] =  cZero  ; K[2][2] = P    ;

                K *= eo;

                // Rotation matrix 
                Matrix<std::complex<double> > T; 

                Rotation_Matrix( T, Bx, By, Bz );

                // TEn = T * K * Tt 
                // Be carefull!!, T.Transpose() changes T permanently
                TEn[i] = T * K;
                TEn[i] = TEn[i] * T.Transpose();
            }
        }
    }

    //********************************************************************************************************************
    //* - Calculates permittivity tensor on Gauss points ( extrapolating TEn with N )
    //********************************************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::PermittivityTensor_OnGaussPoints( Vector<Matrix<std::complex<double> > >& TEgp,  
                                                                           Vector<Matrix<std::complex<double> > >& TEn, 
                                                                           Matrix<double>& N )
    {
        // Number of Gauss points
        int numGaussPoints = N[0].size();

        // Initializing vector of 3x3 permittivity tensors 
        TEgp.resize( numGaussPoints );

        // Complex constants
        std::complex<double> cZero( 0.0, 0.0 );

        // Permittivity tensor on Gauss points  
        for ( int gp=0; gp<numGaussPoints; gp++ )
        {
            TEgp[gp].Resize( 3, 3, cZero );

            for ( int ki=0; ki<3; ki++ )
            {
                for ( int kj=0; kj<3; kj++ )
                {
                    for ( int i=0; i<numNodes; i++ ) 
                    {
                        TEgp[gp][ki][kj] += TEn[i][ki][kj] * N[i][gp];
                    }// end ( int i=0; i<numNodes; i++ ) 
                }// end ( int kj=0; kj<3; kj++ )
            }// end ( int ki=0; ki<3; ki++ )
        }// end ( int gp=0; gp<numGaussPoints; gp++ )
    }

    //***************************************************************************
    //* - Calculate derivatives in nodes
    //***************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::DerivativesInNodes( Matrix<double>& ndNdx, 
												      		 Matrix<double>& ndNdy, 
														     Matrix<double>& ndNdz )
    {
		Vector<double> cX(numNodes); 
		Vector<double> cY(numNodes);
		Vector<double> cZ(numNodes);

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

		LagrangeDerivatives3D_2ndQ( ndNdx, ndNdy, ndNdz, cX, cY, cZ );
	}

    //*****************************************************************************************
    //* - Calculates H = ( 1.0 / j*w*mu ) * rot E
    //*****************************************************************************************
    void VolumeElement_2ndQ_DivD_Plasma::CalculateDerivatives( Vector<std::complex<double> >& ncHx,
		                                                       Vector<std::complex<double> >& ncHy,
												               Vector<std::complex<double> >& ncHz )
    {
        int i,n;

        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mo;

        std::complex<double> cMu     ( mu_real, mu_imag        );
        std::complex<double> cUnit   ( 0.00   , 1.00           );
        std::complex<double> inv_jwmu( 1.00 / (freq*cUnit*cMu) );

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

		std::complex<double> cZero(0.00,0.00);

        ncHx.resize(numNodes);
		ncHy.resize(numNodes);
		ncHz.resize(numNodes);

		for( n=0; n<numNodes; n++ )
		{
			ncHx[n] = cZero;
			ncHy[n] = cZero;
			ncHz[n] = cZero;

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