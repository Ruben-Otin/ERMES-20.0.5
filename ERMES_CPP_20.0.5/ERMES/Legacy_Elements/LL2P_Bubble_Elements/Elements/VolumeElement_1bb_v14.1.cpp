
#include "VolumeElement_1bb.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

#include <unordered_set>

namespace Kratos
{
    //**********************************************************************************************************************
    // - Global index of the nodes
    //**********************************************************************************************************************	
    void VolumeElement_1bb::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( mNumDofs );

        for( int i = 0; i < mNumNodes; ++i ) 
        {  
            EquationId[ i               ] = mNodes[ i ]->pDofcEx()->EquationId();
            EquationId[ i + mNumNodes   ] = mNodes[ i ]->pDofcEy()->EquationId();
			EquationId[ i + mNumNodes*2 ] = mNodes[ i ]->pDofcEz()->EquationId();
        }
    }

    //**********************************************************************************************************************
    // - Calculates the volume of the element
    //**********************************************************************************************************************
    double VolumeElement_1bb::Calculate_Volume()
    {
        double det;

        det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( det / 6.00 );
    }

    //**********************************************************************************************************************
    // - Calculates dN/dx, dN/dy, dN/dz - DN[ X,Y,Z ][ 0,1,2,3 ]
    //**********************************************************************************************************************
    void VolumeElement_1bb::Calculate_Ni_Derivatives( Matrix<double>& DN )
    {
        // [ Coord: X,Y,Z ][ Local_node: 0,1,2,3 ]
        DN.Resize( 3 , 4 );

        double cte = 1.00 / ( 6.00 * mVolume ); 

        // dN/dx
        DN[0][0] = cte * ( Y(3)*Z(2) + Y(4)*Z(3) + Y(2)*Z(4) - Y(3)*Z(4) - Y(2)*Z(3) - Y(4)*Z(2) );
        DN[0][1] = cte * ( Y(3)*Z(4) + Y(4)*Z(1) + Y(1)*Z(3) - Y(3)*Z(1) - Y(1)*Z(4) - Y(4)*Z(3) );
        DN[0][2] = cte * ( Y(2)*Z(1) + Y(1)*Z(4) + Y(4)*Z(2) - Y(2)*Z(4) - Y(4)*Z(1) - Y(1)*Z(2) );
        DN[0][3] = cte * ( Y(2)*Z(3) + Y(3)*Z(1) + Y(1)*Z(2) - Y(2)*Z(1) - Y(3)*Z(2) - Y(1)*Z(3) );
                         
        // dN/dy         
        DN[1][0] = cte * ( X(3)*Z(4) + X(4)*Z(2) + X(2)*Z(3) - X(3)*Z(2) - X(2)*Z(4) - X(4)*Z(3) ); 
        DN[1][1] = cte * ( X(3)*Z(1) + X(1)*Z(4) + X(4)*Z(3) - X(1)*Z(3) - X(3)*Z(4) - X(4)*Z(1) );
        DN[1][2] = cte * ( X(2)*Z(4) + X(4)*Z(1) + X(1)*Z(2) - X(2)*Z(1) - X(4)*Z(2) - X(1)*Z(4) );
        DN[1][3] = cte * ( X(1)*Z(3) + X(3)*Z(2) + X(2)*Z(1) - X(3)*Z(1) - X(1)*Z(2) - X(2)*Z(3) );
                         
        // dN/dz         
        DN[2][0] = cte * ( X(4)*Y(3) + X(2)*Y(4) + X(3)*Y(2) - X(2)*Y(3) - X(4)*Y(2) - X(3)*Y(4) );
        DN[2][1] = cte * ( X(3)*Y(4) + X(4)*Y(1) + X(1)*Y(3) - X(3)*Y(1) - X(4)*Y(3) - X(1)*Y(4) ); 
        DN[2][2] = cte * ( X(2)*Y(1) + X(4)*Y(2) + X(1)*Y(4) - X(4)*Y(1) - X(2)*Y(4) - X(1)*Y(2) ); 
        DN[2][3] = cte * ( X(2)*Y(3) + X(3)*Y(1) + X(1)*Y(2) - X(2)*Y(1) - X(3)*Y(2) - X(1)*Y(3) );    
    }

    //**********************************************************************************************************************
    // - Set singularities
    //**********************************************************************************************************************
    void VolumeElement_1bb::SetPeso( double Peso )
    {
        mPeso = Peso;
    }

	//**********************************************************************************************************************
	// - Distance to an edge
	//**********************************************************************************************************************
	double VolumeElement_1bb::RtoS( double alpha, double beta, double gamma )
    {
		double X1 = mNodes[0]->X(), Y1 = mNodes[0]->Y();
		double X2 = mNodes[1]->X(), Y2 = mNodes[1]->Y();
		double X3 = mNodes[2]->X(), Y3 = mNodes[2]->Y();
		double X4 = mNodes[3]->X(), Y4 = mNodes[3]->Y();
        
		double X = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + gamma * ( X4 - X1 ) + X1;
        double Y = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + gamma * ( Y4 - Y1 ) + Y1;

        double X0 = 20.00;
        double Y0 =  5.00;

        return sqrt( (X-X0)*(X-X0) + (Y-Y0)*(Y-Y0) );
    }

    //**********************************************************************************************************************
    // - Calculation of dN/dx, dN/dy, dN/dz on cX, cY, cZ points
    //**********************************************************************************************************************
    void VolumeElement_1bb::LagrangeDerivatives3D_1bb( Matrix<double>& dNdx,
		                                               Matrix<double>& dNdy,
													   Matrix<double>& dNdz,
													   std::vector<double>& cX, 
													   std::vector<double>& cY, 
													   std::vector<double>& cZ )
    {
		int nPoints = cX.size();

		dNdx.Resize( mNumNodes, nPoints );
		dNdy.Resize( mNumNodes, nPoints );
		dNdz.Resize( mNumNodes, nPoints );
		
		Matrix<double> DN; Calculate_Ni_Derivatives( DN );

		double dL1dx, dL2dx, dL3dx, dL4dx;
        double dL1dy, dL2dy, dL3dy, dL4dy;
		double dL1dz, dL2dz, dL3dz, dL4dz;

        dL1dx = DN[ 0 ][ 0 ], dL2dx = DN[ 0 ][ 1 ], dL3dx = DN[ 0 ][ 2 ], dL4dx = DN[ 0 ][ 3 ];
        dL1dy = DN[ 1 ][ 0 ], dL2dy = DN[ 1 ][ 1 ], dL3dy = DN[ 1 ][ 2 ], dL4dy = DN[ 1 ][ 3 ];
		dL1dz = DN[ 2 ][ 0 ], dL2dz = DN[ 2 ][ 1 ], dL3dz = DN[ 2 ][ 2 ], dL4dz = DN[ 2 ][ 3 ];

		// dNi/dx , dNi/dy, dNi/dz
        for( int p=0; p<nPoints; p++ ) 
        {
			double L2 = cX[ p ];
            double L3 = cY[ p ];
			double L4 = cZ[ p ];
            double L1 = 1.00 - L2 - L3 - L4;			
			
			dNdx[ 0 ][ p ] = dL1dx;
			dNdy[ 0 ][ p ] = dL1dy;
			dNdz[ 0 ][ p ] = dL1dz;
				       	 
			dNdx[ 1 ][ p ] = dL2dx;
			dNdy[ 1 ][ p ] = dL2dy;
			dNdz[ 1 ][ p ] = dL2dz;
				       	 
			dNdx[ 2 ][ p ] = dL3dx;
			dNdy[ 2 ][ p ] = dL3dy;
			dNdz[ 2 ][ p ] = dL3dz;
				       	 
			dNdx[ 3 ][ p ] = dL4dx;
			dNdy[ 3 ][ p ] = dL4dy;
			dNdz[ 3 ][ p ] = dL4dz;
				       	 
			dNdx[ 4 ][ p ] = 256.0 * ( L1*L2*L3*dL4dx + L1*L2*dL3dx*L4 + L1*dL2dx*L3*L4 + dL1dx*L2*L3*L4 );
			dNdy[ 4 ][ p ] = 256.0 * ( L1*L2*L3*dL4dy + L1*L2*dL3dy*L4 + L1*dL2dy*L3*L4 + dL1dy*L2*L3*L4 );
			dNdz[ 4 ][ p ] = 256.0 * ( L1*L2*L3*dL4dz + L1*L2*dL3dz*L4 + L1*dL2dz*L3*L4 + dL1dz*L2*L3*L4 );
        }
	}

    //**********************************************************************************************************************
    // - Calculate derivatives in nodes
    //**********************************************************************************************************************
    void VolumeElement_1bb::DerivativesInNodes( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz )
    {
		std::vector<double> cX( mNumNodes ); 
		std::vector<double>	cY( mNumNodes );
		std::vector<double> cZ( mNumNodes );

		//Points 0,1,2,3,4
		cX[0] = 0.00; cY[0] = 0.00; cZ[0] = 0.00;
		cX[1] = 1.00; cY[1] = 0.00; cZ[1] = 0.00;
		cX[2] = 0.00; cY[2] = 1.00; cZ[2] = 0.00;
		cX[3] = 0.00; cY[3] = 0.00; cZ[3] = 1.00;
        cX[4] = 0.25; cY[4] = 0.25; cZ[4] = 0.25;

		LagrangeDerivatives3D_1bb( ndNdx, ndNdy, ndNdz, cX, cY, cZ );
	}

    //**********************************************************************************************************************
    // - Calculates H = ( 1.0 / i*w*mu0 ) * rot(E)
    //**********************************************************************************************************************
    void VolumeElement_1bb::CalculateDerivatives( Vector<std::complex<double> >& ncHx,
		                                          Vector<std::complex<double> >& ncHy,
												  Vector<std::complex<double> >& ncHz )
    {
        int i,n;

		std::complex<double> cUnit( 0.00, 1.00 );
        std::complex<double> cZero( 0.00, 0.00 );
        
        double mu_real = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
		double mu_imag = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

		std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> inv_jwmu = 1.00 / ( mFreq * cUnit * cMu );

		Matrix<double> ndNdx; 
		Matrix<double> ndNdy; 
		Matrix<double> ndNdz;

		DerivativesInNodes( ndNdx, ndNdy, ndNdz );

		Vector<std::complex<double> > ncEx( mNumNodes );
		Vector<std::complex<double> > ncEy( mNumNodes );
		Vector<std::complex<double> > ncEz( mNumNodes );
		     
		for( n=0; n<mNumNodes; n++ )
		{
			ncEx[n] = (*mProperties)(cEx, *mNodes[n]);
			ncEy[n] = (*mProperties)(cEy, *mNodes[n]);
			ncEz[n] = (*mProperties)(cEz, *mNodes[n]);
		}

        ncHx.resize( mNumNodes );
		ncHy.resize( mNumNodes );
		ncHz.resize( mNumNodes );

		for( n=0; n<mNumNodes; n++ )
		{
			ncHx[n] = cZero;
			ncHy[n] = cZero;
			ncHz[n] = cZero;

			for( i=0; i<mNumNodes; i++ )
			{
			    ncHx[n] += ndNdy[i][n]*ncEz[i] - ndNdz[i][n]*ncEy[i];
				ncHy[n] += ndNdz[i][n]*ncEx[i] - ndNdx[i][n]*ncEz[i];
				ncHz[n] += ndNdx[i][n]*ncEy[i] - ndNdy[i][n]*ncEx[i];
			}						   
		}

		//ncHx[4] = 0.25 * ( ncHx[0] + ncHx[1] + ncHx[2] + ncHx[3] );
		//ncHy[4] = 0.25 * ( ncHy[0] + ncHy[1] + ncHy[2] + ncHy[3] );
		//ncHz[4] = 0.25 * ( ncHz[0] + ncHz[1] + ncHz[2] + ncHz[3] );
			
		ncHx *= inv_jwmu;
		ncHy *= inv_jwmu;
		ncHz *= inv_jwmu;      

        ///////////////////////////////////////////////////////////

		/*
		// Loop counters
		int i, j, n, gp;
		
		Matrix< std::complex<double> > AllE; AllE.Resize( mNumDofs , 1 );
		     
		for( n=0; n<mNumNodes; n++ )
		{
			AllE[ n               ][0] = (*mProperties)(cEx, *mNodes[n]);
			AllE[ n + mNumNodes   ][0] = (*mProperties)(cEy, *mNodes[n]);
			AllE[ n + mNumNodes*2 ][0] = (*mProperties)(cEz, *mNodes[n]);
		}

		std::complex<double> cUnit( 0.00, 1.00 );

		// Material properties
		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
        std::complex<double> cMu( mu_real , mu_imag );
		std::complex<double> inv_jwmu = 1.00 / ( mFreq * cUnit * cMu );
      
		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order8( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_1bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< std::complex<double> > C_cur; C_cur.Resize( 12, mNumDofs, mZero );

		for( i=0; i<4; i++ )
		{
		    for( j=0; j<mNumNodes; j++ )
		    {
		    	double vintg_Ni_dNjdX = 0.00;
		        double vintg_Ni_dNjdY = 0.00;
		        double vintg_Ni_dNjdZ = 0.00;
		    
                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	vintg_Ni_dNjdX += W[gp] * N[i][gp] * dNdx[j][gp]; 
		    		vintg_Ni_dNjdY += W[gp] * N[i][gp] * dNdy[j][gp];
		    		vintg_Ni_dNjdZ += W[gp] * N[i][gp] * dNdz[j][gp];
                }      
		    	
		    	C_cur[ i     ][ j               ] =            0.00; 
                C_cur[ i     ][ j + mNumNodes   ] = -vintg_Ni_dNjdZ; 
                C_cur[ i     ][ j + mNumNodes*2 ] = +vintg_Ni_dNjdY; 
		    								    
		    	C_cur[ i + 4 ][ j               ] = +vintg_Ni_dNjdZ; 
                C_cur[ i + 4 ][ j + mNumNodes   ] =            0.00; 
                C_cur[ i + 4 ][ j + mNumNodes*2 ] = -vintg_Ni_dNjdX; 
		    								    
		    	C_cur[ i + 8 ][ j               ] = -vintg_Ni_dNjdY; 
                C_cur[ i + 8 ][ j + mNumNodes   ] = +vintg_Ni_dNjdX; 
                C_cur[ i + 8 ][ j + mNumNodes*2 ] =            0.00; 
		    }
		}
		
		//Add_C_NxN_Matrix( C_cur ); 
		
		// Projection basis matrices
		Matrix< std::complex<double> > iA_div; iA_div.Resize(  4,  4, mZero );
		Matrix< std::complex<double> > iA_cur; iA_cur.Resize( 12, 12, mZero );
		
		for( i=0; i<4; i++ )
		{
            for( j=0; j<4; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -4.00 / mVolume;
        		else          iA_div[ i ][ j ] = +16.0 / mVolume;                          
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}

		Matrix< std::complex<double> > CurlProj; CurlProj.Resize( 12, 1, mZero ); 
		
		CurlProj += ( iA_cur * C_cur ) * AllE;

		CurlProj *= inv_jwmu;

		ncHx.resize( mNumNodes, mZero );
		ncHy.resize( mNumNodes, mZero );
		ncHz.resize( mNumNodes, mZero );

		for( i=0; i<4; i++ )
		{
			ncHx[i] = CurlProj[i  ][0];
			ncHy[i] = CurlProj[i+4][0];
			ncHz[i] = CurlProj[i+8][0];
		}

		ncHx[4] = 0.25 * ( ncHx[0] + ncHx[1] + ncHx[2] + ncHx[3] );
		ncHy[4] = 0.25 * ( ncHy[0] + ncHy[1] + ncHy[2] + ncHy[3] );
		ncHz[4] = 0.25 * ( ncHz[0] + ncHz[1] + ncHz[2] + ncHz[3] );
		*/
    }

    //**********************************************************************************************************************
    // - Calculates matrix: w2 * volume_intg( Ni * ep * Nj )
    //**********************************************************************************************************************
    void VolumeElement_1bb::Calculate_NiNj_matrix( Matrix< std::complex<double> >& w2_Ni_ep_Nj )
    {
        // Resizing matrix
		w2_Ni_ep_Nj.Resize( mNumNodes , mNumNodes );

		// Loop indexes
        int i , j;

		// Material properties values
		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY)       ;
        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
		double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;

		std::complex<double> cEps( eps_real, eps_imag + ( sigma / mFreq ) );
		std::complex<double> cCte = mFreq * mFreq * cEps;

        // Nodal base
        std::complex<double> cCteNiNj = mVolume * cCte / 20.00;
        std::complex<double> cCteNiNi = mVolume * cCte / 10.00;

        for( i=0; i<4; ++i )
        {
        	for( j=0; j<4; ++j )
        	{
        		if ( i != j ) w2_Ni_ep_Nj[ i ][ j ] = cCteNiNj;
        		else          w2_Ni_ep_Nj[ i ][ j ] = cCteNiNi;		         
        	}
        }

        // Bubble base
        std::complex<double> cCteNiBB = mVolume * cCte * (   256.0 /   3360.0 );
        std::complex<double> cCteBBBB = mVolume * cCte * ( 65536.0 / 415800.0 );
         
        for( i=0; i<4; ++i ) w2_Ni_ep_Nj[ i ][ 4 ] = cCteNiBB;
        for( j=0; j<4; ++j ) w2_Ni_ep_Nj[ 4 ][ j ] = cCteNiBB;

        w2_Ni_ep_Nj[ 4 ][ 4 ] = cCteBBBB;
    }

	//**********************************************************************************************************************
    // - Set connectivity objects
    //**********************************************************************************************************************
    void VolumeElement_1bb::Set_Mesh_Connectivities( std::map< unsigned int, Properties::Pointer >& PropertiesMap    , 
		                                             Vector  < Vector<unsigned int>              >& NodesConnectivity, 
                                                     Vector  < Vector<unsigned int>              >& ElementsInfo     ,
	                                                 std::map< unsigned int, char                >& mType_Of_BCNormal)
    {
        gpPropertiesMap     = &PropertiesMap    ;
        mpNodesConnectivity = &NodesConnectivity; 
        mpElementsInfo      = &ElementsInfo     ; 
		gpNormalType        = &mType_Of_BCNormal;
    }

	//**********************************************************************************************************************
    // - Set global derivatives 
    //**********************************************************************************************************************
    void VolumeElement_1bb::Set_Global_Derivatives( Vector<        double  >& Element_Volumes, 
                                                    Vector< Vector<double> >& Element_dNk_dXs, 
                                                    Vector< Vector<double> >& Element_dNk_dYs, 
		                                            Vector< Vector<double> >& Element_dNk_dZs )
	{
        gpVolumes = &Element_Volumes;
        gp_dNk_dX = &Element_dNk_dXs;
        gp_dNk_dY = &Element_dNk_dYs;
        gp_dNk_dZ = &Element_dNk_dZs; 	
	}

    //**********************************************************************************************************************
    // - Get the list of elements connected to This element 
    //**********************************************************************************************************************
    void VolumeElement_1bb::Get_Elements_List( Vector<unsigned int>& Elements_List )
	{
		std::unordered_set<unsigned int> ElementsSet;
		
        for( int i=0; i<4; ++i ) 
        {   
            ElementsSet.insert( (*mpNodesConnectivity)[ mNodes[i]->Id() - 1 ].begin(), 
				                (*mpNodesConnectivity)[ mNodes[i]->Id() - 1 ].end  () );
        }

		Elements_List.assign( ElementsSet.begin(), ElementsSet.end() );
	}

    //**********************************************************************************************************************
    // - Get the material properites of all the elements on the element list
    //**********************************************************************************************************************
    void VolumeElement_1bb::Get_Elements_Properties( Vector  < unsigned int >& Elements_List, 
                                                     std::map< unsigned int, std::complex<double> >& eEp,
		                                             std::map< unsigned int, std::complex<double> >& eMu )
	{
        Vector<unsigned int>::iterator el_it;

		for( el_it = Elements_List.begin(); el_it != Elements_List.end(); ++el_it )
	    {
			unsigned int ePropertiesId = (*mpElementsInfo)[ *el_it ][ 4 ];

			Properties::Pointer epProperties = (*gpPropertiesMap)[ ePropertiesId ];

			double sigma   = (*epProperties)(IHL_ELECTRIC_CONDUCTIVITY);
            
			double ep_real = (*epProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
		    double ep_imag = (*epProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
			
			double mu_real = (*epProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
		    double mu_imag = (*epProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

		    std::complex<double> cEp( ep_real, ep_imag + ( sigma / mFreq ) );
			std::complex<double> cMu( mu_real, mu_imag );

			eEp[ *el_it ] = cEp;
			eMu[ *el_it ] = cMu;
		}
	}

    //**********************************************************************************************************************
    // - Computing the alpha_kk and tbeta_kk terms
    //**********************************************************************************************************************
    void VolumeElement_1bb::Get_alpha_tbeta_kk( std::map< unsigned int, std::complex<double> >& alpha_kk,
										        std::map< unsigned int, std::complex<double> >& tbeta_kk )
	{
        // Get list of all the lements connected to This element
		Vector<unsigned int> Elements_List;
		
		Get_Elements_List( Elements_List );

		// Get material properites for all the elements on the element list
		std::map< unsigned int, std::complex<double> > eEp;
		std::map< unsigned int, std::complex<double> > eMu;

        Get_Elements_Properties( Elements_List, eEp, eMu );

		// List of all nodes included on the element list
	    Vector<unsigned int> k_Nodes_List;

		std::unordered_set<unsigned int> NodesSet;
		 
		Vector<unsigned int>::iterator el_it;

        for( el_it = Elements_List.begin(); el_it != Elements_List.end(); ++el_it ) 
		{
		    NodesSet.insert( (*mpElementsInfo)[*el_it].begin(), (*mpElementsInfo)[*el_it].end() - 1 );
		}
				
        k_Nodes_List.assign( NodesSet.begin(), NodesSet.end() );

		// Loop over the nodes on the list
		Vector<unsigned int>::iterator knl_it;

		// Computing alpha_kk and tbeta_kk
		for( knl_it = k_Nodes_List.begin(); knl_it != k_Nodes_List.end(); ++knl_it ) 
		{
			Vector<unsigned int> kElements = (*mpNodesConnectivity)[ *knl_it ];

			std::complex<double> inv_alpha( 0.00, 0.00 );
			std::complex<double> inv_tbeta( 0.00, 0.00 );

			Vector<unsigned int>::iterator ke_it;
		
			for( ke_it = kElements.begin(); ke_it != kElements.end(); ++ke_it ) 
            {
				inv_alpha += 0.25 * (*gpVolumes)[ *ke_it ] * eMu[ *ke_it ];
				inv_tbeta += 0.25 * (*gpVolumes)[ *ke_it ] * eMu[ *ke_it ] * eEp[ *ke_it ] * eEp[ *ke_it ];
			}

			alpha_kk[ *knl_it ] = 1.0 / inv_alpha;
			tbeta_kk[ *knl_it ] = 1.0 / inv_tbeta;
		}
	}

    //**********************************************************************************************************************
    // - Computing curl and div matrix elements for the L2 projection method
    //**********************************************************************************************************************
    void VolumeElement_1bb::Calculate_L2P_CurlDiv_Terms( Vector< std::map< unsigned int, double> >& Ci_dNkdX, 
		                                                 Vector< std::map< unsigned int, double> >& Ci_dNkdY,    
		                                                 Vector< std::map< unsigned int, double> >& Ci_dNkdZ )
	{ 
	    // Resizing (Ci,dNk/dX), (Ci,dNk/dY), (Ci,dNk/dZ) terms
		Ci_dNkdX.resize( mNumNodes );
		Ci_dNkdY.resize( mNumNodes );
		Ci_dNkdZ.resize( mNumNodes );

		// Computing ( Ci, dNk/d{X,Y,Z} ) on corner nodes
		for( int i=0; i<4; ++i ) 
		{
			// Elements connected to corner node
			Vector<unsigned int> cnElements = (*mpNodesConnectivity)[ mNodes[i]->Id() - 1 ];

			// Iteration over elements connected to corner node 
			Vector<unsigned int>::iterator ce_it;

			for( ce_it = cnElements.begin(); ce_it != cnElements.end(); ++ce_it )
			{
				Vector<unsigned int> eInfo = (*mpElementsInfo)[ *ce_it ];

				double volInt_Ni = 0.25 * (*gpVolumes)[ *ce_it ];

				// Loop over the nodes of (*ce_it) element
				for( int k=0; k<4; ++k ) 
				{
					unsigned int kId = eInfo[k];

				    double volInt_Ni_dNkdX = volInt_Ni * (*gp_dNk_dX)[*ce_it][k];
					double volInt_Ni_dNkdY = volInt_Ni * (*gp_dNk_dY)[*ce_it][k];
					double volInt_Ni_dNkdZ = volInt_Ni * (*gp_dNk_dZ)[*ce_it][k];

					if( Ci_dNkdX[i].find(kId) == Ci_dNkdX[i].end() )
					{
                        Ci_dNkdX[i][kId] = volInt_Ni_dNkdX;
						Ci_dNkdY[i][kId] = volInt_Ni_dNkdY;
						Ci_dNkdZ[i][kId] = volInt_Ni_dNkdZ;
				    }
					else
					{
						Ci_dNkdX[i][kId] += volInt_Ni_dNkdX;
						Ci_dNkdY[i][kId] += volInt_Ni_dNkdY;
						Ci_dNkdZ[i][kId] += volInt_Ni_dNkdZ;			
					}

			    }// End nodes in connected element loop - for( int k=0; k<4; ++k ) 

			}//End elements connected to corner node loop - for( ce_it = cnElements.begin(); ce_it != cnElements.end(); ++ce_it )

		}// End corner nodes loop - for( int i=0; i<4; ++i ) 

		// Computing ( Ci, dNk/d{X,Y,Z} ) on bubble node
        Matrix<double> DN; 
		
		Calculate_Ni_Derivatives( DN );

		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY)       ;
		double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
		double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;

		std::complex<double> cEps( eps_real, eps_imag + ( sigma / mFreq ) );

        double volInt_Bi = ( 256.0 / 840.0 ) * mVolume;

	    for( int k=0; k<4; ++k ) 
	    {
	    	unsigned int kId = mNodes[k]->Id() - 1;
	    
	        double volInt_Bi_dNkdX = volInt_Bi * DN[0][k];
	    	double volInt_Bi_dNkdY = volInt_Bi * DN[1][k];
	    	double volInt_Bi_dNkdZ = volInt_Bi * DN[2][k];

			Ci_dNkdX[4][kId] = volInt_Bi_dNkdX;
	    	Ci_dNkdY[4][kId] = volInt_Bi_dNkdY;
	    	Ci_dNkdZ[4][kId] = volInt_Bi_dNkdZ;
	    }
	}     

    //**********************************************************************************************************************
    // - Computing curl and ep*div matrix elements for the L2 projection method
    //**********************************************************************************************************************
    void VolumeElement_1bb::Calculate_L2P_Curl_ep_Div_Terms( Vector< std::map< unsigned int, std::complex<double> > >& Ci_ep_dNkdX, 
		                                                     Vector< std::map< unsigned int, std::complex<double> > >& Ci_ep_dNkdY,    
		                                                     Vector< std::map< unsigned int, std::complex<double> > >& Ci_ep_dNkdZ )
	{ 
		// Resizing (Ci,ep*dNk/dX), (Ci,ep*dNk/dY), (Ci,ep*dNk/dZ) terms
		Ci_ep_dNkdX.resize( mNumNodes );
		Ci_ep_dNkdY.resize( mNumNodes );
		Ci_ep_dNkdZ.resize( mNumNodes );

		// Computing ( Ci, ep * dNk/d{X,Y,Z} ) on corner nodes
		for( int i=0; i<4; ++i ) 
		{
			// Elements connected to corner node
			Vector<unsigned int> cnElements = (*mpNodesConnectivity)[ mNodes[i]->Id() - 1 ];

			// Iteration over elements connected to corner node 
			Vector<unsigned int>::iterator ce_it;

			for( ce_it = cnElements.begin(); ce_it != cnElements.end(); ++ce_it )
			{
				Vector<unsigned int> eInfo = (*mpElementsInfo)[ *ce_it ];

				unsigned int ePropertiesId = (*mpElementsInfo)[ *ce_it ][ 4 ];

			    Properties::Pointer epProperties = (*gpPropertiesMap)[ ePropertiesId ];

			    double sigma   = (*epProperties)(IHL_ELECTRIC_CONDUCTIVITY);
			    double ep_real = (*epProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
		        double ep_imag = (*epProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
			    
		        std::complex<double> cEp( ep_real, ep_imag + ( sigma / mFreq ) );

				double volInt_Ni = 0.25 * (*gpVolumes)[ *ce_it ];

				// Loop over the nodes of (*ce_it) element
				for( int k=0; k<4; ++k ) 
				{
					unsigned int kId = eInfo[k];

				    double volInt_Ni_dNkdX = volInt_Ni * (*gp_dNk_dX)[*ce_it][k];
					double volInt_Ni_dNkdY = volInt_Ni * (*gp_dNk_dY)[*ce_it][k];
					double volInt_Ni_dNkdZ = volInt_Ni * (*gp_dNk_dZ)[*ce_it][k];

					if( Ci_ep_dNkdX[i].find(kId) == Ci_ep_dNkdX[i].end() )
					{
                        Ci_ep_dNkdX[i][kId] = cEp * volInt_Ni_dNkdX;
						Ci_ep_dNkdY[i][kId] = cEp * volInt_Ni_dNkdY;
						Ci_ep_dNkdZ[i][kId] = cEp * volInt_Ni_dNkdZ;
				    }
					else
					{
						Ci_ep_dNkdX[i][kId] += cEp * volInt_Ni_dNkdX;
						Ci_ep_dNkdY[i][kId] += cEp * volInt_Ni_dNkdY;
						Ci_ep_dNkdZ[i][kId] += cEp * volInt_Ni_dNkdZ;					
					}

			    }// End nodes in connected element loop - for( int k=0; k<4; ++k ) 

			}//End elements connected to corner node loop - for( ce_it = cnElements.begin(); ce_it != cnElements.end(); ++ce_it )

		}// End corner nodes loop - for( int i=0; i<4; ++i ) 

		// Computing ( Ci, ep * dNk/d{X,Y,Z} ) on bubble node
        Matrix<double> DN; 
		
		Calculate_Ni_Derivatives( DN );

		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY)       ;
		double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
		double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;

		std::complex<double> cEps( eps_real, eps_imag + ( sigma / mFreq ) );

        double volInt_Bi = ( 256.0 / 840.0 ) * mVolume;

	    for( int k=0; k<4; ++k ) 
	    {
	    	unsigned int kId = mNodes[k]->Id() - 1;
	    
	        double volInt_Bi_dNkdX = volInt_Bi * DN[0][k];
	    	double volInt_Bi_dNkdY = volInt_Bi * DN[1][k];
	    	double volInt_Bi_dNkdZ = volInt_Bi * DN[2][k];

			Ci_ep_dNkdX[4][kId] = cEps * volInt_Bi_dNkdX;
	    	Ci_ep_dNkdY[4][kId] = cEps * volInt_Bi_dNkdY;
	    	Ci_ep_dNkdZ[4][kId] = cEps * volInt_Bi_dNkdZ;
	    }
	}    

    //**********************************************************************************************************************
    // - Calculates the number of neighbour elements for the pair of nodes i, j
    //**********************************************************************************************************************
    int VolumeElement_1bb::NumberOfNeigbours( int i, int j )
	{ 
		// If bubble node the return 1
		if ( i == 4 || j == 4 )
		{
		    return 1;
		}
		
		// If corner node then return its connections
		if ( i == j )
		{
		    return (*mpNodesConnectivity)[ mNodes[i]->Id() - 1 ].size();
		}

	    // If nodes i,j on element edge then count the number of common elements 
		Vector<unsigned int> elements_i = (*mpNodesConnectivity)[ mNodes[i]->Id() - 1 ];
		Vector<unsigned int> elements_j = (*mpNodesConnectivity)[ mNodes[j]->Id() - 1 ];
        
		Vector<unsigned int> vIntersection; 

		std::sort( elements_i.begin(), elements_i.end() ); 
        std::sort( elements_j.begin(), elements_j.end() ); 

        std::set_intersection( elements_i.begin(), elements_i.end(), 
                               elements_j.begin(), elements_j.end(), 
                               std::back_inserter( vIntersection ) ); 

		return vIntersection.size();
	}

	//**********************************************************************************************************************
    // - Gets the list of common nodes for the pair of nodes i, j
    //**********************************************************************************************************************
    void VolumeElement_1bb::Get_k_Nodes_List( int i, int j, Vector<unsigned int>& k_Nodes_List )
	{ 
		// Clear vector 
		k_Nodes_List.clear();
		
		// If bubble node then return the list of corner nodes of This element
		if ( i == 4 || j == 4 )
		{
		    k_Nodes_List.push_back( mNodes[0]->Id() - 1 );
			k_Nodes_List.push_back( mNodes[1]->Id() - 1 );
			k_Nodes_List.push_back( mNodes[2]->Id() - 1 );
			k_Nodes_List.push_back( mNodes[3]->Id() - 1 );

			return;
		}

		// If corner node then return the list of nodes connected to the i node 
		if ( i == j )
		{
            Vector<unsigned int> elements_i = (*mpNodesConnectivity)[ mNodes[i]->Id() - 1 ];

			std::unordered_set<unsigned int> NodesSet;

			Vector<unsigned int>::iterator e_it;
            
            for( e_it = elements_i.begin(); e_it != elements_i.end(); ++e_it ) 
			{
				NodesSet.insert( (*mpElementsInfo)[*e_it].begin(), (*mpElementsInfo)[*e_it].end() - 1 );
			}
				
            k_Nodes_List.assign( NodesSet.begin(), NodesSet.end() );

		    return;
		}

		// If nodes i,j on element edge then get list of common nodes
		Vector<unsigned int> elements_i = (*mpNodesConnectivity)[ mNodes[i]->Id() - 1 ];
		Vector<unsigned int> elements_j = (*mpNodesConnectivity)[ mNodes[j]->Id() - 1 ];
        
		Vector<unsigned int> vIntersection; 

		std::sort( elements_i.begin(), elements_i.end() ); 
        std::sort( elements_j.begin(), elements_j.end() ); 

        std::set_intersection( elements_i.begin(), elements_i.end(), 
                               elements_j.begin(), elements_j.end(), 
                               std::back_inserter( vIntersection ) ); 

		 std::unordered_set<unsigned int> NodesSet;
		 
		 Vector<unsigned int>::iterator e_it;

         for( e_it = vIntersection.begin(); e_it != vIntersection.end(); ++e_it ) 
		 {
		     NodesSet.insert( (*mpElementsInfo)[*e_it].begin(), (*mpElementsInfo)[*e_it].end() - 1 );
		 }
				
         k_Nodes_List.assign( NodesSet.begin(), NodesSet.end() );

		 return;
	}

	//**********************************************************************************************************************
    // - Calculates the stiffness matrix with the selected formulation
    //**********************************************************************************************************************
    void VolumeElement_1bb::GetStiffMatrix( Matrix< std::complex<double> >& StiffMatrix )
	{   
		//StiffMatrix.Resize( mNumDofs, mNumDofs, mZero );

		//GetStiffMatrix_Plain_Regular( StiffMatrix );

		//GetStiffMatrix_PLL2P_SimpTau( StiffMatrix );
		//GetStiffMatrix_PLL2P_FullTau( StiffMatrix );

        //GetStiffMatrix_LL2P_VolumAvg( StiffMatrix ); 
		//GetStiffMatrix_LL2P_CurlDivE( StiffMatrix ); 
        //GetStiffMatrix_LL2P_CurlDivQ( StiffMatrix ); 
		GetStiffMatrix_LL2P_EdgeElem( StiffMatrix );

		//Add_hk2CurlDiv_Stabilization( StiffMatrix ); 
		//Add_CurlDivGBk_Stabilization( StiffMatrix );
	}

    //**********************************************************************************************************************
    // - Calculates the stiffness matrix with the L2 projection method using a simplified tau
    //**********************************************************************************************************************
    void VolumeElement_1bb::GetStiffMatrix_PLL2P_SimpTau( Matrix< std::complex<double> >& StiffMatrix )
	{          
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( mNumDofs, mNumDofs );

		// Computing matrix [ w2 * volume_intg( Ni * ep * Nj ) ]
		Matrix< std::complex<double> > w2_Ni_ep_Nj; 

		Calculate_NiNj_matrix( w2_Ni_ep_Nj );

		// Computing the alpha_kk and tbeta_kk terms
        std::map< unsigned int, std::complex<double> > alpha_kk;
		std::map< unsigned int, std::complex<double> > tbeta_kk;

		Get_alpha_tbeta_kk( alpha_kk, tbeta_kk );

		// Computing curl and div matrix elements for the L2 projection method
		Vector< std::map< unsigned int, double> > Ci_dNkdX; 
		Vector< std::map< unsigned int, double> > Ci_dNkdY; 
		Vector< std::map< unsigned int, double> > Ci_dNkdZ; 
        
		Calculate_L2P_CurlDiv_Terms( Ci_dNkdX, Ci_dNkdY, Ci_dNkdZ );
        
		// Building local stiffness matrix
		std::complex<double> Kij;
		std::complex<double> cZero( 0.00, 0.00 );

		Vector<unsigned int> k_Nodes_List; 
		Vector<unsigned int>::iterator kit;

		double NeighbourCoef;

		for(int i=0; i<mNumNodes; i++)
		{
            for(int j=0; j<mNumNodes; j++) 
            {
            	NeighbourCoef = 1.0 / NumberOfNeigbours( i, j );

				Get_k_Nodes_List( i, j, k_Nodes_List );

				/////////////  Bloques: [Kxx], [Kyy], [Kzz]  //////////////////////////
            	Kij = cZero;

				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij += alpha_kk[*kit] * ( Ci_dNkdX[i][*kit] * Ci_dNkdX[j][*kit] + Ci_dNkdY[i][*kit] * Ci_dNkdY[j][*kit] + Ci_dNkdZ[i][*kit] * Ci_dNkdZ[j][*kit] );
			    }

				Kij *= NeighbourCoef;

				Kij -= w2_Ni_ep_Nj[i][j];

            	StiffMatrix[ i               ][ j               ] = Kij;
            	StiffMatrix[ i + mNumNodes   ][ j + mNumNodes   ] = Kij;
            	StiffMatrix[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = Kij;
            
            	/////////////  Bloques: [Kxy], [Kyx]  //////////////////////////
            	Kij = cZero;

				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij += alpha_kk[*kit] * ( Ci_dNkdX[i][*kit] * Ci_dNkdY[j][*kit] - Ci_dNkdY[i][*kit] * Ci_dNkdX[j][*kit] );
			    }

				Kij *= NeighbourCoef;

				StiffMatrix[ i             ][ j + mNumNodes ] =  Kij;
            	StiffMatrix[ i + mNumNodes ][ j             ] = -Kij;
            
            	/////////////  Bloques: [Kxz] , [Kzx]  //////////////////////////
            	Kij = cZero;

				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij += alpha_kk[*kit] * ( Ci_dNkdX[i][*kit] * Ci_dNkdZ[j][*kit] - Ci_dNkdZ[i][*kit] * Ci_dNkdX[j][*kit] );
			    }

				Kij *= NeighbourCoef;

				StiffMatrix[ i               ][ j + mNumNodes*2 ] =  Kij;
            	StiffMatrix[ i + mNumNodes*2 ][ j               ] = -Kij;
            
            	/////////////  Bloques: [Kyz], [Kzy]  //////////////////////////
            	Kij = cZero;

				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij += alpha_kk[*kit] * ( Ci_dNkdY[i][*kit] * Ci_dNkdZ[j][*kit] - Ci_dNkdZ[i][*kit] * Ci_dNkdY[j][*kit] );
			    }

				Kij *= NeighbourCoef;

				StiffMatrix[ i + mNumNodes   ][ j + mNumNodes*2 ] =  Kij;
            	StiffMatrix[ i + mNumNodes*2 ][ j + mNumNodes   ] = -Kij;
            }
		}
	}

    //****************************************************************************************************
    // - Calculates the stiffness matrix with the L2 projection method using a full tau
    //****************************************************************************************************
    void VolumeElement_1bb::GetStiffMatrix_PLL2P_FullTau( Matrix< std::complex<double> >& StiffMatrix )
	{          
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( mNumDofs, mNumDofs );

		// Computing matrix [ w2 * volume_intg( Ni * ep * Nj ) ]
		Matrix< std::complex<double> > w2_Ni_ep_Nj; 

		Calculate_NiNj_matrix( w2_Ni_ep_Nj );

		// Computing the alpha_kk and tbeta_kk terms
        std::map< unsigned int, std::complex<double> > alpha_kk;
		std::map< unsigned int, std::complex<double> > tbeta_kk;

		Get_alpha_tbeta_kk( alpha_kk, tbeta_kk );

		// Computing curl and div matrix elements for the L2 projection method
		Vector< std::map< unsigned int, double> > Ci_dNkdX; 
		Vector< std::map< unsigned int, double> > Ci_dNkdY; 
		Vector< std::map< unsigned int, double> > Ci_dNkdZ; 
        
		Calculate_L2P_CurlDiv_Terms( Ci_dNkdX, Ci_dNkdY, Ci_dNkdZ );

		// Computing curl and ep*div matrix elements for the L2 projection method
        Vector< std::map< unsigned int, std::complex<double> > > Ci_ep_dNkdX; 
		Vector< std::map< unsigned int, std::complex<double> > > Ci_ep_dNkdY; 
		Vector< std::map< unsigned int, std::complex<double> > > Ci_ep_dNkdZ; 

		Calculate_L2P_Curl_ep_Div_Terms( Ci_ep_dNkdX, Ci_ep_dNkdY, Ci_ep_dNkdZ );

		// Building local stiffness matrix
		std::complex<double> Kij;
		std::complex<double> cZero( 0.00, 0.00 );

		Vector<unsigned int> k_Nodes_List; 
		Vector<unsigned int>::iterator kit;

		double NeighbourCoef;

		// Diagonal matrices
	    for( int i=0; i<mNumNodes; i++ )
	    {
	        for( int j=i; j<mNumNodes; j++ ) 
	    	{
                NeighbourCoef = 1.0 / NumberOfNeigbours( i, j );

				Get_k_Nodes_List( i, j, k_Nodes_List );
				
				// [Kxx]
				Kij = cZero;

				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij += alpha_kk[*kit] * ( Ci_dNkdY[i][*kit] * Ci_dNkdY[j][*kit] +    Ci_dNkdZ[i][*kit] * Ci_dNkdZ[j][*kit] );
				    Kij += tbeta_kk[*kit] * (                  Ci_ep_dNkdX[i][*kit] * Ci_ep_dNkdX[j][*kit]                     );
			    }

				Kij *= NeighbourCoef;

	    		Kij -= w2_Ni_ep_Nj[i][j];
 	    
	    		StiffMatrix[ i ][ j ] = Kij;

				// [Kyy]  
				Kij = cZero;

				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij += alpha_kk[*kit] * ( Ci_dNkdX[i][*kit] * Ci_dNkdX[j][*kit] +    Ci_dNkdZ[i][*kit] * Ci_dNkdZ[j][*kit] );
				    Kij += tbeta_kk[*kit] * (                  Ci_ep_dNkdY[i][*kit] * Ci_ep_dNkdY[j][*kit]                     );
			    }

				Kij *= NeighbourCoef;

	    		Kij -= w2_Ni_ep_Nj[i][j];
 	    
	    		StiffMatrix[ i + mNumNodes ][ j + mNumNodes ] = Kij;

				// [Kzz]  
				Kij = cZero;

				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij += alpha_kk[*kit] * ( Ci_dNkdX[i][*kit] * Ci_dNkdX[j][*kit] +    Ci_dNkdY[i][*kit] * Ci_dNkdY[j][*kit] );
				    Kij += tbeta_kk[*kit] * (                  Ci_ep_dNkdZ[i][*kit] * Ci_ep_dNkdZ[j][*kit]                     );
			    }

				Kij *= NeighbourCoef;

	    		Kij -= w2_Ni_ep_Nj[i][j];
 	    
	    		StiffMatrix[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = Kij;
			}
		}

		// Off-diagonal matrices
	    for( int i=0; i<mNumNodes; i++ )
	    {
	        for( int j=0; j<mNumNodes; j++ ) 
	    	{
				NeighbourCoef = 1.0 / NumberOfNeigbours( i, j );

				Get_k_Nodes_List( i, j, k_Nodes_List );

	    		// [Kxy] 
	    		Kij = cZero;
				
				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij -= alpha_kk[*kit] * (    Ci_dNkdY[i][*kit] *    Ci_dNkdX[j][*kit] );
				    Kij += tbeta_kk[*kit] * ( Ci_ep_dNkdX[i][*kit] * Ci_ep_dNkdY[j][*kit] );
			    }

				Kij *= NeighbourCoef;
             
	    		StiffMatrix[ i ][ j + mNumNodes ] = Kij;

				// [Kxz] 
	    		Kij = cZero;
				
				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij -= alpha_kk[*kit] * (    Ci_dNkdZ[i][*kit] *    Ci_dNkdX[j][*kit] );
				    Kij += tbeta_kk[*kit] * ( Ci_ep_dNkdX[i][*kit] * Ci_ep_dNkdZ[j][*kit] );
			    }

				Kij *= NeighbourCoef;
             
	    		StiffMatrix[ i ][ j + mNumNodes*2 ] = Kij;

				// [Kyz] 
	    		Kij = cZero;
				
				for( kit = k_Nodes_List.begin(); kit != k_Nodes_List.end(); ++kit )
				{
				    Kij -= alpha_kk[*kit] * (    Ci_dNkdZ[i][*kit] *    Ci_dNkdY[j][*kit] );
				    Kij += tbeta_kk[*kit] * ( Ci_ep_dNkdY[i][*kit] * Ci_ep_dNkdZ[j][*kit] );
			    }

				Kij *= NeighbourCoef;
             
	    		StiffMatrix[ i + mNumNodes ][ j + mNumNodes*2 ] = Kij;
	    	}
	    }

		// Lower diagonal matrix 
        for( int i=0; i<mNumDofs; i++ )
        {
            for( int j=i+1; j<mNumDofs; j++ )
            {
                StiffMatrix[j][i] = StiffMatrix[i][j];
            }
        }      
	}

    //******************************************************************************************************
    //* - Calculates the local stiffness matrix with a plain regularization
    //******************************************************************************************************	
    void VolumeElement_1bb::GetStiffMatrix_Plain_Regular( Matrix< std::complex<double> >& StiffMatrix ) 
    {  
        // Resizing local stiffness matrix 
        StiffMatrix.Resize( mNumDofs, mNumDofs );

		// Loop counters
		int i, j, gp;

		// Material properties
		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
        std::complex<double> cMu( mu_real , mu_imag );
      
		// Gaussian derivation weights
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions and derivatives
		Matrix<double> N;

		Lagrange3D_Ni_1bb( N, cX, cY, cZ );

		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_1bb( dNdx, dNdy, dNdz, cX, cY, cZ );
		
		if ( mPeso != 0.00 ) 
		{
			double intCurlDiv, intNN;

			std::complex<double> eKij;

			std::complex<double> cteCurlDiv = 1.00 / cMu;
			std::complex<double> cteNN      = mFreq * mFreq * cEp;  

			// [Kxx],[Kyy],[Kzz] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=i; j<mNumNodes; j++ )
				{
					intCurlDiv = 0.00; 
					intNN      = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
						intNN      += W[gp] * (    N[i][gp] *    N[j][gp] );
					}

					eKij = (cteCurlDiv * intCurlDiv) - (cteNN * intNN);

					StiffMatrix[ i               ][ j               ] = eKij;
					StiffMatrix[ i + mNumNodes   ][ j + mNumNodes   ] = eKij;
					StiffMatrix[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = eKij;
				}
			}

			// [Kxy] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=0; j<mNumNodes; j++ )
				{
					intCurlDiv = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
					    intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdy[j][gp] - dNdy[i][gp] * dNdx[j][gp] );
					}
	                
					StiffMatrix[ i ][ j + mNumNodes ] = (cteCurlDiv * intCurlDiv);
				}
			}

			// [Kxz] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=0; j<mNumNodes; j++ )
				{
					intCurlDiv = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdx[j][gp] );
					}
	                
					StiffMatrix[ i ][ j + mNumNodes*2 ] = (cteCurlDiv * intCurlDiv);
				}
			}

			// [Kyz] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=0; j<mNumNodes; j++ )
				{
					intCurlDiv = 0.00; 

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurlDiv += W[gp] * ( dNdy[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdy[j][gp] );
					}
	                
					StiffMatrix[ i + mNumNodes ][ j + mNumNodes*2 ] = (cteCurlDiv * intCurlDiv);
				}
			}
		}
		else
		{
			double intCurl, intNN;

			std::complex<double> cteCurl = 1.00 / cMu;
			std::complex<double> cteNN   = mFreq * mFreq * cEp;  

			// [Kxx]
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=i; j<mNumNodes; j++ )
				{
					intCurl = 0.00;
					intNN   = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
						intNN   += W[gp] * (    N[i][gp] *    N[j][gp] );
					}
	                
					StiffMatrix[ i ][ j ] = (cteCurl * intCurl) - (cteNN * intNN);
				}
			}

			// [Kyy] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=i; j<mNumNodes; j++ )
				{
					intCurl = 0.00;
					intNN   = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
						intNN   += W[gp] * (    N[i][gp] *    N[j][gp] );
					}
	                
					StiffMatrix[ i + mNumNodes ][ j + mNumNodes ] = (cteCurl * intCurl) - (cteNN * intNN);
				}
			}

			// [Kzz] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=i; j<mNumNodes; j++ )
				{
					intCurl = 0.00; 
					intNN   = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] );
						intNN   += W[gp] * (    N[i][gp] *    N[j][gp] );
					}
	                
					StiffMatrix[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = (cteCurl * intCurl) - (cteNN * intNN);
				}
			}

			// [Kxy] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=0; j<mNumNodes; j++ )
				{
					intCurl = 0.00; 

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdy[i][gp] * dNdx[j][gp] );
					}
	                
					StiffMatrix[ i ][ j + mNumNodes ] = - (cteCurl * intCurl);
				}
			}

			// [Kxz] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=0; j<mNumNodes; j++ )
				{
					intCurl = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdz[i][gp] * dNdx[j][gp] );
					}
	                
					StiffMatrix[ i ][ j + mNumNodes*2 ] = - (cteCurl * intCurl);
				}
			}

			// [Kyz] 
			for ( i=0; i<mNumNodes; i++ )
			{
				for ( j=0; j<mNumNodes; j++ )
				{
					intCurl = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurl += W[gp] * ( dNdz[i][gp] * dNdy[j][gp] );
					}
	                
					StiffMatrix[ i + mNumNodes ][ j + mNumNodes*2 ] = - (cteCurl * intCurl);
				}
			}
		}

        // Symmetric block
        for ( i=0; i<mNumDofs; i++ )
        {
            for ( j=i+1; j<mNumDofs; j++ )
            {
                StiffMatrix[ j ][ i ] =  StiffMatrix[ i ][ j ];
            }
        }      
    }

	//*************************************************************************************************
    //* - Calculates element diameter h_k
    //*************************************************************************************************
	double VolumeElement_1bb::Calculate_Element_Diameter() 
	{
		// Nodal coordinates
		double X1 = mNodes[0]->X(), Y1 = mNodes[0]->Y(), Z1 = mNodes[0]->Z();
		double X2 = mNodes[1]->X(), Y2 = mNodes[1]->Y(), Z2 = mNodes[1]->Z();
		double X3 = mNodes[2]->X(), Y3 = mNodes[2]->Y(), Z3 = mNodes[2]->Z();
		double X4 = mNodes[3]->X(), Y4 = mNodes[3]->Y(), Z4 = mNodes[3]->Z();

		// Element edge lenghts 
		double Lenght_12 = sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) + (Z1-Z2)*(Z1-Z2) );
		double Lenght_13 = sqrt( (X1-X3)*(X1-X3) + (Y1-Y3)*(Y1-Y3) + (Z1-Z3)*(Z1-Z3) );
		double Lenght_14 = sqrt( (X1-X4)*(X1-X4) + (Y1-Y4)*(Y1-Y4) + (Z1-Z4)*(Z1-Z4) );
		double Lenght_23 = sqrt( (X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3) + (Z2-Z3)*(Z2-Z3) );
		double Lenght_24 = sqrt( (X2-X4)*(X2-X4) + (Y2-Y4)*(Y2-Y4) + (Z2-Z4)*(Z2-Z4) );
		double Lenght_34 = sqrt( (X3-X4)*(X3-X4) + (Y3-Y4)*(Y3-Y4) + (Z3-Z4)*(Z3-Z4) );

        // Maximum lenght 
		double Max_Lenght = Lenght_12;

		if ( Lenght_13 > Max_Lenght ) Max_Lenght = Lenght_13;
		if ( Lenght_14 > Max_Lenght ) Max_Lenght = Lenght_14;
		if ( Lenght_23 > Max_Lenght ) Max_Lenght = Lenght_23;
		if ( Lenght_24 > Max_Lenght ) Max_Lenght = Lenght_24;
		if ( Lenght_34 > Max_Lenght ) Max_Lenght = Lenght_34;

		return Max_Lenght;
	}

	//**************************************************************************************************
    //* - Adds a stabilization term S_h = hk^2 * mu-1 * ( curl_F, curl_E ) * ( div_F, div_E )
    //**************************************************************************************************	
    void VolumeElement_1bb::Add_hk2CurlDiv_Stabilization( Matrix< std::complex<double> >& StiffMatrix ) 
    {  
        // Loop counters
		int i, j, gp;

		// Material properties
        double mu_real = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cMu( mu_real, mu_imag );
      
		// Gaussian derivation weights
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Derivatives
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_1bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Diameter of the element
		double hk = Calculate_Element_Diameter();
		
		std::complex<double> cteCurlDiv = ( hk * hk ) / cMu;

		// CurlDiv volumetric integral 
		double intCurlDiv;

		// [Kxx], [Kyy], [Kzz] 
		for ( i=0; i<mNumNodes; i++ )
		{
			for ( j=i; j<mNumNodes; j++ )
			{
				intCurlDiv = 0.00; 
				
				for ( gp=0; gp<nGaussPoints; gp++ )
				{
					intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
				}

				StiffMatrix[ i               ][ j               ] += cteCurlDiv * intCurlDiv;
				StiffMatrix[ i + mNumNodes   ][ j + mNumNodes   ] += cteCurlDiv * intCurlDiv;
				StiffMatrix[ i + mNumNodes*2 ][ j + mNumNodes*2 ] += cteCurlDiv * intCurlDiv;
			}
		}

		// [Kxy] 
		for ( i=0; i<mNumNodes; i++ )
		{
			for ( j=0; j<mNumNodes; j++ )
			{
				intCurlDiv = 0.00;

				for ( gp=0; gp<nGaussPoints; gp++ )
				{
				    intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdy[j][gp] - dNdy[i][gp] * dNdx[j][gp] );
				}
	            
				StiffMatrix[ i ][ j + mNumNodes ] += cteCurlDiv * intCurlDiv;
			}
		}

		// [Kxz] 
		for ( i=0; i<mNumNodes; i++ )
		{
			for ( j=0; j<mNumNodes; j++ )
			{
				intCurlDiv = 0.00;

				for ( gp=0; gp<nGaussPoints; gp++ )
				{
					intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdx[j][gp] );
				}
	            
				StiffMatrix[ i ][ j + mNumNodes*2 ] += cteCurlDiv * intCurlDiv;
			}
		}

		// [Kyz] 
		for ( i=0; i<mNumNodes; i++ )
		{
			for ( j=0; j<mNumNodes; j++ )
			{
				intCurlDiv = 0.00; 

				for ( gp=0; gp<nGaussPoints; gp++ )
				{
					intCurlDiv += W[gp] * ( dNdy[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdy[j][gp] );
				}
	            
				StiffMatrix[ i + mNumNodes ][ j + mNumNodes*2 ] += cteCurlDiv * intCurlDiv;
			}
		}
		
        // Symmetric block
        for ( i=0; i<mNumDofs; i++ )
        {
            for ( j=i+1; j<mNumDofs; j++ )
            {
                StiffMatrix[ j ][ i ] +=  StiffMatrix[ i ][ j ];
            }
        }      
    }

	//**************************************************************************************************
    //* - Adds a stabilization term S_h = 
    //**************************************************************************************************	
    void VolumeElement_1bb::Add_CurlDivGBk_Stabilization( Matrix< std::complex<double> >& StiffMatrix ) 
    {  
	}

    //******************************************************************************************************
    //* - Stiffness matrix with local L2 projection method volume average
    //******************************************************************************************************	
	void VolumeElement_1bb::GetStiffMatrix_LL2P_VolumAvg( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( mNumDofs, mNumDofs, mZero );

		// Loop counters
		int i, j, gp;

		// Material properties
		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
        std::complex<double> cMu( mu_real , mu_imag );
      
		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_1bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< std::complex<double> > D_div; D_div.Resize( 1, mNumDofs, mZero );
		Matrix< std::complex<double> > C_cur; C_cur.Resize( 3, mNumDofs, mZero );

		for( j=0; j<mNumNodes; j++ )
		{
			double vintg_dNjdX = 0.00;
		    double vintg_dNjdY = 0.00;
		    double vintg_dNjdZ = 0.00;

            for( gp=0; gp<nGaussPoints; gp++ )
            {
            	vintg_dNjdX += W[gp] * dNdx[j][gp]; 
				vintg_dNjdY += W[gp] * dNdy[j][gp];
				vintg_dNjdZ += W[gp] * dNdz[j][gp];
            }      
			
			D_div[ 0 ][ j               ] = +vintg_dNjdX; 
			D_div[ 0 ][ j + mNumNodes   ] = +vintg_dNjdY; 
			D_div[ 0 ][ j + mNumNodes*2 ] = +vintg_dNjdZ; 

			C_cur[ 0 ][ j               ] =        0.00; 
            C_cur[ 0 ][ j + mNumNodes   ] = -vintg_dNjdZ; 
            C_cur[ 0 ][ j + mNumNodes*2 ] = +vintg_dNjdY; 
										    
			C_cur[ 1 ][ j               ] = +vintg_dNjdZ; 
            C_cur[ 1 ][ j + mNumNodes   ] =         0.00; 
            C_cur[ 1 ][ j + mNumNodes*2 ] = -vintg_dNjdX; 
										    
			C_cur[ 2 ][ j               ] = -vintg_dNjdY; 
            C_cur[ 2 ][ j + mNumNodes   ] = +vintg_dNjdX; 
            C_cur[ 2 ][ j + mNumNodes*2 ] =         0.00; 
		}

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij; N_nij.Resize( mNumDofs, mNumDofs, mZero );
		
		std::complex<double> cteNiNj = mFreq * mFreq * cEp; 

		for( i=0; i<mNumNodes; i++ )
		{
            for( j=0; j<mNumNodes; j++ )
            {
                double int_NiNj = 0.00;

                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	int_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
                }
            
            	N_nij[ i               ][ j               ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes   ][ j + mNumNodes   ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = cteNiNj * int_NiNj; 
            }
		}

		// Projection basis matrices
		Matrix< std::complex<double> > iA_div; iA_div.Resize( 1, 1, mZero );
		Matrix< std::complex<double> > iA_cur; iA_cur.Resize( 3, 3, mZero );
		
		iA_div[0][0] = 1.00 / mVolume;

		iA_cur[0][0] = iA_div[0][0]; 
		iA_cur[1][1] = iA_div[0][0];
		iA_cur[2][2] = iA_div[0][0];

	    // Scaling the projection matrices
		iA_div /= cMu ;
		iA_cur /= cMu ;

		// Final assembly: [ CT * Ac^(-1) * C ] + [ DT * Ad^(-1) * D ] - [ Ni * Nj ]
		StiffMatrix += C_cur.Transpose() * ( iA_cur * C_cur );

		if ( mPeso != 0.00 ) 
		{
			StiffMatrix += D_div.Transpose() * ( iA_div * D_div );
		}

		StiffMatrix -= N_nij;
	}

	//******************************************************************************************************
    //* - Stiffness matrix with local L2 projection method ( div(v), q ) - ( curl(v), w )
    //******************************************************************************************************	
    void VolumeElement_1bb::GetStiffMatrix_LL2P_CurlDivE( Matrix< std::complex<double> >& StiffMatrix ) 
    { 
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( mNumDofs, mNumDofs, mZero );

		// Loop counters
		int i, j, gp;

		// Material properties
		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
        std::complex<double> cMu( mu_real , mu_imag );
      
		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_1bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< std::complex<double> > D_div; D_div.Resize(  4, mNumDofs, mZero );
		Matrix< std::complex<double> > C_cur; C_cur.Resize( 12, mNumDofs, mZero );

		for( i=0; i<4; i++ )
		{
		    for( j=0; j<mNumNodes; j++ )
		    {
		    	double vintg_Ni_dNjdX = 0.00;
		        double vintg_Ni_dNjdY = 0.00;
		        double vintg_Ni_dNjdZ = 0.00;
		    
                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	vintg_Ni_dNjdX += W[gp] * N[i][gp] * dNdx[j][gp]; 
		    		vintg_Ni_dNjdY += W[gp] * N[i][gp] * dNdy[j][gp];
		    		vintg_Ni_dNjdZ += W[gp] * N[i][gp] * dNdz[j][gp];
                }      
		    	
		    	D_div[ i     ][ j               ] = +vintg_Ni_dNjdX; 
		    	D_div[ i     ][ j + mNumNodes   ] = +vintg_Ni_dNjdY; 
		    	D_div[ i     ][ j + mNumNodes*2 ] = +vintg_Ni_dNjdZ; 
		    			  
		    	C_cur[ i     ][ j               ] =            0.00; 
                C_cur[ i     ][ j + mNumNodes   ] = -vintg_Ni_dNjdZ; 
                C_cur[ i     ][ j + mNumNodes*2 ] = +vintg_Ni_dNjdY; 
		    								    
		    	C_cur[ i + 4 ][ j               ] = +vintg_Ni_dNjdZ; 
                C_cur[ i + 4 ][ j + mNumNodes   ] =            0.00; 
                C_cur[ i + 4 ][ j + mNumNodes*2 ] = -vintg_Ni_dNjdX; 
		    								    
		    	C_cur[ i + 8 ][ j               ] = -vintg_Ni_dNjdY; 
                C_cur[ i + 8 ][ j + mNumNodes   ] = +vintg_Ni_dNjdX; 
                C_cur[ i + 8 ][ j + mNumNodes*2 ] =            0.00; 
		    }
		}

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij; N_nij.Resize( mNumDofs, mNumDofs, mZero );
		
		std::complex<double> cteNiNj = mFreq * mFreq * cEp; 

		for( i=0; i<mNumNodes; i++ )
		{
            for( j=0; j<mNumNodes; j++ )
            {
                double int_NiNj = 0.00;

                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	int_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
                }
            
            	N_nij[ i               ][ j               ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes   ][ j + mNumNodes   ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = cteNiNj * int_NiNj; 
            }
		}

		// Projection basis matrices
		Matrix< std::complex<double> > iA_div; iA_div.Resize(  4,  4, mZero );
		Matrix< std::complex<double> > iA_cur; iA_cur.Resize( 12, 12, mZero );
		
		for( i=0; i<4; i++ )
		{
            for( j=0; j<4; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -4.00 / mVolume;
        		else          iA_div[ i ][ j ] = +16.0 / mVolume;                          
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}

		// Scaling the projection matrices
		iA_div /= cMu ;
		iA_cur /= cMu ;

		// Final assembly: [ CT * Ac^(-1) * C ] + [ DT * Ad^(-1) * D ] - [ Ni * Nj ]
		StiffMatrix += C_cur.Transpose() * ( iA_cur * C_cur );

		if ( mPeso != 0.00 ) 
		{
			StiffMatrix += D_div.Transpose() * ( iA_div * D_div );
		}

		StiffMatrix -= N_nij;
    }

    //******************************************************************************************************
    //* - Stiffness matrix with local L2 projection method ( v, -grad(q) ) - ( v, curl(w) )
    //******************************************************************************************************	
	void VolumeElement_1bb::GetStiffMatrix_LL2P_CurlDivQ( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( mNumDofs, mNumDofs, mZero );

		// Loop counters
		int i, j, gp;

		// Material properties
		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
        std::complex<double> cMu( mu_real , mu_imag );
      
		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_1bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< std::complex<double> > D_div; D_div.Resize(  4, mNumDofs, mZero );
		Matrix< std::complex<double> > C_cur; C_cur.Resize( 12, mNumDofs, mZero );

		for( i=0; i<4; i++ )
		{
		    for( j=0; j<mNumNodes; j++ )
		    {
		    	double vintg_dNidX_Nj = 0.00;
		        double vintg_dNidY_Nj = 0.00;
		        double vintg_dNidZ_Nj = 0.00;
		    
                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	vintg_dNidX_Nj += W[gp] * dNdx[i][gp] * N[j][gp]; 
		    		vintg_dNidY_Nj += W[gp] * dNdy[i][gp] * N[j][gp];
		    		vintg_dNidZ_Nj += W[gp] * dNdz[i][gp] * N[j][gp];
                }      
		    	
		    	D_div[ i     ][ j               ] = -vintg_dNidX_Nj; 
		    	D_div[ i     ][ j + mNumNodes   ] = -vintg_dNidY_Nj; 
		    	D_div[ i     ][ j + mNumNodes*2 ] = -vintg_dNidZ_Nj; 
		    			  
		    	C_cur[ i     ][ j               ] =            0.00; 
                C_cur[ i     ][ j + mNumNodes   ] = +vintg_dNidZ_Nj; 
                C_cur[ i     ][ j + mNumNodes*2 ] = -vintg_dNidY_Nj; 
		    								    
		    	C_cur[ i + 4 ][ j               ] = -vintg_dNidZ_Nj; 
                C_cur[ i + 4 ][ j + mNumNodes   ] =            0.00; 
                C_cur[ i + 4 ][ j + mNumNodes*2 ] = +vintg_dNidX_Nj; 
		    								    
		    	C_cur[ i + 8 ][ j               ] = +vintg_dNidY_Nj; 
                C_cur[ i + 8 ][ j + mNumNodes   ] = -vintg_dNidX_Nj; 
                C_cur[ i + 8 ][ j + mNumNodes*2 ] =            0.00; 
		    }
		}

		// Surface integral
		//Add_C_NxN_Matrix( C_cur ); 
	    //Add_D_NxN_Matrix( D_div );

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij; N_nij.Resize( mNumDofs, mNumDofs, mZero );
		
		std::complex<double> cteNiNj = mFreq * mFreq * cEp; 

		for( i=0; i<mNumNodes; i++ )
		{
            for( j=0; j<mNumNodes; j++ )
            {
                double int_NiNj = 0.00;

                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	int_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
                }
            
            	N_nij[ i               ][ j               ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes   ][ j + mNumNodes   ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = cteNiNj * int_NiNj; 
            }
		}

		// Projection basis matrices
		Matrix< std::complex<double> > iA_div; iA_div.Resize(  4,  4, mZero );
		Matrix< std::complex<double> > iA_cur; iA_cur.Resize( 12, 12, mZero );
		
		for( i=0; i<4; i++ )
		{
            for( j=0; j<4; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -4.00 / mVolume;
        		else          iA_div[ i ][ j ] = +16.0 / mVolume;                          
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}

		// Scaling the projection matrices
		iA_div /= cMu ;
		iA_cur /= cMu ;

		// Final assembly: [ CT * Ac^(-1) * C ] + [ DT * Ad^(-1) * D ] - [ Ni * Nj ]
		StiffMatrix += C_cur.Transpose() * ( iA_cur * C_cur );

		if ( mPeso != 0.00 ) 
		{
			StiffMatrix += D_div.Transpose() * ( iA_div * D_div );
		}

		StiffMatrix -= N_nij;
	}

	//*************************************************************************************************************************
    //* - Calculates the exterior normal 
    //*************************************************************************************************************************	
    void VolumeElement_1bb::Calculate_ExtNormal( Vector<double>& extN, Vector<int>& face )
    {
        double v1[3],v2[3];

        v2[0] = mNodes[ face[2] ]->X() - mNodes[ face[0] ]->X();
        v2[1] = mNodes[ face[2] ]->Y() - mNodes[ face[0] ]->Y();
        v2[2] = mNodes[ face[2] ]->Z() - mNodes[ face[0] ]->Z();

        v1[1] = mNodes[ face[1] ]->Y() - mNodes[ face[0] ]->Y();
        v1[2] = mNodes[ face[1] ]->Z() - mNodes[ face[0] ]->Z();
        v1[0] = mNodes[ face[1] ]->X() - mNodes[ face[0] ]->X();

        extN.resize( 3 );

        extN[0] = v2[1]*v1[2] - v2[2]*v1[1];
        extN[1] = v2[2]*v1[0] - v2[0]*v1[2];
        extN[2] = v2[0]*v1[1] - v2[1]*v1[0];

        double inv_extN_norm = 1.00 / sqrt( extN[0]*extN[0] + extN[1]*extN[1] + extN[2]*extN[2] );

		extN *= inv_extN_norm;
    }

    //****************************************************************************************
    //* - Calculates the area of the element
    //****************************************************************************************	
    double VolumeElement_1bb::Calculate_Area( Vector<int>& face ) 
    {
		double n[3],v1[3],v2[3];

        v2[0] = mNodes[ face[2] ]->X() - mNodes[ face[0] ]->X();
        v2[1] = mNodes[ face[2] ]->Y() - mNodes[ face[0] ]->Y();
        v2[2] = mNodes[ face[2] ]->Z() - mNodes[ face[0] ]->Z();

        v1[1] = mNodes[ face[1] ]->Y() - mNodes[ face[0] ]->Y();
        v1[2] = mNodes[ face[1] ]->Z() - mNodes[ face[0] ]->Z();
        v1[0] = mNodes[ face[1] ]->X() - mNodes[ face[0] ]->X();

        //area = 0.5*(v2 x v1)
        n[0] = v2[1]*v1[2] - v2[2]*v1[1];
        n[1] = v2[2]*v1[0] - v2[0]*v1[2];
        n[2] = v2[0]*v1[1] - v2[1]*v1[0];

        double area = 0.5*sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

		return area;
    }

    //******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::Add_C_NxN_Matrix( Matrix< std::complex<double> >& C_cur ) 
	{
		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 3 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; 
		Faces[1].resize( 3 ); Faces[1][0] = 1; Faces[1][1] = 3; Faces[1][2] = 2; 
		Faces[2].resize( 3 ); Faces[2][0] = 2; Faces[2][1] = 3; Faces[2][2] = 0; 
		Faces[3].resize( 3 ); Faces[3][0] = 3; Faces[3][1] = 1; Faces[3][2] = 0; 
		
		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			if ( (*gpNormalType)[ mNodes[ nf[0] ]->Id() ] != 'E' || (*gpNormalType)[ mNodes[ nf[1] ]->Id() ] != 'E' || (*gpNormalType)[ mNodes[ nf[2] ]->Id() ] != 'E' ) continue;

			// Loop indexes
		    int i, j, gp;
		    
            // Gauss points and weights
		    Vector<double> cX; 
		    Vector<double> cY;
		    Vector<double>  W;
		    
		    int nGaussPoints = GaussPoints2D_Order10( cX, cY, W );
		    
            // Jacobian
		    double detJ = 2.00 * Calculate_Area( nf );
		    
		    for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 
		    
            // Lagrangian base on Gauss points
		    Matrix<double> N;
		    
		    Lagrange2D_Ni_1st( N, cX, cY );

		    Vector<double> extN; Calculate_ExtNormal( extN, nf );

			//extN *= -1;
		    
		    for( i=0; i<3; i++ )
		    {
		        for( j=0; j<3; j++ )
		        {
		        	std::complex<double> Sintg_NiNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}
		        		
                    C_cur[ nf[i]     ][ nf[j] + mNumNodes   ] += (extN[2] * Sintg_NiNj); 
                    C_cur[ nf[i]     ][ nf[j] + mNumNodes*2 ] -= (extN[1] * Sintg_NiNj); 
		        					    	 					 
		        	C_cur[ nf[i] + 4 ][ nf[j]               ] -= (extN[2] * Sintg_NiNj); 
                    C_cur[ nf[i] + 4 ][ nf[j] + mNumNodes*2 ] += (extN[0] * Sintg_NiNj);
		        					 	 						
		        	C_cur[ nf[i] + 8 ][ nf[j]               ] += (extN[1] * Sintg_NiNj); 
                    C_cur[ nf[i] + 8 ][ nf[j] + mNumNodes   ] -= (extN[0] * Sintg_NiNj); 
		        }
		    }
		}
	}

	//******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::Add_D_NxN_Matrix( Matrix< std::complex<double> >& D_div ) 
	{
		//return;
		
		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 3 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; 
		Faces[1].resize( 3 ); Faces[1][0] = 1; Faces[1][1] = 3; Faces[1][2] = 2; 
		Faces[2].resize( 3 ); Faces[2][0] = 2; Faces[2][1] = 3; Faces[2][2] = 0; 
		Faces[3].resize( 3 ); Faces[3][0] = 3; Faces[3][1] = 1; Faces[3][2] = 0; 
		
		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			// 'E' = PEC normal | 'D' = Plasma Sheath | 'H' = PMC | 'T' = TEPMC
			if ( (((*gpNormalType)[ mNodes[ nf[0] ]->Id() ] != 'H') && ((*gpNormalType)[ mNodes[ nf[0] ]->Id() ] != 'T')) ||
				 (((*gpNormalType)[ mNodes[ nf[1] ]->Id() ] != 'H') && ((*gpNormalType)[ mNodes[ nf[1] ]->Id() ] != 'T')) ||
				 (((*gpNormalType)[ mNodes[ nf[2] ]->Id() ] != 'H') && ((*gpNormalType)[ mNodes[ nf[2] ]->Id() ] != 'T'))  ) continue;

			// Loop indexes
		    int i, j, gp;
		    
            // Gauss points and weights
		    Vector<double> cX; 
		    Vector<double> cY;
		    Vector<double>  W;
		    
		    int nGaussPoints = GaussPoints2D_Order10( cX, cY, W );
		    
            // Jacobian
		    double detJ = 2.00 * Calculate_Area( nf );
		    
		    for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 
		    
            // Lagrangian base on Gauss points
		    Matrix<double> N;
		    
		    Lagrange2D_Ni_1st( N, cX, cY );
		    
		    Vector<double> extN; Calculate_ExtNormal( extN, nf );

			//extN *= -1;
		    
		    for( i=0; i<3; i++ )
		    {
		        for( j=0; j<3; j++ )
		        {
		        	std::complex<double> Sintg_NiNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}

					D_div[ nf[i] ][ nf[j]               ] += (extN[0] * Sintg_NiNj);
		    	    D_div[ nf[i] ][ nf[j] + mNumNodes   ] += (extN[1] * Sintg_NiNj);
		    	    D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] += (extN[2] * Sintg_NiNj);
		        }
		    }
		}
	}

	//******************************************************************************************************
    //* - Stiffness matrix with local L2 projection method ( v, curl(w) ) - edge elements
    //******************************************************************************************************	
	void VolumeElement_1bb::GetStiffMatrix_LL2P_EdgeElem( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( mNumDofs, mNumDofs, mZero );

		// Loop counters
		int i, j, gp;

		// Material properties
		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
        std::complex<double> cMu( mu_real , mu_imag );
      
		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order8( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_1bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
        Matrix< std::complex<double> > D_div; D_div.Resize( 4, mNumDofs, mZero );
        
        for( i=0; i<4; i++ )
        {
            for( j=0; j<mNumNodes; j++ )
            {
            	double vintg_Ni_dNjdX = 0.00;
                double vintg_Ni_dNjdY = 0.00;
                double vintg_Ni_dNjdZ = 0.00;
            
                for( gp=0; gp<nGaussPoints; gp++ )
                {
                    vintg_Ni_dNjdX += W[gp] * N[i][gp] * dNdx[j][gp]; 
            		vintg_Ni_dNjdY += W[gp] * N[i][gp] * dNdy[j][gp];
            		vintg_Ni_dNjdZ += W[gp] * N[i][gp] * dNdz[j][gp];
                }      
            	
            	D_div[ i ][ j               ] = +vintg_Ni_dNjdX; 
            	D_div[ i ][ j + mNumNodes   ] = +vintg_Ni_dNjdY; 
            	D_div[ i ][ j + mNumNodes*2 ] = +vintg_Ni_dNjdZ; 
            }
        }

		//// ( u , -grad(q) )
		//Matrix< std::complex<double> > D_div; D_div.Resize( 4, mNumDofs, mZero );
		//
		//for( i=0; i<4; i++ )
		//{
		//    for( j=0; j<mNumNodes; j++ )
		//    {
		//    	double vintg_dNidX_Nj = 0.00;
		//        double vintg_dNidY_Nj = 0.00;
		//        double vintg_dNidZ_Nj = 0.00;
		//    
  //              for( gp=0; gp<nGaussPoints; gp++ )
  //              {
  //              	vintg_dNidX_Nj += W[gp] * dNdx[i][gp] * N[j][gp]; 
		//    		vintg_dNidY_Nj += W[gp] * dNdy[i][gp] * N[j][gp];
		//    		vintg_dNidZ_Nj += W[gp] * dNdz[i][gp] * N[j][gp];
  //              }      
		//    	
		//    	D_div[ i ][ j               ] = -vintg_dNidX_Nj; 
		//    	D_div[ i ][ j + mNumNodes   ] = -vintg_dNidY_Nj; 
		//    	D_div[ i ][ j + mNumNodes*2 ] = -vintg_dNidZ_Nj; 
		//    }
		//}

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij; N_nij.Resize( mNumDofs, mNumDofs, mZero );
		
		std::complex<double> cteNiNj = mFreq * mFreq * cEp; 

		for( i=0; i<mNumNodes; i++ )
		{
            for( j=0; j<mNumNodes; j++ )
            {
                double int_NiNj = 0.00;

                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	int_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
                }
            
            	N_nij[ i               ][ j               ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes   ][ j + mNumNodes   ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = cteNiNj * int_NiNj; 
            }
		}

		// Projection basis matrices
		Matrix< std::complex<double> > iA_div; iA_div.Resize( 4, 4, mZero );
		
		for( i=0; i<4; i++ )
		{
            for( j=0; j<4; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -4.00 / mVolume;
        		else          iA_div[ i ][ j ] = +16.0 / mVolume;                          
            }									   
		}

		iA_div /= cMu ;
		
		Matrix< std::complex<double> >  C_cur; 
		Matrix< std::complex<double> > iA_cur; 

		Calculate_EdgeElem_Matrices( C_cur, iA_cur ); 

		StiffMatrix += C_cur.Transpose() * ( iA_cur * C_cur );
		StiffMatrix += D_div.Transpose() * ( iA_div * D_div );
		StiffMatrix -= N_nij;
	}

	//******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::EdgeElementBase_3D( Matrix<double>&   Lx, Matrix<double>&   Ly, Matrix<double>&   Lz, Matrix<double>& N, 
		                                        Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz                    ) 
	{
		int nPoints = N[0].size();
		int nEdges  = 6  ;
		
		Lx.Resize( nEdges, nPoints );
		Ly.Resize( nEdges, nPoints );
		Lz.Resize( nEdges, nPoints );

		double edgeL[6];

		edgeL[0] = sqrt( ( X(2) - X(1) )*( X(2) - X(1) ) + ( Y(2) - Y(1) )*( Y(2) - Y(1) ) + ( Z(2) - Z(1) )*( Z(2) - Z(1) ) );
		edgeL[1] = sqrt( ( X(3) - X(1) )*( X(3) - X(1) ) + ( Y(3) - Y(1) )*( Y(3) - Y(1) ) + ( Z(3) - Z(1) )*( Z(3) - Z(1) ) );
		edgeL[2] = sqrt( ( X(4) - X(1) )*( X(4) - X(1) ) + ( Y(4) - Y(1) )*( Y(4) - Y(1) ) + ( Z(4) - Z(1) )*( Z(4) - Z(1) ) );
		edgeL[3] = sqrt( ( X(3) - X(2) )*( X(3) - X(2) ) + ( Y(3) - Y(2) )*( Y(3) - Y(2) ) + ( Z(3) - Z(2) )*( Z(3) - Z(2) ) );
		edgeL[4] = sqrt( ( X(4) - X(2) )*( X(4) - X(2) ) + ( Y(4) - Y(2) )*( Y(4) - Y(2) ) + ( Z(4) - Z(2) )*( Z(4) - Z(2) ) );
		edgeL[5] = sqrt( ( X(4) - X(3) )*( X(4) - X(3) ) + ( Y(4) - Y(3) )*( Y(4) - Y(3) ) + ( Z(4) - Z(3) )*( Z(4) - Z(3) ) );

        if ( mNodes[1]->Id() > mNodes[0]->Id() ) edgeL[0] *= -1.0;
        if ( mNodes[2]->Id() > mNodes[0]->Id() ) edgeL[1] *= -1.0;
        if ( mNodes[3]->Id() > mNodes[0]->Id() ) edgeL[2] *= -1.0;
        if ( mNodes[2]->Id() > mNodes[1]->Id() ) edgeL[3] *= -1.0;
        if ( mNodes[3]->Id() > mNodes[1]->Id() ) edgeL[4] *= -1.0;
        if ( mNodes[3]->Id() > mNodes[2]->Id() ) edgeL[5] *= -1.0;

		for( int p=0; p<nPoints; p++ ) 
        {
			Lx[ 0 ][ p ] = edgeL[ 0 ] * ( N[ 0 ][ p ] * dNdx[ 1 ][ p ] - N[ 1 ][ p ] * dNdx[ 0 ][ p ] );
			Ly[ 0 ][ p ] = edgeL[ 0 ] * ( N[ 0 ][ p ] * dNdy[ 1 ][ p ] - N[ 1 ][ p ] * dNdy[ 0 ][ p ] );
			Lz[ 0 ][ p ] = edgeL[ 0 ] * ( N[ 0 ][ p ] * dNdz[ 1 ][ p ] - N[ 1 ][ p ] * dNdz[ 0 ][ p ] );
   	 						     	  		  		     	 	 	 	      	   	 
			Lx[ 1 ][ p ] = edgeL[ 1 ] * ( N[ 0 ][ p ] * dNdx[ 2 ][ p ] - N[ 2 ][ p ] * dNdx[ 0 ][ p ] );
			Ly[ 1 ][ p ] = edgeL[ 1 ] * ( N[ 0 ][ p ] * dNdy[ 2 ][ p ] - N[ 2 ][ p ] * dNdy[ 0 ][ p ] );
			Lz[ 1 ][ p ] = edgeL[ 1 ] * ( N[ 0 ][ p ] * dNdz[ 2 ][ p ] - N[ 2 ][ p ] * dNdz[ 0 ][ p ] );
							     	  		  		     	 	 	 	      	   	 
			Lx[ 2 ][ p ] = edgeL[ 2 ] * ( N[ 0 ][ p ] * dNdx[ 3 ][ p ] - N[ 3 ][ p ] * dNdx[ 0 ][ p ] );
			Ly[ 2 ][ p ] = edgeL[ 2 ] * ( N[ 0 ][ p ] * dNdy[ 3 ][ p ] - N[ 3 ][ p ] * dNdy[ 0 ][ p ] );
			Lz[ 2 ][ p ] = edgeL[ 2 ] * ( N[ 0 ][ p ] * dNdz[ 3 ][ p ] - N[ 3 ][ p ] * dNdz[ 0 ][ p ] );
 	 						     	  		  		     	 	 	 	      	   	 
			Lx[ 3 ][ p ] = edgeL[ 3 ] * ( N[ 1 ][ p ] * dNdx[ 2 ][ p ] - N[ 2 ][ p ] * dNdx[ 1 ][ p ] );
			Ly[ 3 ][ p ] = edgeL[ 3 ] * ( N[ 1 ][ p ] * dNdy[ 2 ][ p ] - N[ 2 ][ p ] * dNdy[ 1 ][ p ] );
			Lz[ 3 ][ p ] = edgeL[ 3 ] * ( N[ 1 ][ p ] * dNdz[ 2 ][ p ] - N[ 2 ][ p ] * dNdz[ 1 ][ p ] );
	 						   	  	 	  		     	 	 	 	      	   	 
			Lx[ 4 ][ p ] = edgeL[ 4 ] * ( N[ 1 ][ p ] * dNdx[ 3 ][ p ] - N[ 3 ][ p ] * dNdx[ 1 ][ p ] );
			Ly[ 4 ][ p ] = edgeL[ 4 ] * ( N[ 1 ][ p ] * dNdy[ 3 ][ p ] - N[ 3 ][ p ] * dNdy[ 1 ][ p ] );
			Lz[ 4 ][ p ] = edgeL[ 4 ] * ( N[ 1 ][ p ] * dNdz[ 3 ][ p ] - N[ 3 ][ p ] * dNdz[ 1 ][ p ] );
							     	  		  		     	 	 	 	      	   	 
			Lx[ 5 ][ p ] = edgeL[ 5 ] * ( N[ 2 ][ p ] * dNdx[ 3 ][ p ] - N[ 3 ][ p ] * dNdx[ 2 ][ p ] );
			Ly[ 5 ][ p ] = edgeL[ 5 ] * ( N[ 2 ][ p ] * dNdy[ 3 ][ p ] - N[ 3 ][ p ] * dNdy[ 2 ][ p ] );
			Lz[ 5 ][ p ] = edgeL[ 5 ] * ( N[ 2 ][ p ] * dNdz[ 3 ][ p ] - N[ 3 ][ p ] * dNdz[ 2 ][ p ] );
        }	
	}

	//******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::EdgeElementBase_Derivatives_3D( Matrix<double>& dLxdx, Matrix<double>& dLxdy, Matrix<double>& dLxdz,
		                                                    Matrix<double>& dLydx, Matrix<double>& dLydy, Matrix<double>& dLydz,
		                                                    Matrix<double>& dLzdx, Matrix<double>& dLzdy, Matrix<double>& dLzdz,
		                                                    Matrix<double>& N, 
		                                                    Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz   ) 
	{
		int nPoints = N[0].size();
		int nEdges  = 6  ;
		
		dLxdx.Resize( nEdges, nPoints );
		dLxdy.Resize( nEdges, nPoints );
		dLxdz.Resize( nEdges, nPoints );

		dLydx.Resize( nEdges, nPoints );
		dLydy.Resize( nEdges, nPoints );
		dLydz.Resize( nEdges, nPoints );

		dLzdx.Resize( nEdges, nPoints );
		dLzdy.Resize( nEdges, nPoints );
		dLzdz.Resize( nEdges, nPoints );

		double edgeL[6];

		edgeL[0] = sqrt( ( X(2) - X(1) )*( X(2) - X(1) ) + ( Y(2) - Y(1) )*( Y(2) - Y(1) ) + ( Z(2) - Z(1) )*( Z(2) - Z(1) ) );
		edgeL[1] = sqrt( ( X(3) - X(1) )*( X(3) - X(1) ) + ( Y(3) - Y(1) )*( Y(3) - Y(1) ) + ( Z(3) - Z(1) )*( Z(3) - Z(1) ) );
		edgeL[2] = sqrt( ( X(4) - X(1) )*( X(4) - X(1) ) + ( Y(4) - Y(1) )*( Y(4) - Y(1) ) + ( Z(4) - Z(1) )*( Z(4) - Z(1) ) );
		edgeL[3] = sqrt( ( X(3) - X(2) )*( X(3) - X(2) ) + ( Y(3) - Y(2) )*( Y(3) - Y(2) ) + ( Z(3) - Z(2) )*( Z(3) - Z(2) ) );
		edgeL[4] = sqrt( ( X(4) - X(2) )*( X(4) - X(2) ) + ( Y(4) - Y(2) )*( Y(4) - Y(2) ) + ( Z(4) - Z(2) )*( Z(4) - Z(2) ) );
		edgeL[5] = sqrt( ( X(4) - X(3) )*( X(4) - X(3) ) + ( Y(4) - Y(3) )*( Y(4) - Y(3) ) + ( Z(4) - Z(3) )*( Z(4) - Z(3) ) );

        if ( mNodes[1]->Id() > mNodes[0]->Id() ) edgeL[0] *= -1.0;
        if ( mNodes[2]->Id() > mNodes[0]->Id() ) edgeL[1] *= -1.0;
        if ( mNodes[3]->Id() > mNodes[0]->Id() ) edgeL[2] *= -1.0;
        if ( mNodes[2]->Id() > mNodes[1]->Id() ) edgeL[3] *= -1.0;
        if ( mNodes[3]->Id() > mNodes[1]->Id() ) edgeL[4] *= -1.0;
        if ( mNodes[3]->Id() > mNodes[2]->Id() ) edgeL[5] *= -1.0;

		for( int p=0; p<nPoints; p++ ) 
        {
			dLxdx[ 0 ][ p ] = edgeL[ 0 ] * ( dNdx[ 0 ][ p ] * dNdx[ 1 ][ p ] - dNdx[ 1 ][ p ] * dNdx[ 0 ][ p ] );
			dLydx[ 0 ][ p ] = edgeL[ 0 ] * ( dNdx[ 0 ][ p ] * dNdy[ 1 ][ p ] - dNdx[ 1 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdx[ 0 ][ p ] = edgeL[ 0 ] * ( dNdx[ 0 ][ p ] * dNdz[ 1 ][ p ] - dNdx[ 1 ][ p ] * dNdz[ 0 ][ p ] );			
			
			dLxdy[ 0 ][ p ] = edgeL[ 0 ] * ( dNdy[ 0 ][ p ] * dNdx[ 1 ][ p ] - dNdy[ 1 ][ p ] * dNdx[ 0 ][ p ] );
			dLydy[ 0 ][ p ] = edgeL[ 0 ] * ( dNdy[ 0 ][ p ] * dNdy[ 1 ][ p ] - dNdy[ 1 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdy[ 0 ][ p ] = edgeL[ 0 ] * ( dNdy[ 0 ][ p ] * dNdz[ 1 ][ p ] - dNdy[ 1 ][ p ] * dNdz[ 0 ][ p ] );		

			dLxdz[ 0 ][ p ] = edgeL[ 0 ] * ( dNdz[ 0 ][ p ] * dNdx[ 1 ][ p ] - dNdz[ 1 ][ p ] * dNdx[ 0 ][ p ] );
			dLydz[ 0 ][ p ] = edgeL[ 0 ] * ( dNdz[ 0 ][ p ] * dNdy[ 1 ][ p ] - dNdz[ 1 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdz[ 0 ][ p ] = edgeL[ 0 ] * ( dNdz[ 0 ][ p ] * dNdz[ 1 ][ p ] - dNdz[ 1 ][ p ] * dNdz[ 0 ][ p ] );		
			


			dLxdx[ 1 ][ p ] = edgeL[ 1 ] * ( dNdx[ 0 ][ p ] * dNdx[ 2 ][ p ] - dNdx[ 2 ][ p ] * dNdx[ 0 ][ p ] );
			dLydx[ 1 ][ p ] = edgeL[ 1 ] * ( dNdx[ 0 ][ p ] * dNdy[ 2 ][ p ] - dNdx[ 2 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdx[ 1 ][ p ] = edgeL[ 1 ] * ( dNdx[ 0 ][ p ] * dNdz[ 2 ][ p ] - dNdx[ 2 ][ p ] * dNdz[ 0 ][ p ] );			
				   					 								
			dLxdy[ 1 ][ p ] = edgeL[ 1 ] * ( dNdy[ 0 ][ p ] * dNdx[ 2 ][ p ] - dNdy[ 2 ][ p ] * dNdx[ 0 ][ p ] );
			dLydy[ 1 ][ p ] = edgeL[ 1 ] * ( dNdy[ 0 ][ p ] * dNdy[ 2 ][ p ] - dNdy[ 2 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdy[ 1 ][ p ] = edgeL[ 1 ] * ( dNdy[ 0 ][ p ] * dNdz[ 2 ][ p ] - dNdy[ 2 ][ p ] * dNdz[ 0 ][ p ] );		
				  					 								
			dLxdz[ 1 ][ p ] = edgeL[ 1 ] * ( dNdz[ 0 ][ p ] * dNdx[ 2 ][ p ] - dNdz[ 2 ][ p ] * dNdx[ 0 ][ p ] );
			dLydz[ 1 ][ p ] = edgeL[ 1 ] * ( dNdz[ 0 ][ p ] * dNdy[ 2 ][ p ] - dNdz[ 2 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdz[ 1 ][ p ] = edgeL[ 1 ] * ( dNdz[ 0 ][ p ] * dNdz[ 2 ][ p ] - dNdz[ 2 ][ p ] * dNdz[ 0 ][ p ] );					
			


			dLxdx[ 2 ][ p ] = edgeL[ 2 ] * ( dNdx[ 0 ][ p ] * dNdx[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdx[ 0 ][ p ] );
			dLydx[ 2 ][ p ] = edgeL[ 2 ] * ( dNdx[ 0 ][ p ] * dNdy[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdx[ 2 ][ p ] = edgeL[ 2 ] * ( dNdx[ 0 ][ p ] * dNdz[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdz[ 0 ][ p ] );			
				  					 
			dLxdy[ 2 ][ p ] = edgeL[ 2 ] * ( dNdy[ 0 ][ p ] * dNdx[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdx[ 0 ][ p ] );
			dLydy[ 2 ][ p ] = edgeL[ 2 ] * ( dNdy[ 0 ][ p ] * dNdy[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdy[ 2 ][ p ] = edgeL[ 2 ] * ( dNdy[ 0 ][ p ] * dNdz[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdz[ 0 ][ p ] );		
				  					 
			dLxdz[ 2 ][ p ] = edgeL[ 2 ] * ( dNdz[ 0 ][ p ] * dNdx[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdx[ 0 ][ p ] );
			dLydz[ 2 ][ p ] = edgeL[ 2 ] * ( dNdz[ 0 ][ p ] * dNdy[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdy[ 0 ][ p ] );
			dLzdz[ 2 ][ p ] = edgeL[ 2 ] * ( dNdz[ 0 ][ p ] * dNdz[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdz[ 0 ][ p ] );				
   	 						     	  		  		     	 	 	 	      	   	 


			dLxdx[ 3 ][ p ] = edgeL[ 3 ] * ( dNdx[ 1 ][ p ] * dNdx[ 2 ][ p ] - dNdx[ 2 ][ p ] * dNdx[ 1 ][ p ] );
			dLydx[ 3 ][ p ] = edgeL[ 3 ] * ( dNdx[ 1 ][ p ] * dNdy[ 2 ][ p ] - dNdx[ 2 ][ p ] * dNdy[ 1 ][ p ] );
			dLzdx[ 3 ][ p ] = edgeL[ 3 ] * ( dNdx[ 1 ][ p ] * dNdz[ 2 ][ p ] - dNdx[ 2 ][ p ] * dNdz[ 1 ][ p ] );			
				  
			dLxdy[ 3 ][ p ] = edgeL[ 3 ] * ( dNdy[ 1 ][ p ] * dNdx[ 2 ][ p ] - dNdy[ 2 ][ p ] * dNdx[ 1 ][ p ] );
			dLydy[ 3 ][ p ] = edgeL[ 3 ] * ( dNdy[ 1 ][ p ] * dNdy[ 2 ][ p ] - dNdy[ 2 ][ p ] * dNdy[ 1 ][ p ] );
			dLzdy[ 3 ][ p ] = edgeL[ 3 ] * ( dNdy[ 1 ][ p ] * dNdz[ 2 ][ p ] - dNdy[ 2 ][ p ] * dNdz[ 1 ][ p ] );		
				  	
			dLxdz[ 3 ][ p ] = edgeL[ 3 ] * ( dNdz[ 1 ][ p ] * dNdx[ 2 ][ p ] - dNdz[ 2 ][ p ] * dNdx[ 1 ][ p ] );
			dLydz[ 3 ][ p ] = edgeL[ 3 ] * ( dNdz[ 1 ][ p ] * dNdy[ 2 ][ p ] - dNdz[ 2 ][ p ] * dNdy[ 1 ][ p ] );
			dLzdz[ 3 ][ p ] = edgeL[ 3 ] * ( dNdz[ 1 ][ p ] * dNdz[ 2 ][ p ] - dNdz[ 2 ][ p ] * dNdz[ 1 ][ p ] );				


	
 			dLxdx[ 4 ][ p ] = edgeL[ 4 ] * ( dNdx[ 1 ][ p ] * dNdx[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdx[ 1 ][ p ] );
			dLydx[ 4 ][ p ] = edgeL[ 4 ] * ( dNdx[ 1 ][ p ] * dNdy[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdy[ 1 ][ p ] );
			dLzdx[ 4 ][ p ] = edgeL[ 4 ] * ( dNdx[ 1 ][ p ] * dNdz[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdz[ 1 ][ p ] );			
				   
			dLxdy[ 4 ][ p ] = edgeL[ 4 ] * ( dNdy[ 1 ][ p ] * dNdx[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdx[ 1 ][ p ] );
			dLydy[ 4 ][ p ] = edgeL[ 4 ] * ( dNdy[ 1 ][ p ] * dNdy[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdy[ 1 ][ p ] );
			dLzdy[ 4 ][ p ] = edgeL[ 4 ] * ( dNdy[ 1 ][ p ] * dNdz[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdz[ 1 ][ p ] );		
				  
			dLxdz[ 4 ][ p ] = edgeL[ 4 ] * ( dNdz[ 1 ][ p ] * dNdx[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdx[ 1 ][ p ] );
			dLydz[ 4 ][ p ] = edgeL[ 4 ] * ( dNdz[ 1 ][ p ] * dNdy[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdy[ 1 ][ p ] );
			dLzdz[ 4 ][ p ] = edgeL[ 4 ] * ( dNdz[ 1 ][ p ] * dNdz[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdz[ 1 ][ p ] );				
	 						     	  		  		     	 	 	 	      	   	 
	
			dLxdx[ 5 ][ p ] = edgeL[ 5 ] * ( dNdx[ 2 ][ p ] * dNdx[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdx[ 2 ][ p ] );
			dLydx[ 5 ][ p ] = edgeL[ 5 ] * ( dNdx[ 2 ][ p ] * dNdy[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdy[ 2 ][ p ] );
			dLzdx[ 5 ][ p ] = edgeL[ 5 ] * ( dNdx[ 2 ][ p ] * dNdz[ 3 ][ p ] - dNdx[ 3 ][ p ] * dNdz[ 2 ][ p ] );			
				  
			dLxdy[ 5 ][ p ] = edgeL[ 5 ] * ( dNdy[ 2 ][ p ] * dNdx[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdx[ 2 ][ p ] );
			dLydy[ 5 ][ p ] = edgeL[ 5 ] * ( dNdy[ 2 ][ p ] * dNdy[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdy[ 2 ][ p ] );
			dLzdy[ 5 ][ p ] = edgeL[ 5 ] * ( dNdy[ 2 ][ p ] * dNdz[ 3 ][ p ] - dNdy[ 3 ][ p ] * dNdz[ 2 ][ p ] );		
				   
			dLxdz[ 5 ][ p ] = edgeL[ 5 ] * ( dNdz[ 2 ][ p ] * dNdx[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdx[ 2 ][ p ] );
			dLydz[ 5 ][ p ] = edgeL[ 5 ] * ( dNdz[ 2 ][ p ] * dNdy[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdy[ 2 ][ p ] );
			dLzdz[ 5 ][ p ] = edgeL[ 5 ] * ( dNdz[ 2 ][ p ] * dNdz[ 3 ][ p ] - dNdz[ 3 ][ p ] * dNdz[ 2 ][ p ] );		

        }	
	}

	//******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::Calculate_EdgeElem_Matrices( Matrix< std::complex<double> >& C_cur, Matrix< std::complex<double> >& iA_cur ) 
	{
		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order8( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( int gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_1bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		Matrix<double> Lx, Ly, Lz;

		EdgeElementBase_3D( Lx, Ly, Lz, N, dNdx, dNdy, dNdz );

		Matrix<double> dLxdx, dLxdy, dLxdz;
		Matrix<double> dLydx, dLydy, dLydz;
		Matrix<double> dLzdx, dLzdy, dLzdz;

		EdgeElementBase_Derivatives_3D( dLxdx, dLxdy, dLxdz, 
			                            dLydx, dLydy, dLydz, 
			                            dLzdx, dLzdy, dLzdz, 
			                            N, 
			                            dNdx, dNdy, dNdz   );

		int numEdges = 6;

		// C cur ( rot(u), w )
        C_cur.Resize( numEdges, mNumDofs, mZero );
        
        for( int i=0; i<numEdges; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                double vintg_Lxi_dNjdY = 0.00;
                double vintg_Lxi_dNjdZ = 0.00;
        
        		double vintg_Lyi_dNjdX = 0.00; 
                double vintg_Lyi_dNjdZ = 0.00;
        
        		double vintg_Lzi_dNjdX = 0.00; 
                double vintg_Lzi_dNjdY = 0.00;
            
                for( int gp=0; gp<nGaussPoints; gp++ )
                {
            		vintg_Lxi_dNjdY += W[gp] * Lx[i][gp] * dNdy[j][gp];
            		vintg_Lxi_dNjdZ += W[gp] * Lx[i][gp] * dNdz[j][gp];
        
        			vintg_Lyi_dNjdX += W[gp] * Ly[i][gp] * dNdx[j][gp];
            		vintg_Lyi_dNjdZ += W[gp] * Ly[i][gp] * dNdz[j][gp];
                    
        			vintg_Lzi_dNjdX += W[gp] * Lz[i][gp] * dNdx[j][gp];
            		vintg_Lzi_dNjdY += W[gp] * Lz[i][gp] * dNdy[j][gp];
                }      
            	
            	C_cur[ i ][ j               ] = vintg_Lyi_dNjdZ - vintg_Lzi_dNjdY; 
                C_cur[ i ][ j + mNumNodes   ] = vintg_Lzi_dNjdX - vintg_Lxi_dNjdZ; 
                C_cur[ i ][ j + mNumNodes*2 ] = vintg_Lxi_dNjdY - vintg_Lyi_dNjdX; 
            }
        }

		//// C cur ( u, rot(w) )
		//C_cur.Resize( numEdges, mNumDofs, mZero );

		//for( int i=0; i<numEdges; i++ )
		//{
		//    for( int j=0; j<mNumNodes; j++ )
		//    {
		//        double vintg_rotLix_Nj = 0.00;
		//        double vintg_rotLiy_Nj = 0.00;
		//		double vintg_rotLiz_Nj = 0.00; 
		//        
  //              for( int gp=0; gp<nGaussPoints; gp++ )
  //              {
		//    		vintg_rotLix_Nj += W[gp] * ( dLzdy[i][gp] - dLydz[i][gp] ) * N[j][gp];
		//    		vintg_rotLiy_Nj += W[gp] * ( dLxdz[i][gp] - dLzdx[i][gp] ) * N[j][gp];
		//			vintg_rotLiz_Nj += W[gp] * ( dLydx[i][gp] - dLxdy[i][gp] ) * N[j][gp];
  //              }      
		//    	
		//    	C_cur[ i ][ j               ] = vintg_rotLix_Nj; 
  //              C_cur[ i ][ j + mNumNodes   ] = vintg_rotLiy_Nj; 
  //              C_cur[ i ][ j + mNumNodes*2 ] = vintg_rotLiz_Nj; 
		//    }
		//}

		Matrix<double> WiWj( numEdges, numEdges );

		for( int i=0; i<numEdges; i++ )
		{
            for( int j=0; j<numEdges; j++ )
            {
				double vintg_Li_Lj = 0.00;
		    
                for( int gp=0; gp<nGaussPoints; gp++ )
                {
		    		vintg_Li_Lj += W[gp] * ( Lx[i][gp]*Lx[j][gp] + Ly[i][gp]*Ly[j][gp] + Lz[i][gp]*Lz[j][gp] );
				}
        		
				WiWj[ i ][ j ] = vintg_Li_Lj;
            }									   
		}

		Invert_Matrix( WiWj );
		
		iA_cur.Resize( numEdges, numEdges, mZero );

		for( int i=0; i<numEdges; i++ )
		{
            for( int j=0; j<numEdges; j++ )
            {
				iA_cur[ i ][ j ] = WiWj[ i ][ j ];
            }									   
		}

		double mu_real = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cMu( mu_real, mu_imag );

		iA_cur /= cMu ;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

    //===========================================================================
    void LUFactor(Matrix<double>& a, Vector<int>& ipivot, int n, double &det)
    //---------------------------------------------------------------------------
    // Performs LU factorization of (n x n) matrix a (diag(L) = 1). On exit,
    // replaces upper triangle and diagonal with U, and lower triangle, with L.
    // Uses partial pivoting on columns.
    // a - coefficient matrix (n x n); LU decomposition on exit
    // ipivot - array of pivot row indexes (output)
    // det - determinant of coefficient matrix (output).
    //---------------------------------------------------------------------------
    {
        double amax, sum, t;
        int i, imax, j, k;
        det = 1e0;
        
		for (j=0; j<n; j++) 
		{ 
			// loop over columns
            for (i=0; i<j; i++) 
		    { 
				// elements of matrix U
                sum = a[i][j];
                for (k=0; k<i; k++) sum -= a[i][k]*a[k][j];
                a[i][j] = sum;
            }

            amax = 0e0;

            for (i=j; i<n; i++) 
			{ 
				// elements of matrix L
                sum = a[i][j]; // undivided by pivot
			    
                for (k=0; k<j; k++) sum -= a[i][k]*a[k][j];
			    
                a[i][j] = sum;
			    
                // determine pivot
                if (amax < fabs(a[i][j])) 
			    {
			    	amax = fabs(a[i][j]); 
					imax = i;
			    }
            }

            if (amax == 0e0)
            { 
				std::cout<<"LUFactor: singular matrix !"<<std::endl;
		    	det = 0e0; 
		    	return; 
		    }
		    
            ipivot[j] = imax; // store pivot row index
            
			// interchange rows imax and j
            if (imax != j) 
		    { 
		    	// to put pivot on diagonal
                det = -det;
                for (k=0; k<n; k++)
                { 
		    		t = a[imax][k]; a[imax][k] = a[j][k]; a[j][k] = t; 
		    	}
            }

			det *= a[j][j]; // multiply determinant with pivot
            t = 1e0/a[j][j]; // divide elements of L by pivot

            for (i=j+1; i<n; i++) a[i][j] *= t;
        }
    }

    //===========================================================================
    void LUSystem( Matrix<double>& a, Vector<int>& ipivot, Vector<double>& b, int n )
    //---------------------------------------------------------------------------
    // Solves linear system a x = b of order n by LU factorization.
    // a - LU decomposition of coefficient matrix (returned by LUFactor)
    // ipivot - array of pivot row indexes (input)
    // b - vector of constant terms (input); solution x (on exit)
    //---------------------------------------------------------------------------
    {
        double sum;
        int i, j;

        for (i=0; i<n; i++) 
		{ 
			// solves Ly = b
            sum = b[ipivot[i]];
            b[ipivot[i]] = b[i];
			for (j=0; j<i; j++) sum -= a[i][j]*b[j];
            b[i] = sum;
        }

        for (i=n-1; i>=0; i--) 
		{ 
			// solves Ux = y
            sum = b[i];
            for (j=i+1; j<n; j++) sum -= a[i][j]*b[j];
            b[i] = sum/a[i][i];
        }
    }

    void VolumeElement_1bb::Invert_Matrix( Matrix<double>& a )
    //---------------------------------------------------------------------------
    // Calculates inverse of real matrix by LU factorization.
    // a - (n x n) matrix (input); a^(-1) (output)
    // det - determinant of coefficient matrix (output).
    // Calls: LUFactor, LUSystem.
    //---------------------------------------------------------------------------
    {
		int n = a[0].size();

        int i, j;

		double det;

        Matrix<double> ainv  (n,n);
        Vector<double> b     (n)  ;
        Vector< int  > ipivot(n)  ;

        LUFactor( a, ipivot, n, det ); // LU factorization of a

        if (det == 0e0) // singular matrix
        { 
			std::cout<<"MatInv: singular matrix !"<<std::endl;
			return; 
		}

        for (j=0; j<n; j++) 
		{ 
            for (i=0; i<n; i++) b[i] = 0e0; 
			b[j] = 1e0;
			LUSystem( a, ipivot, b, n );
			for (i=0; i<n; i++) ainv[i][j] = b[i]; 
        }
        
		for (j=0; j<n; j++) 
		{
			// copy inverse in a
            for (i=0; i<n; i++) 
			{
				a[i][j] = ainv[i][j];
			}
		}
    }

} /* end namespace Kratos */ 