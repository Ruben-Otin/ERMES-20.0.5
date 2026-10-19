
#include "VolumeElement_1eg_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{ 
    //******************************************************************************************************************************************
    // - Global index of the DOFs
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( mNumDofs );

		if( mPotentialsOn )
		{
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                EquationId[ i               ] = mNodes[ i ]->pDofcAx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcAy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcAz()->EquationId();
            }	

            for( int i=0; i<4; i++ ) 
            { 
                EquationId[ i + 3*mNumNodes ] = mNodes[ i ]->pDofcVs()->EquationId();
            }	
		}											   
		else										   
		{											   
            for( int i=0; i<mNumNodes; i++ ) 		   
            { 										   
                EquationId[ i               ] = mNodes[ i ]->pDofcEx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcEy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcEz()->EquationId();
            }		
		}
    }

    //******************************************************************************************************************************************
    // - Volume of the element
    //******************************************************************************************************************************************
    double VolumeElement_1eg_FullWave::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //******************************************************************************************************************************************
    // - 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::Calculate_DN( double DN[ 3 ][ 4 ] )
    {
        double Cte = 1.0 / ( 6.0 * mVolume ); 

        // dNi/dx
        DN[ 0 ][ 0 ] = Cte * ( Y(3)*Z(2) + Y(4)*Z(3) + Y(2)*Z(4) - Y(3)*Z(4) - Y(2)*Z(3) - Y(4)*Z(2) );
        DN[ 0 ][ 1 ] = Cte * ( Y(3)*Z(4) + Y(4)*Z(1) + Y(1)*Z(3) - Y(3)*Z(1) - Y(1)*Z(4) - Y(4)*Z(3) );
        DN[ 0 ][ 2 ] = Cte * ( Y(2)*Z(1) + Y(1)*Z(4) + Y(4)*Z(2) - Y(2)*Z(4) - Y(4)*Z(1) - Y(1)*Z(2) );
        DN[ 0 ][ 3 ] = Cte * ( Y(2)*Z(3) + Y(3)*Z(1) + Y(1)*Z(2) - Y(2)*Z(1) - Y(3)*Z(2) - Y(1)*Z(3) );
        
        // dNi/dy
        DN[ 1 ][ 0 ] = Cte * ( X(3)*Z(4) + X(4)*Z(2) + X(2)*Z(3) - X(3)*Z(2) - X(2)*Z(4) - X(4)*Z(3) ); 
        DN[ 1 ][ 1 ] = Cte * ( X(3)*Z(1) + X(1)*Z(4) + X(4)*Z(3) - X(1)*Z(3) - X(3)*Z(4) - X(4)*Z(1) );
        DN[ 1 ][ 2 ] = Cte * ( X(2)*Z(4) + X(4)*Z(1) + X(1)*Z(2) - X(2)*Z(1) - X(4)*Z(2) - X(1)*Z(4) );
        DN[ 1 ][ 3 ] = Cte * ( X(1)*Z(3) + X(3)*Z(2) + X(2)*Z(1) - X(3)*Z(1) - X(1)*Z(2) - X(2)*Z(3) );
        
        // dNi/dz
        DN[ 2 ][ 0 ] = Cte * ( X(4)*Y(3) + X(2)*Y(4) + X(3)*Y(2) - X(2)*Y(3) - X(4)*Y(2) - X(3)*Y(4) );
        DN[ 2 ][ 1 ] = Cte * ( X(3)*Y(4) + X(4)*Y(1) + X(1)*Y(3) - X(3)*Y(1) - X(4)*Y(3) - X(1)*Y(4) ); 
        DN[ 2 ][ 2 ] = Cte * ( X(2)*Y(1) + X(4)*Y(2) + X(1)*Y(4) - X(4)*Y(1) - X(2)*Y(4) - X(1)*Y(2) ); 
        DN[ 2 ][ 3 ] = Cte * ( X(2)*Y(3) + X(3)*Y(1) + X(1)*Y(2) - X(2)*Y(1) - X(3)*Y(2) - X(1)*Y(3) );    
    }

    //******************************************************************************************************************************************
    // - Element edge lengths with sign
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::Calculate_Edge_Lengths( Vector<double>& eLength )
    {
		// Edge coordinates vector
		Vector< Vector<double> > R( mNumEdges ); 
		
		for( int i=0; i<mNumEdges; i++ ) 
		{
			R[ i ].resize( 3 );
		}

		// Edge: [0 1]
		R[ 0 ][ 0 ] = X(2) - X(1); 
		R[ 0 ][ 1 ] = Y(2) - Y(1); 
		R[ 0 ][ 2 ] = Z(2) - Z(1);

		// Edge: [1 2]
		R[ 1 ][ 0 ] = X(3) - X(2); 
		R[ 1 ][ 1 ] = Y(3) - Y(2); 
		R[ 1 ][ 2 ] = Z(3) - Z(2);

		// Edge: [2 0]
		R[ 2 ][ 0 ] = X(1) - X(3); 
		R[ 2 ][ 1 ] = Y(1) - Y(3); 
		R[ 2 ][ 2 ] = Z(1) - Z(3);

		// Edge: [0 3]
		R[ 3 ][ 0 ] = X(4) - X(1); 
		R[ 3 ][ 1 ] = Y(4) - Y(1); 
		R[ 3 ][ 2 ] = Z(4) - Z(1);

		// Edge: [1 3]
		R[ 4 ][ 0 ] = X(4) - X(2); 
		R[ 4 ][ 1 ] = Y(4) - Y(2); 
		R[ 4 ][ 2 ] = Z(4) - Z(2);

		// Edge: [2 3]
		R[ 5 ][ 0 ] = X(4) - X(3); 
		R[ 5 ][ 1 ] = Y(4) - Y(3); 
		R[ 5 ][ 2 ] = Z(4) - Z(3);

		// Edge lengths
		for( int i=0; i<mNumEdges; i++ )
		{
			eLength[ i ] = sqrt( R[ i ][ 0 ]*R[ i ][ 0 ] + R[ i ][ 1 ]*R[ i ][ 1 ] + R[ i ][ 2 ]*R[ i ][ 2 ] );
		}
		
		// Check directions
		if( mNodes[ 0 ]->Id() > mNodes[ 1 ]->Id() ) eLength[ 0 ] *= -1.0;
		if( mNodes[ 1 ]->Id() > mNodes[ 2 ]->Id() ) eLength[ 1 ] *= -1.0;
		if( mNodes[ 2 ]->Id() > mNodes[ 0 ]->Id() ) eLength[ 2 ] *= -1.0;
		if( mNodes[ 0 ]->Id() > mNodes[ 3 ]->Id() ) eLength[ 3 ] *= -1.0;
		if( mNodes[ 1 ]->Id() > mNodes[ 3 ]->Id() ) eLength[ 4 ] *= -1.0;
		if( mNodes[ 2 ]->Id() > mNodes[ 3 ]->Id() ) eLength[ 5 ] *= -1.0;
	}

    //******************************************************************************************************************************************
    // - 1st order Whitney base for 3D tetrahedrons
    //******************************************************************************************************************************************
	void VolumeElement_1eg_FullWave::Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz,
		                                               Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ )
	{
		double L1, L2, L3, L4;
		
		int numPoints = cX.size();

		Nx.Resize( mNumEdges, numPoints );
		Ny.Resize( mNumEdges, numPoints );
		Nz.Resize( mNumEdges, numPoints );
		
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		double dL1_dX = DN[ 0 ][ 0 ], dL1_dY = DN[ 1 ][ 0 ], dL1_dZ = DN[ 2 ][ 0 ];
		double dL2_dX = DN[ 0 ][ 1 ], dL2_dY = DN[ 1 ][ 1 ], dL2_dZ = DN[ 2 ][ 1 ];
		double dL3_dX = DN[ 0 ][ 2 ], dL3_dY = DN[ 1 ][ 2 ], dL3_dZ = DN[ 2 ][ 2 ];
		double dL4_dX = DN[ 0 ][ 3 ], dL4_dY = DN[ 1 ][ 3 ], dL4_dZ = DN[ 2 ][ 3 ];

		Vector<double> eLength( mNumEdges ); Calculate_Edge_Lengths( eLength );

		for( int p=0; p<numPoints; p++ ) 
        {
            L2 = cX[ p ];
            L3 = cY[ p ];
			L4 = cZ[ p ];
            L1 = 1.0 - L2 - L3 - L4;

			// Edge: [0 1]
			Nx[ 0 ][ p ] = ( L1 * dL2_dX - L2 * dL1_dX ) * eLength[ 0 ];
            Ny[ 0 ][ p ] = ( L1 * dL2_dY - L2 * dL1_dY ) * eLength[ 0 ];
		    Nz[ 0 ][ p ] = ( L1 * dL2_dZ - L2 * dL1_dZ ) * eLength[ 0 ];

			// Edge: [1 2]
			Nx[ 1 ][ p ] = ( L2 * dL3_dX - L3 * dL2_dX ) * eLength[ 1 ];
            Ny[ 1 ][ p ] = ( L2 * dL3_dY - L3 * dL2_dY ) * eLength[ 1 ];
		    Nz[ 1 ][ p ] = ( L2 * dL3_dZ - L3 * dL2_dZ ) * eLength[ 1 ];

			// Edge: [2 0]
			Nx[ 2 ][ p ] = ( L3 * dL1_dX - L1 * dL3_dX ) * eLength[ 2 ];
            Ny[ 2 ][ p ] = ( L3 * dL1_dY - L1 * dL3_dY ) * eLength[ 2 ];
		    Nz[ 2 ][ p ] = ( L3 * dL1_dZ - L1 * dL3_dZ ) * eLength[ 2 ];

			// Edge: [0 3]
			Nx[ 3 ][ p ] = ( L1 * dL4_dX - L4 * dL1_dX ) * eLength[ 3 ];
            Ny[ 3 ][ p ] = ( L1 * dL4_dY - L4 * dL1_dY ) * eLength[ 3 ];
		    Nz[ 3 ][ p ] = ( L1 * dL4_dZ - L4 * dL1_dZ ) * eLength[ 3 ];

			// Edge: [1 3]
			Nx[ 4 ][ p ] = ( L2 * dL4_dX - L4 * dL2_dX ) * eLength[ 4 ];
            Ny[ 4 ][ p ] = ( L2 * dL4_dY - L4 * dL2_dY ) * eLength[ 4 ];
		    Nz[ 4 ][ p ] = ( L2 * dL4_dZ - L4 * dL2_dZ ) * eLength[ 4 ];

			// Edge: [2 3]
			Nx[ 5 ][ p ] = ( L3 * dL4_dX - L4 * dL3_dX ) * eLength[ 5 ];
            Ny[ 5 ][ p ] = ( L3 * dL4_dY - L4 * dL3_dY ) * eLength[ 5 ];
		    Nz[ 5 ][ p ] = ( L3 * dL4_dZ - L4 * dL3_dZ ) * eLength[ 5 ];
        }
	}

    //******************************************************************************************************************************************
    // - dNx/dx, dNx/dy, dNx/dz on cX, cY, cZ points
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::Derivatives_On_Points( Vector<double>& dNx_dX, Vector<double>& dNx_dY, Vector<double>& dNx_dZ,
		                                                    Vector<double>& dNy_dX, Vector<double>& dNy_dY, Vector<double>& dNy_dZ,
		                                                    Vector<double>& dNz_dX, Vector<double>& dNz_dY, Vector<double>& dNz_dZ )
    {
		// Li derivatives
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		double dL1_dX = DN[ 0 ][ 0 ], dL1_dY = DN[ 1 ][ 0 ], dL1_dZ = DN[ 2 ][ 0 ];
		double dL2_dX = DN[ 0 ][ 1 ], dL2_dY = DN[ 1 ][ 1 ], dL2_dZ = DN[ 2 ][ 1 ];
		double dL3_dX = DN[ 0 ][ 2 ], dL3_dY = DN[ 1 ][ 2 ], dL3_dZ = DN[ 2 ][ 2 ];
		double dL4_dX = DN[ 0 ][ 3 ], dL4_dY = DN[ 1 ][ 3 ], dL4_dZ = DN[ 2 ][ 3 ];

		// Element edges lengths
		Vector<double> eLength( mNumEdges ); Calculate_Edge_Lengths( eLength );

		// Whitney base derivatives
		dNx_dX.resize( mNumEdges ); 
		dNx_dY.resize( mNumEdges ); 
		dNx_dZ.resize( mNumEdges ); 
				   			
		dNy_dX.resize( mNumEdges ); 
		dNy_dY.resize( mNumEdges ); 
		dNy_dZ.resize( mNumEdges );
				  			
		dNz_dX.resize( mNumEdges ); 
		dNz_dY.resize( mNumEdges ); 
		dNz_dZ.resize( mNumEdges );

        // Edge: [0 1]
		dNx_dX[ 0 ] = ( dL1_dX * dL2_dX - dL2_dX * dL1_dX ) * eLength[ 0 ];
		dNx_dY[ 0 ] = ( dL1_dY * dL2_dX - dL2_dY * dL1_dX ) * eLength[ 0 ];
		dNx_dZ[ 0 ] = ( dL1_dZ * dL2_dX - dL2_dZ * dL1_dX ) * eLength[ 0 ];
			    													   
  		dNy_dX[ 0 ] = ( dL1_dX * dL2_dY - dL2_dX * dL1_dY ) * eLength[ 0 ];
		dNy_dY[ 0 ] = ( dL1_dY * dL2_dY - dL2_dY * dL1_dY ) * eLength[ 0 ];
		dNy_dZ[ 0 ] = ( dL1_dZ * dL2_dY - dL2_dZ * dL1_dY ) * eLength[ 0 ];
			    													   
  		dNz_dX[ 0 ] = ( dL1_dX * dL2_dZ - dL2_dX * dL1_dZ ) * eLength[ 0 ];
		dNz_dY[ 0 ] = ( dL1_dY * dL2_dZ - dL2_dY * dL1_dZ ) * eLength[ 0 ];
		dNz_dZ[ 0 ] = ( dL1_dZ * dL2_dZ - dL2_dZ * dL1_dZ ) * eLength[ 0 ];

		// Edge: [1 2]
		dNx_dX[ 1 ] = ( dL2_dX * dL3_dX - dL3_dX * dL2_dX ) * eLength[ 1 ];
		dNx_dY[ 1 ] = ( dL2_dY * dL3_dX - dL3_dY * dL2_dX ) * eLength[ 1 ];
		dNx_dZ[ 1 ] = ( dL2_dZ * dL3_dX - dL3_dZ * dL2_dX ) * eLength[ 1 ];
			    	 												   
 		dNy_dX[ 1 ] = ( dL2_dX * dL3_dY - dL3_dX * dL2_dY ) * eLength[ 1 ];
		dNy_dY[ 1 ] = ( dL2_dY * dL3_dY - dL3_dY * dL2_dY ) * eLength[ 1 ];
		dNy_dZ[ 1 ] = ( dL2_dZ * dL3_dY - dL3_dZ * dL2_dY ) * eLength[ 1 ];
			    	 										   		   
 		dNz_dX[ 1 ] = ( dL2_dX * dL3_dZ - dL3_dX * dL2_dZ ) * eLength[ 1 ];
		dNz_dY[ 1 ] = ( dL2_dY * dL3_dZ - dL3_dY * dL2_dZ ) * eLength[ 1 ];
		dNz_dZ[ 1 ] = ( dL2_dZ * dL3_dZ - dL3_dZ * dL2_dZ ) * eLength[ 1 ];

		// Edge: [2 0]
		dNx_dX[ 2 ] = ( dL3_dX * dL1_dX - dL1_dX * dL3_dX ) * eLength[ 2 ];
		dNx_dY[ 2 ] = ( dL3_dY * dL1_dX - dL1_dY * dL3_dX ) * eLength[ 2 ];
		dNx_dZ[ 2 ] = ( dL3_dZ * dL1_dX - dL1_dZ * dL3_dX ) * eLength[ 2 ];
			    	 												     
		dNy_dX[ 2 ] = ( dL3_dX * dL1_dY - dL1_dX * dL3_dY ) * eLength[ 2 ];
		dNy_dY[ 2 ] = ( dL3_dY * dL1_dY - dL1_dY * dL3_dY ) * eLength[ 2 ];
		dNy_dZ[ 2 ] = ( dL3_dZ * dL1_dY - dL1_dZ * dL3_dY ) * eLength[ 2 ];
			    	 										   		   
 		dNz_dX[ 2 ] = ( dL3_dX * dL1_dZ - dL1_dX * dL3_dZ ) * eLength[ 2 ];
		dNz_dY[ 2 ] = ( dL3_dY * dL1_dZ - dL1_dY * dL3_dZ ) * eLength[ 2 ];
		dNz_dZ[ 2 ] = ( dL3_dZ * dL1_dZ - dL1_dZ * dL3_dZ ) * eLength[ 2 ];

		// Edge: [0 3]
		dNx_dX[ 3 ] = ( dL1_dX * dL4_dX - dL4_dX * dL1_dX ) * eLength[ 3 ];
		dNx_dY[ 3 ] = ( dL1_dY * dL4_dX - dL4_dY * dL1_dX ) * eLength[ 3 ];
		dNx_dZ[ 3 ] = ( dL1_dZ * dL4_dX - dL4_dZ * dL1_dX ) * eLength[ 3 ];
			  	  												   	   
 		dNy_dX[ 3 ] = ( dL1_dX * dL4_dY - dL4_dX * dL1_dY ) * eLength[ 3 ];
		dNy_dY[ 3 ] = ( dL1_dY * dL4_dY - dL4_dY * dL1_dY ) * eLength[ 3 ];
		dNy_dZ[ 3 ] = ( dL1_dZ * dL4_dY - dL4_dZ * dL1_dY ) * eLength[ 3 ];
			    	 												     
		dNz_dX[ 3 ] = ( dL1_dX * dL4_dZ - dL4_dX * dL1_dZ ) * eLength[ 3 ];
		dNz_dY[ 3 ] = ( dL1_dY * dL4_dZ - dL4_dY * dL1_dZ ) * eLength[ 3 ];
		dNz_dZ[ 3 ] = ( dL1_dZ * dL4_dZ - dL4_dZ * dL1_dZ ) * eLength[ 3 ];

		// Edge: [1 3]
		dNx_dX[ 4 ] = ( dL2_dX * dL4_dX - dL4_dX * dL2_dX ) * eLength[ 4 ];
		dNx_dY[ 4 ] = ( dL2_dY * dL4_dX - dL4_dY * dL2_dX ) * eLength[ 4 ];
		dNx_dZ[ 4 ] = ( dL2_dZ * dL4_dX - dL4_dZ * dL2_dX ) * eLength[ 4 ];
			  	  												   	   
 		dNy_dX[ 4 ] = ( dL2_dX * dL4_dY - dL4_dX * dL2_dY ) * eLength[ 4 ];
		dNy_dY[ 4 ] = ( dL2_dY * dL4_dY - dL4_dY * dL2_dY ) * eLength[ 4 ];
		dNy_dZ[ 4 ] = ( dL2_dZ * dL4_dY - dL4_dZ * dL2_dY ) * eLength[ 4 ];
			    	 												     
		dNz_dX[ 4 ] = ( dL2_dX * dL4_dZ - dL4_dX * dL2_dZ ) * eLength[ 4 ];
		dNz_dY[ 4 ] = ( dL2_dY * dL4_dZ - dL4_dY * dL2_dZ ) * eLength[ 4 ];
		dNz_dZ[ 4 ] = ( dL2_dZ * dL4_dZ - dL4_dZ * dL2_dZ ) * eLength[ 4 ];

		// Edge: [2 3]
		dNx_dX[ 5 ] = ( dL3_dX * dL4_dX - dL4_dX * dL3_dX ) * eLength[ 5 ];
		dNx_dY[ 5 ] = ( dL3_dY * dL4_dX - dL4_dY * dL3_dX ) * eLength[ 5 ];
		dNx_dZ[ 5 ] = ( dL3_dZ * dL4_dX - dL4_dZ * dL3_dX ) * eLength[ 5 ];
			  	  												   	   
 		dNy_dX[ 5 ] = ( dL3_dX * dL4_dY - dL4_dX * dL3_dY ) * eLength[ 5 ];
		dNy_dY[ 5 ] = ( dL3_dY * dL4_dY - dL4_dY * dL3_dY ) * eLength[ 5 ];
		dNy_dZ[ 5 ] = ( dL3_dZ * dL4_dY - dL4_dZ * dL3_dY ) * eLength[ 5 ];
			  	  												   	   
 		dNz_dX[ 5 ] = ( dL3_dX * dL4_dZ - dL4_dX * dL3_dZ ) * eLength[ 5 ];
		dNz_dY[ 5 ] = ( dL3_dY * dL4_dZ - dL4_dY * dL3_dZ ) * eLength[ 5 ];
		dNz_dZ[ 5 ] = ( dL3_dZ * dL4_dZ - dL4_dZ * dL3_dZ ) * eLength[ 5 ];		
	}

    //******************************************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
                                                                  Matrix< std::complex<double> >& EleStiffMatrix, 
                                                                  Vector< std::complex<double> >& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		if( mPotentialsOn ) 
		{
            for( int j=0; j<mNumVerts; j++ ) 
            { 
                if( FixedValue.find( mNodes[ j ]->Id() ) != FixedValue.end() )
                {    
                    std::complex<double> cVoltage = FixedValue[ mNodes[ j ]->Id() ];
		    
                    for( int i=0; i<mNumDofs; i++ ) 
                    {
                        ResidualVector[ i ] -= EleStiffMatrix[ i ][ j + 3*mNumNodes ] * cVoltage;
                    }
                }
            }
		}
    }

    //******************************************************************************************************************************************
    // - Element stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );	

		Get_CurlCurl_NiNj_Matrix( StiffMatrix ); 

		if( mStabilizerOn )
		{
		    Add_Stabilization_Matrix( StiffMatrix );
		}
		
		if( mPotentialsOn ) 
		{
		    Add_Potentials_AV_Matrix( StiffMatrix ); 
		}
	}

    //******************************************************************************************************************************************
    // - Cur-Curl + NiNj contribution to the element stiffness matrix 
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::Get_CurlCurl_NiNj_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{
        // Element material properties
		double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;

		double pi = 3.14159265358979;
        double eo = 8.8541878176e-12;
	    double mo = pi * 4.0e-7;    
		
		double sg_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
		double sg_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        double ep_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double ep_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
		double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

		std::complex<double> w2cEps;

		if( mComplxFreqOn )
		{
		    std::complex<double> jUnit( 0.0, 1.0 );
		    
		    std::complex<double> sLaplace  = mComplxFreq;
		    std::complex<double> sLaplace2 = sLaplace * sLaplace;
		    
            w2cEps = sLaplace * ( sg_real + jUnit * sg_imag ) - sLaplace2 * ( ep_real + jUnit * ep_imag );		 
		}
		else
		{
            w2cEps = std::complex<double>( w2 * ep_real - wf * sg_imag, w2 * ep_imag + wf * sg_real );
		}

		std::complex<double> cMuInv = 1.0 / std::complex<double>( mu_real, mu_imag );

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_1st( Nx, Ny, Nz, cX, cY, cZ );

		// Derivatives on Gauss points
		Vector<double> dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ;

		Derivatives_On_Points( dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ );

		// Edge element matrix
		Matrix< std::complex<double> > EdgeBaseMatrix( mNumEdges, mNumEdges );

		// Edge element matrix upper diagonal 
		for( int i=0; i<mNumEdges; i++ )
		{
			for( int j=i; j<mNumEdges; j++ )
			{
				// Curl-Curl volumetric integral
				double vIntg_Curli_Curlj = mVolume * ( ( dNz_dY[ i ] - dNy_dZ[ i ] ) * ( dNz_dY[ j ] - dNy_dZ[ j ] ) + 
					                                   ( dNx_dZ[ i ] - dNz_dX[ i ] ) * ( dNx_dZ[ j ] - dNz_dX[ j ] ) +
					                                   ( dNy_dX[ i ] - dNx_dY[ i ] ) * ( dNy_dX[ j ] - dNx_dY[ j ] ) ); 

				// Ni*Nj volumetric integral
				double vIntg_Ni_Nj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Ni_Nj += W[ gp ] * ( Nx[ i ][ gp ] * Nx[ j ][ gp ] + Ny[ i ][ gp ] * Ny[ j ][ gp ] + Nz[ i ][ gp ] * Nz[ j ][ gp ] );
				}

				EdgeBaseMatrix[ i ][ j ] = ( cMuInv * vIntg_Curli_Curlj ) - ( w2cEps * vIntg_Ni_Nj );
			}
		}

		// Edge element matrix lower diagonal 
        for( int i=0; i<mNumEdges; i++ )
        {
            for( int j=i+1; j<mNumEdges; j++ )
            {
                EdgeBaseMatrix[ j ][ i ] = EdgeBaseMatrix[ i ][ j ];
            }
        }        

		// Edge element matrix inclusion into stiffness matrix
		for( int i=0; i<mNumEdges; i++ )
        {
            for( int j=0; j<mNumEdges; j++ )
            {
                StiffMatrix[ i + mNumNodes + mNumVerts ][ j + mNumNodes + mNumVerts ] = EdgeBaseMatrix[ i ][ j ];
            }
        }        
	}
    
	//******************************************************************************************************************************************
    // - Edge element stabilization term with Lagrange multipliers
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::Add_Stabilization_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
        // Stabilization constant
		std::complex<double> cStabCte;

		if( mPotentialsOn )
		{
		    cStabCte = 1.0;
		}
		else 
		{
            double eo      = 8.8541878176e-12;
            double freqc   = (*mProperties)( FREQUENCY );
		    
			double sg_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
			double sg_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

            double ep_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		    double ep_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;
		    
		    std::complex<double> cEp( ep_real - ( sg_imag / freqc ), ep_imag + ( sg_real / freqc ) );

		    cStabCte = cEp / eo;
		}

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_1st( Nx, Ny, Nz, cX, cY, cZ );

		// 1st order Lagrange base derivatives 
		double DN[ 3 ][ 4 ]; 
		
		Calculate_DN( DN );

		// [ Wi * grad( Nj ) ] - [ grad( Ni ) * Wj ]
		for( int i=0; i<mNumEdges; i++ )
		{
			for( int j=0; j<mNumVerts; j++ )
			{
				double vIntg_Wi_DNj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Wi_DNj += W[ gp ] * ( Nx[ i ][ gp ] * DN[ 0 ][ j ] + Ny[ i ][ gp ] * DN[ 1 ][ j ] + Nz[ i ][ gp ] * DN[ 2 ][ j ] );
				}

				StiffMatrix[ i + mNumNodes + mNumVerts ][ j + 2*mNumNodes ] += cStabCte * vIntg_Wi_DNj;
				StiffMatrix[ j + 2*mNumNodes ][ i + mNumNodes + mNumVerts ] += cStabCte * vIntg_Wi_DNj;
			}
		}

		// [ DNi * DNj ] 
		for( int i=0; i<mNumVerts; i++ )
        {
            for( int j=0; j<mNumVerts; j++ )
            {
				double vIntg_DNi_DNj = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );

				StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] -= cStabCte * vIntg_DNi_DNj;
            }
        } 
	}

    //******************************************************************************************************************************************
    // - AV-potentials contribution [AV], [VA], [VV] to stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::Add_Potentials_AV_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 		
		// Element material properties
		double wf = (*mProperties)( FREQUENCY );

        double w2 = wf * wf;
        double eo = 8.8541878176e-12;
	    
		double sgm_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
	    double sgm_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        double eps_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double eps_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		std::complex<double> nw2cEps;

		if( mComplxFreqOn )
		{
		    std::complex<double> jUnit( 0.0, 1.0 );
		    
		    std::complex<double> sLaplace  = mComplxFreq;
		    std::complex<double> sLaplace2 = sLaplace * sLaplace;
		    
            nw2cEps = sLaplace2 * ( eps_real + jUnit * eps_imag ) - sLaplace * ( sgm_real + jUnit * sgm_imag );         	 
		}
		else
		{
            nw2cEps = std::complex<double>( - w2 * eps_real + wf * sgm_imag, - w2 * eps_imag - wf * sgm_real );
		}

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_1st( Nx, Ny, Nz, cX, cY, cZ );

		// 1st order Lagrange base derivatives 
		double DN[ 3 ][ 4 ]; 
		
		Calculate_DN( DN );

		// [ Wi * grad( Nj ) ] - [ grad( Ni ) * Wj ] 
		for( int i=0; i<mNumEdges; i++ )
		{
			for( int j=0; j<mNumVerts; j++ )
			{
				double vIntg_Wi_DNj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Wi_DNj += W[ gp ] * ( Nx[ i ][ gp ] * DN[ 0 ][ j ] + Ny[ i ][ gp ] * DN[ 1 ][ j ] + Nz[ i ][ gp ] * DN[ 2 ][ j ] );
				}

				StiffMatrix[ i + mNumNodes + mNumVerts ][ j + 3*mNumNodes ] = nw2cEps * vIntg_Wi_DNj;
				StiffMatrix[ j + 3*mNumNodes ][ i + mNumNodes + mNumVerts ] = nw2cEps * vIntg_Wi_DNj;
			}
		}

		// [ grad( Ni ) * grad( Nj ) ]
		for( int i=0; i<mNumVerts; i++ )
        {
            for( int j=0; j<mNumVerts; j++ ) 
        	{
				StiffMatrix[ i + 3*mNumNodes ][ j + 3*mNumNodes ] = nw2cEps * mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );
        	}
        }
	}

    //******************************************************************************************************************************************
    // - Rotational on nodes
    //******************************************************************************************************************************************
	void VolumeElement_1eg_FullWave::Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes )
	{
        // Derivatives on element
		Vector<double> dNx_dX, dNx_dY, dNx_dZ; 
		Vector<double> dNy_dX, dNy_dY, dNy_dZ; 
		Vector<double> dNz_dX, dNz_dY, dNz_dZ;

		Derivatives_On_Points( dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ );

        // Ft on edges
		Vector< std::complex<double> > cFt( mNumEdges );

		if( mPotentialsOn )
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }
		}
		else
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }
		}

		// Rotational 
		Vector< std::complex<double> > cRot( 3, std::complex<double>( 0.0, 0.0 ) );

		for( int e=0; e<mNumEdges; e++ )
		{
		    cRot[ 0 ] += cFt[ e ] * ( dNz_dY[ e ] - dNy_dZ[ e ] );
		    cRot[ 1 ] += cFt[ e ] * ( dNx_dZ[ e ] - dNz_dX[ e ] );
		    cRot[ 2 ] += cFt[ e ] * ( dNy_dX[ e ] - dNx_dY[ e ] );
		}	

		// Rotational on nodes
		rot_OnNodes.resize( mNumNodes );

		for( int n=0; n<mNumNodes; n++ )
		{
		    rot_OnNodes[ n ] = cRot;
		}
	}
	  
    //******************************************************************************************************************************************
    // - Rotational on Gauss points
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs )
	{
		// GiD internal Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );		

        // Derivatives on element
		Vector<double> dNx_dX, dNx_dY, dNx_dZ; 
		Vector<double> dNy_dX, dNy_dY, dNy_dZ; 
		Vector<double> dNz_dX, dNz_dY, dNz_dZ;

		Derivatives_On_Points( dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ );

        // Ft on edges
		Vector< std::complex<double> > cFt( mNumEdges );

		if( mPotentialsOn )
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }
		}
		else
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }
		}

		// Rotational on element
		Vector< std::complex<double> > cRot( 3, std::complex<double>( 0.0, 0.0 ) );

		for ( int e=0; e<mNumEdges; e++ )
		{
		    cRot[ 0 ] += cFt[ e ] * ( dNz_dY[ e ] - dNy_dZ[ e ] );
		    cRot[ 1 ] += cFt[ e ] * ( dNx_dZ[ e ] - dNz_dX[ e ] );
		    cRot[ 2 ] += cFt[ e ] * ( dNy_dX[ e ] - dNx_dY[ e ] );
		}	

		// Rotational on Gauss points
		rot_OnGPs.resize( numGaussPoints );

		for ( int gp=0; gp<numGaussPoints; gp++ )
		{
		    rot_OnGPs[ gp ] = cRot;
		}
	}

    //******************************************************************************************************************************************
    // - Electric field on nodes
    //******************************************************************************************************************************************
	void VolumeElement_1eg_FullWave::Calculate_E_field_OnNodes( cVector2Type& Ef_OnNodes )
	{
		// Natural coordinates of nodes
		Vector<double> nX( mNumNodes );
		Vector<double> nY( mNumNodes );
		Vector<double> nZ( mNumNodes );

		// Corners 
		nX[ 0 ] = 0.0; nY[ 0 ] = 0.0; nZ[ 0 ] = 0.0;
		nX[ 1 ] = 1.0; nY[ 1 ] = 0.0; nZ[ 1 ] = 0.0;
		nX[ 2 ] = 0.0; nY[ 2 ] = 1.0; nZ[ 2 ] = 0.0;
		nX[ 3 ] = 0.0; nY[ 3 ] = 0.0; nZ[ 3 ] = 1.0;
		
		// Edges 
		nX[ 4 ] = 0.5; nY[ 4 ] = 0.0; nZ[ 4 ] = 0.0; 
		nX[ 5 ] = 0.5; nY[ 5 ] = 0.5; nZ[ 5 ] = 0.0; 
		nX[ 6 ] = 0.0; nY[ 6 ] = 0.5; nZ[ 6 ] = 0.0; 
		nX[ 7 ] = 0.0; nY[ 7 ] = 0.0; nZ[ 7 ] = 0.5; 
		nX[ 8 ] = 0.5; nY[ 8 ] = 0.0; nZ[ 8 ] = 0.5;
		nX[ 9 ] = 0.0; nY[ 9 ] = 0.5; nZ[ 9 ] = 0.5; 

        // 1st order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_1st( Nx, Ny, Nz, nX, nY, nZ );
	   
        // Ft on edges
		Vector< std::complex<double> > cFt( mNumEdges );

		if( mPotentialsOn )
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }
		}
		else
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }
		}

		// Field on nodes
		cVector2Type cF_OnNodes( mNumNodes );

		for( int n=0; n<mNumNodes; n++ )
		{
		    Vector< std::complex<double> > cF_OnNode( 3, std::complex<double>( 0.0, 0.0 ) );

			for( int e=0; e<mNumEdges; e++ )
			{
		        cF_OnNode[ 0 ] += cFt[ e ] * Nx[ e ][ n ];
		        cF_OnNode[ 1 ] += cFt[ e ] * Ny[ e ][ n ];
		        cF_OnNode[ 2 ] += cFt[ e ] * Nz[ e ][ n ];
			}

			cF_OnNodes[ n ] = cF_OnNode;
		}

		// Electric field on nodes
		if( mPotentialsOn )
		{
            double DN[ 3 ][ 4 ]; Calculate_DN( DN );

			std::complex<double> jw;

			if( mComplxFreqOn )
			{
                jw = mComplxFreq;
			}
			else
			{
			    jw = std::complex<double>( 0.0, (*mProperties)( FREQUENCY ) );
			}

			Vector< std::complex<double> > cGradVs( 3, std::complex<double>( 0.0, 0.0 ) ); 

			for( int i=0; i<mNumVerts; i++ )
            {
                cGradVs[ 0 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 0 ][ i ];
                cGradVs[ 1 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 1 ][ i ];
                cGradVs[ 2 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 2 ][ i ];
            }

            // E = jw * ( A + grad( V ) )
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        cF_OnNodes[ n ][ 0 ] = jw * ( cF_OnNodes[ n ][ 0 ] + cGradVs[ 0 ] );
		    	cF_OnNodes[ n ][ 1 ] = jw * ( cF_OnNodes[ n ][ 1 ] + cGradVs[ 1 ] );
		    	cF_OnNodes[ n ][ 2 ] = jw * ( cF_OnNodes[ n ][ 2 ] + cGradVs[ 2 ] );
		    }			
		}

		Ef_OnNodes = cF_OnNodes;
	}

    //******************************************************************************************************************************************
    // - Electric field on Gauss points
    //******************************************************************************************************************************************
    void VolumeElement_1eg_FullWave::Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs )
	{
		// GiD internal Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );

        // 1st order Whitney base functions
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_1st( Nx, Ny, Nz, cX, cY, cZ );
	   
        // Ft on edges
		Vector< std::complex<double> > cFt( mNumEdges );

		if( mPotentialsOn )
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cAy, *mNodes[ e + mNumVerts ] );
		    }
		}
		else
		{
		    for( int e=0; e<mNumEdges; e++ )
		    {
		        cFt[ e ] = (*mProperties)( cEy, *mNodes[ e + mNumVerts ] );
		    }
		}

		// Field on Gauss points
		cVector2Type cF_OnGPs( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
		    Vector< std::complex<double> > cF_OnGP( 3, std::complex<double>( 0.0, 0.0 ) );

			for( int e=0; e<mNumEdges; e++ )
			{
		        cF_OnGP[ 0 ] += cFt[ e ] * Nx[ e ][ gp ];
		        cF_OnGP[ 1 ] += cFt[ e ] * Ny[ e ][ gp ];
		        cF_OnGP[ 2 ] += cFt[ e ] * Nz[ e ][ gp ];
			}

			cF_OnGPs[ gp ] = cF_OnGP;
		}

		// Electric field on Gauss points
		if( mPotentialsOn )
		{
            double DN[ 3 ][ 4 ]; Calculate_DN( DN );

			std::complex<double> jw;

			if( mComplxFreqOn )
			{
                jw = mComplxFreq;
			}
			else
			{
			    jw = std::complex<double>( 0.0, (*mProperties)( FREQUENCY ) );
			}

			Vector< std::complex<double> > cGradVs( 3, std::complex<double>( 0.0, 0.0 ) ); 

			for( int i=0; i<mNumVerts; i++ )
            {
                cGradVs[ 0 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 0 ][ i ];
                cGradVs[ 1 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 1 ][ i ];
                cGradVs[ 2 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 2 ][ i ];
            }

            // E = jw * ( A + grad( V ) )
		    for( int gp=0; gp<numGaussPoints; gp++ )
		    {
		        cF_OnGPs[ gp ][ 0 ] = jw * ( cF_OnGPs[ gp ][ 0 ] + cGradVs[ 0 ] );
		    	cF_OnGPs[ gp ][ 1 ] = jw * ( cF_OnGPs[ gp ][ 1 ] + cGradVs[ 1 ] );
		    	cF_OnGPs[ gp ][ 2 ] = jw * ( cF_OnGPs[ gp ][ 2 ] + cGradVs[ 2 ] );
		    }			
		}

		Ef_OnGPs = cF_OnGPs;
	}

    //******************************************************************************************************************************************
    // - Inner GiD Gauss points 
    //******************************************************************************************************************************************
    int VolumeElement_1eg_FullWave::GetInnerGiDGaussPoints( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, int numResultsOnGPs )
	{
		if( numResultsOnGPs <= 1 ) 
		{
			GaussPoints3D_InnerGiD_0rd( cX, cY, cZ );
		}
		else    
		{
			GaussPoints3D_InnerGiD_1st( cX, cY, cZ );
		}
	
        return cX.size();
	}
} 