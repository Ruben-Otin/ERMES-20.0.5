
#include "VolumeElement_1eg_ColdPlasma.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
    //******************************************************************************************************************************************
    // - Global indexes of the DOFs
    //******************************************************************************************************************************************
    void VolumeElement_1eg_ColdPlasma::GetEquationIdVector( Vector<int>& EquationId )
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
    double VolumeElement_1eg_ColdPlasma::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //******************************************************************************************************************************************
    //* - Set element as IHL or plasma 
    //******************************************************************************************************************************************	
    void VolumeElement_1eg_ColdPlasma::SetPlasmaIHL() 
    { 
        // Material properties
        double sg_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
        double sg_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        double ep_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY );
        double ep_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY );

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY );
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY );

        // If the material is vaccum then it is considered as plasma. 
        // We only consider the first one written in the plasma file *-11.dat
        // This allows the use of different plasmas volumes (for visualization). 
        if( ( sg_real == 0.0 ) && ( sg_imag == 0.0 ) &&
            ( ep_real == 1.0 ) && ( ep_imag == 0.0 ) && 
            ( mu_real == 1.0 ) && ( mu_imag == 0.0 )  )
        {
            mIs_IHL = false;        
        }
        else
        {
            mIs_IHL = true;        
        }
    }

    //******************************************************************************************************************************************
    // - Lagrangian 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //******************************************************************************************************************************************
    void VolumeElement_1eg_ColdPlasma::Calculate_DN( double DN[ 3 ][ 4 ] )
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
    void VolumeElement_1eg_ColdPlasma::Calculate_Edge_Lengths( Vector<double>& eLength )
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
	void VolumeElement_1eg_ColdPlasma::Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz,
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
    void VolumeElement_1eg_ColdPlasma::Derivatives_On_Points( Vector<double>& dNx_dX, Vector<double>& dNx_dY, Vector<double>& dNx_dZ,
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
    //* - Permittivity tensor on nodes 
    //******************************************************************************************************************************************
    void VolumeElement_1eg_ColdPlasma::PermittivityTensor_OnNodes( Vector< Matrix< std::complex<double> > >& TEn )
    {
        // Permittivity tensor on vertice nodes
        TEn.resize( mNumVerts );
        
        // Constants
        double wf = (*mProperties)( FREQUENCY );
        double eo = 8.8541878176e-12;
        
        std::complex<double> cZero( 0.0, 0.0 );
        std::complex<double> cUnit( 0.0, 1.0 );

        if( mIs_IHL )
        {
            for( int i=0; i<TEn.size(); i++ ) 
            { 
                TEn[ i ].Resize( 3, 3 );

                double sgm_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
                double sgm_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

                double eps_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
                double eps_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

                std::complex<double> cEps( eps_real - ( sgm_imag / wf ), eps_imag + ( sgm_real / wf ) );

                // IHL permittivity tensor 
                Matrix< std::complex<double> > K( 3, 3 );

                K[ 0 ][ 0 ] = cEps ; K[ 0 ][ 1 ] = cZero; K[ 0 ][ 2 ] = cZero; 
                K[ 1 ][ 0 ] = cZero; K[ 1 ][ 1 ] = cEps ; K[ 1 ][ 2 ] = cZero; 
                K[ 2 ][ 0 ] = cZero; K[ 2 ][ 1 ] = cZero; K[ 2 ][ 2 ] = cEps ;

                TEn[ i ] = K;
            }
        }
        else
        {
            for( int i=0; i<TEn.size(); i++ ) 
            { 
                mpColdPlasma->Get_PermittivityTensor_OnNode( TEn[ i ], mNodes[ i ], wf );                
            }
        }
    }

    //******************************************************************************************************************************************
    //* - Permittivity tensor on Gauss points ( extrapolating TEn to GPs with N )
    //******************************************************************************************************************************************
    void VolumeElement_1eg_ColdPlasma::PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, Matrix<double>& N )
    {
        // Permittivity tensor on vertice nodes
        Vector< Matrix< std::complex<double> > > TEn;

        PermittivityTensor_OnNodes( TEn );
        
        // Number of Gauss points
        int numGaussPoints = N[ 0 ].size();

        // Initializing vector of 3x3 permittivity tensors 
        TEgp.resize( numGaussPoints );

        // Permittivity tensor on Gauss points  
        for( int gp=0; gp<numGaussPoints; gp++ )
        {
            TEgp[ gp ].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );

            for( int ki=0; ki<3; ki++ )
            {
                for( int kj=0; kj<3; kj++ )
                {
                    for( int i=0; i<TEn.size(); i++ ) 
                    {
                        TEgp[ gp ][ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * N[ i ][ gp ];
                    }
                }
            }
        }
    }

    //******************************************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //******************************************************************************************************************************************
    void VolumeElement_1eg_ColdPlasma::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
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
	//* - Stiffness matrix
	//******************************************************************************************************************************************
	void VolumeElement_1eg_ColdPlasma::GetStiffnessMatrix( Matrix<std::complex<double> >& StiffMatrix )
	{
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );	

		Get_Edge_Formulation_Matrix_EA( StiffMatrix ); 
	}

    //******************************************************************************************************************************************
    //* - Stiffness matrix ( Ef and AV formulation )
    //******************************************************************************************************************************************	
    void VolumeElement_1eg_ColdPlasma::Get_Edge_Formulation_Matrix_EA( Matrix< std::complex<double> >& StiffMatrix )
    {  
        // Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        double pi = 3.14159265358979;
        double mo = pi * 4.0e-7;

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> cInvMu = 1.0 / cMu;
		
		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order Whitney base functions on Gauss points
		Matrix<double> Nx, Ny, Nz;

		Whitney3D_Ni_1st( Nx, Ny, Nz, cX, cY, cZ );

		// Derivatives on Gauss points
		Vector<double> dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ;

		Derivatives_On_Points( dNx_dX, dNx_dY, dNx_dZ, dNy_dX, dNy_dY, dNy_dZ, dNz_dX, dNz_dY, dNz_dZ );

        // Lagrange 1st order basis functions for the permittivity tensor 
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );     

        // 1st order Lagrange base derivatives 
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        // Permittivity tensor on Gauss points
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, N );

		// Edge element matrix
		Matrix< std::complex<double> > EdgeBaseMatrix( mNumEdges, mNumEdges );

		for( int i=0; i<mNumEdges; i++ )
		{
			for( int j=0; j<mNumEdges; j++ )
			{
				// Curl-Curl volumetric integral
				double vIntg_Curli_Curlj = mVolume * ( ( dNz_dY[ i ] - dNy_dZ[ i ] ) * ( dNz_dY[ j ] - dNy_dZ[ j ] ) + 
					                                   ( dNx_dZ[ i ] - dNz_dX[ i ] ) * ( dNx_dZ[ j ] - dNz_dX[ j ] ) +
					                                   ( dNy_dX[ i ] - dNx_dY[ i ] ) * ( dNy_dX[ j ] - dNx_dY[ j ] ) ); 

				// Ni*TEp*Nj volumetric integral
				std::complex<double> vIntg_Ni_TEp_Nj( 0.0, 0.0 );

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Ni_TEp_Nj += W[ gp ] * Nx[ i ][ gp ] * ( TEgp[ gp ][ 0 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 0 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 0 ][ 2 ] * Nz[ j ][ gp ] );
                    vIntg_Ni_TEp_Nj += W[ gp ] * Ny[ i ][ gp ] * ( TEgp[ gp ][ 1 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 1 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 1 ][ 2 ] * Nz[ j ][ gp ] );
                    vIntg_Ni_TEp_Nj += W[ gp ] * Nz[ i ][ gp ] * ( TEgp[ gp ][ 2 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 2 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 2 ][ 2 ] * Nz[ j ][ gp ] );
				}

				EdgeBaseMatrix[ i ][ j ] = ( cInvMu * vIntg_Curli_Curlj ) - ( w2 * vIntg_Ni_TEp_Nj );
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

        // Stabilization matrix
        if( mStabilizerOn )
        {
            Add_Edge_Stabilizatn_Matrix_EA( StiffMatrix, TEgp, W, Nx, Ny, Nz );
        }
        
        // [AV], [VA], [VV]
        if( mPotentialsOn )
        {
            Add_Edge_Potentials_Matrix_AVV( StiffMatrix, TEgp, W, Nx, Ny, Nz );            
        }        
    }

    //******************************************************************************************************************************************
    //* - Edge element stabilization term with Lagrange multipliers
    //******************************************************************************************************************************************	
    void VolumeElement_1eg_ColdPlasma::Add_Edge_Stabilizatn_Matrix_EA( Matrix< std::complex<double> >& StiffMatrix, Vector< Matrix< std::complex<double> > >& TEgp, Vector<double>& W, Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz )
    {
        // Number of Gauss points
        int numGaussPoints = W.size();

        // 1st order Lagrange base derivatives 
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        if( mPotentialsOn )
        {
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
		    
		    		StiffMatrix[ i + mNumNodes + mNumVerts ][ j + 2*mNumNodes ] += vIntg_Wi_DNj;
		    		StiffMatrix[ j + 2*mNumNodes ][ i + mNumNodes + mNumVerts ] += vIntg_Wi_DNj;
		    	}
		    }

		    // [ DNi * DNj ] 
		    for( int i=0; i<mNumVerts; i++ )
            {
                for( int j=0; j<mNumVerts; j++ )
                {
		    		double vIntg_DNi_DNj = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );
		    
		    		StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] -= vIntg_DNi_DNj;
                }
            } 
        }
        else
        {
		    // Stabilization constant
            double StabCte = 1.0 / 8.8541878176e-12;
            
            // [ TEp * Wi * grad( Nj ) ]
		    for( int i=0; i<mNumEdges; i++ )
		    {
		    	for( int j=0; j<mNumVerts; j++ )
		    	{
		    		std::complex<double> vIntg_Wi_TEp_DNj( 0.0, 0.0 );
            
		    		for( int gp=0; gp<numGaussPoints; gp++ )
		    		{
		    			vIntg_Wi_TEp_DNj += W[ gp ] * DN[ 0 ][ j ] * ( TEgp[ gp ][ 0 ][ 0 ] * Nx[ i ][ gp ] + TEgp[ gp ][ 0 ][ 1 ] * Ny[ i ][ gp ] + TEgp[ gp ][ 0 ][ 2 ] * Nz[ i ][ gp ] );
                        vIntg_Wi_TEp_DNj += W[ gp ] * DN[ 1 ][ j ] * ( TEgp[ gp ][ 1 ][ 0 ] * Nx[ i ][ gp ] + TEgp[ gp ][ 1 ][ 1 ] * Ny[ i ][ gp ] + TEgp[ gp ][ 1 ][ 2 ] * Nz[ i ][ gp ] );
                        vIntg_Wi_TEp_DNj += W[ gp ] * DN[ 2 ][ j ] * ( TEgp[ gp ][ 2 ][ 0 ] * Nx[ i ][ gp ] + TEgp[ gp ][ 2 ][ 1 ] * Ny[ i ][ gp ] + TEgp[ gp ][ 2 ][ 2 ] * Nz[ i ][ gp ] );
		    		}
            
		    		StiffMatrix[ i + mNumNodes + mNumVerts ][ j + 2*mNumNodes ] += StabCte * vIntg_Wi_TEp_DNj;
		    	}
		    }
            
		    // [ grad( Ni ) * TEp * Wj ]
		    for( int i=0; i<mNumVerts; i++ )
		    {
		    	for( int j=0; j<mNumEdges; j++ )
		    	{
		    		std::complex<double> vIntg_DNi_TEp_Wj( 0.0, 0.0 );
            
		    		for( int gp=0; gp<numGaussPoints; gp++ )
		    		{
		    			vIntg_DNi_TEp_Wj += W[ gp ] * DN[ 0 ][ i ] * ( TEgp[ gp ][ 0 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 0 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 0 ][ 2 ] * Nz[ j ][ gp ] );
                        vIntg_DNi_TEp_Wj += W[ gp ] * DN[ 1 ][ i ] * ( TEgp[ gp ][ 1 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 1 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 1 ][ 2 ] * Nz[ j ][ gp ] );
                        vIntg_DNi_TEp_Wj += W[ gp ] * DN[ 2 ][ i ] * ( TEgp[ gp ][ 2 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 2 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 2 ][ 2 ] * Nz[ j ][ gp ] );
		    		}
            
		    		StiffMatrix[ i + 2*mNumNodes ][ j + mNumNodes + mNumVerts ] += StabCte * vIntg_DNi_TEp_Wj;
		    	}
		    }

		    // [ DNi * DNj ] 
		    for( int i=0; i<mNumVerts; i++ )
            {
                for( int j=0; j<mNumVerts; j++ )
                {
		    		double vIntg_DNi_DNj = mVolume * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );
		    
		    		StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] -= StabCte * vIntg_DNi_DNj;
                }
            } 
        }
    }

    //******************************************************************************************************************************************
    // - AV-potentials contribution [AV], [VA], [VV] to stiffness matrix
    //******************************************************************************************************************************************
    void VolumeElement_1eg_ColdPlasma::Add_Edge_Potentials_Matrix_AVV( Matrix< std::complex<double> >& StiffMatrix, Vector< Matrix< std::complex<double> > >& TEgp, Vector<double>& W, Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz ) 
	{ 		
	    // Frequency
		double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;

        // Number of Gauss points
        int numGaussPoints = W.size();

        // 1st order Lagrange base derivatives 
		double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        // [ Wi * TEp * grad( Nj ) ]
		for( int i=0; i<mNumEdges; i++ )
		{
			for( int j=0; j<mNumVerts; j++ )
			{
				std::complex<double> vIntg_Wi_TEp_DNj( 0.0, 0.0 );
        
				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_Wi_TEp_DNj += W[ gp ] * Nx[ i ][ gp ] * ( TEgp[ gp ][ 0 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 0 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 0 ][ 2 ] * DN[ 2 ][ j ] );
                    vIntg_Wi_TEp_DNj += W[ gp ] * Ny[ i ][ gp ] * ( TEgp[ gp ][ 1 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 1 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 1 ][ 2 ] * DN[ 2 ][ j ] );
                    vIntg_Wi_TEp_DNj += W[ gp ] * Nz[ i ][ gp ] * ( TEgp[ gp ][ 2 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 2 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 2 ][ 2 ] * DN[ 2 ][ j ] );
				}
        
				StiffMatrix[ i + mNumNodes + mNumVerts ][ j + 3*mNumNodes ] -= w2 * vIntg_Wi_TEp_DNj;
			}
		}
        
		// [ grad( Ni ) * TEp * Wj ]
		for( int i=0; i<mNumVerts; i++ )
		{
			for( int j=0; j<mNumEdges; j++ )
			{
				std::complex<double> vIntg_DNi_TEp_Wj( 0.0, 0.0 );
        
				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_DNi_TEp_Wj += W[ gp ] * DN[ 0 ][ i ] * ( TEgp[ gp ][ 0 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 0 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 0 ][ 2 ] * Nz[ j ][ gp ] );
                    vIntg_DNi_TEp_Wj += W[ gp ] * DN[ 1 ][ i ] * ( TEgp[ gp ][ 1 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 1 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 1 ][ 2 ] * Nz[ j ][ gp ] );
                    vIntg_DNi_TEp_Wj += W[ gp ] * DN[ 2 ][ i ] * ( TEgp[ gp ][ 2 ][ 0 ] * Nx[ j ][ gp ] + TEgp[ gp ][ 2 ][ 1 ] * Ny[ j ][ gp ] + TEgp[ gp ][ 2 ][ 2 ] * Nz[ j ][ gp ] );
				}
        
				StiffMatrix[ i + 3*mNumNodes ][ j + mNumNodes + mNumVerts ] -= w2 * vIntg_DNi_TEp_Wj;
			}
		}

		// [ grad( Ni ) * TEp * grad( Nj ) ]
		for( int i=0; i<mNumVerts; i++ )
        {
            for( int j=0; j<mNumVerts; j++ ) 
        	{
                std::complex<double> vIntg_DNi_TEp_DNj( 0.0, 0.0 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					vIntg_DNi_TEp_DNj += W[ gp ] * DN[ 0 ][ i ] * ( TEgp[ gp ][ 0 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 0 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 0 ][ 2 ] * DN[ 2 ][ j ] );
                    vIntg_DNi_TEp_DNj += W[ gp ] * DN[ 1 ][ i ] * ( TEgp[ gp ][ 1 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 1 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 1 ][ 2 ] * DN[ 2 ][ j ] );
                    vIntg_DNi_TEp_DNj += W[ gp ] * DN[ 2 ][ i ] * ( TEgp[ gp ][ 2 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 2 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 2 ][ 2 ] * DN[ 2 ][ j ] );
				}				

                StiffMatrix[ i + 3*mNumNodes ][ j + 3*mNumNodes ] -= w2 * vIntg_DNi_TEp_DNj;
        	}
        }
	}
} 