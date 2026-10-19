
#include "VolumeElement_3sb_ColdPlasma.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
    //********************************************************************************************************************************
    // - Global index of the DOFs
    //********************************************************************************************************************************	
    void VolumeElement_3sb_ColdPlasma::GetEquationIdVector( Vector<int>& EquationId )
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

    //********************************************************************************************************************************
    // - Volume of the element
    //********************************************************************************************************************************
    double VolumeElement_3sb_ColdPlasma::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //********************************************************************************************************************************
    //* - Set element as IHL or plasma 
    //********************************************************************************************************************************	
    void VolumeElement_3sb_ColdPlasma::SetPlasmaIHL() 
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

    //********************************************************************************************************************************
    // - 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //********************************************************************************************************************************
    void VolumeElement_3sb_ColdPlasma::Calculate_DN( double DN[ 3 ][ 4 ] )
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

    //********************************************************************************************************************************
    // - dN/dx, dN/dy, dN/dz on cX, cY, cZ points
    //********************************************************************************************************************************
    void VolumeElement_3sb_ColdPlasma::Derivatives_On_Points( Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz,
													          Vector<double>&   cX, Vector<double>&   cY, Vector<double>&   cZ )
    {
		int nPoints = cX.size();

		dNdx.Resize( mNumNodes, nPoints );
		dNdy.Resize( mNumNodes, nPoints );
		dNdz.Resize( mNumNodes, nPoints );
		
		double DN[ 3 ][ 4 ];

        Calculate_DN( DN );

		double dL1dx, dL2dx, dL3dx, dL4dx;
        double dL1dy, dL2dy, dL3dy, dL4dy;
		double dL1dz, dL2dz, dL3dz, dL4dz;

        dL1dx = DN[ 0 ][ 0 ]; dL2dx = DN[ 0 ][ 1 ]; dL3dx = DN[ 0 ][ 2 ]; dL4dx = DN[ 0 ][ 3 ];
        dL1dy = DN[ 1 ][ 0 ]; dL2dy = DN[ 1 ][ 1 ]; dL3dy = DN[ 1 ][ 2 ]; dL4dy = DN[ 1 ][ 3 ];
		dL1dz = DN[ 2 ][ 0 ]; dL2dz = DN[ 2 ][ 1 ]; dL3dz = DN[ 2 ][ 2 ]; dL4dz = DN[ 2 ][ 3 ];

		// dNi/dx ,dNi/dy, dNi/dz
        for( int p=0; p<nPoints; p++ ) 
        {
			// Natural coordinates
			double L2 = cX[ p ];
            double L3 = cY[ p ];
			double L4 = cZ[ p ];
            double L1 = 1.0 - L2 - L3 - L4;

			// Nodes [ 0 1 2 3 ] 
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

			// Face bubbles [ 0 1 2 ]
			dNdx[ 4 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL1dx ) ) + ( L1 * L2 * dL3dx * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL2dx * L3 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dx * L2 * L3 * ( 4.0 * L1 - 1.0 ) ) );
			dNdy[ 4 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL1dy ) ) + ( L1 * L2 * dL3dy * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL2dy * L3 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dy * L2 * L3 * ( 4.0 * L1 - 1.0 ) ) );
			dNdz[ 4 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL1dz ) ) + ( L1 * L2 * dL3dz * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL2dz * L3 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dz * L2 * L3 * ( 4.0 * L1 - 1.0 ) ) );
				 	 			                                                                                                                                                         
			dNdx[ 5 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL2dx ) ) + ( L1 * L2 * dL3dx * ( 4.0 * L2 - 1.0 ) ) + ( L1 * dL2dx * L3 * ( 4.0 * L2 - 1.0 ) ) + ( dL1dx * L2 * L3 * ( 4.0 * L2 - 1.0 ) ) );
			dNdy[ 5 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL2dy ) ) + ( L1 * L2 * dL3dy * ( 4.0 * L2 - 1.0 ) ) + ( L1 * dL2dy * L3 * ( 4.0 * L2 - 1.0 ) ) + ( dL1dy * L2 * L3 * ( 4.0 * L2 - 1.0 ) ) );
			dNdz[ 5 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL2dz ) ) + ( L1 * L2 * dL3dz * ( 4.0 * L2 - 1.0 ) ) + ( L1 * dL2dz * L3 * ( 4.0 * L2 - 1.0 ) ) + ( dL1dz * L2 * L3 * ( 4.0 * L2 - 1.0 ) ) );
				 	 			                                                                                                                                                         
			dNdx[ 6 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL3dx ) ) + ( L1 * L2 * dL3dx * ( 4.0 * L3 - 1.0 ) ) + ( L1 * dL2dx * L3 * ( 4.0 * L3 - 1.0 ) ) + ( dL1dx * L2 * L3 * ( 4.0 * L3 - 1.0 ) ) );
			dNdy[ 6 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL3dy ) ) + ( L1 * L2 * dL3dy * ( 4.0 * L3 - 1.0 ) ) + ( L1 * dL2dy * L3 * ( 4.0 * L3 - 1.0 ) ) + ( dL1dy * L2 * L3 * ( 4.0 * L3 - 1.0 ) ) );
			dNdz[ 6 ][ p ] = 32.0 * ( ( L1 * L2 * L3 * ( 4.0 * dL3dz ) ) + ( L1 * L2 * dL3dz * ( 4.0 * L3 - 1.0 ) ) + ( L1 * dL2dz * L3 * ( 4.0 * L3 - 1.0 ) ) + ( dL1dz * L2 * L3 * ( 4.0 * L3 - 1.0 ) ) );
                                                                                                                                                                                          
			// Face bubbles [ 0 1 3 ]                                                                                                                                                     
			dNdx[ 7 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL1dx ) ) + ( L1 * L2 * dL4dx * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL2dx * L4 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dx * L2 * L4 * ( 4.0 * L1 - 1.0 ) ) );
			dNdy[ 7 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL1dy ) ) + ( L1 * L2 * dL4dy * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL2dy * L4 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dy * L2 * L4 * ( 4.0 * L1 - 1.0 ) ) );
			dNdz[ 7 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL1dz ) ) + ( L1 * L2 * dL4dz * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL2dz * L4 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dz * L2 * L4 * ( 4.0 * L1 - 1.0 ) ) );
				 	 			 		                                                                                                                                                  
			dNdx[ 8 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL2dx ) ) + ( L1 * L2 * dL4dx * ( 4.0 * L2 - 1.0 ) ) + ( L1 * dL2dx * L4 * ( 4.0 * L2 - 1.0 ) ) + ( dL1dx * L2 * L4 * ( 4.0 * L2 - 1.0 ) ) );
			dNdy[ 8 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL2dy ) ) + ( L1 * L2 * dL4dy * ( 4.0 * L2 - 1.0 ) ) + ( L1 * dL2dy * L4 * ( 4.0 * L2 - 1.0 ) ) + ( dL1dy * L2 * L4 * ( 4.0 * L2 - 1.0 ) ) );
			dNdz[ 8 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL2dz ) ) + ( L1 * L2 * dL4dz * ( 4.0 * L2 - 1.0 ) ) + ( L1 * dL2dz * L4 * ( 4.0 * L2 - 1.0 ) ) + ( dL1dz * L2 * L4 * ( 4.0 * L2 - 1.0 ) ) );
				 	 			 			                                                                                                                                               
			dNdx[ 9 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL4dx ) ) + ( L1 * L2 * dL4dx * ( 4.0 * L4 - 1.0 ) ) + ( L1 * dL2dx * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL1dx * L2 * L4 * ( 4.0 * L4 - 1.0 ) ) );
			dNdy[ 9 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL4dy ) ) + ( L1 * L2 * dL4dy * ( 4.0 * L4 - 1.0 ) ) + ( L1 * dL2dy * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL1dy * L2 * L4 * ( 4.0 * L4 - 1.0 ) ) );
			dNdz[ 9 ][ p ] = 32.0 * ( ( L1 * L2 * L4 * ( 4.0 * dL4dz ) ) + ( L1 * L2 * dL4dz * ( 4.0 * L4 - 1.0 ) ) + ( L1 * dL2dz * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL1dz * L2 * L4 * ( 4.0 * L4 - 1.0 ) ) );

			// Face bubbles [ 1 2 3 ]
			dNdx[ 10 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL2dx ) ) + ( L2 * L3 * dL4dx * ( 4.0 * L2 - 1.0 ) ) + ( L2 * dL3dx * L4 * ( 4.0 * L2 - 1.0 ) ) + ( dL2dx * L3 * L4 * ( 4.0 * L2 - 1.0 ) ) );
			dNdy[ 10 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL2dy ) ) + ( L2 * L3 * dL4dy * ( 4.0 * L2 - 1.0 ) ) + ( L2 * dL3dy * L4 * ( 4.0 * L2 - 1.0 ) ) + ( dL2dy * L3 * L4 * ( 4.0 * L2 - 1.0 ) ) );
			dNdz[ 10 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL2dz ) ) + ( L2 * L3 * dL4dz * ( 4.0 * L2 - 1.0 ) ) + ( L2 * dL3dz * L4 * ( 4.0 * L2 - 1.0 ) ) + ( dL2dz * L3 * L4 * ( 4.0 * L2 - 1.0 ) ) );
							 		                                                                                                                                                         
			dNdx[ 11 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL3dx ) ) + ( L2 * L3 * dL4dx * ( 4.0 * L3 - 1.0 ) ) + ( L2 * dL3dx * L4 * ( 4.0 * L3 - 1.0 ) ) + ( dL2dx * L3 * L4 * ( 4.0 * L3 - 1.0 ) ) );
			dNdy[ 11 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL3dy ) ) + ( L2 * L3 * dL4dy * ( 4.0 * L3 - 1.0 ) ) + ( L2 * dL3dy * L4 * ( 4.0 * L3 - 1.0 ) ) + ( dL2dy * L3 * L4 * ( 4.0 * L3 - 1.0 ) ) );
			dNdz[ 11 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL3dz ) ) + ( L2 * L3 * dL4dz * ( 4.0 * L3 - 1.0 ) ) + ( L2 * dL3dz * L4 * ( 4.0 * L3 - 1.0 ) ) + ( dL2dz * L3 * L4 * ( 4.0 * L3 - 1.0 ) ) );
							 		                                                                                                                                                         
			dNdx[ 12 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL4dx ) ) + ( L2 * L3 * dL4dx * ( 4.0 * L4 - 1.0 ) ) + ( L2 * dL3dx * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL2dx * L3 * L4 * ( 4.0 * L4 - 1.0 ) ) );
			dNdy[ 12 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL4dy ) ) + ( L2 * L3 * dL4dy * ( 4.0 * L4 - 1.0 ) ) + ( L2 * dL3dy * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL2dy * L3 * L4 * ( 4.0 * L4 - 1.0 ) ) );
			dNdz[ 12 ][ p ] = 32.0 * ( ( L2 * L3 * L4 * ( 4.0 * dL4dz ) ) + ( L2 * L3 * dL4dz * ( 4.0 * L4 - 1.0 ) ) + ( L2 * dL3dz * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL2dz * L3 * L4 * ( 4.0 * L4 - 1.0 ) ) );
                                                                                                                                                                                            
			// Face bubbles [ 0 2 3 ]                                                                                                                                                        
			dNdx[ 13 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL1dx ) ) + ( L1 * L3 * dL4dx * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL3dx * L4 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dx * L3 * L4 * ( 4.0 * L1 - 1.0 ) ) );
			dNdy[ 13 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL1dy ) ) + ( L1 * L3 * dL4dy * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL3dy * L4 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dy * L3 * L4 * ( 4.0 * L1 - 1.0 ) ) );
			dNdz[ 13 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL1dz ) ) + ( L1 * L3 * dL4dz * ( 4.0 * L1 - 1.0 ) ) + ( L1 * dL3dz * L4 * ( 4.0 * L1 - 1.0 ) ) + ( dL1dz * L3 * L4 * ( 4.0 * L1 - 1.0 ) ) );
							 		                                                                                                                                                         
			dNdx[ 14 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL3dx ) ) + ( L1 * L3 * dL4dx * ( 4.0 * L3 - 1.0 ) ) + ( L1 * dL3dx * L4 * ( 4.0 * L3 - 1.0 ) ) + ( dL1dx * L3 * L4 * ( 4.0 * L3 - 1.0 ) ) );
			dNdy[ 14 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL3dy ) ) + ( L1 * L3 * dL4dy * ( 4.0 * L3 - 1.0 ) ) + ( L1 * dL3dy * L4 * ( 4.0 * L3 - 1.0 ) ) + ( dL1dy * L3 * L4 * ( 4.0 * L3 - 1.0 ) ) );
			dNdz[ 14 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL3dz ) ) + ( L1 * L3 * dL4dz * ( 4.0 * L3 - 1.0 ) ) + ( L1 * dL3dz * L4 * ( 4.0 * L3 - 1.0 ) ) + ( dL1dz * L3 * L4 * ( 4.0 * L3 - 1.0 ) ) );
							 		                                                                                                                                                         
			dNdx[ 15 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL4dx ) ) + ( L1 * L3 * dL4dx * ( 4.0 * L4 - 1.0 ) ) + ( L1 * dL3dx * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL1dx * L3 * L4 * ( 4.0 * L4 - 1.0 ) ) );
			dNdy[ 15 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL4dy ) ) + ( L1 * L3 * dL4dy * ( 4.0 * L4 - 1.0 ) ) + ( L1 * dL3dy * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL1dy * L3 * L4 * ( 4.0 * L4 - 1.0 ) ) );
			dNdz[ 15 ][ p ] = 32.0 * ( ( L1 * L3 * L4 * ( 4.0 * dL4dz ) ) + ( L1 * L3 * dL4dz * ( 4.0 * L4 - 1.0 ) ) + ( L1 * dL3dz * L4 * ( 4.0 * L4 - 1.0 ) ) + ( dL1dz * L3 * L4 * ( 4.0 * L4 - 1.0 ) ) );

			// Inner bubble
			dNdx[ 16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dx ) + ( L1 * L2 * dL3dx * L4 ) + ( L1 * dL2dx * L3 * L4 ) + ( dL1dx * L2 * L3 * L4 ) );
			dNdy[ 16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dy ) + ( L1 * L2 * dL3dy * L4 ) + ( L1 * dL2dy * L3 * L4 ) + ( dL1dy * L2 * L3 * L4 ) );
			dNdz[ 16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dz ) + ( L1 * L2 * dL3dz * L4 ) + ( L1 * dL2dz * L3 * L4 ) + ( dL1dz * L2 * L3 * L4 ) );
        }
	}

	//********************************************************************************************************************************
    //* - Exterior normal 
    //********************************************************************************************************************************	
    void VolumeElement_3sb_ColdPlasma::Calculate_ExtNormal( Vector<double>& extN, Vector<int>& face )
    {
        double v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ face[2] ]->X() - mNodes[ face[0] ]->X();
        v2[ 1 ] = mNodes[ face[2] ]->Y() - mNodes[ face[0] ]->Y();
        v2[ 2 ] = mNodes[ face[2] ]->Z() - mNodes[ face[0] ]->Z();

        v1[ 1 ] = mNodes[ face[1] ]->Y() - mNodes[ face[0] ]->Y();
        v1[ 2 ] = mNodes[ face[1] ]->Z() - mNodes[ face[0] ]->Z();
        v1[ 0 ] = mNodes[ face[1] ]->X() - mNodes[ face[0] ]->X();

        extN.resize( 3 );

        extN[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        extN[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        extN[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

        extN *= ( 1.0 / sqrt( extN[ 0 ]*extN[ 0 ] + extN[ 1 ]*extN[ 1 ] + extN[ 2 ]*extN[ 2 ] ) );
    }

    //********************************************************************************************************************************
    //* - Area of the element
    //********************************************************************************************************************************	
    double VolumeElement_3sb_ColdPlasma::Calculate_Area( Vector<int>& face ) 
    {
		double n[ 3 ], v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ face[2] ]->X() - mNodes[ face[0] ]->X();
        v2[ 1 ] = mNodes[ face[2] ]->Y() - mNodes[ face[0] ]->Y();
        v2[ 2 ] = mNodes[ face[2] ]->Z() - mNodes[ face[0] ]->Z();

        v1[ 1 ] = mNodes[ face[1] ]->Y() - mNodes[ face[0] ]->Y();
        v1[ 2 ] = mNodes[ face[1] ]->Z() - mNodes[ face[0] ]->Z();
        v1[ 0 ] = mNodes[ face[1] ]->X() - mNodes[ face[0] ]->X();

        // area = 0.5*(v2 x v1)
        n[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        n[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        n[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

        return ( 0.5 * sqrt( n[ 0 ]*n[ 0 ] + n[ 1 ]*n[ 1 ] + n[ 2 ]*n[ 2 ] ) );
    }

    //********************************************************************************************************************************
    //* - Permittivity tensor on nodes 
    //********************************************************************************************************************************
    void VolumeElement_3sb_ColdPlasma::PermittivityTensor_OnNodes( Vector< Matrix< std::complex<double> > >& TEn )
    {
        // Permittivity tensor on main nodes
        TEn.resize( 4 );
        
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

    //********************************************************************************************************************************
    //* - Permittivity tensor on Gauss points ( extrapolating TEn to GPs with N )
    //********************************************************************************************************************************
    void VolumeElement_3sb_ColdPlasma::PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp,  
                                                                         Vector< Matrix< std::complex<double> > >& TEn , 
                                                                         Matrix< double >& N )
    {
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

    //********************************************************************************************************************************
    //* - Derivatives of the permittivity tensor on Gauss points
    //********************************************************************************************************************************
    void VolumeElement_3sb_ColdPlasma::PermittivityTensor_Derivatives( Vector< Matrix< std::complex<double> > >& dTEdx, 
                                                                       Vector< Matrix< std::complex<double> > >& dTEdy,
                                                                       Vector< Matrix< std::complex<double> > >& dTEdz,
                                                                       Vector< Matrix< std::complex<double> > >& TEn  ,
                                                                       Matrix< double >& dNdx, 
                                                                       Matrix< double >& dNdy, 
                                                                       Matrix< double >& dNdz )
    {
        // Number of Gauss points
        int numGaussPoints = dNdx[ 0 ].size();

        // Initializing vector of 3x3 permittivity tensors derivatives
        dTEdx.resize( numGaussPoints );
        dTEdy.resize( numGaussPoints );
        dTEdz.resize( numGaussPoints );

        // TE derivatives: (dTE/dx)|gp = Sum_i{ TEni * (dNi/dx)|gp } 
        for( int gp=0; gp<numGaussPoints; gp++ )
        { 
            dTEdx[gp].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );
            dTEdy[gp].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );
            dTEdz[gp].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );

            for( int ki=0; ki<3; ki++ )
            {
                for( int kj=0; kj<3; kj++ )
                {
                    for( int i=0; i<TEn.size(); i++ )
                    {
                        dTEdx[ gp ][ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * dNdx[ i ][ gp ];
                        dTEdy[ gp ][ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * dNdy[ i ][ gp ];
                        dTEdz[ gp ][ ki ][ kj ] += TEn[ i ][ ki ][ kj ] * dNdz[ i ][ gp ];
                    }
                }
            }
        }
    }

    //******************************************************************************************************************************** 
    //* - Inverse of the regularization weight tau
    //********************************************************************************************************************************
    void VolumeElement_3sb_ColdPlasma::Get_invTau_OnGaussPoints( Vector< std::complex<double> >& invTau, Vector< Matrix< std::complex<double> > >& TEgp )
    {
        int numGaussPoints = TEgp.size();

        invTau.resize( numGaussPoints );
        
        for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            std::complex<double> Tr_EcE;
            
            Tr_EcE = ( TEgp[ gp ][ 0 ][ 0 ] * std::conj( TEgp[ gp ][ 0 ][ 0 ] ) + 
                       TEgp[ gp ][ 0 ][ 1 ] * std::conj( TEgp[ gp ][ 0 ][ 1 ] ) + 
                       TEgp[ gp ][ 0 ][ 2 ] * std::conj( TEgp[ gp ][ 0 ][ 2 ] ) + 
                             
                       TEgp[ gp ][ 1 ][ 0 ] * std::conj( TEgp[ gp ][ 1 ][ 0 ] ) + 
                       TEgp[ gp ][ 1 ][ 1 ] * std::conj( TEgp[ gp ][ 1 ][ 1 ] ) + 
                       TEgp[ gp ][ 1 ][ 2 ] * std::conj( TEgp[ gp ][ 1 ][ 2 ] ) + 
                             
                       TEgp[ gp ][ 2 ][ 0 ] * std::conj( TEgp[ gp ][ 2 ][ 0 ] ) + 
                       TEgp[ gp ][ 2 ][ 1 ] * std::conj( TEgp[ gp ][ 2 ][ 1 ] ) + 
                       TEgp[ gp ][ 2 ][ 2 ] * std::conj( TEgp[ gp ][ 2 ][ 2 ] ) );             
            
            invTau[ gp ] = 3.0 / Tr_EcE;
        }    
    }

    //********************************************************************************************************************************
    //* - Surface integral (nxu)*q in PEC faces to C_cur matrix
    //********************************************************************************************************************************	
	void VolumeElement_3sb_ColdPlasma::Add_C_NxN_Matrix( Matrix< std::complex<double> >& C_cur ) 
	{
		Vector< Vector<int> > Faces( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;

		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			// Face vector
			Vector<int> nf = *fs_it;

			// Check if the surface is not PEC 
			if( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) 
			{
				continue;
			}

			// Gauss points and weights
		    Vector<double> cX, cY, W;
		    
		    int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );
		    
            // Jacobian
		    double detJ = 2.0 * Calculate_Area( nf );
		    
		    for( int gp=0; gp<numGaussPoints; gp++ ) 
			{
				W[ gp ] *= detJ; 
			}
		    
            // 3sb base on Gauss points
		    Matrix<double> N;
		    
		    Lagrange2D_Ni_3sb( N, cX, cY );

		    Vector<double> extN; 
			
			Calculate_ExtNormal( extN, nf );

			for( int i=0; i<3; i++ )
		    {
		        for( int j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_NiNj( 0.0, 0.0 );
		            
                    for( int gp=0; gp<numGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[ gp ] * N[ i ][ gp ] * N[ j ][ gp ]; 
					}
   		
                    C_cur[ nf[i]     ][ nf[j] +   mNumNodes ] += ( extN[ 2 ] * Sintg_NiNj ); 
                    C_cur[ nf[i]     ][ nf[j] + 2*mNumNodes ] -= ( extN[ 1 ] * Sintg_NiNj ); 
		        					    	 					   		 			    
		        	C_cur[ nf[i] + 4 ][ nf[j]               ] -= ( extN[ 2 ] * Sintg_NiNj ); 
                    C_cur[ nf[i] + 4 ][ nf[j] + 2*mNumNodes ] += ( extN[ 0 ] * Sintg_NiNj );
		        					 	 						   		 			    
		        	C_cur[ nf[i] + 8 ][ nf[j]               ] += ( extN[ 1 ] * Sintg_NiNj ); 
                    C_cur[ nf[i] + 8 ][ nf[j] +   mNumNodes ] -= ( extN[ 0 ] * Sintg_NiNj ); 
		        }
		    }
		}
	}

	//********************************************************************************************************************************
    //* - Stabilization term: ( curl_F, curl_E ) + ( div_F, div_E ) 
    //********************************************************************************************************************************		
    void VolumeElement_3sb_ColdPlasma::Add_CurlDiv_Stab_Matrix_EA( Matrix< std::complex<double> >& StiffMatrix ) 
    {  
		// Matrix with only [E|A] DOFs
		Matrix< std::complex<double> > Stab_Matrix( 51, 51 ); 
		
		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		Derivatives_On_Points( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Upper diagonal [Kxx], [Kyy], [Kzz]
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=i; j<mNumNodes; j++ )
			{
				double intg_CurlDiv = 0.0; 
				
				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_CurlDiv += W[ gp ] * ( dNdx[ i ][ gp ] * dNdx[ j ][ gp ] + dNdy[ i ][ gp ] * dNdy[ j ][ gp ] + dNdz[ i ][ gp ] * dNdz[ j ][ gp ] );
				}

				Stab_Matrix[ i               ][ j               ] = intg_CurlDiv;
				Stab_Matrix[ i +   mNumNodes ][ j +   mNumNodes ] = intg_CurlDiv;
				Stab_Matrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = intg_CurlDiv;
			}
		}

		// [Kxy], [Kxz], [Kyz] 
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=0; j<mNumNodes; j++ )
			{
				double intg_CurlDiv_XY = 0.0;
				double intg_CurlDiv_XZ = 0.0;
				double intg_CurlDiv_YZ = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
				    intg_CurlDiv_XY += W[ gp ] * ( dNdx[ i ][ gp ] * dNdy[ j ][ gp ] - dNdy[ i ][ gp ] * dNdx[ j ][ gp ] );
                    intg_CurlDiv_XZ += W[ gp ] * ( dNdx[ i ][ gp ] * dNdz[ j ][ gp ] - dNdz[ i ][ gp ] * dNdx[ j ][ gp ] );
                    intg_CurlDiv_YZ += W[ gp ] * ( dNdy[ i ][ gp ] * dNdz[ j ][ gp ] - dNdz[ i ][ gp ] * dNdy[ j ][ gp ] );
				}
	            
				Stab_Matrix[ i             ][ j +   mNumNodes ] = intg_CurlDiv_XY;
				Stab_Matrix[ i             ][ j + 2*mNumNodes ] = intg_CurlDiv_XZ;
				Stab_Matrix[ i + mNumNodes ][ j + 2*mNumNodes ] = intg_CurlDiv_YZ;
			}
		}

		// Lower diagonal ( [EE] | [AA] )
        for( int i=0; i<51; i++ )
        {
            for( int j=i+1; j<51; j++ )
            {
                Stab_Matrix[ j ][ i ] = Stab_Matrix[ i ][ j ];
            }
        }  
	    
		for( int i=0; i<51; i++ )
        {
            for( int j=0; j<51; j++ )
            {
                StiffMatrix[ i ][ j ] += Stab_Matrix[ i ][ j ];
            }
        } 
	}

    //********************************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //********************************************************************************************************************************
    void VolumeElement_3sb_ColdPlasma::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
                                                                    Matrix< std::complex<double> >& EleStiffMatrix, 
                                                                    Vector< std::complex<double> >& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		if( mPotentialsOn ) 
		{
            for( int j=0; j<4; j++ ) 
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

	//********************************************************************************************************************************
    // - Stiffness matrix
    //********************************************************************************************************************************
    void VolumeElement_3sb_ColdPlasma::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix )
	{   
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		if( mPotentialsOn ) 
	    {
	        Get_CurlDiv_LL2P_Matrix_AV( StiffMatrix ); 
	    }
        else
        {
            Get_CurlDiv_LL2P_Matrix_Ef( StiffMatrix ); 
        }	

		if( mhkFactor != 0.0 )
		{
		    Add_CurlDiv_Stab_Matrix_EA( StiffMatrix );
		}
	}

	//********************************************************************************************************************************
    //* - Local L2 projection method stiffness matrix ( Ef formulation )
    //********************************************************************************************************************************	
    void VolumeElement_3sb_ColdPlasma::Get_CurlDiv_LL2P_Matrix_Ef( Matrix< std::complex<double> >& StiffMatrix ) 
    { 
        // Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        double pi = 3.14159265358979;
        double mo = pi * 4.0e-7;

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> invMu = 1.0 / cMu;

		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[gp] *= detJ; 
		}

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_3sb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		Derivatives_On_Points( dNdx, dNdy, dNdz, cX, cY, cZ );

        // Permitivitty tensor on nodes
        Vector< Matrix< std::complex<double> > > TEn;

        PermittivityTensor_OnNodes( TEn );

        // Permittivity tensor on Gauss points
        Vector< Matrix< std::complex<double> > > TEgp;

        PermittivityTensor_OnGaussPoints( TEgp, TEn, N );

        // Permittivity tensor derivatives
        Vector< Matrix< std::complex<double> > > dTEdx, dTEdy, dTEdz;

        PermittivityTensor_Derivatives( dTEdx, dTEdy, dTEdz, TEn, dNdx, dNdy, dNdz );

        // Regularization constant
        Vector< std::complex<double> > invTau;

        Get_invTau_OnGaussPoints( invTau, TEgp );

		// Inverted P1 projection matrices
		Matrix< std::complex<double> > iA_div(  4,  4 );
		Matrix< std::complex<double> > iA_cur( 12, 12 );
		
		for( int i=0; i<4; i++ )
		{
            for( int j=0; j<4; j++ )
            {
        		if( i != j ) 
				{
					iA_div[ i ][ j ] = -4.00 / mVolume;
				}
        		else
				{
					iA_div[ i ][ j ] = +16.0 / mVolume;                          
				}
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}
		
		// Curl matrix 
		Matrix< std::complex<double> > C_cur( 12, 51 );

		for( int i=0; i<4; i++ )
		{
		    for( int j=0; j<mNumNodes; j++ )
		    {
		    	double intg_Ni_dNjdX = 0.0;
		        double intg_Ni_dNjdY = 0.0;
		        double intg_Ni_dNjdZ = 0.0;
		    
                for( int gp=0; gp<numGaussPoints; gp++ )
                {
                	intg_Ni_dNjdX += W[ gp ] * N[ i ][ gp ] * dNdx[ j ][ gp ]; 
		    		intg_Ni_dNjdY += W[ gp ] * N[ i ][ gp ] * dNdy[ j ][ gp ];
		    		intg_Ni_dNjdZ += W[ gp ] * N[ i ][ gp ] * dNdz[ j ][ gp ];
                }      
		    	
		    	C_cur[ i     ][ j               ] =            0.0; 
                C_cur[ i     ][ j +   mNumNodes ] = -intg_Ni_dNjdZ; 
                C_cur[ i     ][ j + 2*mNumNodes ] = +intg_Ni_dNjdY; 
		    								    
		    	C_cur[ i + 4 ][ j               ] = +intg_Ni_dNjdZ; 
                C_cur[ i + 4 ][ j +   mNumNodes ] =            0.0; 
                C_cur[ i + 4 ][ j + 2*mNumNodes ] = -intg_Ni_dNjdX; 
		    								    
		    	C_cur[ i + 8 ][ j               ] = -intg_Ni_dNjdY; 
                C_cur[ i + 8 ][ j +   mNumNodes ] = +intg_Ni_dNjdX; 
                C_cur[ i + 8 ][ j + 2*mNumNodes ] =            0.0; 
		    }
		}
  
		Add_C_NxN_Matrix( C_cur ); 

        iA_cur *= invMu;

		Matrix< std::complex<double> > iAc_x_Ccur = iA_cur * C_cur;

		Matrix< std::complex<double> > Ccur_tran  = C_cur.Transpose();

		StiffMatrix += ( Ccur_tran * iAc_x_Ccur );

        // Regularization weight matrix
        if( mPeso == 1.0 ) 
        {
            Matrix< std::complex<double> > A_reg( 4, 4 );
            
            for( int i=0; i<4; i++ )
		    {
                for( int j=0; j<4; j++ )
                {
                    std::complex<double> intg_Ni_invTau_Nj( 0.0, 0.0 );
                    
                    for( int gp=0; gp<numGaussPoints; gp++ )
                    {
                        intg_Ni_invTau_Nj += W[ gp ] * N[ i ][ gp ] * invTau[ gp ] * N[ j ][ gp ];
                    }
            		
		    		A_reg[ i ][ j ] = invMu * intg_Ni_invTau_Nj; 
                }									   
		    }
            
            // Div matrix
		    Matrix< std::complex<double> > D_div( 4, 51 );
            
		    for( int i=0; i<4; i++ )
		    {
		        for( int j=0; j<mNumNodes; j++ )
		        {
		        	std::complex<double> intg_Ni_Div_eEX( 0.0, 0.0 );
		            std::complex<double> intg_Ni_Div_eEY( 0.0, 0.0 );
		            std::complex<double> intg_Ni_Div_eEZ( 0.0, 0.0 );
		        
                    for( int gp=0; gp<numGaussPoints; gp++ )
                    {                    
                        intg_Ni_Div_eEX += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 0 ][ 0 ] * dNdx[ j ][ gp ] + N[ j ][ gp ] * dTEdx[ gp ][ 0 ][ 0 ] + 
                                                                      TEgp[ gp ][ 1 ][ 0 ] * dNdy[ j ][ gp ] + N[ j ][ gp ] * dTEdy[ gp ][ 1 ][ 0 ] +
                                                                      TEgp[ gp ][ 2 ][ 0 ] * dNdz[ j ][ gp ] + N[ j ][ gp ] * dTEdz[ gp ][ 2 ][ 0 ] );
                        
                    	intg_Ni_Div_eEY += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 0 ][ 1 ] * dNdx[ j ][ gp ] + N[ j ][ gp ] * dTEdx[ gp ][ 0 ][ 1 ] +
                                                                      TEgp[ gp ][ 1 ][ 1 ] * dNdy[ j ][ gp ] + N[ j ][ gp ] * dTEdy[ gp ][ 1 ][ 1 ] +
                                                                      TEgp[ gp ][ 2 ][ 1 ] * dNdz[ j ][ gp ] + N[ j ][ gp ] * dTEdz[ gp ][ 2 ][ 1 ] );
                        
                    	intg_Ni_Div_eEZ += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 0 ][ 2 ] * dNdx[ j ][ gp ] + N[ j ][ gp ] * dTEdx[ gp ][ 0 ][ 2 ] +
                                                                      TEgp[ gp ][ 1 ][ 2 ] * dNdy[ j ][ gp ] + N[ j ][ gp ] * dTEdy[ gp ][ 1 ][ 2 ] +
                                                                      TEgp[ gp ][ 2 ][ 2 ] * dNdz[ j ][ gp ] + N[ j ][ gp ] * dTEdz[ gp ][ 2 ][ 2 ] );
                    }      
            
		        	D_div[ i ][ j               ] = intg_Ni_Div_eEX; 
		        	D_div[ i ][ j +   mNumNodes ] = intg_Ni_Div_eEY; 
		        	D_div[ i ][ j + 2*mNumNodes ] = intg_Ni_Div_eEZ; 
		        }
		    }
            
		    Matrix< std::complex<double> > iAd_x_Ddiv = iA_div * D_div;
            
            Matrix< std::complex<double> > Ar_x_iAdxDdiv = A_reg * iAd_x_Ddiv;
            
            Matrix< std::complex<double> > DdivT_x_iAdT = D_div.Transpose() * iA_div.Transpose();
            
		    StiffMatrix += ( DdivT_x_iAdT * Ar_x_iAdxDdiv );
        }

		// NiNj matrix
        for( int i=0; i<mNumNodes; i++ )
		{
            for( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Ni_TEp_Nj_XX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZZ( 0.0, 0.0 ); 

                std::complex<double> intg_Ni_TEp_Nj_XY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_XZ( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YZ( 0.0, 0.0 ); 

                std::complex<double> intg_Ni_TEp_Nj_YX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZY( 0.0, 0.0 ); 

                for( int gp=0; gp<numGaussPoints; gp++ )
                {
                	intg_Ni_TEp_Nj_XX += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 0 ] * N[ j ][ gp ];
                	intg_Ni_TEp_Nj_XY += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 1 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_XZ += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 2 ] * N[ j ][ gp ];

                	intg_Ni_TEp_Nj_YX += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 0 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_YY += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 1 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_YZ += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 2 ] * N[ j ][ gp ];

                    intg_Ni_TEp_Nj_ZX += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 0 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_ZY += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 1 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_ZZ += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 2 ] * N[ j ][ gp ];
                }
            
            	StiffMatrix[ i               ][ j               ] -= ( w2 * intg_Ni_TEp_Nj_XX ); 
            	StiffMatrix[ i               ][ j +   mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_XY ); 
            	StiffMatrix[ i               ][ j + 2*mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_XZ ); 

            	StiffMatrix[ i +   mNumNodes ][ j               ] -= ( w2 * intg_Ni_TEp_Nj_YX ); 
            	StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_YY ); 
            	StiffMatrix[ i +   mNumNodes ][ j + 2*mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_YZ ); 

            	StiffMatrix[ i + 2*mNumNodes ][ j               ] -= ( w2 * intg_Ni_TEp_Nj_ZX ); 
            	StiffMatrix[ i + 2*mNumNodes ][ j +   mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_ZY ); 
            	StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_ZZ ); 
            }
		}
    }

	//********************************************************************************************************************************
    //* - Local L2 projection method stiffness matrix ( AV formulation )
    //********************************************************************************************************************************	
    void VolumeElement_3sb_ColdPlasma::Get_CurlDiv_LL2P_Matrix_AV( Matrix< std::complex<double> >& StiffMatrix ) 
    { 
		// Matrix with only A DOFs
		Matrix< std::complex<double> > LL2PCurlDiv_Matrix;

		LL2PCurlDiv_Matrix.Resize( 51, 51, std::complex<double>( 0.0, 0.0 ) );

		// Material properties
        double wf = (*mProperties)( FREQUENCY );
        double w2 = wf * wf;   

        double pi = 3.14159265358979;
        double mo = pi * 4.0e-7;

        double mu_real = (*mProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
        double mu_imag = (*mProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu( mu_real, mu_imag );
		std::complex<double> invMu = 1.0 / cMu;

		// Gauss points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_3sb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		Derivatives_On_Points( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< std::complex<double> > D_div(  4, 51 );
		Matrix< std::complex<double> > C_cur( 12, 51 );

		for( int i=0; i<4; i++ )
		{
		    for( int j=0; j<mNumNodes; j++ )
		    {
		    	double intg_Ni_dNjdX = 0.0;
		        double intg_Ni_dNjdY = 0.0;
		        double intg_Ni_dNjdZ = 0.0;
		    
                for( int gp=0; gp<numGaussPoints; gp++ )
                {
                	intg_Ni_dNjdX += W[ gp ] * N[ i ][ gp ] * dNdx[ j ][ gp ]; 
		    		intg_Ni_dNjdY += W[ gp ] * N[ i ][ gp ] * dNdy[ j ][ gp ];
		    		intg_Ni_dNjdZ += W[ gp ] * N[ i ][ gp ] * dNdz[ j ][ gp ];
                } 

				// D_div
				D_div[ i     ][ j               ] = intg_Ni_dNjdX; 
		    	D_div[ i     ][ j +   mNumNodes ] = intg_Ni_dNjdY; 
		    	D_div[ i     ][ j + 2*mNumNodes ] = intg_Ni_dNjdZ; 

				// C_curl
		    	C_cur[ i     ][ j               ] =            0.0; 
                C_cur[ i     ][ j +   mNumNodes ] = -intg_Ni_dNjdZ; 
                C_cur[ i     ][ j + 2*mNumNodes ] = +intg_Ni_dNjdY; 
		    								    
		    	C_cur[ i + 4 ][ j               ] = +intg_Ni_dNjdZ; 
                C_cur[ i + 4 ][ j +   mNumNodes ] =            0.0; 
                C_cur[ i + 4 ][ j + 2*mNumNodes ] = -intg_Ni_dNjdX; 
		    								    
		    	C_cur[ i + 8 ][ j               ] = -intg_Ni_dNjdY; 
                C_cur[ i + 8 ][ j +   mNumNodes ] = +intg_Ni_dNjdX; 
                C_cur[ i + 8 ][ j + 2*mNumNodes ] =            0.0; 
		    }
		}
  
		// Adding surface integral (nxu)*q
		Add_C_NxN_Matrix( C_cur ); 

		// Projection in q basis matrices
		Matrix< std::complex<double> > iA_div(  4,  4 );
		Matrix< std::complex<double> > iA_cur( 12, 12 );
		
		for( int i=0; i<4; i++ )
		{
            for( int j=0; j<4; j++ )
            {
        		if( i != j ) 
				{
					iA_div[ i ][ j ] = -4.00 / mVolume;
				}
        		else
				{
					iA_div[ i ][ j ] = +16.0 / mVolume;                          
				}
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}
		
		// Scaling the projection in q matrices
		iA_div /= cMu;
		iA_cur /= cMu;
		
		// Final assembly: [ CT * Ac^(-1) * C ] + [ DT * Ad^(-1) * D ] - [ Ni * Nj ]
		Matrix< std::complex<double> > iAc_x_Ccur = iA_cur * C_cur;
		Matrix< std::complex<double> > Ccur_tran  = C_cur.Transpose();

		LL2PCurlDiv_Matrix += ( Ccur_tran * iAc_x_Ccur );

        if( mPeso == 1.0 ) 
        {
		    Matrix< std::complex<double> > iAd_x_Ddiv = iA_div * D_div;
		    Matrix< std::complex<double> > Ddiv_tran  = D_div.Transpose();

		    LL2PCurlDiv_Matrix += ( Ddiv_tran * iAd_x_Ddiv );
        }

		// Add LL2PCurlDiv_Matrix to StiffMatrix
		for( int i=0; i<51; i++ )
        {
            for( int j=0; j<51; j++ )
            {
                StiffMatrix[ i ][ j ] = LL2PCurlDiv_Matrix[ i ][ j ];
            }
        }  

		// NiNj matrix
        Vector< Matrix< std::complex<double> > > TEn, TEgp; 

        PermittivityTensor_OnNodes( TEn ); 
        PermittivityTensor_OnGaussPoints( TEgp, TEn, N );

        for( int i=0; i<mNumNodes; i++ )
		{
            for( int j=0; j<mNumNodes; j++ )
            {
                std::complex<double> intg_Ni_TEp_Nj_XX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZZ( 0.0, 0.0 ); 

                std::complex<double> intg_Ni_TEp_Nj_XY( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_XZ( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_YZ( 0.0, 0.0 ); 

                std::complex<double> intg_Ni_TEp_Nj_YX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZX( 0.0, 0.0 ); 
                std::complex<double> intg_Ni_TEp_Nj_ZY( 0.0, 0.0 ); 

                for( int gp=0; gp<numGaussPoints; gp++ )
                {
                	intg_Ni_TEp_Nj_XX += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 0 ] * N[ j ][ gp ];
                	intg_Ni_TEp_Nj_XY += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 1 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_XZ += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 0 ][ 2 ] * N[ j ][ gp ];

                	intg_Ni_TEp_Nj_YX += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 0 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_YY += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 1 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_YZ += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 1 ][ 2 ] * N[ j ][ gp ];

                    intg_Ni_TEp_Nj_ZX += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 0 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_ZY += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 1 ] * N[ j ][ gp ];
                    intg_Ni_TEp_Nj_ZZ += W[ gp ] * N[ i ][ gp ] * TEgp[ gp ][ 2 ][ 2 ] * N[ j ][ gp ];
                }
            
            	StiffMatrix[ i               ][ j               ] -= ( w2 * intg_Ni_TEp_Nj_XX ); 
            	StiffMatrix[ i               ][ j +   mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_XY ); 
            	StiffMatrix[ i               ][ j + 2*mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_XZ ); 

            	StiffMatrix[ i +   mNumNodes ][ j               ] -= ( w2 * intg_Ni_TEp_Nj_YX ); 
            	StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_YY ); 
            	StiffMatrix[ i +   mNumNodes ][ j + 2*mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_YZ ); 

            	StiffMatrix[ i + 2*mNumNodes ][ j               ] -= ( w2 * intg_Ni_TEp_Nj_ZX ); 
            	StiffMatrix[ i + 2*mNumNodes ][ j +   mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_ZY ); 
            	StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] -= ( w2 * intg_Ni_TEp_Nj_ZZ ); 
            }
        }

        // [AV], [VA], [VV] contribution
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

        // [Kvv]
		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
                std::complex<double> intg_dNi_TEp_dNj( 0.0, 0.0 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_dNi_TEp_dNj += W[ gp ] * DN[ 0 ][ i ] * ( TEgp[ gp ][ 0 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 0 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 0 ][ 2 ] * DN[ 2 ][ j ] );
                    intg_dNi_TEp_dNj += W[ gp ] * DN[ 1 ][ i ] * ( TEgp[ gp ][ 1 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 1 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 1 ][ 2 ] * DN[ 2 ][ j ] );
                    intg_dNi_TEp_dNj += W[ gp ] * DN[ 2 ][ i ] * ( TEgp[ gp ][ 2 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 2 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 2 ][ 2 ] * DN[ 2 ][ j ] );
				}
	                
				StiffMatrix[ i + 3*mNumNodes ][ j + 3*mNumNodes ] -= ( w2 * intg_dNi_TEp_dNj );
        	}
        }

		// [Kxv], [Kyv], [Kzv]
		for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
                std::complex<double> intg_NiX_TEp_dNj( 0.0, 0.0 ); 
                std::complex<double> intg_NiY_TEp_dNj( 0.0, 0.0 ); 
                std::complex<double> intg_NiZ_TEp_dNj( 0.0, 0.0 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_NiX_TEp_dNj += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 0 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 0 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 0 ][ 2 ] * DN[ 2 ][ j ] );
                    intg_NiY_TEp_dNj += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 1 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 1 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 1 ][ 2 ] * DN[ 2 ][ j ] );
                    intg_NiZ_TEp_dNj += W[ gp ] * N[ i ][ gp ] * ( TEgp[ gp ][ 2 ][ 0 ] * DN[ 0 ][ j ] + TEgp[ gp ][ 2 ][ 1 ] * DN[ 1 ][ j ] + TEgp[ gp ][ 2 ][ 2 ] * DN[ 2 ][ j ] );
				}
	                
				StiffMatrix[ i               ][ j + 3*mNumNodes ] -= ( w2 * intg_NiX_TEp_dNj );
                StiffMatrix[ i +   mNumNodes ][ j + 3*mNumNodes ] -= ( w2 * intg_NiY_TEp_dNj );
                StiffMatrix[ i + 2*mNumNodes ][ j + 3*mNumNodes ] -= ( w2 * intg_NiZ_TEp_dNj );
        	}
        }

		// [Kvx], [Kvy], [Kvz]
		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<mNumNodes; j++ ) 
        	{
                std::complex<double> intg_dNi_TEp_NjX( 0.0, 0.0 ); 
                std::complex<double> intg_dNi_TEp_NjY( 0.0, 0.0 ); 
                std::complex<double> intg_dNi_TEp_NjZ( 0.0, 0.0 ); 

				for ( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_dNi_TEp_NjX += W[ gp ] * N[ j ][ gp ] * ( TEgp[ gp ][ 0 ][ 0 ] * DN[ 0 ][ i ] + TEgp[ gp ][ 1 ][ 0 ] * DN[ 1 ][ i ] + TEgp[ gp ][ 2 ][ 0 ] * DN[ 2 ][ i ] );
                    intg_dNi_TEp_NjY += W[ gp ] * N[ j ][ gp ] * ( TEgp[ gp ][ 0 ][ 1 ] * DN[ 0 ][ i ] + TEgp[ gp ][ 1 ][ 1 ] * DN[ 1 ][ i ] + TEgp[ gp ][ 2 ][ 1 ] * DN[ 2 ][ i ] );
                    intg_dNi_TEp_NjZ += W[ gp ] * N[ j ][ gp ] * ( TEgp[ gp ][ 0 ][ 2 ] * DN[ 0 ][ i ] + TEgp[ gp ][ 1 ][ 2 ] * DN[ 1 ][ i ] + TEgp[ gp ][ 2 ][ 2 ] * DN[ 2 ][ i ] );
				}
	                
				StiffMatrix[ i + 3*mNumNodes ][ j               ] -= ( w2 * intg_dNi_TEp_NjX );
                StiffMatrix[ i + 3*mNumNodes ][ j +   mNumNodes ] -= ( w2 * intg_dNi_TEp_NjY );
                StiffMatrix[ i + 3*mNumNodes ][ j + 2*mNumNodes ] -= ( w2 * intg_dNi_TEp_NjZ );
        	}
        }
    }
} 