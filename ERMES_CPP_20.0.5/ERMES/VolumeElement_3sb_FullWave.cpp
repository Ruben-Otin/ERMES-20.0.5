
#include "VolumeElement_3sb_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{
    //********************************************************************************************************************************
    // - Global index of the DOFs
    //********************************************************************************************************************************	
    void VolumeElement_3sb_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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
    double VolumeElement_3sb_FullWave::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //********************************************************************************************************************************
    // - 1st order base derivatives - DN[ X,Y,Z ][ 0,1,2,3 ]
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::Calculate_DN( double DN[ 3 ][ 4 ] )
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
    void VolumeElement_3sb_FullWave::Derivatives_On_Points( Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz,
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
    //* - Dirichlet contribution to the residual vector
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
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
    void VolumeElement_3sb_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix )
	{   
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		Add_CurlDiv_LL2P_Matrix( StiffMatrix ); 

		if( mhkFactor != 0.0 )
		{
		    Add_CurlDiv_Stab_Matrix( StiffMatrix );
		}

		if( mPotentialsOn ) 
		{
		    Add_AV_VA_VV_Matrix( StiffMatrix ); 
		}
	}

	//********************************************************************************************************************************
    //* - Stabilization term: ( curl_F, curl_E ) + ( div_F, div_E ) 
    //********************************************************************************************************************************	
    void VolumeElement_3sb_FullWave::Add_CurlDiv_Stab_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
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
    //* - Exterior normal 
    //********************************************************************************************************************************	
    void VolumeElement_3sb_FullWave::Calculate_ExtNormal( Vector<double>& extN, Vector<int>& face )
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
    double VolumeElement_3sb_FullWave::Calculate_Area( Vector<int>& face ) 
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
    //* - Surface integral (nxu)*q in PEC faces to C_cur matrix
    //********************************************************************************************************************************	
	void VolumeElement_3sb_FullWave::Add_C_NxN_Matrix( Matrix< std::complex<double> >& C_cur ) 
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
    //* - Local L2 projection method stiffness matrix ( div(v), q ) , ( curl(v), w )
    //********************************************************************************************************************************	
    void VolumeElement_3sb_FullWave::Add_CurlDiv_LL2P_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
    { 
		// Matrix with only [ E | A ] DOFs
		Matrix< std::complex<double> > LL2PCurlDiv_Matrix;

		LL2PCurlDiv_Matrix.Resize( 51, 51, std::complex<double>( 0.0, 0.0 ) );

		// Material properties
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
		std::complex<double> cEps;

		if( mComplxFreqOn )
		{
		    std::complex<double> jUnit( 0.0, 1.0 );

		    std::complex<double> sLaplace  = mComplxFreq;
		    std::complex<double> sLaplace2 = sLaplace * sLaplace;
		    
            w2cEps = ( + ( sg_real + jUnit * sg_imag ) * sLaplace ) - ( ( ep_real + jUnit * ep_imag ) * sLaplace2 );
			cEps   = ( - ( sg_real + jUnit * sg_imag ) / sLaplace ) + ( ( ep_real + jUnit * ep_imag )             );
		}
		else
		{
            w2cEps = std::complex<double>( ep_real * w2 - ( sg_imag * wf ), ep_imag * w2 + ( sg_real * wf ) );
			cEps   = std::complex<double>( ep_real      - ( sg_imag / wf ), ep_imag      + ( sg_real / wf ) );
		}

		std::complex<double> cMu( mu_real, mu_imag );

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
				if( mPotentialsOn )
				{
				    D_div[ i ][ j               ] = intg_Ni_dNjdX; 
		    	    D_div[ i ][ j +   mNumNodes ] = intg_Ni_dNjdY; 
		    	    D_div[ i ][ j + 2*mNumNodes ] = intg_Ni_dNjdZ; 
				}
				else
				{
					D_div[ i ][ j               ] = cEps * intg_Ni_dNjdX; 
		    	    D_div[ i ][ j +   mNumNodes ] = cEps * intg_Ni_dNjdY; 
		    	    D_div[ i ][ j + 2*mNumNodes ] = cEps * intg_Ni_dNjdZ; 
				}
		    	
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

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij( 51, 51 );
		
		for( int i=0; i<mNumNodes; i++ )
		{
            for( int j=0; j<mNumNodes; j++ )
            {
                double intg_NiNj = 0.0;

                for( int gp=0; gp<numGaussPoints; gp++ )
                {
                	intg_NiNj += W[ gp ] * N[ i ][ gp ] * N[ j ][ gp ]; 
                }
            
            	N_nij[ i               ][ j               ] = w2cEps * intg_NiNj; 
            	N_nij[ i +   mNumNodes ][ j +   mNumNodes ] = w2cEps * intg_NiNj; 
            	N_nij[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = w2cEps * intg_NiNj; 
            }
		}

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
		if( mPotentialsOn )
		{
		    iA_div /= cMu;
		}
		else
		{
		    iA_div /= ( cMu * cEps * cEps );
		}

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

		LL2PCurlDiv_Matrix -= N_nij;

		// Add LL2PCurlDiv_Matrix to StiffMatrix
		for( int i=0; i<51; i++ )
        {
            for( int j=0; j<51; j++ )
            {
                StiffMatrix[ i ][ j ] = LL2PCurlDiv_Matrix[ i ][ j ];
            }
        }  
    }

    //**********************************************************************************************************************
    // - AV-potentials contribution [AV], [VA], [VV] to stiffness matrix
    //**********************************************************************************************************************
    void VolumeElement_3sb_FullWave::Add_AV_VA_VV_Matrix( Matrix< std::complex<double> >& StiffMatrix ) 
	{ 
        // Material properties
		double wf = (*mProperties)( FREQUENCY );

        double w2 = wf * wf;
        double eo = 8.8541878176e-12;
	    
		double sgm_real = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
		double sgm_imag = (*mProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        double eps_real = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double eps_imag = (*mProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		std::complex<double> w2cEps;

		if( mComplxFreqOn )
		{
		    std::complex<double> jUnit( 0.0, 1.0 );
		    
		    std::complex<double> sLaplace  = mComplxFreq;
		    std::complex<double> sLaplace2 = sLaplace * sLaplace;
		    
            w2cEps = sLaplace * ( sgm_real + jUnit * sgm_imag ) - sLaplace2 * ( eps_real + jUnit * eps_imag );		 
		}
		else
		{
            w2cEps = std::complex<double>( w2 * eps_real - wf * sgm_imag, w2 * eps_imag + wf * sgm_real );
		}

		std::complex<double> nw2cEpsVol = -mVolume * w2cEps;

		// Volume integral of Ni 
		Vector< std::complex<double> > nw2cEps_intgVolNi( mNumNodes );

		for( int i=0; i< 4; i++ ) nw2cEps_intgVolNi[ i ] = nw2cEpsVol / 4.0;

		for( int i=4; i<16; i++ ) nw2cEps_intgVolNi[ i ] = nw2cEpsVol * ( ( 128.0 / 420.0 ) - ( 32.0 / 120.0 ) );

		nw2cEps_intgVolNi[ 16 ] = nw2cEpsVol * 1536.0 / 5040.0;

		// Derivatives 1st order basis
        double DN[ 3 ][ 4 ]; Calculate_DN( DN );

		// [Kvv]
		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
				StiffMatrix[ i + 51 ][ j + 51 ] = nw2cEpsVol * ( DN[ 0 ][ i ]*DN[ 0 ][ j ] + DN[ 1 ][ i ]*DN[ 1 ][ j ] + DN[ 2 ][ i ]*DN[ 2 ][ j ] );
        	}
        }

		// [Kxv], [Kyv], [Kzv]
		for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<4; j++ ) 
        	{
        		StiffMatrix[ i               ][ j + 3*mNumNodes ] = nw2cEps_intgVolNi[ i ] * DN[ 0 ][ j ];
        		StiffMatrix[ i +   mNumNodes ][ j + 3*mNumNodes ] = nw2cEps_intgVolNi[ i ] * DN[ 1 ][ j ];
        		StiffMatrix[ i + 2*mNumNodes ][ j + 3*mNumNodes ] = nw2cEps_intgVolNi[ i ] * DN[ 2 ][ j ];
        	}
        }

        // [Kvx], [Kvy], [Kvz]
		for( int i=0; i<4; i++ )
        {
            for( int j=0; j<mNumNodes; j++ ) 
        	{
        		StiffMatrix[ i + 3*mNumNodes ][ j               ] = nw2cEps_intgVolNi[ j ] * DN[ 0 ][ i ];
        		StiffMatrix[ i + 3*mNumNodes ][ j +   mNumNodes ] = nw2cEps_intgVolNi[ j ] * DN[ 1 ][ i ];
        		StiffMatrix[ i + 3*mNumNodes ][ j + 2*mNumNodes ] = nw2cEps_intgVolNi[ j ] * DN[ 2 ][ i ];
        	}
        }
	}

    //********************************************************************************************************************************
    // - Rotational P1 projection on vertice nodes
    //********************************************************************************************************************************
	void VolumeElement_3sb_FullWave::Calculate_Rot_P1_Projection( cVectorType& cRotP1x, cVectorType& cRotP1y, cVectorType& cRotP1z )
	{
		// Field on nodes
		Matrix< std::complex<double> > cF( 51, 1 );
		     
		if( mPotentialsOn )
		{
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        cF[ n               ][ 0 ] = (*mProperties)( cAx, *mNodes[ n ] );
		        cF[ n +   mNumNodes ][ 0 ] = (*mProperties)( cAy, *mNodes[ n ] );
		        cF[ n + 2*mNumNodes ][ 0 ] = (*mProperties)( cAz, *mNodes[ n ] );
		    }															   
 		}																   
 		else															   
 		{																   
 		    for( int n=0; n<mNumNodes; n++ )							   
 		    {															   
 		        cF[ n               ][ 0 ] = (*mProperties)( cEx, *mNodes[ n ] );
		        cF[ n +   mNumNodes ][ 0 ] = (*mProperties)( cEy, *mNodes[ n ] );
		        cF[ n + 2*mNumNodes ][ 0 ] = (*mProperties)( cEz, *mNodes[ n ] );
		    }		
		}

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

		// Building projection matrices (curl contribution)
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
                C_cur[ i     ][ j + mNumNodes   ] = -intg_Ni_dNjdZ; 
                C_cur[ i     ][ j + mNumNodes*2 ] = +intg_Ni_dNjdY; 
		    								    
		    	C_cur[ i + 4 ][ j               ] = +intg_Ni_dNjdZ; 
                C_cur[ i + 4 ][ j + mNumNodes   ] =            0.0; 
                C_cur[ i + 4 ][ j + mNumNodes*2 ] = -intg_Ni_dNjdX; 
		    								    
		    	C_cur[ i + 8 ][ j               ] = -intg_Ni_dNjdY; 
                C_cur[ i + 8 ][ j + mNumNodes   ] = +intg_Ni_dNjdX; 
                C_cur[ i + 8 ][ j + mNumNodes*2 ] =            0.0; 
		    }
		}

		// Adding surface integral nxu*q
		Add_C_NxN_Matrix( C_cur ); 

		// Projection q basis matrices
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

		// Rotational P1 projection on vertices
		Matrix< std::complex<double> > cRotP1Proj( 12, 1 ); 
		
		cRotP1Proj = ( iA_cur * C_cur ) * cF;

		cRotP1x.resize( 4 );
		cRotP1y.resize( 4 );
		cRotP1z.resize( 4 );

        for( int i=0; i<4; i++ )
		{
			cRotP1x[ i ] = cRotP1Proj[ i     ][ 0 ];
			cRotP1y[ i ] = cRotP1Proj[ i + 4 ][ 0 ];
			cRotP1z[ i ] = cRotP1Proj[ i + 8 ][ 0 ];
		}		
	}

    //********************************************************************************************************************************
    // - Rotational on nodes
    //********************************************************************************************************************************
	void VolumeElement_3sb_FullWave::Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes )
	{
		// Rotational on vertices (P1 projection)
		cVectorType cRotP1x, cRotP1y, cRotP1z;

		Calculate_Rot_P1_Projection( cRotP1x, cRotP1y, cRotP1z );
		
		// P1 basis functions on nodal points
		Vector<double> nX( mNumNodes );
		Vector<double> nY( mNumNodes );
		Vector<double> nZ( mNumNodes );

		// Vertices 0,1,2,3
		nX[ 0 ] = 0.00; nY[ 0 ] = 0.00; nZ[ 0 ] = 0.00;
		nX[ 1 ] = 1.00; nY[ 1 ] = 0.00; nZ[ 1 ] = 0.00;
		nX[ 2 ] = 0.00; nY[ 2 ] = 1.00; nZ[ 2 ] = 0.00;
		nX[ 3 ] = 0.00; nY[ 3 ] = 0.00; nZ[ 3 ] = 1.00;
							  	 	 	  		 
		// Face [0 1 2]		  	 	 	  		 
		nX[ 4 ] = 0.25; nY[ 4 ] = 0.25; nZ[ 4 ] = 0.00;
		nX[ 5 ] = 0.50; nY[ 5 ] = 0.25; nZ[ 5 ] = 0.00;
		nX[ 6 ] = 0.25; nY[ 6 ] = 0.50; nZ[ 6 ] = 0.00;
							  	 	 	  		 
		// Face [0 1 3]		  	 	 	  		 
		nX[ 7 ] = 0.25; nY[ 7 ] = 0.00; nZ[ 7 ] = 0.25;
		nX[ 8 ] = 0.50; nY[ 8 ] = 0.00; nZ[ 8 ] = 0.25;
		nX[ 9 ] = 0.25; nY[ 9 ] = 0.00; nZ[ 9 ] = 0.50;
						   
		// Face [1 2 3]	   
		nX[ 10 ] = 0.50; nY[ 10 ] = 0.25; nZ[ 10 ] = 0.25;
		nX[ 11 ] = 0.25; nY[ 11 ] = 0.50; nZ[ 11 ] = 0.25;
		nX[ 12 ] = 0.25; nY[ 12 ] = 0.25; nZ[ 12 ] = 0.50;
						   
		// Face [0 2 3]	   
		nX[ 13 ] = 0.00; nY[ 13 ] = 0.25; nZ[ 13 ] = 0.25;
		nX[ 14 ] = 0.00; nY[ 14 ] = 0.50; nZ[ 14 ] = 0.25;
		nX[ 15 ] = 0.00; nY[ 15 ] = 0.25; nZ[ 15 ] = 0.50;

		// Volume [0 1 2 3]
		nX[ 16 ] = 0.25; nY[ 16 ] = 0.25; nZ[ 16 ] = 0.25; 

		// P1 basis on nodes
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, nX, nY, nZ );

		// Rotational in nodes
		rot_OnNodes.resize( mNumNodes );

		for( int n=0; n<mNumNodes; n++ )
		{
		    Vector< std::complex<double> > ncRot( 3, std::complex<double>( 0.0, 0.0 ) );
		 
		    for( int i=0; i<4; i++ )
		    {
		        ncRot[ 0 ] += cRotP1x[ i ] * N[ i ][ n ];
			    ncRot[ 1 ] += cRotP1y[ i ] * N[ i ][ n ];
			    ncRot[ 2 ] += cRotP1z[ i ] * N[ i ][ n ];
		    }

		    rot_OnNodes[ n ] = ncRot;
		}		
	}
	  
    //********************************************************************************************************************************
    // - Rotational on Gauss points
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs )
	{
		// Rotational on vertices (P1 projection)
		cVectorType cRotP1x, cRotP1y, cRotP1z;

		Calculate_Rot_P1_Projection( cRotP1x, cRotP1y, cRotP1z );

		// Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGaussPoints = GetInnerGiDGaussPoints( cX, cY, cZ, numResultsOnGPs );

		// P1 basis on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

		// Rotational on Gauss points
		rot_OnGPs.resize( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
		    Vector< std::complex<double> > gpcRot( 3, std::complex<double>( 0.0, 0.0 ) );
		 
		    for( int i=0; i<4; i++ )
		    {
		        gpcRot[ 0 ] += cRotP1x[ i ] * N[ i ][ gp ];
			    gpcRot[ 1 ] += cRotP1y[ i ] * N[ i ][ gp ];
			    gpcRot[ 2 ] += cRotP1z[ i ] * N[ i ][ gp ];
		    }

		    rot_OnGPs[ gp ] = gpcRot;
		}
	}

    //********************************************************************************************************************************
    // - Electric field on nodes
    //********************************************************************************************************************************
	void VolumeElement_3sb_FullWave::Calculate_E_field_OnNodes( cVector2Type& Ef_OnNodes )
	{
		Ef_OnNodes.resize( mNumNodes );

		for( int n=0; n<mNumNodes; n++ )
		{
		    Ef_OnNodes[ n ].resize( 3 );
		}

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

            for( int i=0; i<4; i++ )
            {
                cGradVs[ 0 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 0 ][ i ];
                cGradVs[ 1 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 1 ][ i ];
                cGradVs[ 2 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 2 ][ i ];
            }

		    for( int n=0; n<mNumNodes; n++ )
		    {
		        Ef_OnNodes[ n ][ 0 ] = jw * ( (*mProperties)( cAx, *mNodes[ n ] ) + cGradVs[ 0 ] );
		    	Ef_OnNodes[ n ][ 1 ] = jw * ( (*mProperties)( cAy, *mNodes[ n ] ) + cGradVs[ 1 ] );
		    	Ef_OnNodes[ n ][ 2 ] = jw * ( (*mProperties)( cAz, *mNodes[ n ] ) + cGradVs[ 2 ] );
		    }	
		}
		else
		{
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        Ef_OnNodes[ n ][ 0 ] = (*mProperties)( cEx, *mNodes[ n ] );
		    	Ef_OnNodes[ n ][ 1 ] = (*mProperties)( cEy, *mNodes[ n ] );
		    	Ef_OnNodes[ n ][ 2 ] = (*mProperties)( cEz, *mNodes[ n ] );
		    }			
		}
	}

    //********************************************************************************************************************************
    // - Electric field on Gauss points
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs )
	{
        // Vector magnitude on nodes [ E | A ]
		Matrix< std::complex<double> > cF( 51, 1 );

		if( mPotentialsOn )
		{
		    for( int n=0; n<mNumNodes; n++ )
		    {
		        cF[ n               ][ 0 ] = (*mProperties)( cAx, *mNodes[ n ] );
		        cF[ n +   mNumNodes ][ 0 ] = (*mProperties)( cAy, *mNodes[ n ] );
		        cF[ n + 2*mNumNodes ][ 0 ] = (*mProperties)( cAz, *mNodes[ n ] );
		    }															   
 		}																   
 		else															   
 		{																   
 		    for( int n=0; n<mNumNodes; n++ )							   
 		    {															   
 		        cF[ n               ][ 0 ] = (*mProperties)( cEx, *mNodes[ n ] );
		        cF[ n +   mNumNodes ][ 0 ] = (*mProperties)( cEy, *mNodes[ n ] );
		        cF[ n + 2*mNumNodes ][ 0 ] = (*mProperties)( cEz, *mNodes[ n ] );
		    }		
		}

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

		// Building projection matrices 
		Matrix< std::complex<double> > matrix_NiNj( 12, 51 );

		for( int i=0; i<4; i++ )
		{
		    for( int j=0; j<mNumNodes; j++ )
		    {
		    	double intg_NiNj = 0.0;
		    
                for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_NiNj += W[ gp ] * N[ i ][ gp ] * N[ j ][ gp ]; 
                }      
		    	
				matrix_NiNj[ i     ][ j               ] = intg_NiNj; 
                matrix_NiNj[ i + 4 ][ j +   mNumNodes ] = intg_NiNj;  						 
                matrix_NiNj[ i + 8 ][ j + 2*mNumNodes ] = intg_NiNj; 
		    }
		}

		// Projection basis matrices
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

		// L2 projection on nodes 
		Matrix< std::complex<double> > FProj( 12, 1 ); 
		
		FProj += ( iA_cur * matrix_NiNj ) * cF;

		Vector< std::complex<double> > pnFx( 4 );
		Vector< std::complex<double> > pnFy( 4 );
		Vector< std::complex<double> > pnFz( 4 );

		for( int i=0; i<4; i++ )
		{
			pnFx[ i ] = FProj[ i     ][ 0 ];
			pnFy[ i ] = FProj[ i + 4 ][ 0 ];
			pnFz[ i ] = FProj[ i + 8 ][ 0 ];
		}

		// GiD internal Gauss points
		Vector<double> cXig, cYig, cZig;
	
		int numGaussPoints_ig = GetInnerGiDGaussPoints( cXig, cYig, cZig, numResultsOnGPs );

        // Basis functions on Gauss points
		Matrix<double> Nig;

		Lagrange3D_Ni_1st( Nig, cXig, cYig, cZig );

	    // Electric field on Gauss points
		Ef_OnGPs.resize( numGaussPoints_ig );

		for( int ig=0; ig<numGaussPoints_ig; ig++ )
		{
			Vector< std::complex<double> > cF_ig( 3, std::complex<double>( 0.0, 0.0 ) );

			for( int n=0; n<4; n++ )
		    {       
				cF_ig[ 0 ] += pnFx[ n ] * Nig[ n ][ ig ];
			    cF_ig[ 1 ] += pnFy[ n ] * Nig[ n ][ ig ];
			    cF_ig[ 2 ] += pnFz[ n ] * Nig[ n ][ ig ];
			}

			Ef_OnGPs[ ig ] = cF_ig;
		}

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
		    
            for( int i=0; i<4; i++ )
            {
                cGradVs[ 0 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 0 ][ i ];
                cGradVs[ 1 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 1 ][ i ];
                cGradVs[ 2 ] += (*mProperties)( cVs, *mNodes[ i ] ) * DN[ 2 ][ i ];
            }
		    
		    for( int ig=0; ig<numGaussPoints_ig; ig++ )
		    {
		        Ef_OnGPs[ ig ] += cGradVs;
		    	Ef_OnGPs[ ig ] *= jw;
		    }
		}
	}

    //********************************************************************************************************************************
    // - Inner GiD Gauss points 
    //********************************************************************************************************************************
    int VolumeElement_3sb_FullWave::GetInnerGiDGaussPoints( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, int numResultsOnGPs )
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