
#include "VolumeElement_3sb_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{
    //********************************************************************************************************************************
    // - Global index of the nodes
    //********************************************************************************************************************************	
    void VolumeElement_3sb_FullWave::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( mNumDofs );

        for( int i=0; i<mNumNodes; i++ ) 
        {  
            EquationId[ i               ] = mNodes[ i ]->pDofcEx()->EquationId();
            EquationId[ i + mNumNodes   ] = mNodes[ i ]->pDofcEy()->EquationId();
			EquationId[ i + mNumNodes*2 ] = mNodes[ i ]->pDofcEz()->EquationId();
        }
    }

    //********************************************************************************************************************************
    // - Calculates the volume of the element
    //********************************************************************************************************************************
    double VolumeElement_3sb_FullWave::Calculate_Volume()
    {
        double det;

        det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( det / 6.00 );
    }

    //********************************************************************************************************************************
    // - Calculates dN/dx, dN/dy, dN/dz - DN[ X,Y,Z ][ 0,1,2,3 ]
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::Calculate_Ni_Derivatives( Matrix<double>& DN )
    {
        // [ Coord: X,Y,Z ][ Local_node: 0,1,2,3 ]
        DN.Resize( 3, 4 );

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

    //********************************************************************************************************************************
    // - Calculation of dN/dx, dN/dy, dN/dz on cX, cY, cZ points
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::LagrangeDerivatives3D_3sb( Matrix<double>&    dNdx,
		                                                        Matrix<double>&    dNdy,
													            Matrix<double>&    dNdz,
													            std::vector<double>& cX, 
													            std::vector<double>& cY, 
													            std::vector<double>& cZ )
    {
		int nPoints = cX.size();

		dNdx.Resize( mNumNodes, nPoints );
		dNdy.Resize( mNumNodes, nPoints );
		dNdz.Resize( mNumNodes, nPoints );
		
		Matrix<double> DN; 
		
		Calculate_Ni_Derivatives( DN );

		double dL1dx, dL2dx, dL3dx, dL4dx;
        double dL1dy, dL2dy, dL3dy, dL4dy;
		double dL1dz, dL2dz, dL3dz, dL4dz;

        dL1dx = DN[0][0]; dL2dx = DN[0][1]; dL3dx = DN[0][2]; dL4dx = DN[0][3];
        dL1dy = DN[1][0]; dL2dy = DN[1][1]; dL3dy = DN[1][2]; dL4dy = DN[1][3];
		dL1dz = DN[2][0]; dL2dz = DN[2][1]; dL3dz = DN[2][2]; dL4dz = DN[2][3];

		// dNi/dx ,dNi/dy, dNi/dz
        for( int p=0; p<nPoints; p++ ) 
        {
			// Natural coordinates
			double L2 = cX[p];
            double L3 = cY[p];
			double L4 = cZ[p];
            double L1 = 1.00 - L2 - L3 - L4;

			// Node [ 0 1 2 3 ] 
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

			// Face bubble [ 0 1 2 ]
			dNdx[ 4 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL1dx ) ) + ( L1 * L2 * dL3dx * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL2dx * L3 * ( 4.00*L1 - 1.00 ) ) + ( dL1dx * L2 * L3 * ( 4.00*L1 - 1.00 ) ) );
			dNdy[ 4 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL1dy ) ) + ( L1 * L2 * dL3dy * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL2dy * L3 * ( 4.00*L1 - 1.00 ) ) + ( dL1dy * L2 * L3 * ( 4.00*L1 - 1.00 ) ) );
			dNdz[ 4 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL1dz ) ) + ( L1 * L2 * dL3dz * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL2dz * L3 * ( 4.00*L1 - 1.00 ) ) + ( dL1dz * L2 * L3 * ( 4.00*L1 - 1.00 ) ) );
						 			  
			dNdx[ 5 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL2dx ) ) + ( L1 * L2 * dL3dx * ( 4.00*L2 - 1.00 ) ) + ( L1 * dL2dx * L3 * ( 4.00*L2 - 1.00 ) ) + ( dL1dx * L2 * L3 * ( 4.00*L2 - 1.00 ) ) );
			dNdy[ 5 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL2dy ) ) + ( L1 * L2 * dL3dy * ( 4.00*L2 - 1.00 ) ) + ( L1 * dL2dy * L3 * ( 4.00*L2 - 1.00 ) ) + ( dL1dy * L2 * L3 * ( 4.00*L2 - 1.00 ) ) );
			dNdz[ 5 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL2dz ) ) + ( L1 * L2 * dL3dz * ( 4.00*L2 - 1.00 ) ) + ( L1 * dL2dz * L3 * ( 4.00*L2 - 1.00 ) ) + ( dL1dz * L2 * L3 * ( 4.00*L2 - 1.00 ) ) );
						 			  
			dNdx[ 6 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL3dx ) ) + ( L1 * L2 * dL3dx * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL2dx * L3 * ( 4.00*L3 - 1.00 ) ) + ( dL1dx * L2 * L3 * ( 4.00*L3 - 1.00 ) ) );
			dNdy[ 6 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL3dy ) ) + ( L1 * L2 * dL3dy * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL2dy * L3 * ( 4.00*L3 - 1.00 ) ) + ( dL1dy * L2 * L3 * ( 4.00*L3 - 1.00 ) ) );
			dNdz[ 6 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 4.00 * dL3dz ) ) + ( L1 * L2 * dL3dz * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL2dz * L3 * ( 4.00*L3 - 1.00 ) ) + ( dL1dz * L2 * L3 * ( 4.00*L3 - 1.00 ) ) );

			// Face bubble [ 0 1 3 ]
			dNdx[ 7 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL1dx ) ) + ( L1 * L2 * dL4dx * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL2dx * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dx * L2 * L4 * ( 4.00*L1 - 1.00 ) ) );
			dNdy[ 7 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL1dy ) ) + ( L1 * L2 * dL4dy * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL2dy * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dy * L2 * L4 * ( 4.00*L1 - 1.00 ) ) );
			dNdz[ 7 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL1dz ) ) + ( L1 * L2 * dL4dz * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL2dz * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dz * L2 * L4 * ( 4.00*L1 - 1.00 ) ) );
						 			  		
			dNdx[ 8 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL2dx ) ) + ( L1 * L2 * dL4dx * ( 4.00*L2 - 1.00 ) ) + ( L1 * dL2dx * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL1dx * L2 * L4 * ( 4.00*L2 - 1.00 ) ) );
			dNdy[ 8 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL2dy ) ) + ( L1 * L2 * dL4dy * ( 4.00*L2 - 1.00 ) ) + ( L1 * dL2dy * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL1dy * L2 * L4 * ( 4.00*L2 - 1.00 ) ) );
			dNdz[ 8 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL2dz ) ) + ( L1 * L2 * dL4dz * ( 4.00*L2 - 1.00 ) ) + ( L1 * dL2dz * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL1dz * L2 * L4 * ( 4.00*L2 - 1.00 ) ) );
						 			  			
			dNdx[ 9 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL4dx ) ) + ( L1 * L2 * dL4dx * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL2dx * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dx * L2 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdy[ 9 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL4dy ) ) + ( L1 * L2 * dL4dy * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL2dy * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dy * L2 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdz[ 9 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 4.00 * dL4dz ) ) + ( L1 * L2 * dL4dz * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL2dz * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dz * L2 * L4 * ( 4.00*L4 - 1.00 ) ) );

			// Face bubble [ 1 2 3 ]
			dNdx[ 10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL2dx ) ) + ( L2 * L3 * dL4dx * ( 4.00*L2 - 1.00 ) ) + ( L2 * dL3dx * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL2dx * L3 * L4 * ( 4.00*L2 - 1.00 ) ) );
			dNdy[ 10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL2dy ) ) + ( L2 * L3 * dL4dy * ( 4.00*L2 - 1.00 ) ) + ( L2 * dL3dy * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL2dy * L3 * L4 * ( 4.00*L2 - 1.00 ) ) );
			dNdz[ 10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL2dz ) ) + ( L2 * L3 * dL4dz * ( 4.00*L2 - 1.00 ) ) + ( L2 * dL3dz * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL2dz * L3 * L4 * ( 4.00*L2 - 1.00 ) ) );
							 			
			dNdx[ 11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL3dx ) ) + ( L2 * L3 * dL4dx * ( 4.00*L3 - 1.00 ) ) + ( L2 * dL3dx * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL2dx * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
			dNdy[ 11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL3dy ) ) + ( L2 * L3 * dL4dy * ( 4.00*L3 - 1.00 ) ) + ( L2 * dL3dy * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL2dy * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
			dNdz[ 11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL3dz ) ) + ( L2 * L3 * dL4dz * ( 4.00*L3 - 1.00 ) ) + ( L2 * dL3dz * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL2dz * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
							 			  
			dNdx[ 12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL4dx ) ) + ( L2 * L3 * dL4dx * ( 4.00*L4 - 1.00 ) ) + ( L2 * dL3dx * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL2dx * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdy[ 12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL4dy ) ) + ( L2 * L3 * dL4dy * ( 4.00*L4 - 1.00 ) ) + ( L2 * dL3dy * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL2dy * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdz[ 12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL4dz ) ) + ( L2 * L3 * dL4dz * ( 4.00*L4 - 1.00 ) ) + ( L2 * dL3dz * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL2dz * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );

			// Face bubble [ 0 2 3 ]
			dNdx[ 13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL1dx ) ) + ( L1 * L3 * dL4dx * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL3dx * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dx * L3 * L4 * ( 4.00*L1 - 1.00 ) ) );
			dNdy[ 13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL1dy ) ) + ( L1 * L3 * dL4dy * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL3dy * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dy * L3 * L4 * ( 4.00*L1 - 1.00 ) ) );
			dNdz[ 13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL1dz ) ) + ( L1 * L3 * dL4dz * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL3dz * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dz * L3 * L4 * ( 4.00*L1 - 1.00 ) ) );
							 			  
			dNdx[ 14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL3dx ) ) + ( L1 * L3 * dL4dx * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL3dx * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL1dx * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
			dNdy[ 14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL3dy ) ) + ( L1 * L3 * dL4dy * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL3dy * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL1dy * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
			dNdz[ 14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL3dz ) ) + ( L1 * L3 * dL4dz * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL3dz * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL1dz * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
							 			 
			dNdx[ 15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL4dx ) ) + ( L1 * L3 * dL4dx * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL3dx * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dx * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdy[ 15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL4dy ) ) + ( L1 * L3 * dL4dy * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL3dy * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dy * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdz[ 15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL4dz ) ) + ( L1 * L3 * dL4dz * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL3dz * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dz * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );

			// Inner bubble
			dNdx[ 16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dx ) + ( L1 * L2 * dL3dx * L4 ) + ( L1 * dL2dx * L3 * L4 ) + ( dL1dx * L2 * L3 * L4 ) );
			dNdy[ 16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dy ) + ( L1 * L2 * dL3dy * L4 ) + ( L1 * dL2dy * L3 * L4 ) + ( dL1dy * L2 * L3 * L4 ) );
			dNdz[ 16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dz ) + ( L1 * L2 * dL3dz * L4 ) + ( L1 * dL2dz * L3 * L4 ) + ( dL1dz * L2 * L3 * L4 ) );
        }
	}

	//********************************************************************************************************************************
    // - Calculates the stiffness matrix 
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::GetStiffMatrix( Matrix< std::complex<double> >& StiffMatrix )
	{   
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.00, 0.00 ) );

		Matrix< std::complex<double> > StiffMatrix_CDE;

		GetStiffMatrix_LL2P_CurlDivE( StiffMatrix_CDE ); 

		StiffMatrix += StiffMatrix_CDE;

		if ( mhkFactor != 0.0 )
		{
		    Matrix< std::complex<double> > StiffMatrix_hk2;

		    Add_hk2CurlDiv_Stabilization( StiffMatrix_hk2 ); 
		    
		    if ( mhkFactor != 1.0 ) StiffMatrix_hk2 *= mhkFactor;
		    
		    StiffMatrix += StiffMatrix_hk2;
		}
	}

	//********************************************************************************************************************************
    //* - Calculates element diameter h_k
    //********************************************************************************************************************************
	double VolumeElement_3sb_FullWave::Calculate_Element_Diameter() 
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

	//********************************************************************************************************************************
    //* - Adds a stabilization term: hk^2 * mu-1 * ( curl_F, curl_E ) * ( div_F, div_E ) 
    //********************************************************************************************************************************	
    void VolumeElement_3sb_FullWave::Add_hk2CurlDiv_Stabilization( Matrix< std::complex<double> >& StiffMatrix ) 
    {  
		// Resize matrix
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.00, 0.00 ) );
		
		// Loop counters
		int i, j, gp;

		// Material properties
        double freq  = (*mProperties)(FREQUENCY);
        double sigma = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);

		double eo = 8.8541878176e-12;
		double pi = 3.1415926535897932384626433832795;
		double mo = pi * 4.00e-7;

        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / freq ) );
        std::complex<double> cMu( mu_real , mu_imag );

		// Gaussian derivation weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order5( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Derivatives
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_3sb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Diameter of the element
		double hk = Calculate_Element_Diameter();
		
		std::complex<double> cteCurlDiv = ( hk * hk ) / cMu;

		// CurlDiv volumetric integral 
		double intCurlDiv;

		// [Kxx], [Kyy], [Kzz] 
		for( i=0; i<mNumNodes; i++ )
		{
			for( j=i; j<mNumNodes; j++ )
			{
				intCurlDiv = 0.00; 
				
				for( gp=0; gp<nGaussPoints; gp++ )
				{
					intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
				}

				StiffMatrix[ i               ][ j               ] += cteCurlDiv * intCurlDiv;
				StiffMatrix[ i + mNumNodes   ][ j + mNumNodes   ] += cteCurlDiv * intCurlDiv;
				StiffMatrix[ i + mNumNodes*2 ][ j + mNumNodes*2 ] += cteCurlDiv * intCurlDiv;
			}
		}

		// [Kxy] 
		for( i=0; i<mNumNodes; i++ )
		{
			for( j=0; j<mNumNodes; j++ )
			{
				intCurlDiv = 0.00;

				for( gp=0; gp<nGaussPoints; gp++ )
				{
				    intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdy[j][gp] - dNdy[i][gp] * dNdx[j][gp] );
				}
	            
				StiffMatrix[ i ][ j + mNumNodes ] += cteCurlDiv * intCurlDiv;
			}
		}

		// [Kxz] 
		for( i=0; i<mNumNodes; i++ )
		{
			for( j=0; j<mNumNodes; j++ )
			{
				intCurlDiv = 0.00;

				for( gp=0; gp<nGaussPoints; gp++ )
				{
					intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdx[j][gp] );
				}
	            
				StiffMatrix[ i ][ j + mNumNodes*2 ] += cteCurlDiv * intCurlDiv;
			}
		}

		// [Kyz] 
		for( i=0; i<mNumNodes; i++ )
		{
			for( j=0; j<mNumNodes; j++ )
			{
				intCurlDiv = 0.00; 

				for( gp=0; gp<nGaussPoints; gp++ )
				{
					intCurlDiv += W[gp] * ( dNdy[i][gp] * dNdz[j][gp] - dNdz[i][gp] * dNdy[j][gp] );
				}
	            
				StiffMatrix[ i + mNumNodes ][ j + mNumNodes*2 ] += cteCurlDiv * intCurlDiv;
			}
		}
		
        // Symmetric block
        for( i=0; i<mNumDofs; i++ )
        {
            for( j=i+1; j<mNumDofs; j++ )
            {
                StiffMatrix[ j ][ i ] +=  StiffMatrix[ i ][ j ];
            }
        }  
		
		//Adding surface integral
		Vector< Vector<int> > Faces; Faces.resize( 4 );

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

			// Check if face is not PEC
			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) 
			{
				continue;
			}

			// Gauss points and weights
		    Vector<double> cXs, cYs, Ws;
		    
		    int nGaussPoints = GaussPoints2D_Order12( cXs, cYs, Ws );
		    
            // Jacobian
		    double detJs = 2.00 * Calculate_Area( nf );
		    
		    for( int gp=0; gp<nGaussPoints; gp++ ) Ws[gp] *= detJs; 
		    
            // Lagrangian base on Gauss points
		    Matrix<double> Ns;
		    
		    Lagrange2D_Ni_3sb( Ns, cXs, cYs );

		    Vector<double> extN; 
			
			Calculate_ExtNormal( extN, nf );

			double X1 = mNodes[nf[0]]->X(), Y1 = mNodes[nf[0]]->Y(), Z1 = mNodes[nf[0]]->Z();
		    double X2 = mNodes[nf[1]]->X(), Y2 = mNodes[nf[1]]->Y(), Z2 = mNodes[nf[1]]->Z();
		    double X3 = mNodes[nf[2]]->X(), Y3 = mNodes[nf[2]]->Y(), Z3 = mNodes[nf[2]]->Z();
		    
		    // Element edge lenghts 
		    double Lenght_12 = sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) + (Z1-Z2)*(Z1-Z2) );
		    double Lenght_13 = sqrt( (X1-X3)*(X1-X3) + (Y1-Y3)*(Y1-Y3) + (Z1-Z3)*(Z1-Z3) );
		    double Lenght_23 = sqrt( (X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3) + (Z2-Z3)*(Z2-Z3) );
		    
            // Maximum lenght 
		    double hF = Lenght_12;
		    
		    if ( Lenght_13 > hF ) hF = Lenght_13;
		    if ( Lenght_23 > hF ) hF = Lenght_23;

			std::complex<double> cteNxN = ( hF / cMu );

			for( int i=0; i<6; i++ )
			{
				for( int j=0; j<6; j++ )
				{
		       		std::complex<double> Sintg_NsiNsj = 0.00;
		           
					for( int gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_NsiNsj += Ws[gp] * Ns[i][gp] * Ns[j][gp]; 
					}

					Sintg_NsiNsj *= cteNxN;
		       		
					StiffMatrix[ nf[i]      ][ nf[j]      ] += ( 1.0 - extN[0] * extN[0] ) * Sintg_NsiNsj; 
					StiffMatrix[ nf[i] +  6 ][ nf[j]      ] += (     - extN[1] * extN[0] ) * Sintg_NsiNsj; 
					StiffMatrix[ nf[i] + 12 ][ nf[j]      ] += (     - extN[2] * extN[0] ) * Sintg_NsiNsj; 

					StiffMatrix[ nf[i]      ][ nf[j] +  6 ] += (     - extN[0] * extN[1] ) * Sintg_NsiNsj; 
					StiffMatrix[ nf[i] +  6 ][ nf[j] +  6 ] += ( 1.0 - extN[1] * extN[1] ) * Sintg_NsiNsj; 
					StiffMatrix[ nf[i] + 12 ][ nf[j] +  6 ] += (     - extN[2] * extN[1] ) * Sintg_NsiNsj; 

					StiffMatrix[ nf[i]      ][ nf[j] + 12 ] += (     - extN[0] * extN[2] ) * Sintg_NsiNsj; 
					StiffMatrix[ nf[i] +  6 ][ nf[j] + 12 ] += (     - extN[1] * extN[2] ) * Sintg_NsiNsj; 
					StiffMatrix[ nf[i] + 12 ][ nf[j] + 12 ] += ( 1.0 - extN[2] * extN[2] ) * Sintg_NsiNsj; 
				}
			}
		}
	}
	
	//********************************************************************************************************************************
    //* - Calculates the exterior normal 
    //********************************************************************************************************************************	
    void VolumeElement_3sb_FullWave::Calculate_ExtNormal( Vector<double>& extN, Vector<int>& face )
    {
        double v1[3], v2[3];

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

    //********************************************************************************************************************************
    //* - Calculates the area of the element
    //********************************************************************************************************************************	
    double VolumeElement_3sb_FullWave::Calculate_Area( Vector<int>& face ) 
    {
		double n[3], v1[3], v2[3];

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

    //********************************************************************************************************************************
    //* - Adds surface integral (nxu)*q in PEC faces to C_cur matrix
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
			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) 
			{
				continue;
			}

			// Gauss points and weights
		    Vector<double> cX, cY, W;
		    
		    int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );
		    
            // Jacobian
		    double detJ = 2.00 * Calculate_Area( nf );
		    
		    for( int gp=0; gp<numGaussPoints; gp++ ) W[gp] *= detJ; 
		    
            // Lagrangian base on Gauss points
		    Matrix<double> N;
		    
		    Lagrange2D_Ni_3sb( N, cX, cY );

		    Vector<double> extN; 
			
			Calculate_ExtNormal( extN, nf );

			for( int i=0; i<3; i++ )
		    {
		        for( int j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_NiNj( 0.00, 0.00 );
		            
                    for( int gp=0; gp<numGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}
   		
                    C_cur[ nf[i]     ][ nf[j] + mNumNodes   ] += ( extN[2] * Sintg_NiNj ); 
                    C_cur[ nf[i]     ][ nf[j] + mNumNodes*2 ] -= ( extN[1] * Sintg_NiNj ); 
		        					    	 					   					    
		        	C_cur[ nf[i] + 4 ][ nf[j]               ] -= ( extN[2] * Sintg_NiNj ); 
                    C_cur[ nf[i] + 4 ][ nf[j] + mNumNodes*2 ] += ( extN[0] * Sintg_NiNj );
		        					 	 						   					    
		        	C_cur[ nf[i] + 8 ][ nf[j]               ] += ( extN[1] * Sintg_NiNj ); 
                    C_cur[ nf[i] + 8 ][ nf[j] + mNumNodes   ] -= ( extN[0] * Sintg_NiNj ); 
		        }
		    }
		}
	}

	//********************************************************************************************************************************
    //* - Adds surface integral (n*u)*q in PMC faces to D_div matrix
    //********************************************************************************************************************************	
	void VolumeElement_3sb_FullWave::Add_D_NxN_Matrix( Matrix< std::complex<double> >& D_div ) 
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

			// Check if the surface is not PMC ( 'H' = PMC || 'T' = TEPMC || 'R' = R_node-CtC || 'L' = L_node-CtC )
			if ( ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'H' ) && ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'T' ) )
			{
				continue;
			}

			// Gauss points and weights
		    Vector<double> cX, cY, W;
		    
		    int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );
		    
            // Jacobian
		    double detJ = 2.00 * Calculate_Area( nf );
		    
		    for( int gp=0; gp<numGaussPoints; gp++ ) W[gp] *= detJ; 
		    
            // Lagrangian base on Gauss points
		    Matrix<double> N;
		    
		    Lagrange2D_Ni_3sb( N, cX, cY );

		    Vector<double> extN; 
			
			Calculate_ExtNormal( extN, nf );

			for( int i=0; i<3; i++ )
		    {
		        for( int j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_NiNj( 0.00, 0.00 );
		            
                    for( int gp=0; gp<numGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}
   		
					D_div[ nf[i] ][ nf[j]               ] += ( extN[0] * Sintg_NiNj );
		    	    D_div[ nf[i] ][ nf[j] + mNumNodes   ] += ( extN[1] * Sintg_NiNj );
		    	    D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] += ( extN[2] * Sintg_NiNj );
		        }
		    }
		}
	}

	//********************************************************************************************************************************
    //* - Stiffness matrix with local L2 projection method ( div(v), q ) - ( curl(v), w )
    //********************************************************************************************************************************	
    void VolumeElement_3sb_FullWave::GetStiffMatrix_LL2P_CurlDivE( Matrix< std::complex<double> >& StiffMatrix ) 
    { 
		// Resizing local stiffness matrix 
		StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.00, 0.00 ) );

		// Material properties
        double freq  = (*mProperties)(FREQUENCY);
        double sigma = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);

		double eo = 8.8541878176e-12;
		double pi = 3.1415926535897932384626433832795;
		double mo = pi * 4.00e-7;

        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / freq ) );
        std::complex<double> cMu( mu_real , mu_imag );

		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order5( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( int gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_3sb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_3sb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< std::complex<double> > D_div(  4, mNumDofs );
		Matrix< std::complex<double> > C_cur( 12, mNumDofs );

		for( int i=0; i<4; i++ )
		{
		    for( int j=0; j<mNumNodes; j++ )
		    {
		    	double vintg_Ni_dNjdX = 0.00;
		        double vintg_Ni_dNjdY = 0.00;
		        double vintg_Ni_dNjdZ = 0.00;
		    
                for( int gp=0; gp<nGaussPoints; gp++ )
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
	
		// Adding surface integral (nxu)*q
		Add_C_NxN_Matrix( C_cur ); 

		// Adding surface integral (n*u)*q
		// Add_D_NxN_Matrix( D_div ); 

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij( mNumDofs, mNumDofs );
		
		std::complex<double> cteNiNj = freq * freq * cEp; 

		for( int i=0; i<mNumNodes; i++ )
		{
            for( int j=0; j<mNumNodes; j++ )
            {
                double int_NiNj = 0.00;

                for( int gp=0; gp<nGaussPoints; gp++ )
                {
                	int_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
                }
            
            	N_nij[ i               ][ j               ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes   ][ j + mNumNodes   ] = cteNiNj * int_NiNj; 
            	N_nij[ i + mNumNodes*2 ][ j + mNumNodes*2 ] = cteNiNj * int_NiNj; 
            }
		}
		
		// Projection in q basis matrices
		Matrix< std::complex<double> > iA_div(  4,  4 );
		Matrix< std::complex<double> > iA_cur( 12, 12 );
		
		for( int i=0; i<4; i++ )
		{
            for( int j=0; j<4; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -4.00 / mVolume;
        		else          iA_div[ i ][ j ] = +16.0 / mVolume;                          
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}
		
		// Scaling the projection in q matrices
		iA_div /= cMu ;
		iA_cur /= cMu ;
		
		// Final assembly: [ CT * Ac^(-1) * C ] + [ DT * Ad^(-1) * D ] - [ Ni * Nj ]
		Matrix< std::complex<double> > iAc_x_Ccur = iA_cur * C_cur;
		Matrix< std::complex<double> > Ccur_tran  = C_cur.Transpose();
		StiffMatrix += ( Ccur_tran * iAc_x_Ccur );

		Matrix< std::complex<double> > iAd_x_Ddiv = iA_div * D_div;
		Matrix< std::complex<double> > Ddiv_tran  = D_div.Transpose();
		StiffMatrix += ( Ddiv_tran * iAd_x_Ddiv );

		StiffMatrix -= N_nij;
    }

    //********************************************************************************************************************************
    // - Calculates rotational P1 projection on vertice nodes
    //********************************************************************************************************************************
	void VolumeElement_3sb_FullWave::Calculate_Rot_P1_Projection( cVectorType& cRotP1x, cVectorType& cRotP1y, cVectorType& cRotP1z )
	{
		// Electric field on nodes
		Matrix< std::complex<double> > ncE( mNumDofs, 1 );
		     
		for( int n=0; n<mNumNodes; n++ )
		{
			ncE[ n               ][ 0 ] = (*mProperties)( cEx, *mNodes[n] );
			ncE[ n + mNumNodes   ][ 0 ] = (*mProperties)( cEy, *mNodes[n] );
			ncE[ n + mNumNodes*2 ][ 0 ] = (*mProperties)( cEz, *mNodes[n] );
		}

		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order5( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( int gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_3sb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_3sb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl contribution)
		Matrix< std::complex<double> > C_cur( 12, mNumDofs );

		for( int i=0; i<4; i++ )
		{
		    for( int j=0; j<mNumNodes; j++ )
		    {
		    	double vintg_Ni_dNjdX = 0.00;
		        double vintg_Ni_dNjdY = 0.00;
		        double vintg_Ni_dNjdZ = 0.00;
		    
                for( int gp=0; gp<nGaussPoints; gp++ )
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

		// Adding surface integral nxu*q
		Add_C_NxN_Matrix( C_cur ); 

		// Projection q basis matrices
		Matrix< std::complex<double> > iA_div(  4,  4 );
		Matrix< std::complex<double> > iA_cur( 12, 12 );
		
		for( int i=0; i<4; i++ )
		{
            for( int j=0; j<4; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -4.00 / mVolume;
        		else          iA_div[ i ][ j ] = +16.0 / mVolume;                          
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}

		// Rotational P1 projection on vertices
		Matrix< std::complex<double> > cRotP1Proj( 12, 1 ); 
		
		cRotP1Proj = ( iA_cur * C_cur ) * ncE;

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
    // - Calculates rotational on nodes
    //********************************************************************************************************************************
	void VolumeElement_3sb_FullWave::Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes )
	{
		// Rotational on vertices (P1 projection)
		cVectorType cRotP1x, cRotP1y, cRotP1z;

		Calculate_Rot_P1_Projection( cRotP1x, cRotP1y, cRotP1z );
		
		// P1 basis functions on nodal points
		std::vector<double> nX( mNumNodes ), nY( mNumNodes ), nZ( mNumNodes );

		// Vertices 0,1,2,3
		nX[0] = 0.00; nY[0] = 0.00; nZ[0] = 0.00;
		nX[1] = 1.00; nY[1] = 0.00; nZ[1] = 0.00;
		nX[2] = 0.00; nY[2] = 1.00; nZ[2] = 0.00;
		nX[3] = 0.00; nY[3] = 0.00; nZ[3] = 1.00;
		
		// Face [0 1 2]
		nX[4] = 0.25; nY[4] = 0.25; nZ[4] = 0.00;
		nX[5] = 0.50; nY[5] = 0.25; nZ[5] = 0.00;
		nX[6] = 0.25; nY[6] = 0.50; nZ[6] = 0.00;

		// Face [0 1 3]
		nX[7] = 0.25; nY[7] = 0.00; nZ[7] = 0.25;
		nX[8] = 0.50; nY[8] = 0.00; nZ[8] = 0.25;
		nX[9] = 0.25; nY[9] = 0.00; nZ[9] = 0.50;

		// Face [1 2 3]
		nX[10] = 0.50; nY[10] = 0.25; nZ[10] = 0.25;
		nX[11] = 0.25; nY[11] = 0.50; nZ[11] = 0.25;
		nX[12] = 0.25; nY[12] = 0.25; nZ[12] = 0.50;

		// Face [0 2 3]
		nX[13] = 0.00; nY[13] = 0.25; nZ[13] = 0.25;
		nX[14] = 0.00; nY[14] = 0.50; nZ[14] = 0.25;
		nX[15] = 0.00; nY[15] = 0.25; nZ[15] = 0.50;

		// Volume [0 1 2 3]
		nX[16] = 0.25; nY[16] = 0.25; nZ[16] = 0.25; 

		// P1 basis on nodes
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, nX, nY, nZ );

		// Rotational in nodes
		rot_OnNodes.resize( mNumNodes );

		Vector< std::complex<double> > ncRot( 3 );

		std::complex<double> cZero( 0.00, 0.00 );

		for( int n=0; n<mNumNodes; n++ )
		{
		    ncRot[0] = cZero; 
		    ncRot[1] = cZero; 
		    ncRot[2] = cZero;
		 
		    for( int i=0; i<4; i++ )
		    {
		        ncRot[0] += cRotP1x[i] * N[i][n];
			    ncRot[1] += cRotP1y[i] * N[i][n];
			    ncRot[2] += cRotP1z[i] * N[i][n];
		    }

		    rot_OnNodes[n] = ncRot;
		}		
	}
	  
    //********************************************************************************************************************************
    // - Calculates rotational on Gauss points
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs )
	{
		// Rotational on vertices (P1 projection)
		cVectorType cRotP1x, cRotP1y, cRotP1z;

		Calculate_Rot_P1_Projection( cRotP1x, cRotP1y, cRotP1z );

		// Gauss points
		std::vector<double> cX, cY, cZ;
	
		int numGaussPoints = Get_InnerGiD_GaussPoints( cX, cY, cZ, numResultsOnGPs );

		// P1 basis on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_1st( N, cX, cY, cZ );

		// Rotational on Gauss points
		rot_OnGPs.resize( numGaussPoints );

		Vector< std::complex<double> > gpcRot( 3 );

		std::complex<double> cZero( 0.00, 0.00 );
        
		for( int gp=0; gp<numGaussPoints; gp++ )
		{
		    gpcRot[0] = cZero; 
		    gpcRot[1] = cZero; 
		    gpcRot[2] = cZero;
		 
		    for( int i=0; i<4; i++ )
		    {
		        gpcRot[0] += cRotP1x[i] * N[i][gp];
			    gpcRot[1] += cRotP1y[i] * N[i][gp];
			    gpcRot[2] += cRotP1z[i] * N[i][gp];
		    }

		    rot_OnGPs[gp] = gpcRot;
		}
	}

    //********************************************************************************************************************************
    // - Calculates Electric field on Gauss points
    //********************************************************************************************************************************
    void VolumeElement_3sb_FullWave::Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs )
	{
		// Electric field on nodes
        Matrix< std::complex<double> > cE_n( mNumDofs, 1 );
		     
		for( int n=0; n<mNumNodes; n++ )
		{
			cE_n[ n               ][0] = (*mProperties)(cEx, *mNodes[n]);
			cE_n[ n + mNumNodes   ][0] = (*mProperties)(cEy, *mNodes[n]);
			cE_n[ n + mNumNodes*2 ][0] = (*mProperties)(cEz, *mNodes[n]);
		}

		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order5( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_3sb( N, cX, cY, cZ );

		// Building projection matrices 
		Matrix< std::complex<double> > matrix_NiNj( 12, mNumDofs );

		for( int i=0; i<4; i++ )
		{
		    for( int j=0; j<mNumNodes; j++ )
		    {
		    	double intg_NiNj = 0.00;
		    
                for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
                }      
		    	
				matrix_NiNj[ i     ][ j               ] = intg_NiNj; 
                matrix_NiNj[ i + 4 ][ j + mNumNodes   ] = intg_NiNj;  						 
                matrix_NiNj[ i + 8 ][ j + mNumNodes*2 ] = intg_NiNj; 
		    }
		}

		// Projection basis matrices
		Matrix< std::complex<double> > iA_div(  4,  4 );
		Matrix< std::complex<double> > iA_cur( 12, 12 );
		
		for( int i=0; i<4; i++ )
		{
            for( int j=0; j<4; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -4.00 / mVolume;
        		else          iA_div[ i ][ j ] = +16.0 / mVolume;                          
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}

		// L2 projection on nodes 
		Matrix< std::complex<double> > EProj( 12, 1 ); 
		
		EProj += ( iA_cur * matrix_NiNj ) * cE_n;

		Vector< std::complex<double> > pnEx( 4 );
		Vector< std::complex<double> > pnEy( 4 );
		Vector< std::complex<double> > pnEz( 4 );

		for( int i=0; i<4; i++ )
		{
			pnEx[ i ] = EProj[ i     ][ 0 ];
			pnEy[ i ] = EProj[ i + 4 ][ 0 ];
			pnEz[ i ] = EProj[ i + 8 ][ 0 ];
		}

		// GiD internal Gauss points
		std::vector<double> cXig, cYig, cZig;
	
		int numGaussPoints_ig = Get_InnerGiD_GaussPoints( cXig, cYig, cZig, numResultsOnGPs );

        // Basis functions on Gauss points
		Matrix<double> Nig;

		Lagrange3D_Ni_1st( Nig, cXig, cYig, cZig );

	    // Electric field on Gauss points
		Ef_OnGPs.resize( numGaussPoints_ig );

		Vector< std::complex<double> > cE_ig( 3 );

		std::complex<double> cZero( 0.00, 0.00 );

	    for( int ig=0; ig<numGaussPoints_ig; ig++ )
		{
			cE_ig[ 0 ] = cZero; 
			cE_ig[ 1 ] = cZero; 
			cE_ig[ 2 ] = cZero;

			for( int n=0; n<4; n++ )
		    {       
				cE_ig[ 0 ] += pnEx[ n ] * Nig[ n ][ ig ];
			    cE_ig[ 1 ] += pnEy[ n ] * Nig[ n ][ ig ];
			    cE_ig[ 2 ] += pnEz[ n ] * Nig[ n ][ ig ];
			}

			Ef_OnGPs[ ig ] = cE_ig;
		}	
	}

    //********************************************************************************************************************************
    // - Get inner GiD Gauss points 
    //********************************************************************************************************************************
    int VolumeElement_3sb_FullWave::Get_InnerGiD_GaussPoints( std::vector<double>& cX, 
		                                                      std::vector<double>& cY, 
		                                                      std::vector<double>& cZ, 
		                                                      int numResultsOnGPs    )
	{
		if ( numResultsOnGPs <= 1 ) 
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