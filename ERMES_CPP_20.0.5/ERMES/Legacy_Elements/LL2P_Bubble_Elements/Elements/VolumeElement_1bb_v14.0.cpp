
#include "VolumeElement_1bb.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

#include <unordered_set>

//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//
//#define NR_END 1
//
//#define FREE_ARG char*
//
//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
//
//static double dmaxarg1,dmaxarg2;
//#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
//
//static int iminarg1,iminarg2;
//#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : - fabs(a))

static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1 = (a),maxarg2 = (b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1 = (a),iminarg2 = (b),(iminarg1 < (iminarg2) ? (iminarg1) : iminarg2))

static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)


namespace Kratos
{
    //**********************************************************************************************************************
    // - Global index of the nodes
    //**********************************************************************************************************************	
    void VolumeElement_1bb::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( ctcNumDofs );

        for( int i = 0; i < ctcNumNodes; ++i ) 
        {  
            EquationId[ i                 ] = mNodes[ i ]->pDofcEx()->EquationId();
            EquationId[ i + ctcNumNodes   ] = mNodes[ i ]->pDofcEy()->EquationId();
			EquationId[ i + ctcNumNodes*2 ] = mNodes[ i ]->pDofcEz()->EquationId();
        }

		//EquationId.resize( 12 );

  //      for( int i = 0; i < 4; ++i ) 
  //      {  
  //          EquationId[ i       ] = mNodes[ i ]->pDofcEx()->EquationId();
  //          EquationId[ i + 4   ] = mNodes[ i ]->pDofcEy()->EquationId();
		//	EquationId[ i + 4*2 ] = mNodes[ i ]->pDofcEz()->EquationId();
  //      }
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

			// Node 0 
			dNdx[ 0 ][ p ] = dL1dx;
			dNdy[ 0 ][ p ] = dL1dy;
			dNdz[ 0 ][ p ] = dL1dz;

			// Node 1  
			dNdx[ 1 ][ p ] = dL2dx;
			dNdy[ 1 ][ p ] = dL2dy;
			dNdz[ 1 ][ p ] = dL2dz;

			// Node 2  
			dNdx[ 2 ][ p ] = dL3dx;
			dNdy[ 2 ][ p ] = dL3dy;
			dNdz[ 2 ][ p ] = dL3dz;

			// Node 3  
			dNdx[ 3 ][ p ] = dL4dx;
			dNdy[ 3 ][ p ] = dL4dy;
			dNdz[ 3 ][ p ] = dL4dz;

			// Face bubble [0 1 2]
			dNdx[ 4 ][ p ] = 27.0 * ( ( L1 * L2 * dL3dx ) + ( L1 * dL2dx * L3 ) + ( dL1dx * L2 * L3 ) );
			dNdy[ 4 ][ p ] = 27.0 * ( ( L1 * L2 * dL3dy ) + ( L1 * dL2dy * L3 ) + ( dL1dy * L2 * L3 ) );
			dNdz[ 4 ][ p ] = 27.0 * ( ( L1 * L2 * dL3dz ) + ( L1 * dL2dz * L3 ) + ( dL1dz * L2 * L3 ) );

			// Face bubble [0 1 3]
			dNdx[ 5 ][ p ] = 27.0 * ( ( L1 * L2 * dL4dx ) + ( L1 * dL2dx * L4 ) + ( dL1dx * L2 * L4 ) );
			dNdy[ 5 ][ p ] = 27.0 * ( ( L1 * L2 * dL4dy ) + ( L1 * dL2dy * L4 ) + ( dL1dy * L2 * L4 ) );
			dNdz[ 5 ][ p ] = 27.0 * ( ( L1 * L2 * dL4dz ) + ( L1 * dL2dz * L4 ) + ( dL1dz * L2 * L4 ) );

			// Face bubble [1 2 3]
			dNdx[ 6 ][ p ] = 27.0 * ( ( L2 * L3 * dL4dx ) + ( L2 * dL3dx * L4 ) + ( dL2dx * L3 * L4 ) );
			dNdy[ 6 ][ p ] = 27.0 * ( ( L2 * L3 * dL4dy ) + ( L2 * dL3dy * L4 ) + ( dL2dy * L3 * L4 ) );
			dNdz[ 6 ][ p ] = 27.0 * ( ( L2 * L3 * dL4dz ) + ( L2 * dL3dz * L4 ) + ( dL2dz * L3 * L4 ) );

			// Face bubble [0 2 3]
			dNdx[ 7 ][ p ] = 27.0 * ( ( L1 * L3 * dL4dx ) + ( L1 * dL3dx * L4 ) + ( dL1dx * L3 * L4 ) );
			dNdy[ 7 ][ p ] = 27.0 * ( ( L1 * L3 * dL4dy ) + ( L1 * dL3dy * L4 ) + ( dL1dy * L3 * L4 ) );
			dNdz[ 7 ][ p ] = 27.0 * ( ( L1 * L3 * dL4dz ) + ( L1 * dL3dz * L4 ) + ( dL1dz * L3 * L4 ) );

			// Inner bubble
			dNdx[ 8 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dx ) + ( L1 * L2 * dL3dx * L4 ) + ( L1 * dL2dx * L3 * L4 ) + ( dL1dx * L2 * L3 * L4 ) );
			dNdy[ 8 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dy ) + ( L1 * L2 * dL3dy * L4 ) + ( L1 * dL2dy * L3 * L4 ) + ( dL1dy * L2 * L3 * L4 ) );
			dNdz[ 8 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dz ) + ( L1 * L2 * dL3dz * L4 ) + ( L1 * dL2dz * L3 * L4 ) + ( dL1dz * L2 * L3 * L4 ) );
        }
	}

    //**********************************************************************************************************************
    // - Calculation of dN/dx, dN/dy, dN/dz on cX, cY, cZ points
    //**********************************************************************************************************************
    void VolumeElement_1bb::LagrangeDerivatives3D_2bb( Matrix<double>& dNdx,
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

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
			dNdx[10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL2dx ) ) + ( L2 * L3 * dL4dx * ( 4.00*L2 - 1.00 ) ) + ( L2 * dL3dx * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL2dx * L3 * L4 * ( 4.00*L2 - 1.00 ) ) );
			dNdy[10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL2dy ) ) + ( L2 * L3 * dL4dy * ( 4.00*L2 - 1.00 ) ) + ( L2 * dL3dy * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL2dy * L3 * L4 * ( 4.00*L2 - 1.00 ) ) );
			dNdz[10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL2dz ) ) + ( L2 * L3 * dL4dz * ( 4.00*L2 - 1.00 ) ) + ( L2 * dL3dz * L4 * ( 4.00*L2 - 1.00 ) ) + ( dL2dz * L3 * L4 * ( 4.00*L2 - 1.00 ) ) );
							 			
			dNdx[11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL3dx ) ) + ( L2 * L3 * dL4dx * ( 4.00*L3 - 1.00 ) ) + ( L2 * dL3dx * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL2dx * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
			dNdy[11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL3dy ) ) + ( L2 * L3 * dL4dy * ( 4.00*L3 - 1.00 ) ) + ( L2 * dL3dy * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL2dy * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
			dNdz[11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL3dz ) ) + ( L2 * L3 * dL4dz * ( 4.00*L3 - 1.00 ) ) + ( L2 * dL3dz * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL2dz * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
							 			  
			dNdx[12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL4dx ) ) + ( L2 * L3 * dL4dx * ( 4.00*L4 - 1.00 ) ) + ( L2 * dL3dx * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL2dx * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdy[12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL4dy ) ) + ( L2 * L3 * dL4dy * ( 4.00*L4 - 1.00 ) ) + ( L2 * dL3dy * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL2dy * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdz[12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 4.00 * dL4dz ) ) + ( L2 * L3 * dL4dz * ( 4.00*L4 - 1.00 ) ) + ( L2 * dL3dz * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL2dz * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );

			// Face bubble [ 0 2 3 ]
			dNdx[13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL1dx ) ) + ( L1 * L3 * dL4dx * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL3dx * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dx * L3 * L4 * ( 4.00*L1 - 1.00 ) ) );
			dNdy[13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL1dy ) ) + ( L1 * L3 * dL4dy * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL3dy * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dy * L3 * L4 * ( 4.00*L1 - 1.00 ) ) );
			dNdz[13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL1dz ) ) + ( L1 * L3 * dL4dz * ( 4.00*L1 - 1.00 ) ) + ( L1 * dL3dz * L4 * ( 4.00*L1 - 1.00 ) ) + ( dL1dz * L3 * L4 * ( 4.00*L1 - 1.00 ) ) );
							 			  
			dNdx[14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL3dx ) ) + ( L1 * L3 * dL4dx * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL3dx * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL1dx * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
			dNdy[14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL3dy ) ) + ( L1 * L3 * dL4dy * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL3dy * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL1dy * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
			dNdz[14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL3dz ) ) + ( L1 * L3 * dL4dz * ( 4.00*L3 - 1.00 ) ) + ( L1 * dL3dz * L4 * ( 4.00*L3 - 1.00 ) ) + ( dL1dz * L3 * L4 * ( 4.00*L3 - 1.00 ) ) );
							 			 
			dNdx[15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL4dx ) ) + ( L1 * L3 * dL4dx * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL3dx * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dx * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdy[15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL4dy ) ) + ( L1 * L3 * dL4dy * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL3dy * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dy * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );
			dNdz[15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 4.00 * dL4dz ) ) + ( L1 * L3 * dL4dz * ( 4.00*L4 - 1.00 ) ) + ( L1 * dL3dz * L4 * ( 4.00*L4 - 1.00 ) ) + ( dL1dz * L3 * L4 * ( 4.00*L4 - 1.00 ) ) );

///////////////////////////////////

			//// Face bubble [ 0 1 2 ]
			//dNdx[ 4 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL1dx ) ) + ( L1 * L2 * dL3dx * ( 2.00*L1 ) ) + ( L1 * dL2dx * L3 * ( 2.00*L1 ) ) + ( dL1dx * L2 * L3 * ( 2.00*L1 ) ) );
			//dNdy[ 4 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL1dy ) ) + ( L1 * L2 * dL3dy * ( 2.00*L1 ) ) + ( L1 * dL2dy * L3 * ( 2.00*L1 ) ) + ( dL1dy * L2 * L3 * ( 2.00*L1 ) ) );
			//dNdz[ 4 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL1dz ) ) + ( L1 * L2 * dL3dz * ( 2.00*L1 ) ) + ( L1 * dL2dz * L3 * ( 2.00*L1 ) ) + ( dL1dz * L2 * L3 * ( 2.00*L1 ) ) );
			//				 			  				  
			//dNdx[ 5 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL2dx ) ) + ( L1 * L2 * dL3dx * ( 2.00*L2 ) ) + ( L1 * dL2dx * L3 * ( 2.00*L2 ) ) + ( dL1dx * L2 * L3 * ( 2.00*L2 ) ) );
			//dNdy[ 5 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL2dy ) ) + ( L1 * L2 * dL3dy * ( 2.00*L2 ) ) + ( L1 * dL2dy * L3 * ( 2.00*L2 ) ) + ( dL1dy * L2 * L3 * ( 2.00*L2 ) ) );
			//dNdz[ 5 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL2dz ) ) + ( L1 * L2 * dL3dz * ( 2.00*L2 ) ) + ( L1 * dL2dz * L3 * ( 2.00*L2 ) ) + ( dL1dz * L2 * L3 * ( 2.00*L2 ) ) );
			//				 			  				 
			//dNdx[ 6 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL3dx ) ) + ( L1 * L2 * dL3dx * ( 2.00*L3 ) ) + ( L1 * dL2dx * L3 * ( 2.00*L3 ) ) + ( dL1dx * L2 * L3 * ( 2.00*L3 ) ) );
			//dNdy[ 6 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL3dy ) ) + ( L1 * L2 * dL3dy * ( 2.00*L3 ) ) + ( L1 * dL2dy * L3 * ( 2.00*L3 ) ) + ( dL1dy * L2 * L3 * ( 2.00*L3 ) ) );
			//dNdz[ 6 ][ p ] = 32.00 * ( ( L1 * L2 * L3 * ( 2.00 * dL3dz ) ) + ( L1 * L2 * dL3dz * ( 2.00*L3 ) ) + ( L1 * dL2dz * L3 * ( 2.00*L3 ) ) + ( dL1dz * L2 * L3 * ( 2.00*L3 ) ) );
			//											  
			//// Face bubble [ 0 1 3 ]					  
			//dNdx[ 7 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL1dx ) ) + ( L1 * L2 * dL4dx * ( 2.00*L1 ) ) + ( L1 * dL2dx * L4 * ( 2.00*L1 ) ) + ( dL1dx * L2 * L4 * ( 2.00*L1 ) ) );
			//dNdy[ 7 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL1dy ) ) + ( L1 * L2 * dL4dy * ( 2.00*L1 ) ) + ( L1 * dL2dy * L4 * ( 2.00*L1 ) ) + ( dL1dy * L2 * L4 * ( 2.00*L1 ) ) );
			//dNdz[ 7 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL1dz ) ) + ( L1 * L2 * dL4dz * ( 2.00*L1 ) ) + ( L1 * dL2dz * L4 * ( 2.00*L1 ) ) + ( dL1dz * L2 * L4 * ( 2.00*L1 ) ) );
			//				 			  				
			//dNdx[ 8 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL2dx ) ) + ( L1 * L2 * dL4dx * ( 2.00*L2 ) ) + ( L1 * dL2dx * L4 * ( 2.00*L2 ) ) + ( dL1dx * L2 * L4 * ( 2.00*L2 ) ) );
			//dNdy[ 8 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL2dy ) ) + ( L1 * L2 * dL4dy * ( 2.00*L2 ) ) + ( L1 * dL2dy * L4 * ( 2.00*L2 ) ) + ( dL1dy * L2 * L4 * ( 2.00*L2 ) ) );
			//dNdz[ 8 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL2dz ) ) + ( L1 * L2 * dL4dz * ( 2.00*L2 ) ) + ( L1 * dL2dz * L4 * ( 2.00*L2 ) ) + ( dL1dz * L2 * L4 * ( 2.00*L2 ) ) );
			//				 			  				 
			//dNdx[ 9 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL4dx ) ) + ( L1 * L2 * dL4dx * ( 2.00*L4 ) ) + ( L1 * dL2dx * L4 * ( 2.00*L4 ) ) + ( dL1dx * L2 * L4 * ( 2.00*L4 ) ) );
			//dNdy[ 9 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL4dy ) ) + ( L1 * L2 * dL4dy * ( 2.00*L4 ) ) + ( L1 * dL2dy * L4 * ( 2.00*L4 ) ) + ( dL1dy * L2 * L4 * ( 2.00*L4 ) ) );
			//dNdz[ 9 ][ p ] = 32.00 * ( ( L1 * L2 * L4 * ( 2.00 * dL4dz ) ) + ( L1 * L2 * dL4dz * ( 2.00*L4 ) ) + ( L1 * dL2dz * L4 * ( 2.00*L4 ) ) + ( dL1dz * L2 * L4 * ( 2.00*L4 ) ) );
			//											  
			//// Face bubble [ 1 2 3 ]					  
			//dNdx[10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL2dx ) ) + ( L2 * L3 * dL4dx * ( 2.00*L2 ) ) + ( L2 * dL3dx * L4 * ( 2.00*L2 ) ) + ( dL2dx * L3 * L4 * ( 2.00*L2 ) ) );
			//dNdy[10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL2dy ) ) + ( L2 * L3 * dL4dy * ( 2.00*L2 ) ) + ( L2 * dL3dy * L4 * ( 2.00*L2 ) ) + ( dL2dy * L3 * L4 * ( 2.00*L2 ) ) );
			//dNdz[10 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL2dz ) ) + ( L2 * L3 * dL4dz * ( 2.00*L2 ) ) + ( L2 * dL3dz * L4 * ( 2.00*L2 ) ) + ( dL2dz * L3 * L4 * ( 2.00*L2 ) ) );
			//				 							 
			//dNdx[11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL3dx ) ) + ( L2 * L3 * dL4dx * ( 2.00*L3 ) ) + ( L2 * dL3dx * L4 * ( 2.00*L3 ) ) + ( dL2dx * L3 * L4 * ( 2.00*L3 ) ) );
			//dNdy[11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL3dy ) ) + ( L2 * L3 * dL4dy * ( 2.00*L3 ) ) + ( L2 * dL3dy * L4 * ( 2.00*L3 ) ) + ( dL2dy * L3 * L4 * ( 2.00*L3 ) ) );
			//dNdz[11 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL3dz ) ) + ( L2 * L3 * dL4dz * ( 2.00*L3 ) ) + ( L2 * dL3dz * L4 * ( 2.00*L3 ) ) + ( dL2dz * L3 * L4 * ( 2.00*L3 ) ) );
			//				 			  				 
			//dNdx[12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL4dx ) ) + ( L2 * L3 * dL4dx * ( 2.00*L4 ) ) + ( L2 * dL3dx * L4 * ( 2.00*L4 ) ) + ( dL2dx * L3 * L4 * ( 2.00*L4 ) ) );
			//dNdy[12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL4dy ) ) + ( L2 * L3 * dL4dy * ( 2.00*L4 ) ) + ( L2 * dL3dy * L4 * ( 2.00*L4 ) ) + ( dL2dy * L3 * L4 * ( 2.00*L4 ) ) );
			//dNdz[12 ][ p ] = 32.00 * ( ( L2 * L3 * L4 * ( 2.00 * dL4dz ) ) + ( L2 * L3 * dL4dz * ( 2.00*L4 ) ) + ( L2 * dL3dz * L4 * ( 2.00*L4 ) ) + ( dL2dz * L3 * L4 * ( 2.00*L4 ) ) );
			//											 
			//// Face bubble [ 0 2 3 ]					 
			//dNdx[13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL1dx ) ) + ( L1 * L3 * dL4dx * ( 2.00*L1 ) ) + ( L1 * dL3dx * L4 * ( 2.00*L1 ) ) + ( dL1dx * L3 * L4 * ( 2.00*L1 ) ) );
			//dNdy[13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL1dy ) ) + ( L1 * L3 * dL4dy * ( 2.00*L1 ) ) + ( L1 * dL3dy * L4 * ( 2.00*L1 ) ) + ( dL1dy * L3 * L4 * ( 2.00*L1 ) ) );
			//dNdz[13 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL1dz ) ) + ( L1 * L3 * dL4dz * ( 2.00*L1 ) ) + ( L1 * dL3dz * L4 * ( 2.00*L1 ) ) + ( dL1dz * L3 * L4 * ( 2.00*L1 ) ) );
			//				 			  				 
			//dNdx[14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL3dx ) ) + ( L1 * L3 * dL4dx * ( 2.00*L3 ) ) + ( L1 * dL3dx * L4 * ( 2.00*L3 ) ) + ( dL1dx * L3 * L4 * ( 2.00*L3 ) ) );
			//dNdy[14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL3dy ) ) + ( L1 * L3 * dL4dy * ( 2.00*L3 ) ) + ( L1 * dL3dy * L4 * ( 2.00*L3 ) ) + ( dL1dy * L3 * L4 * ( 2.00*L3 ) ) );
			//dNdz[14 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL3dz ) ) + ( L1 * L3 * dL4dz * ( 2.00*L3 ) ) + ( L1 * dL3dz * L4 * ( 2.00*L3 ) ) + ( dL1dz * L3 * L4 * ( 2.00*L3 ) ) );
			//				 			 				  
			//dNdx[15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL4dx ) ) + ( L1 * L3 * dL4dx * ( 2.00*L4 ) ) + ( L1 * dL3dx * L4 * ( 2.00*L4 ) ) + ( dL1dx * L3 * L4 * ( 2.00*L4 ) ) );
			//dNdy[15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL4dy ) ) + ( L1 * L3 * dL4dy * ( 2.00*L4 ) ) + ( L1 * dL3dy * L4 * ( 2.00*L4 ) ) + ( dL1dy * L3 * L4 * ( 2.00*L4 ) ) );
			//dNdz[15 ][ p ] = 32.00 * ( ( L1 * L3 * L4 * ( 2.00 * dL4dz ) ) + ( L1 * L3 * dL4dz * ( 2.00*L4 ) ) + ( L1 * dL3dz * L4 * ( 2.00*L4 ) ) + ( dL1dz * L3 * L4 * ( 2.00*L4 ) ) );

////////////////////////////////////////////

			// Inner bubble
			dNdx[16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dx ) + ( L1 * L2 * dL3dx * L4 ) + ( L1 * dL2dx * L3 * L4 ) + ( dL1dx * L2 * L3 * L4 ) );
			dNdy[16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dy ) + ( L1 * L2 * dL3dy * L4 ) + ( L1 * dL2dy * L3 * L4 ) + ( dL1dy * L2 * L3 * L4 ) );
			dNdz[16 ][ p ] = 256.0 * ( ( L1 * L2 * L3 * dL4dz ) + ( L1 * L2 * dL3dz * L4 ) + ( L1 * dL2dz * L3 * L4 ) + ( dL1dz * L2 * L3 * L4 ) );
        }
	}

	//**************************************************************************************************************
    // - Calculation of dN/dx, dN/dy, dN/dz on cX,cY,cZ points
    //**************************************************************************************************************
    void VolumeElement_1bb::LagrangeDerivatives3D_1st( Matrix<double>& dNdx,
		                                               Matrix<double>& dNdy,
													   Matrix<double>& dNdz,
													   std::vector<double>& cX, 
													   std::vector<double>& cY, 
													   std::vector<double>& cZ )
    {
		int nNodes  = 4;
		int nPoints = cX.size();

		dNdx.Resize(nNodes, nPoints);
		dNdy.Resize(nNodes, nPoints);
		dNdz.Resize(nNodes, nPoints);
		
		Matrix<double> DN; 
		
		Calculate_Ni_Derivatives( DN );

		double dL1dx,dL2dx,dL3dx,dL4dx;
        double dL1dy,dL2dy,dL3dy,dL4dy;
		double dL1dz,dL2dz,dL3dz,dL4dz;

        dL1dx = DN[0][0], dL2dx = DN[0][1], dL3dx = DN[0][2], dL4dx = DN[0][3];
        dL1dy = DN[1][0], dL2dy = DN[1][1], dL3dy = DN[1][2], dL4dy = DN[1][3];
		dL1dz = DN[2][0], dL2dz = DN[2][1], dL3dz = DN[2][2], dL4dz = DN[2][3];

		//dNi/dx ,dNi/dy, dNi/dz
        for (int p=0; p<nPoints; p++ ) 
        {
			int nn;

			/////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////   POINTS   /////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////////////////////////////
			
			//////////////////////////////////////// Points: 0  /////////////////////////////////////////////////
			nn = 0;
			
			dNdx[nn][p] = dL1dx;
			dNdy[nn][p] = dL1dy;
			dNdz[nn][p] = dL1dz;

			//////////////////////////////////////// Points: 1  /////////////////////////////////////////////////
			nn = 1;
			
			dNdx[nn][p] = dL2dx;
			dNdy[nn][p] = dL2dy;
			dNdz[nn][p] = dL2dz;

			//////////////////////////////////////// Points: 2  /////////////////////////////////////////////////
			nn = 2;
			
			dNdx[nn][p] = dL3dx;
			dNdy[nn][p] = dL3dy;
			dNdz[nn][p] = dL3dz;

			//////////////////////////////////////// Points: 3  /////////////////////////////////////////////////
			nn = 3;
			
			dNdx[nn][p] = dL4dx;
			dNdy[nn][p] = dL4dy;
			dNdz[nn][p] = dL4dz;
        }
	}

	//***********************************************************************************
    //* - Calculation of dN/dx, dN/dy, dN/dz on cX,cY,cZ points
    //***********************************************************************************
    void VolumeElement_1bb::LagrangeDerivatives3D_3th(Matrix<double>& dNdx,
		                                              Matrix<double>& dNdy,
													  Matrix<double>& dNdz,
													  std::vector<double>& cX, 
													  std::vector<double>& cY, 
													  std::vector<double>& cZ)
    {
		int nNodes  = 20;
		int nPoints = cX.size();

		dNdx.Resize(nNodes, nPoints);
		dNdy.Resize(nNodes, nPoints);
		dNdz.Resize(nNodes, nPoints);
		
		Matrix<double> DN;
        
		Calculate_Ni_Derivatives( DN );

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
			dNdL = 13.50 * (L*L) - 9.00 * L + 1.00;
			
			dNdx[nn][p] = dNdL * dL1dx;
			dNdy[nn][p] = dNdL * dL1dy;
			dNdz[nn][p] = dNdL * dL1dz;

			//////////////////////////////////////// Points: 1  /////////////////////////////////////////////////
			nn   = 1;
			L    = L2;
			dNdL = 13.50 * (L*L) - 9.00 * L + 1.00;
			
			dNdx[nn][p] = dNdL * dL2dx;
			dNdy[nn][p] = dNdL * dL2dy;
			dNdz[nn][p] = dNdL * dL2dz;

			//////////////////////////////////////// Points: 2  /////////////////////////////////////////////////
			nn   = 2;
			L    = L3;
			dNdL = 13.50 * (L*L) - 9.00 * L + 1.00;
			
			dNdx[nn][p] = dNdL * dL3dx;
			dNdy[nn][p] = dNdL * dL3dy;
			dNdz[nn][p] = dNdL * dL3dz;

			//////////////////////////////////////// Points: 3  /////////////////////////////////////////////////
			nn   = 3;
			L    = L4;
			dNdL = 13.50 * (L*L) - 9.00 * L + 1.00;
			
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
			dNdLi = 27.00 * Li * Lj - 4.50 * Lj;
			dNdLj = 13.50 * (Li*Li) - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 5
			nn    = 5;
			dNdLi = 13.50 * (Lj*Lj) - 4.50 * Lj;
		    dNdLj = 27.00 * Lj * Li - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [1 2] /////////////////////////////////////////////////
			Li = L2; dLidx = dL2dx; dLidy = dL2dy; dLidz = dL2dz;
			Lj = L3; dLjdx = dL3dx; dLjdy = dL3dy; dLjdz = dL3dz;

            //nodo: 6
            nn    = 6;
			dNdLi = 27.00 * Li * Lj - 4.50 * Lj;
			dNdLj = 13.50 * (Li*Li) - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 7
			nn    = 7;
			dNdLi = 13.50 * (Lj*Lj) - 4.50 * Lj;
		    dNdLj = 27.00 * Lj * Li - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [2 0] /////////////////////////////////////////////////
			Li = L3; dLidx = dL3dx; dLidy = dL3dy; dLidz = dL3dz;
			Lj = L1; dLjdx = dL1dx; dLjdy = dL1dy; dLjdz = dL1dz;

			//nodo: 8
            nn    = 8;
			dNdLi = 27.00 * Li * Lj - 4.50 * Lj;
			dNdLj = 13.50 * (Li*Li) - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 9
			nn    = 9;
			dNdLi = 13.50 * (Lj*Lj) - 4.50 * Lj;
		    dNdLj = 27.00 * Lj * Li - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [0 3] /////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 10
            nn    = 10;
			dNdLi = 27.00 * Li * Lj - 4.50 * Lj;
			dNdLj = 13.50 * (Li*Li) - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 11
			nn    = 11;
			dNdLi = 13.50 * (Lj*Lj) - 4.50 * Lj;
		    dNdLj = 27.00 * Lj * Li - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [1 3] /////////////////////////////////////////////////
			Li = L2; dLidx = dL2dx; dLidy = dL2dy; dLidz = dL2dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 12
            nn    = 12;
			dNdLi = 27.00 * Li * Lj - 4.50 * Lj;
			dNdLj = 13.50 * (Li*Li) - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 13
			nn    = 13;
			dNdLi = 13.50 * (Lj*Lj) - 4.50 * Lj;
		    dNdLj = 27.00 * Lj * Li - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [2 3] /////////////////////////////////////////////////
			Li = L3; dLidx = dL3dx; dLidy = dL3dy; dLidz = dL3dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 14
            nn    = 14;
			dNdLi = 27.00 * Li * Lj - 4.50 * Lj;
			dNdLj = 13.50 * (Li*Li) - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 15
			nn    = 15;
			dNdLi = 13.50 * (Lj*Lj) - 4.50 * Lj;
		    dNdLj = 27.00 * Lj * Li - 4.50 * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////   FACES   /////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
			double Lk;
            double dNdLk;
			double dLkdx, dLkdy, dLkdz;
			
			///////////////////////////////////// Face: [0 1 2] ////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L2; dLjdx = dL2dx; dLjdy = dL2dy; dLjdz = dL2dz;
			Lk = L3; dLkdx = dL3dx; dLkdy = dL3dy; dLkdz = dL3dz;

			//nodo: 16
            nn    = 16;
			dNdLi = 27.00 * Lj * Lk;
			dNdLj = 27.00 * Li * Lk ;
			dNdLk = 27.00 * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			///////////////////////////////////// Face: [0 1 3] ////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L2; dLjdx = dL2dx; dLjdy = dL2dy; dLjdz = dL2dz;
			Lk = L4; dLkdx = dL4dx; dLkdy = dL4dy; dLkdz = dL4dz;

			//nodo: 17
            nn    = 17;
			dNdLi = 27.00 * Lj * Lk;
			dNdLj = 27.00 * Li * Lk ;
			dNdLk = 27.00 * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			///////////////////////////////////// Face: [1 2 3] ////////////////////////////////////////////////
			Li = L2; dLidx = dL2dx; dLidy = dL2dy; dLidz = dL2dz;
			Lj = L3; dLjdx = dL3dx; dLjdy = dL3dy; dLjdz = dL3dz;
			Lk = L4; dLkdx = dL4dx; dLkdy = dL4dy; dLkdz = dL4dz;

			//nodo: 18
            nn    = 18;
			dNdLi = 27.00 * Lj * Lk;
			dNdLj = 27.00 * Li * Lk ;
			dNdLk = 27.00 * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			///////////////////////////////////// Face: [0 2 3] ////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L3; dLjdx = dL3dx; dLjdy = dL3dy; dLjdz = dL3dz;
			Lk = L4; dLkdx = dL4dx; dLkdy = dL4dy; dLkdz = dL4dz;

			//nodo: 19
            nn    = 19;
			dNdLi = 27.00 * Lj * Lk;
			dNdLj = 27.00 * Li * Lk ;
			dNdLk = 27.00 * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;		
        }
	}

	//***********************************************************************************
    //* - Calculation of dN/dx, dN/dy, dN/dz on cX,cY,cZ points
    //***********************************************************************************
    void VolumeElement_1bb::LagrangeDerivatives3D_2nd(Matrix<double>& dNdx,
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
		
		Matrix<double> DN;
        
		Calculate_Ni_Derivatives( DN );

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

	//***********************************************************************************
    //* - Calculation of dN/dx, dN/dy, dN/dz on cX,cY,cZ points
    //***********************************************************************************
    void VolumeElement_1bb::LagrangeDerivatives3D_4th(Matrix<double>& dNdx,
		                                              Matrix<double>& dNdy,
													  Matrix<double>& dNdz,
													  std::vector<double>& cX, 
													  std::vector<double>& cY, 
													  std::vector<double>& cZ)
    {
		int nNodes  = 35;
		int nPoints = cX.size();

		dNdx.Resize(nNodes, nPoints);
		dNdy.Resize(nNodes, nPoints);
		dNdz.Resize(nNodes, nPoints);
		
		Matrix<double> DN;
        
		Calculate_Ni_Derivatives( DN );

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
			dNdL = (128.00/3.00) * (L*L*L) - (48.00) * (L*L) + (44.00/3.00) * L - 1.00;
			
			dNdx[nn][p] = dNdL * dL1dx;
			dNdy[nn][p] = dNdL * dL1dy;
			dNdz[nn][p] = dNdL * dL1dz;

			//////////////////////////////////////// Points: 1  /////////////////////////////////////////////////
			nn   = 1;
			L    = L2;
			dNdL = (128.00/3.00) * (L*L*L) - (48.00) * (L*L) + (44.00/3.00) * L - 1.00;
			
			dNdx[nn][p] = dNdL * dL2dx;
			dNdy[nn][p] = dNdL * dL2dy;
			dNdz[nn][p] = dNdL * dL2dz;

			//////////////////////////////////////// Points: 2  /////////////////////////////////////////////////
			nn   = 2;
			L    = L3;
			dNdL = (128.00/3.00) * (L*L*L) - (48.00) * (L*L) + (44.00/3.00) * L - 1.00;
			
			dNdx[nn][p] = dNdL * dL3dx;
			dNdy[nn][p] = dNdL * dL3dy;
			dNdz[nn][p] = dNdL * dL3dz;

			//////////////////////////////////////// Points: 3  /////////////////////////////////////////////////
			nn   = 3;
			L    = L4;
			dNdL = (128.00/3.00) * (L*L*L) - (48.00) * (L*L) + (44.00/3.00) * L - 1.00;
			
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
			dNdLi = (128.00)      * (Li*Li) * Lj - (64.00) * Li * Lj + (16.00/3.00) * Lj;
			dNdLj = (128.00/3.00) * (Li*Li*Li)   - (32.00) * (Li*Li) + (16.00/3.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 5
			nn    = 5;
			dNdLi = (128.00) * Li * (Lj*Lj) - (32.00) * Li * Lj - (16.00) * (Lj*Lj) + (4.00) * Lj;
		    dNdLj = (128.00) * Lj * (Li*Li) - (32.00) * Li * Lj - (16.00) * (Li*Li) + (4.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

            //nodo: 6
			nn    = 6;
			dNdLi = (128.00/3.00) * (Lj*Lj*Lj)   - (32.00) * (Lj*Lj) + (16.00/3.00) * Lj;
		    dNdLj = (128.00)      * (Lj*Lj) * Li - (64.00) * Lj * Li + (16.00/3.00) * Li; 

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [1 2] /////////////////////////////////////////////////
			Li = L2; dLidx = dL2dx; dLidy = dL2dy; dLidz = dL2dz;
			Lj = L3; dLjdx = dL3dx; dLjdy = dL3dy; dLjdz = dL3dz;

            //nodo: 7
            nn    = 7;
			dNdLi = (128.00)      * (Li*Li) * Lj - (64.00) * Li * Lj + (16.00/3.00) * Lj;
			dNdLj = (128.00/3.00) * (Li*Li*Li)   - (32.00) * (Li*Li) + (16.00/3.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 8
			nn    = 8;
			dNdLi = (128.00) * Li * (Lj*Lj) - (32.00) * Li * Lj - (16.00) * (Lj*Lj) + (4.00) * Lj;
		    dNdLj = (128.00) * Lj * (Li*Li) - (32.00) * Li * Lj - (16.00) * (Li*Li) + (4.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

            //nodo: 9
			nn    = 9;
			dNdLi = (128.00/3.00) * (Lj*Lj*Lj)   - (32.00) * (Lj*Lj) + (16.00/3.00) * Lj;
		    dNdLj = (128.00)      * (Lj*Lj) * Li - (64.00) * Lj * Li + (16.00/3.00) * Li; 

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [2 0] /////////////////////////////////////////////////
			Li = L3; dLidx = dL3dx; dLidy = dL3dy; dLidz = dL3dz;
			Lj = L1; dLjdx = dL1dx; dLjdy = dL1dy; dLjdz = dL1dz;

			//nodo: 10
            nn    = 10;
			dNdLi = (128.00)      * (Li*Li) * Lj - (64.00) * Li * Lj + (16.00/3.00) * Lj;
			dNdLj = (128.00/3.00) * (Li*Li*Li)   - (32.00) * (Li*Li) + (16.00/3.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 11
			nn    = 11;
			dNdLi = (128.00) * Li * (Lj*Lj) - (32.00) * Li * Lj - (16.00) * (Lj*Lj) + (4.00) * Lj;
		    dNdLj = (128.00) * Lj * (Li*Li) - (32.00) * Li * Lj - (16.00) * (Li*Li) + (4.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

            //nodo: 12
			nn    = 12;
			dNdLi = (128.00/3.00) * (Lj*Lj*Lj)   - (32.00) * (Lj*Lj) + (16.00/3.00) * Lj;
		    dNdLj = (128.00)      * (Lj*Lj) * Li - (64.00) * Lj * Li + (16.00/3.00) * Li; 

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [0 3] /////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 13
            nn    = 13;
			dNdLi = (128.00)      * (Li*Li) * Lj - (64.00) * Li * Lj + (16.00/3.00) * Lj;
			dNdLj = (128.00/3.00) * (Li*Li*Li)   - (32.00) * (Li*Li) + (16.00/3.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 14
			nn    = 14;
			dNdLi = (128.00) * Li * (Lj*Lj) - (32.00) * Li * Lj - (16.00) * (Lj*Lj) + (4.00) * Lj;
		    dNdLj = (128.00) * Lj * (Li*Li) - (32.00) * Li * Lj - (16.00) * (Li*Li) + (4.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

            //nodo: 15
			nn    = 15;
			dNdLi = (128.00/3.00) * (Lj*Lj*Lj)   - (32.00) * (Lj*Lj) + (16.00/3.00) * Lj;
		    dNdLj = (128.00)      * (Lj*Lj) * Li - (64.00) * Lj * Li + (16.00/3.00) * Li; 

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [1 3] /////////////////////////////////////////////////
			Li = L2; dLidx = dL2dx; dLidy = dL2dy; dLidz = dL2dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 16
            nn    = 16;
			dNdLi = (128.00)      * (Li*Li) * Lj - (64.00) * Li * Lj + (16.00/3.00) * Lj;
			dNdLj = (128.00/3.00) * (Li*Li*Li)   - (32.00) * (Li*Li) + (16.00/3.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 17
			nn    = 17;
			dNdLi = (128.00) * Li * (Lj*Lj) - (32.00) * Li * Lj - (16.00) * (Lj*Lj) + (4.00) * Lj;
		    dNdLj = (128.00) * Lj * (Li*Li) - (32.00) * Li * Lj - (16.00) * (Li*Li) + (4.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

            //nodo: 18
			nn    = 18;
			dNdLi = (128.00/3.00) * (Lj*Lj*Lj)   - (32.00) * (Lj*Lj) + (16.00/3.00) * Lj;
		    dNdLj = (128.00)      * (Lj*Lj) * Li - (64.00) * Lj * Li + (16.00/3.00) * Li; 

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////// Edge: [2 3] /////////////////////////////////////////////////
			Li = L3; dLidx = dL3dx; dLidy = dL3dy; dLidz = dL3dz;
			Lj = L4; dLjdx = dL4dx; dLjdy = dL4dy; dLjdz = dL4dz;

            //nodo: 19
            nn    = 19;
			dNdLi = (128.00)      * (Li*Li) * Lj - (64.00) * Li * Lj + (16.00/3.00) * Lj;
			dNdLj = (128.00/3.00) * (Li*Li*Li)   - (32.00) * (Li*Li) + (16.00/3.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			//nodo: 20
			nn    = 20;
			dNdLi = (128.00) * Li * (Lj*Lj) - (32.00) * Li * Lj - (16.00) * (Lj*Lj) + (4.00) * Lj;
		    dNdLj = (128.00) * Lj * (Li*Li) - (32.00) * Li * Lj - (16.00) * (Li*Li) + (4.00) * Li;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

            //nodo: 21
			nn    = 21;
			dNdLi = (128.00/3.00) * (Lj*Lj*Lj)   - (32.00) * (Lj*Lj) + (16.00/3.00) * Lj;
		    dNdLj = (128.00)      * (Lj*Lj) * Li - (64.00) * Lj * Li + (16.00/3.00) * Li; 

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz;

			////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////   FACES   /////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
			double Lk;
            double dNdLk;
			double dLkdx, dLkdy, dLkdz;
			
			///////////////////////////////////// Face: [0 1 2] ////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L2; dLjdx = dL2dx; dLjdy = dL2dy; dLjdz = dL2dz;
			Lk = L3; dLkdx = dL3dx; dLkdy = dL3dy; dLkdz = dL3dz;

			//nodo: 22
            nn    = 22;
			dNdLi = (256.00) * Li * Lj * Lk - (32.00) * Lj * Lk;
			dNdLj = (128.00) * (Li*Li) * Lk - (32.00) * Li * Lk;
			dNdLk = (128.00) * (Li*Li) * Lj - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			//nodo: 23
            nn    = 23;
			dNdLi = (128.00) * (Lj*Lj) * Lk - (32.00) * Lj * Lk;
			dNdLj = (256.00) * Li * Lj * Lk - (32.00) * Li * Lk;
			dNdLk = (128.00) * (Lj*Lj) * Li - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			//nodo: 24
            nn    = 24;
			dNdLi = (128.00) * (Lk*Lk) * Lj - (32.00) * Lj * Lk;
			dNdLj = (128.00) * (Lk*Lk) * Li - (32.00) * Li * Lk;
			dNdLk = (256.00) * Li * Lj * Lk - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			///////////////////////////////////// Face: [0 1 3] ////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L2; dLjdx = dL2dx; dLjdy = dL2dy; dLjdz = dL2dz;
			Lk = L4; dLkdx = dL4dx; dLkdy = dL4dy; dLkdz = dL4dz;

			//nodo: 25
            nn    = 25;
			dNdLi = (256.00) * Li * Lj * Lk - (32.00) * Lj * Lk;
			dNdLj = (128.00) * (Li*Li) * Lk - (32.00) * Li * Lk;
			dNdLk = (128.00) * (Li*Li) * Lj - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			//nodo: 26
            nn    = 26;
			dNdLi = (128.00) * (Lj*Lj) * Lk - (32.00) * Lj * Lk;
			dNdLj = (256.00) * Li * Lj * Lk - (32.00) * Li * Lk;
			dNdLk = (128.00) * (Lj*Lj) * Li - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			//nodo: 27
            nn    = 27;
			dNdLi = (128.00) * (Lk*Lk) * Lj - (32.00) * Lj * Lk;
			dNdLj = (128.00) * (Lk*Lk) * Li - (32.00) * Li * Lk;
			dNdLk = (256.00) * Li * Lj * Lk - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			///////////////////////////////////// Face: [1 2 3] ////////////////////////////////////////////////
			Li = L2; dLidx = dL2dx; dLidy = dL2dy; dLidz = dL2dz;
			Lj = L3; dLjdx = dL3dx; dLjdy = dL3dy; dLjdz = dL3dz;
			Lk = L4; dLkdx = dL4dx; dLkdy = dL4dy; dLkdz = dL4dz;

			//nodo: 28
            nn    = 28;
			dNdLi = (256.00) * Li * Lj * Lk - (32.00) * Lj * Lk;
			dNdLj = (128.00) * (Li*Li) * Lk - (32.00) * Li * Lk;
			dNdLk = (128.00) * (Li*Li) * Lj - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			//nodo: 29
            nn    = 29;
			dNdLi = (128.00) * (Lj*Lj) * Lk - (32.00) * Lj * Lk;
			dNdLj = (256.00) * Li * Lj * Lk - (32.00) * Li * Lk;
			dNdLk = (128.00) * (Lj*Lj) * Li - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			//nodo: 30
            nn    = 30;
			dNdLi = (128.00) * (Lk*Lk) * Lj - (32.00) * Lj * Lk;
			dNdLj = (128.00) * (Lk*Lk) * Li - (32.00) * Li * Lk;
			dNdLk = (256.00) * Li * Lj * Lk - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			///////////////////////////////////// Face: [0 2 3] ////////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L3; dLjdx = dL3dx; dLjdy = dL3dy; dLjdz = dL3dz;
			Lk = L4; dLkdx = dL4dx; dLkdy = dL4dy; dLkdz = dL4dz;

			//nodo: 31
            nn    = 31;
			dNdLi = (256.00) * Li * Lj * Lk - (32.00) * Lj * Lk;
			dNdLj = (128.00) * (Li*Li) * Lk - (32.00) * Li * Lk;
			dNdLk = (128.00) * (Li*Li) * Lj - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			//nodo: 32
            nn    = 32;
			dNdLi = (128.00) * (Lj*Lj) * Lk - (32.00) * Lj * Lk;
			dNdLj = (256.00) * Li * Lj * Lk - (32.00) * Li * Lk;
			dNdLk = (128.00) * (Lj*Lj) * Li - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			//nodo: 33
            nn    = 33;
			dNdLi = (128.00) * (Lk*Lk) * Lj - (32.00) * Lj * Lk;
			dNdLj = (128.00) * (Lk*Lk) * Li - (32.00) * Li * Lk;
			dNdLk = (256.00) * Li * Lj * Lk - (32.00) * Li * Lj;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz;

			////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////   VOLUME   ////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////////////////
			double Lm;
            double dNdLm;
			double dLmdx, dLmdy, dLmdz;

			////////////////////////////////// Volume: [0 1 2 3 4] /////////////////////////////////////////////
			Li = L1; dLidx = dL1dx; dLidy = dL1dy; dLidz = dL1dz;
			Lj = L2; dLjdx = dL2dx; dLjdy = dL2dy; dLjdz = dL2dz;
			Lk = L3; dLkdx = dL3dx; dLkdy = dL3dy; dLkdz = dL3dz;
			Lm = L4; dLmdx = dL4dx; dLmdy = dL4dy; dLmdz = dL4dz;

			//nodo: 34
            nn    = 34;
			dNdLi = (256.00) * Lj * Lk * Lm;
			dNdLj = (256.00) * Li * Lk * Lm;
			dNdLk = (256.00) * Li * Lj * Lm;
			dNdLm = (256.00) * Li * Lj * Lk;

			dNdx[nn][p] = dNdLi * dLidx + dNdLj * dLjdx + dNdLk * dLkdx + dNdLm * dLmdx;
			dNdy[nn][p] = dNdLi * dLidy + dNdLj * dLjdy + dNdLk * dLkdy + dNdLm * dLmdy;
			dNdz[nn][p] = dNdLi * dLidz + dNdLj * dLjdz + dNdLk * dLkdz + dNdLm * dLmdz;		
        }
	}

    //**********************************************************************************************************************
    // - Calculate derivatives in nodes
    //**********************************************************************************************************************
    void VolumeElement_1bb::DerivativesInNodes1bb( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz )
    {
		std::vector<double> cX( mNumNodes ); 
		std::vector<double>	cY( mNumNodes );
		std::vector<double> cZ( mNumNodes );

		// Corner nodes 
		cX[0] = 0.00; cY[0] = 0.00; cZ[0] = 0.00;
		cX[1] = 1.00; cY[1] = 0.00; cZ[1] = 0.00;
		cX[2] = 0.00; cY[2] = 1.00; cZ[2] = 0.00;
		cX[3] = 0.00; cY[3] = 0.00; cZ[3] = 1.00;

		// Face bubbles
		double  inv3 = 1.0/3.0;
		cX[4] = inv3; cY[4] = inv3; cZ[4] = 0.00;
		cX[5] = inv3; cY[5] = 0.00; cZ[5] = inv3;
		cX[6] = inv3; cY[6] = inv3; cZ[6] = inv3;
		cX[7] = 0.00; cY[7] = inv3; cZ[7] = inv3;

		// Inner bubble
        cX[8] = 0.25; cY[8] = 0.25; cZ[8] = 0.25;

		LagrangeDerivatives3D_1bb( ndNdx, ndNdy, ndNdz, cX, cY, cZ );
	}

	//**********************************************************************************************************************
    // - Calculate derivatives in nodes
    //**********************************************************************************************************************
	void VolumeElement_1bb::DerivativesInNodes2bb( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz )
    {
		std::vector<double> cX( mNumNodes ); 
		std::vector<double>	cY( mNumNodes );
		std::vector<double> cZ( mNumNodes );

		// Points 0,1,2,3
		cX[0] = 0.00; cY[0] = 0.00; cZ[0] = 0.00;
		cX[1] = 1.00; cY[1] = 0.00; cZ[1] = 0.00;
		cX[2] = 0.00; cY[2] = 1.00; cZ[2] = 0.00;
		cX[3] = 0.00; cY[3] = 0.00; cZ[3] = 1.00;
		
		// Face [0 1 2]
		cX[4 ] = 0.25; cY[4 ] = 0.25; cZ[4 ] = 0.00;
		cX[5 ] = 0.50; cY[5 ] = 0.25; cZ[5 ] = 0.00;
		cX[6 ] = 0.25; cY[6 ] = 0.50; cZ[6 ] = 0.00;

		// Face [0 1 3]
		cX[7 ] = 0.25; cY[7 ] = 0.00; cZ[7 ] = 0.25;
		cX[8 ] = 0.50; cY[8 ] = 0.00; cZ[8 ] = 0.25;
		cX[9 ] = 0.25; cY[9 ] = 0.00; cZ[9 ] = 0.50;

		// Face [1 2 3]
		cX[10] = 0.50; cY[10] = 0.25; cZ[10] = 0.25;
		cX[11] = 0.25; cY[11] = 0.50; cZ[11] = 0.25;
		cX[12] = 0.25; cY[12] = 0.25; cZ[12] = 0.50;

		// Face [0 2 3]
		cX[13] = 0.00; cY[13] = 0.25; cZ[13] = 0.25;
		cX[14] = 0.00; cY[14] = 0.50; cZ[14] = 0.25;
		cX[15] = 0.00; cY[15] = 0.25; cZ[15] = 0.50;

		// Volume [0 1 2 3]
		cX[16] = 0.25; cY[16] = 0.25; cZ[16] = 0.25; 

		LagrangeDerivatives3D_2bb( ndNdx, ndNdy, ndNdz, cX, cY, cZ );
	}

    //**********************************************************************************************************************
    // - Calculates H = ( 1.0 / i*w*mu0 ) * rot(E)
    //**********************************************************************************************************************
    void VolumeElement_1bb::CalculateDerivatives( Vector<std::complex<double> >& ncHx,
		                                          Vector<std::complex<double> >& ncHy,
												  Vector<std::complex<double> >& ncHz )
    {
		/*
        //int i,n;
		//
		//std::complex<double> cUnit( 0.00, 1.00 );
        //std::complex<double> cZero( 0.00, 0.00 );
        //
        //double mu_real = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
		//double mu_imag = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;
		//
		//std::complex<double> cMu( mu_real, mu_imag );
		//std::complex<double> inv_jwmu = 1.00 / ( mFreq * cUnit * cMu );
		//
		//Matrix<double> ndNdx; 
		//Matrix<double> ndNdy; 
		//Matrix<double> ndNdz;
		//
		//DerivativesInNodes1bb( ndNdx, ndNdy, ndNdz );
		//
		//Vector<std::complex<double> > ncEx(mNumNodes);
		//Vector<std::complex<double> > ncEy(mNumNodes);
		//Vector<std::complex<double> > ncEz(mNumNodes);
		//     
		//for(n=0; n<mNumNodes; n++)
		//{
		//	ncEx[n] = (*mProperties)(cEx, *mNodes[n]);
		//	ncEy[n] = (*mProperties)(cEy, *mNodes[n]);
		//	ncEz[n] = (*mProperties)(cEz, *mNodes[n]);
		//}
		//
        //ncHx.resize( mNumNodes );
		//ncHy.resize( mNumNodes );
		//ncHz.resize( mNumNodes );
        //      
		//for( n=0; n<mNumNodes; n++ )
		//{
		//	ncHx[n] = cZero;
		//	ncHy[n] = cZero;
		//	ncHz[n] = cZero;
		//
		//	for( i=0; i<mNumNodes; i++ )
		//	{
		//	    ncHx[n] += ndNdy[i][n]*ncEz[i] - ndNdz[i][n]*ncEy[i];
		//		ncHy[n] += ndNdz[i][n]*ncEx[i] - ndNdx[i][n]*ncEz[i];
		//		ncHz[n] += ndNdx[i][n]*ncEy[i] - ndNdy[i][n]*ncEx[i];
		//	}
		//}

		//ncHx *= inv_jwmu;
		//ncHy *= inv_jwmu;
		//ncHz *= inv_jwmu;  

		//double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        //double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        //double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        
        //std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
		//std::complex<double> inv_jwep = 1.00 / ( mFreq * cUnit * cEp );

		//ncHx *= inv_jwep;
		//ncHy *= inv_jwep;
		//ncHz *= inv_jwep;  
		*/

        /////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////

	
        //Loop counters
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
		//std::complex<double> cEp( sigma, mFreq * eps_real );

        std::complex<double> cMu( mu_real , mu_imag );
		std::complex<double> inv_jwmu = 1.00 / ( mFreq * cUnit * cMu );
		std::complex<double> inv_jwep = 1.00 / ( mFreq * cUnit * cEp );

      
		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order8( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

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
                	/*vintg_Ni_dNjdX += W[gp] * N[i][gp] * dNdx[j][gp]; 
		    		vintg_Ni_dNjdY += W[gp] * N[i][gp] * dNdy[j][gp];
		    		vintg_Ni_dNjdZ += W[gp] * N[i][gp] * dNdz[j][gp];*/

					vintg_Ni_dNjdX += W[gp] * N[i][gp] * N[j][gp]; 
		    		
                }      
		    	
		    	//C_cur[ i     ][ j               ] =            0.00; 
       //         C_cur[ i     ][ j + mNumNodes   ] = -vintg_Ni_dNjdZ; 
       //         C_cur[ i     ][ j + mNumNodes*2 ] = +vintg_Ni_dNjdY; 
		    	//							    
		    	//C_cur[ i + 4 ][ j               ] = +vintg_Ni_dNjdZ; 
       //         C_cur[ i + 4 ][ j + mNumNodes   ] =            0.00; 
       //         C_cur[ i + 4 ][ j + mNumNodes*2 ] = -vintg_Ni_dNjdX; 
		    	//							    
		    	//C_cur[ i + 8 ][ j               ] = -vintg_Ni_dNjdY; 
       //         C_cur[ i + 8 ][ j + mNumNodes   ] = +vintg_Ni_dNjdX; 
       //         C_cur[ i + 8 ][ j + mNumNodes*2 ] =            0.00; 

				C_cur[ i     ][ j               ] = vintg_Ni_dNjdX; 
                C_cur[ i + 4 ][ j + mNumNodes   ] = vintg_Ni_dNjdX;  						 
                C_cur[ i + 8 ][ j + mNumNodes*2 ] = vintg_Ni_dNjdX; 
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

		//CurlProj *= inv_jwmu;
		//CurlProj *=  sigma * inv_jwep;

		ncHx.resize( mNumNodes, mZero );
		ncHy.resize( mNumNodes, mZero );
		ncHz.resize( mNumNodes, mZero );

		for( i=0; i<4; i++ )
		{
			ncHx[i] = CurlProj[i  ][0];
			ncHy[i] = CurlProj[i+4][0];
			ncHz[i] = CurlProj[i+8][0];
		}

        /////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////
		//ncHx.resize(35, mZero );
		//ncHy.resize(35, mZero );
		//ncHz.resize(35, mZero );
  //      return;

  //      //Loop counters
		//int i, j, n, gp;

		//int p4NumNodes = 35;
		//int p4NumDofs  = p4NumNodes*3;

		//// Gauss points and weights
		//std::vector<double> cX, cY, cZ, W;

		//int nGaussPoints = GaussPoints3D_Order10( cX, cY, cZ, W );

		//double detJ = 6.00 * mVolume;

		//for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

  //      // Basis functions on Gauss points
		//Matrix<double> N;

		//Lagrange3D_Ni_4th( N, cX, cY, cZ );

		////std::cout<<N<<std::endl;

		//// Basis functions on Gauss points
		//Matrix<double> M;

		//Lagrange3D_Ni_2bb( M, cX, cY, cZ );

		////std::cout<<M<<std::endl;

		//// Building projection matrices (curl and div contributions)
		//Matrix< std::complex<double> > C_cur; C_cur.Resize( p4NumDofs, mNumDofs, mZero );

		//for( i=0; i<p4NumNodes; i++ )
		//{
		//    for( j=0; j<mNumNodes; j++ )
		//    {
		//    	double vintg_Ni_Mj = 0.00;

  //              for( gp=0; gp<nGaussPoints; gp++ )
  //              {
		//			vintg_Ni_Mj += W[gp] * N[i][gp] * M[j][gp];
  //              }      
		//    	
		//		C_cur[ i                ][ j               ] += vintg_Ni_Mj; 
  //              C_cur[ i + p4NumNodes   ][ j + mNumNodes   ] += vintg_Ni_Mj;  						 
  //              C_cur[ i + p4NumNodes*2 ][ j + mNumNodes*2 ] += vintg_Ni_Mj; 
		//    }
		//}

		////std::cout<<C_cur<<std::endl;

		///////////////////////////////////////////////////////////////////////////////////////

		//Matrix< double > p4A; p4A.Resize( p4NumNodes, p4NumNodes, 0.00 );

		//for( i=0; i<p4NumNodes; i++ )
		//{
  //          for( j=0; j<p4NumNodes; j++ )
  //          {
		//		double vintg_Ni_Nj = 0.00;

  //              for( gp=0; gp<nGaussPoints; gp++ )
  //              {
		//			vintg_Ni_Nj += W[gp] * N[i][gp] * N[j][gp];
  //              }      
  //           
		//		p4A[ i ][ j ] = vintg_Ni_Nj/mVolume; 
  //          }									   
		//}

		////std::cout<<p4A<<std::endl;

		//Matrix< double > ip4A = p4A;

		//Invert_Matrix( ip4A );

  //      //std::cout<<"2"<<std::endl;

  //      Matrix< std::complex<double> > invp4A; invp4A.Resize( p4NumDofs, p4NumDofs, mZero );
		//
		//for( i=0; i<p4NumNodes; i++ )
		//{
  //          for( j=0; j<p4NumNodes; j++ )
  //          {
		//		invp4A[ i                ][ j                ] = mVolume*ip4A[ i ][ j ]; 
  //              invp4A[ i + p4NumNodes   ][ j + p4NumNodes   ] = mVolume*ip4A[ i ][ j ];
  //              invp4A[ i + p4NumNodes*2 ][ j + p4NumNodes*2 ] = mVolume*ip4A[ i ][ j ];
  //          }									   
		//}

  //      //std::cout<<invp4A<<std::endl;

		/////////////////////////////////////////////////////////////////////

		////std::cout<<( invp4A * C_cur )<<std::endl;

		//Matrix< std::complex<double> > L2Proj = ( invp4A * C_cur );

		//Matrix< std::complex<double> > AllE; AllE.Resize( mNumDofs, 1 );
		//     
		//for( n=0; n<4; n++ )
		//{
		//	AllE[ n               ][0] = (*mProperties)(cEx, *mNodes[n]);
		//	AllE[ n + mNumNodes   ][0] = (*mProperties)(cEy, *mNodes[n]);
		//	AllE[ n + mNumNodes*2 ][0] = (*mProperties)(cEz, *mNodes[n]);
		//}

		//for( n=22; n<35; n++ )
		//{
		//	int m = n-18;
		//	AllE[ m               ][0] = (*mProperties)(cEx, *mNodes[n]);
		//	AllE[ m + mNumNodes   ][0] = (*mProperties)(cEy, *mNodes[n]);
		//	AllE[ m + mNumNodes*2 ][0] = (*mProperties)(cEz, *mNodes[n]);
		//}

  //      //std::cout << AllE << std::endl;

		//Matrix< std::complex<double> > CurlProj = ( L2Proj * AllE );

		////std::cout << CurlProj << std::endl;

		//ncHx.resize( p4NumNodes, mZero );
		//ncHy.resize( p4NumNodes, mZero );
		//ncHz.resize( p4NumNodes, mZero );

		//for( i=0; i<p4NumNodes; i++ )
		//{
		//	ncHx[i] = CurlProj[i             ][0];
		//	ncHy[i] = CurlProj[i+p4NumNodes  ][0];
		//	ncHz[i] = CurlProj[i+p4NumNodes*2][0];
		//}
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

        // Nodal basis x Nodal basis
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

		// Face bubble x Nodal basis
        std::complex<double> cCteFBNNNN = mVolume * cCte * ( 27.0 / 840.0 );
		std::complex<double> cCteFBN2NN = mVolume * cCte * ( 27.0 / 420.0 );
		
		std::complex<double> cCteFBNi;

		// Positions where we have N1*N2*N3*N4
		std::vector<int> N2index{ 6, 7, 5, 4 };

		for( i=0; i<4; ++i ) 
		{
			for( j=4; j<8; ++j ) 
			{
				if ( j == N2index[i] ) cCteFBNi = cCteFBNNNN;
				else                   cCteFBNi = cCteFBN2NN;
				
				w2_Ni_ep_Nj[ i ][ j ] = cCteFBNi;
			    w2_Ni_ep_Nj[ j ][ i ] = cCteFBNi;
			}
	    }

		// Face bubble x Face bubble
		std::complex<double> cCteFBiFBi = mVolume * cCte * ( 729.0 /  7560.0 );
		std::complex<double> cCteFBiFBj = mVolume * cCte * ( 729.0 / 15120.0 );

		for( i=4; i<8; ++i ) 
		{
			for( j=4; j<8; ++j ) 
			{
        		if ( i != j ) w2_Ni_ep_Nj[ i ][ j ] = cCteFBiFBj;
        		else          w2_Ni_ep_Nj[ i ][ j ] = cCteFBiFBi;		      
			}
	    }

        // Inner bubble x Nodal basis
        std::complex<double> cCteNiBB = mVolume * cCte * ( 256.0 / 3360.0 );

		for( i=0; i<4; ++i ) 
		{
			w2_Ni_ep_Nj[ i ][ 8 ] = cCteNiBB;
			w2_Ni_ep_Nj[ 8 ][ i ] = cCteNiBB;
	    }

		// Inner bubble x Face bubble
		std::complex<double> cCteFBBB = mVolume * cCte * ( 6912.0 / 75600.0 );

        for( i=4; i<8; ++i ) 
		{
		    w2_Ni_ep_Nj[ i ][ 8 ] = cCteFBBB;
		    w2_Ni_ep_Nj[ 8 ][ i ] = cCteFBBB;
		}
			
		// Inner bubble x Inner bubble
		std::complex<double> cCteBBBB = mVolume * cCte * ( 65536.0 / 415800.0 );

        w2_Ni_ep_Nj[ 8 ][ 8 ] = cCteBBBB;
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
			Vector<unsigned int> cnElements = (*mpNodesConnectivity)[ mNodes[ i ]->Id() - 1 ];

			// Iteration over elements connected to corner node 
			Vector<unsigned int>::iterator ce_it;

			for( ce_it = cnElements.begin(); ce_it != cnElements.end(); ++ce_it )
			{
				Vector<unsigned int> eInfo = (*mpElementsInfo)[ *ce_it ];

				double volInt_Ni = 0.25 * (*gpVolumes)[ *ce_it ];

				// Loop over the nodes of (*ce_it) element
				for( int k=0; k<4; ++k ) 
				{
					unsigned int kId = eInfo[ k ];

				    double volInt_Ni_dNkdX = volInt_Ni * (*gp_dNk_dX)[ *ce_it ][ k ];
					double volInt_Ni_dNkdY = volInt_Ni * (*gp_dNk_dY)[ *ce_it ][ k ];
					double volInt_Ni_dNkdZ = volInt_Ni * (*gp_dNk_dZ)[ *ce_it ][ k ];

					if( Ci_dNkdX[ i ].find( kId ) != Ci_dNkdX[ i ].end() )
					{
                        Ci_dNkdX[ i ][ kId ] += volInt_Ni_dNkdX;
						Ci_dNkdY[ i ][ kId ] += volInt_Ni_dNkdY;
						Ci_dNkdZ[ i ][ kId ] += volInt_Ni_dNkdZ;
				    }
					else
					{
						Ci_dNkdX[ i ][ kId ]  = volInt_Ni_dNkdX;
						Ci_dNkdY[ i ][ kId ]  = volInt_Ni_dNkdY;
						Ci_dNkdZ[ i ][ kId ]  = volInt_Ni_dNkdZ;			
					}

			    }// End nodes in connected element loop - for( int k=0; k<4; ++k ) 

			}//End elements connected to corner node loop - for( ce_it = cnElements.begin(); ce_it != cnElements.end(); ++ce_it )

		}// End corner nodes loop - for( int i=0; i<4; ++i ) 

		// Computing ( Ci, dNk/d{X,Y,Z} ) on face bubble nodes
		for( int i=4; i<8; ++i ) 
		{
			// Elements connected to this face bubble node
			Vector<unsigned int> fbElements; Get_Node_Connectivity( i, fbElements );

			// Iteration over elements connected to face bubble node 
			Vector<unsigned int>::iterator fb_it;

			for( fb_it = fbElements.begin(); fb_it != fbElements.end(); ++fb_it )
			{
				Vector<unsigned int> eInfo = (*mpElementsInfo)[ *fb_it ];

				double volInt_FBi = ( 27.0 / 120.0 ) * (*gpVolumes)[ *fb_it ];

				// Loop over the nodes of (*fb_it) element
				for( int k=0; k<4; ++k ) 
				{
					unsigned int kId = eInfo[ k ];

				    double volInt_FBi_dNkdX = volInt_FBi * (*gp_dNk_dX)[ *fb_it ][ k ];
					double volInt_FBi_dNkdY = volInt_FBi * (*gp_dNk_dY)[ *fb_it ][ k ];
					double volInt_FBi_dNkdZ = volInt_FBi * (*gp_dNk_dZ)[ *fb_it ][ k ];

					if( Ci_dNkdX[ i ].find( kId ) != Ci_dNkdX[ i ].end() )
					{
                        Ci_dNkdX[ i ][ kId ] += volInt_FBi_dNkdX;
						Ci_dNkdY[ i ][ kId ] += volInt_FBi_dNkdY;
						Ci_dNkdZ[ i ][ kId ] += volInt_FBi_dNkdZ;
				    }
					else
					{
						Ci_dNkdX[ i ][ kId ]  = volInt_FBi_dNkdX;
						Ci_dNkdY[ i ][ kId ]  = volInt_FBi_dNkdY;
						Ci_dNkdZ[ i ][ kId ]  = volInt_FBi_dNkdZ;			
					}

			    }// End nodes in connected element loop - for( int k=0; k<4; ++k ) 

			}// End elements connected to face bubble node loop - for( fb_it = fbElements.begin(); fb_it != fbElements.end(); ++fb_it )

		}// End face bubble nodes loop - for( int i=4; i<8; ++i ) 

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
	    	unsigned int kId = mNodes[ k ]->Id() - 1;
	    
	        double volInt_Bi_dNkdX = volInt_Bi * DN[ 0 ][ k ];
	    	double volInt_Bi_dNkdY = volInt_Bi * DN[ 1 ][ k ];
	    	double volInt_Bi_dNkdZ = volInt_Bi * DN[ 2 ][ k ];

			Ci_dNkdX[ 8 ][ kId ] = volInt_Bi_dNkdX;
	    	Ci_dNkdY[ 8 ][ kId ] = volInt_Bi_dNkdY;
	    	Ci_dNkdZ[ 8 ][ kId ] = volInt_Bi_dNkdZ;
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


		// Computing ( Ci, dNk/d{X,Y,Z} ) on corner nodes
		for( int i=0; i<4; ++i ) 
		{
			// Elements connected to corner node
			Vector<unsigned int> cnElements = (*mpNodesConnectivity)[ mNodes[ i ]->Id() - 1 ];

			// Iteration over elements connected to corner node 
			Vector<unsigned int>::iterator ce_it;

			for( ce_it = cnElements.begin(); ce_it != cnElements.end(); ++ce_it )
			{
				Vector<unsigned int> eInfo = (*mpElementsInfo)[ *ce_it ];

				Properties::Pointer epProperties = (*gpPropertiesMap)[ eInfo[4] ];

			    double sigma   = (*epProperties)(IHL_ELECTRIC_CONDUCTIVITY);
			    double ep_real = (*epProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
		        double ep_imag = (*epProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;

		        std::complex<double> cEp( ep_real, ep_imag + ( sigma / mFreq ) );

				double volInt_Ni = 0.25 * (*gpVolumes)[ *ce_it ];

				// Loop over the nodes of (*ce_it) element
				for( int k=0; k<4; ++k ) 
				{
					unsigned int kId = eInfo[ k ];

				    double volInt_Ni_dNkdX = volInt_Ni * (*gp_dNk_dX)[ *ce_it ][ k ];
					double volInt_Ni_dNkdY = volInt_Ni * (*gp_dNk_dY)[ *ce_it ][ k ];
					double volInt_Ni_dNkdZ = volInt_Ni * (*gp_dNk_dZ)[ *ce_it ][ k ];

					if( Ci_ep_dNkdX[ i ].find( kId ) != Ci_ep_dNkdX[ i ].end() )
					{
                        Ci_ep_dNkdX[ i ][ kId ] += cEp * volInt_Ni_dNkdX;
						Ci_ep_dNkdY[ i ][ kId ] += cEp * volInt_Ni_dNkdY;
						Ci_ep_dNkdZ[ i ][ kId ] += cEp * volInt_Ni_dNkdZ;
				    }
					else
					{
						Ci_ep_dNkdX[ i ][ kId ]  = cEp * volInt_Ni_dNkdX;
						Ci_ep_dNkdY[ i ][ kId ]  = cEp * volInt_Ni_dNkdY;
						Ci_ep_dNkdZ[ i ][ kId ]  = cEp * volInt_Ni_dNkdZ;			
					}

			    }// End nodes in connected element loop - for( int k=0; k<4; ++k ) 

			}//End elements connected to corner node loop - for( ce_it = cnElements.begin(); ce_it != cnElements.end(); ++ce_it )

		}// End corner nodes loop - for( int i=0; i<4; ++i ) 

		// Computing ( Ci, dNk/d{X,Y,Z} ) on face bubble nodes
		for( int i=4; i<8; ++i ) 
		{
			// Elements connected to this face bubble node
			Vector<unsigned int> fbElements; Get_Node_Connectivity( i, fbElements );

			// Iteration over elements connected to face bubble node 
			Vector<unsigned int>::iterator fb_it;

			for( fb_it = fbElements.begin(); fb_it != fbElements.end(); ++fb_it )
			{
				Vector<unsigned int> eInfo = (*mpElementsInfo)[ *fb_it ];

				Properties::Pointer epProperties = (*gpPropertiesMap)[ eInfo[4] ];

			    double sigma   = (*epProperties)(IHL_ELECTRIC_CONDUCTIVITY);
			    double ep_real = (*epProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
		        double ep_imag = (*epProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;

		        std::complex<double> cEp( ep_real, ep_imag + ( sigma / mFreq ) );

				double volInt_FBi = ( 27.0 / 120.0 ) * (*gpVolumes)[ *fb_it ];

				// Loop over the nodes of (*fb_it) element
				for( int k=0; k<4; ++k ) 
				{
					unsigned int kId = eInfo[ k ];

				    double volInt_FBi_dNkdX = volInt_FBi * (*gp_dNk_dX)[ *fb_it ][ k ];
					double volInt_FBi_dNkdY = volInt_FBi * (*gp_dNk_dY)[ *fb_it ][ k ];
					double volInt_FBi_dNkdZ = volInt_FBi * (*gp_dNk_dZ)[ *fb_it ][ k ];

					if( Ci_ep_dNkdX[ i ].find( kId ) != Ci_ep_dNkdX[ i ].end() )
					{
                        Ci_ep_dNkdX[ i ][ kId ] += cEp * volInt_FBi_dNkdX;
						Ci_ep_dNkdY[ i ][ kId ] += cEp * volInt_FBi_dNkdY;
						Ci_ep_dNkdZ[ i ][ kId ] += cEp * volInt_FBi_dNkdZ;
				    }
					else
					{
						Ci_ep_dNkdX[ i ][ kId ]  = cEp * volInt_FBi_dNkdX;
						Ci_ep_dNkdY[ i ][ kId ]  = cEp * volInt_FBi_dNkdY;
						Ci_ep_dNkdZ[ i ][ kId ]  = cEp * volInt_FBi_dNkdZ;			
					}

			    }// End nodes in connected element loop - for( int k=0; k<4; ++k ) 

			}// End elements connected to face bubble node loop - for( fb_it = fbElements.begin(); fb_it != fbElements.end(); ++fb_it )

		}// End face bubble nodes loop - for( int i=4; i<8; ++i ) 

		// Computing ( Ci, dNk/d{X,Y,Z} ) on bubble node
        Matrix<double> DN; 
		
		Calculate_Ni_Derivatives( DN );

		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY)       ;
		double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
		double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;

		std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );

        double volInt_Bi = ( 256.0 / 840.0 ) * mVolume;

	    for( int k=0; k<4; ++k ) 
	    {
	    	unsigned int kId = mNodes[ k ]->Id() - 1;
	    
	        double volInt_Bi_dNkdX = volInt_Bi * DN[ 0 ][ k ];
	    	double volInt_Bi_dNkdY = volInt_Bi * DN[ 1 ][ k ];
	    	double volInt_Bi_dNkdZ = volInt_Bi * DN[ 2 ][ k ];

			Ci_ep_dNkdX[ 8 ][ kId ] = cEp * volInt_Bi_dNkdX;
	    	Ci_ep_dNkdY[ 8 ][ kId ] = cEp * volInt_Bi_dNkdY;
	    	Ci_ep_dNkdZ[ 8 ][ kId ] = cEp * volInt_Bi_dNkdZ;
	    }
	}    

    //**********************************************************************************************************************
    // - Calculates the number of neighbour elements for the pair of local nodes i, j
    //**********************************************************************************************************************
    int VolumeElement_1bb::NumberOfNeigbours( int locNodeId_i, int locNodeId_j )
	{ 
		// If a node in the pair is an inner bubble then return 1
		if( locNodeId_i == 8 || locNodeId_j == 8 )
		{
		    return 1;
		}
		
		// If the nodes in the pair are equal then return the elements connected to the first one
		if( locNodeId_i == locNodeId_j )
		{
		    Vector<unsigned int> elements_ii; 
			
			Get_Node_Connectivity( locNodeId_i, elements_ii );
			
			return elements_ii.size();
		}

	    // If nodes i,j are different then return intersection
		Vector<unsigned int> elements_i; Get_Node_Connectivity( locNodeId_i, elements_i );
		Vector<unsigned int> elements_j; Get_Node_Connectivity( locNodeId_j, elements_j );

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
    void VolumeElement_1bb::Get_k_Nodes_List( int locNodeId_i, int locNodeId_j, Vector<unsigned int>& k_Nodes_List )
	{ 
		// Clear vector 
		k_Nodes_List.clear();
		
		// If inner bubble node then return the list of corner nodes of This element
		if ( locNodeId_i == 8 || locNodeId_j == 8 )
		{
		    k_Nodes_List.push_back( mNodes[0]->Id() - 1 );
			k_Nodes_List.push_back( mNodes[1]->Id() - 1 );
			k_Nodes_List.push_back( mNodes[2]->Id() - 1 );
			k_Nodes_List.push_back( mNodes[3]->Id() - 1 );

			return;
		}

		// If the nodes in the pair are equal then return the list of nodes connected to the first one
		if( locNodeId_i == locNodeId_j )
		{
            Vector<unsigned int> elements_i; Get_Node_Connectivity( locNodeId_i, elements_i );

			std::unordered_set<unsigned int> NodesSet;

			Vector<unsigned int>::iterator e_it;
            
            for( e_it = elements_i.begin(); e_it != elements_i.end(); ++e_it ) 
			{
				NodesSet.insert( (*mpElementsInfo)[*e_it].begin(), (*mpElementsInfo)[*e_it].end() - 1 );
			}
				
            k_Nodes_List.assign( NodesSet.begin(), NodesSet.end() );

		    return;
		}

		// If nodes i,j are different then return the list of common nodes
		Vector<unsigned int> elements_i; Get_Node_Connectivity( locNodeId_i, elements_i );
		Vector<unsigned int> elements_j; Get_Node_Connectivity( locNodeId_j, elements_j );
        
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
    // - Gives the elements connected to the node with the local Id == locNodeId ( < 8 )
    //**********************************************************************************************************************
	void VolumeElement_1bb::Get_Node_Connectivity( int locNodeId, Vector<unsigned int>& Node_Connectivity )
	{
		// Corener nodes
		if( locNodeId < 4 )
		{
		    Node_Connectivity = (*mpNodesConnectivity)[ mNodes[ locNodeId ]->Id() - 1 ];
			return;
		}

		// Node connectivity for face bubble nodes
		int lcId_i, lcId_j, lcId_k;

		if      ( locNodeId == 4 ) { lcId_i = 0; lcId_j = 1; lcId_k = 2; }
		else if ( locNodeId == 5 ) { lcId_i = 0; lcId_j = 1; lcId_k = 3; }
		else if ( locNodeId == 6 ) { lcId_i = 1; lcId_j = 2; lcId_k = 3; }
		else if ( locNodeId == 7 ) { lcId_i = 0; lcId_j = 2; lcId_k = 3; }
       
		Vector<unsigned int> elements_i = (*mpNodesConnectivity)[ mNodes[ lcId_i ]->Id() - 1 ];
		Vector<unsigned int> elements_j = (*mpNodesConnectivity)[ mNodes[ lcId_j ]->Id() - 1 ];
		Vector<unsigned int> elements_k = (*mpNodesConnectivity)[ mNodes[ lcId_k ]->Id() - 1 ];
        
		Vector<unsigned int> vIntsct_ij; 

		std::sort( elements_i.begin(), elements_i.end() ); 
        std::sort( elements_j.begin(), elements_j.end() ); 

        std::set_intersection( elements_i.begin(), elements_i.end(), 
                               elements_j.begin(), elements_j.end(), 
                               std::back_inserter( vIntsct_ij )   ); 

		Vector<unsigned int> vIntsct_kij; 

		std::sort( vIntsct_ij.begin(), vIntsct_ij.end() ); 
        std::sort( elements_k.begin(), elements_k.end() ); 

        std::set_intersection( vIntsct_ij.begin(), vIntsct_ij.end(), 
                               elements_k.begin(), elements_k.end(), 
                               std::back_inserter( vIntsct_kij )  ); 

		Node_Connectivity = vIntsct_kij;
	}

	//**********************************************************************************************************************
    // - 
    //**********************************************************************************************************************
	void VolumeElement_1bb::Set_Contact_Info( Vector< Vector<Node::Pointer> >& CtC_pNodes, Vector< std::complex<double>  >& CtC_Properties )
	{
		return;
		
		mCtC_pNodes     = CtC_pNodes;
        mCtC_Properties = CtC_Properties;

		ctcNumNodes = mNumNodes;
       
		for( int i=0; i<4; i++ )
		{
		    if ( mCtC_pNodes[i].size() > 0 ) 
			{
				for( int fn=0; fn<6; fn++ ) mNodes.push_back( mCtC_pNodes[i][fn] );
				
				ctcNumNodes += 6;
			}
		}

        ctcNumDofs = ctcNumNodes*3; 

		/*
		int ThereIsACtCFace = 0;
		int FaceId;

		for( int iface=0; iface<4; iface++ ) if ( mCtC_pNodes[iface].size()>0 ) { ThereIsACtCFace++; FaceId = iface; }

		if ( ThereIsACtCFace > 0 )
		{
			for ( int i=0; i<6; i++ ) 
			{
				mNodes.push_back( mCtC_pNodes[FaceId][i] );
			}
		     
		    ctcNumNodes = mNumNodes   + 6;
		    ctcNumDofs  = ctcNumNodes * 3;
		}
*/
		/*
		//mCTC_FaceUse.resize(4,false);

		//Vector< Vector<int> > Faces; Faces.resize( 4 );

		//Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
  //      Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		//Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		//Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;

		//mCTC_DuplicateFace.resize(4);
		//

		//if ( mCtC_pNodes[0].size() > 0 )
		//{
		//    mCTC_FaceUse[0] = true;
		//	mCTC_DuplicateFace[0].resize( 6 );

		//	for ( int i=0; i<6; i++ ) 
		//	{
		//		mNodes.push_back( mCtC_pNodes[0][i] );
		//		ctcNumNodes++;
		//        ctcNumDofs += 3;
		//		mCTC_DuplicateFace[0][i] = 17+i;
		//	}
		//}

		//if ( mCtC_pNodes[1].size() > 0 ) 
		//{
		//     mCTC_FaceUse[1] = true;
		//	 mCTC_DuplicateFace[1].resize( 6 );
		//	
		//	if ( mCTC_FaceUse[0] == false )
		//	{
		//	    for ( int i=0; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[1][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[1][i] = 17+i;
		//	    }
		//    }
		//	else
		//	{
		//		mNodes.push_back( mCtC_pNodes[1][0] );
		//	    ctcNumNodes++;
		//        ctcNumDofs += 3;

		//		mCTC_DuplicateFace[1][0] = 23;
		//		mCTC_DuplicateFace[1][1] = 18;
		//		mCTC_DuplicateFace[1][2] = 17;

		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[1][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[1][i] = 21+i;
		//	    }
		//	}
		//}

		//if ( mCtC_pNodes[2].size() > 0 ) 
		//{
		//     mCTC_FaceUse[2] = true;
		//	 mCTC_DuplicateFace[2].resize( 6 );
		//	
		//	if ( (mCTC_FaceUse[0] == false) && (mCTC_FaceUse[1] == false) )
		//	{
		//	    for ( int i=0; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[2][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[2][i] = 17+i;
		//	    }
		//    }
		//	else if ( (mCTC_FaceUse[0] == true) && (mCTC_FaceUse[1] == false) )
		//	{
		//		mNodes.push_back( mCtC_pNodes[2][1] );
		//	    ctcNumNodes++;
		//        ctcNumDofs += 3;

		//		mCTC_DuplicateFace[2][0] = 18;
		//		mCTC_DuplicateFace[2][1] = 23;
		//		mCTC_DuplicateFace[2][2] = 19;

		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[2][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[2][i] = 21+i;
		//	    }
		//	}
		//	else if ( (mCTC_FaceUse[0] == false) && (mCTC_FaceUse[1] == true) )
		//	{
		//		mNodes.push_back( mCtC_pNodes[2][2] );
		//	    ctcNumNodes++;
		//        ctcNumDofs += 3;

		//		mCTC_DuplicateFace[2][0] = 18;
		//		mCTC_DuplicateFace[2][1] = 17;
		//		mCTC_DuplicateFace[2][2] = 23;

		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[2][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[2][i] = 21+i;
		//	    }
		//	}
		//	else if ( (mCTC_FaceUse[0] == true) && (mCTC_FaceUse[1] == true) )
		//	{
		//		mCTC_DuplicateFace[2][0] = 19;
		//		mCTC_DuplicateFace[2][1] = 23;
		//		mCTC_DuplicateFace[2][2] = 17;
		//		
		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[2][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[2][i] = 24+i;
		//	    }
		//	}
		//}

		//if ( mCtC_pNodes[3].size() > 0 ) 
		//{
		//     mCTC_FaceUse[3] = true;
		//	 mCTC_DuplicateFace[3].resize( 6 );
		//	
		//	if ( (mCTC_FaceUse[0] == false) && (mCTC_FaceUse[1] == false) && (mCTC_FaceUse[2] == false) )
		//	{
		//	    for ( int i=0; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[3][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[3][i] = 17+i;
		//	    }
		//    }
		//	else if  ( (mCTC_FaceUse[0] == false) && (mCTC_FaceUse[1] == false) && (mCTC_FaceUse[2] == true) )
		//	{
		//		mNodes.push_back( mCtC_pNodes[3][2] );
		//	    ctcNumNodes++;
		//        ctcNumDofs += 3;

		//		mCTC_DuplicateFace[3][0] = 19;
		//		mCTC_DuplicateFace[3][1] = 18;
		//		mCTC_DuplicateFace[3][2] = 23;

		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[3][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[3][i] = 21+i;
		//	    }
		//	}
		//	else if  ( (mCTC_FaceUse[0] == false) && (mCTC_FaceUse[1] == true) && (mCTC_FaceUse[2] == false) )
		//	{
		//		mNodes.push_back( mCtC_pNodes[3][0] );
		//	    ctcNumNodes++;
		//        ctcNumDofs += 3;

		//		mCTC_DuplicateFace[3][0] = 23;
		//		mCTC_DuplicateFace[3][1] = 17;
		//		mCTC_DuplicateFace[3][2] = 19;

		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[3][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;
		//			mCTC_DuplicateFace[3][i] = 21+i;
		//	    }
		//	}
		//	else if ( (mCTC_FaceUse[0] == true) && (mCTC_FaceUse[1] == false) && (mCTC_FaceUse[2] == false) )
		//	{
		//		mNodes.push_back( mCtC_pNodes[3][1] );
		//	    ctcNumNodes++;
		//        ctcNumDofs += 3;

		//		mCTC_DuplicateFace[3][0] = 19;
		//		mCTC_DuplicateFace[3][1] = 23;
		//		mCTC_DuplicateFace[3][2] = 17;
		//		
		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[3][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;

		//			mCTC_DuplicateFace[3][i] = 21+i;
		//	    }
		//	}
		//	else if ( (mCTC_FaceUse[0] == true) && (mCTC_FaceUse[1] == true) && (mCTC_FaceUse[2] == false) )
		//	{
		//		mCTC_DuplicateFace[3][0] = 19;
		//		mCTC_DuplicateFace[3][1] = 23;
		//		mCTC_DuplicateFace[3][2] = 17;
		//		
		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[3][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;

		//			mCTC_DuplicateFace[3][i] = 24+i;
		//	    }
		//	}
		//	else if ( (mCTC_FaceUse[0] == true) && (mCTC_FaceUse[1] == false) && (mCTC_FaceUse[2] == true) )
		//	{
		//		mCTC_DuplicateFace[3][0] = 19;
		//		mCTC_DuplicateFace[3][1] = 23;
		//		mCTC_DuplicateFace[3][2] = 17;
		//		
		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[3][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;

		//			mCTC_DuplicateFace[3][i] = 24+i;
		//	    }
		//	}
		//	else if ( (mCTC_FaceUse[0] == false) && (mCTC_FaceUse[1] == true) && (mCTC_FaceUse[2] == true) )
		//	{
		//		mCTC_DuplicateFace[3][0] = 19;
		//		mCTC_DuplicateFace[3][1] = 17;
		//		mCTC_DuplicateFace[3][2] = 23;
		//		
		//		for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[3][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;

		//			mCTC_DuplicateFace[3][i] = 24+i;
		//	    }
		//	}
		//	else if ( (mCTC_FaceUse[0] == true) && (mCTC_FaceUse[1] == true) && (mCTC_FaceUse[2] == true) )
		//	{
		//		mCTC_DuplicateFace[3][0] = 19;
		//		mCTC_DuplicateFace[3][1] = 23;
		//		mCTC_DuplicateFace[3][2] = 17;
		//	
		//	    for ( int i=3; i<6; i++ ) 
		//	    {
		//	    	mNodes.push_back( mCtC_pNodes[3][i] );
		//	    	ctcNumNodes++;
		//            ctcNumDofs += 3;

		//			mCTC_DuplicateFace[3][i] = 27+i;
		//	    }
		//	}
		//}
        */
	}

	//**********************************************************************************************************************
    // - Calculates the stiffness matrix with the selected formulation
    //**********************************************************************************************************************
    void VolumeElement_1bb::GetStiffMatrix( Matrix< std::complex<double> >& StiffMatrix )
	{   
		StiffMatrix.Resize( ctcNumDofs, ctcNumDofs, mZero );

		////Matrix< std::complex<double> > StiffMatrix_PRE;
		////GetStiffMatrix_Plain_Regular( StiffMatrix_PRE );
		////StiffMatrix += StiffMatrix_PRE;
		//
		////GetStiffMatrix_PLL2P_SimpTau( StiffMatrix );
		////GetStiffMatrix_PLL2P_FullTau( StiffMatrix );
		//
		////Matrix< std::complex<double> > StiffMatrix_AV;
  //      //GetStiffMatrix_LL2P_VolumAvg ( StiffMatrix_AV ); 
		////StiffMatrix += StiffMatrix_AV;
		//
		Matrix< std::complex<double> > StiffMatrix_CDE;
		GetStiffMatrix_LL2P_CurlDivE ( StiffMatrix_CDE ); 
		StiffMatrix += StiffMatrix_CDE;
	 //   
		////Matrix< std::complex<double> > StiffMatrix_CDE_CtC;
		////GetStiffMatrix_LL2P_CurlDivE_CtC ( StiffMatrix_CDE_CtC ); 
		////StiffMatrix += StiffMatrix_CDE_CtC;
		//
		////GetStiffMatrix_LL2P_CurlDivE_ctc( StiffMatrix_CDE ); 
		////GetStiffMatrix_LL2P_EdgeElem( StiffMatrix );
		//
		////Matrix< std::complex<double> > StiffMatrix_CDQ;
		////GetStiffMatrix_LL2P_CurlDivQ( StiffMatrix_CDQ ); 
		////StiffMatrix += StiffMatrix_CDQ;
		//
		////Matrix< std::complex<double> > StiffMatrix_CtC;
		////GetStiffMatrix_LL2P_CtC      ( StiffMatrix_CtC ); 
		////StiffMatrix += StiffMatrix_CtC;
		//
		Matrix< std::complex<double> > StiffMatrix_hk2;
		Add_hk2CurlDiv_Stabilization ( StiffMatrix_hk2 ); 
		StiffMatrix += StiffMatrix_hk2;

		//if ( mPeso != 0.00 ) GetStiffMatrix_Plain_Regular( StiffMatrix );
		//else                 GetStiffMatrix_LL2P_CurlDivE( StiffMatrix ); 

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//// Material properties
        //double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        //double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        //double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        //double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        //double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;
        
        //std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
        //std::complex<double> cMu( mu_real , mu_imag );

		//double hk = Calculate_Element_Diameter();
		//double D  = 1e-3;

		//std::complex<double> NormCte = std::exp( -mFreq * std::sqrt( std::abs(cEp) * std::abs(cMu) ) * D     );
		//std::complex<double> NormCte = std::exp( -mFreq * mFreq * std::abs(cEp) * std::abs(cMu) * D * D );
		//std::complex<double> NormCte = ( 0.0 / ( cMu * D * D  ) );
        //std::complex<double> NormCte = 0.0;
		//std::complex<double> NormCte = ( 1.00 / ( cMu * hk * hk ) );
		//std::complex<double> NormCte = ( 1.00e-6 / ( cMu ) );

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Matrix< std::complex<double> > Sh_div, Sh_cur  ; 
		//Add_DivCurGBk_Stabilization  ( Sh_div, Sh_cur ); 
		//Sh_div *= NormCte;
		//Sh_cur *= NormCte;
		//StiffMatrix += Sh_div;
		//StiffMatrix += Sh_cur;

	    //Matrix< std::complex<double> > CtC_Sh_div; 
		//Add_CtCDiv_Stabilization     ( CtC_Sh_div ); 
		//CtC_Sh_div *= NormCte;
		//StiffMatrix += CtC_Sh_div;

		//Matrix< std::complex<double> > nxn_Sh_cur; 
		//Add_nxnCur_Stabilization     ( nxn_Sh_cur ); 
		//nxn_Sh_cur *= NormCte;
		//StiffMatrix += nxn_Sh_cur;

		////////////////////////////////////////////////////////////////////////

		return;

		//StiffMatrix.Resize( 12, 12, mZero );

		//Matrix< std::complex<double> > Full_StiffMatrix;
		//Full_StiffMatrix.Resize( ctcNumDofs, ctcNumDofs, mZero );

		//Matrix< std::complex<double> > StiffMatrix_CDE;
		//GetStiffMatrix_LL2P_CurlDivE ( StiffMatrix_CDE ); 
		//Full_StiffMatrix += StiffMatrix_CDE;
	 //   
		//Matrix< std::complex<double> > StiffMatrix_hk2;
		//Add_hk2CurlDiv_Stabilization ( StiffMatrix_hk2 ); 
		//Full_StiffMatrix += StiffMatrix_hk2;

		////////////////////////////////////////////////////////////////////////

		 //Loop counters
		int i, j, n, gp;

		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order8( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< double > C_cur; C_cur.Resize( 12, mNumDofs, 0.00 );

		for( int i=0; i<4; i++ )
		{
		    for( j=0; j<mNumNodes; j++ )
		    {
		    	double vintg_Ni_Nj = 0.00;
		        
                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	vintg_Ni_Nj += W[gp] * N[i][gp] * N[j][gp]; 
                }      
		    	
				C_cur[ i     ][ j               ] = vintg_Ni_Nj; 
                C_cur[ i + 4 ][ j + mNumNodes   ] = vintg_Ni_Nj;  						 
                C_cur[ i + 8 ][ j + mNumNodes*2 ] = vintg_Ni_Nj; 
		    }
		}

		// Projection basis matrices
		Matrix< double > iA_div; iA_div.Resize(  4,  4, 0.00 );
		Matrix< double > iA_cur; iA_cur.Resize( 12, 12, 0.00 );
		
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

        // 12 x 51
		Matrix< double > iA_x_C = ( iA_cur * C_cur );
		
		////////////////////////////////////////////////////////////////////////////////////////

		int nRows = 12;
		int nCols = 51;
		
		double **A_svd = (double**) malloc(sizeof(double*)*nRows);
        for ( i=0; i<nRows; i++ ) A_svd[i] = (double*) malloc(sizeof(double)*nCols);
		for ( i=0; i<nRows; i++ ) { for ( j=0; j<nCols; j++ ) { A_svd[i][j] = iA_x_C[i][j]; } }

	    double **V_svd = (double**) malloc(sizeof(double*)*nCols);
        for ( i=0; i<nCols; i++ ) V_svd[i] = (double*) malloc(sizeof(double)*nCols);
		for ( i=0; i<nCols; i++ ) { for ( j=0; j<nCols; j++ ) { V_svd[i][j] = 0.00; } }

        double *W_svd = (double*) malloc(sizeof(double*)*nCols);
		for ( i=0; i<nCols; i++ ) { W_svd[i] = 0.00; }

        /**************************************************************************
        Modified from Numerical Recipes in C
        Given a matrix a[nRows][nCols], svdcmp() computes its singular value
        decomposition, A = U * W * Vt.  A is replaced by U when svdcmp
        returns.  The diagonal matrix W is output as a vector w[nCols].
        V (not V transpose) is output as the matrix V[nCols][nCols].
        ***************************************************************************/
	    svdcmp( A_svd, nRows, nCols, W_svd, V_svd );

		Matrix< double > _U_ ( nRows, nRows, 0.00 );
		for ( i=0; i<nRows; i++ ) { for ( j=0; j<nRows; j++ ) { _U_[i][j] = A_svd[i][j]; } }
        
		//std::cout<<_U_<<std::endl;
		
		Matrix< double > _V_( nCols, nCols, 0.00 );
		for ( i=0; i<nCols; i++ ) { for ( j=0; j<nCols; j++ ) { _V_[i][j] = V_svd[i][j]; } }

		//std::cout<<_V_<<std::endl;

        Matrix< double > _W_( nRows, nCols, 0.00);
		for ( i=0; i<nRows; i++ ) { for ( j=0; j<nCols; j++ ) { if ( i==j ) _W_[i][j] = W_svd[j]; } }

		//std::cout<<_W_<<std::endl;

        Matrix< double > _invW_( nCols, nRows, 0.00 );
		for ( i=0; i<nCols; i++ ) { for ( j=0; j<nRows; j++ ) { if ( ( i==j ) && ( abs( W_svd[j] ) > 0.0 ) ) _invW_[i][j] = (1.00/W_svd[j]); } }

		//std::cout<<_invW_<<std::endl;

  //      printMatrix(A_svd, nRows, nRows);
		//std::cout<<_U_<<std::endl;

  //      printVector(W_svd, nCols);
		//std::cout<<_W_<<std::endl;

		//printMatrix(V_svd, nCols, nCols);
		//std::cout<<_V_<<std::endl;


		Matrix< double > V_invW    = ( _V_ * _invW_    );
        Matrix< double > U_W       = ( _U_ *    _W_    );

		Matrix< double > Ut        = ( _U_.Transpose() );
        Matrix< double > Vt        = ( _V_.Transpose() );

        Matrix< double > V_invW_Ut = ( V_invW * Ut     );
        Matrix< double > U_W_Vt    = ( U_W    * Vt     );

        
		//std::cout<<U_<<std::endl;
		//std::cout<<_W_<<std::endl;
		//std::cout<<_Vt<<std::endl;
        
		

		std::cout<<iA_x_C<<std::endl;
		std::cout<<U_W_Vt<<std::endl;
		std::cout<<V_invW_Ut*iA_x_C<<std::endl;


		exit(EXIT_SUCCESS);
	
		//////////////////////////////////////////////////////////////////////////////////////////

		Matrix< double > real_pinv_iAxC = V_invW_Ut;// real_pinv_iAxC.Resize( 51, 12, 0.00 );

        // 51 x 12
        Matrix< std::complex<double> > pinv_iAxC; pinv_iAxC.Resize( 51, 12, mZero );
		
        for ( i=0; i<51; i++ )
		{
			for ( j=0; j<12; j++ ) 
			{
				std::complex<double> aij( real_pinv_iAxC[i][j], 0.00 ); pinv_iAxC[i][j] = aij;
			}
		}

		/*Matrix< std::complex<double> > Full_x_pinv = ( Full_StiffMatrix * pinv_iAxC );

		Matrix< std::complex<double> > pinv_iAxC_tran = pinv_iAxC.Transpose();

		StiffMatrix += ( pinv_iAxC_tran * Full_x_pinv );*/
	}

    //**********************************************************************************************************************
    // - Calculation of the stiffness matrix with the RME L2 projection method - Homogeneous media
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

    //**********************************************************************************************************************
    // - Calculation of the stiffness matrix with the RME L2 projection method - Inhomogeneous media
    //**********************************************************************************************************************
    void VolumeElement_1bb::GetStiffMatrix_InHomogeneous( Matrix< std::complex<double> >& StiffMatrix )
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
                StiffMatrix[ j ][ i ] = StiffMatrix[ i ][ j ];
            }
        }      
	}

	//***********************************************************************************************
    //* - Calculation of the stiffness matrix with the RME method
    //***********************************************************************************************	
    void VolumeElement_1bb::GetStiffMatrix_Plain_Regular( Matrix<std::complex<double> >& StiffMatrix ) 
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

		//std::complex<double> cEp( mu_real , mu_imag );
  //      std::complex<double> cMu( eps_real, eps_imag + ( sigma / mFreq ) );

		// Gaussian derivation weights
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order6( cX, cY, cZ, W );

		double detJ = 6.00 * Calculate_Volume();

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions and derivatives
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		//mPeso = 0.0;

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

					//double intCurlDivX = 0.00;
					//double intCurlDivY = 0.00;
					//double intCurlDivZ = 0.00;

					for ( gp=0; gp<nGaussPoints; gp++ )
					{
						intCurlDiv += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp] );
						//intCurlDivX += W[gp] * ( dNdx[i][gp] * dNdx[j][gp] );
						//intCurlDivY += W[gp] * ( dNdy[i][gp] * dNdy[j][gp] );
						//intCurlDivZ += W[gp] * ( dNdz[i][gp] * dNdz[j][gp] );
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

	//**************************************************************************************************
    //* - Adds a stabilization term S_h = 
    //**************************************************************************************************	
	void VolumeElement_1bb::Add_DivCurGBk_Stabilization( Matrix< std::complex<double> >& Sh_div, Matrix< std::complex<double> >& Sh_cur ) 
	{
		Sh_div.Resize( ctcNumDofs, ctcNumDofs, mZero );
		Sh_cur.Resize( ctcNumDofs, ctcNumDofs, mZero );

	    // Loop counters
		int i, j, l, gp;
      
		// Gaussian derivation weights
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order10( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Derivatives
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Basis functions T
		Matrix<double> Tk;
        Matrix<double> dTkdx, dTkdy, dTkdz;

		///////////////////////////////////////////////////////////////////////////////////////

		//int TkSize = 4; 
		//Lagrange3D_Ni_1st        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_1st( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		//int TkSize = 10; 
		//Lagrange3D_Ni_2nd        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_2nd( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

        //int TkSize = 20; 
		//Lagrange3D_Ni_3th        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_3th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		int TkSize = 35; 
		Lagrange3D_Ni_4th        ( Tk,                  cX, cY, cZ );
		LagrangeDerivatives3D_4th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		////////////////////////////////////////////////////////////////////////////////////////

		// Norm
		double gradX, gradY, gradZ;
		double norm_gradTkBk = 0.00;
		
		for ( l=0; l<TkSize; l++ ) 
		{
			double norm_gradTkBk_l = 0.00;

			for ( gp=0; gp<nGaussPoints; gp++ ) 
			{	
				gradX = N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp];
				gradY = N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp];
				gradZ = N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp];
				
				norm_gradTkBk_l += W[gp] * ( gradX*gradX + gradY*gradY + gradZ*gradZ );
			}

			norm_gradTkBk += norm_gradTkBk_l;
		}
		
		// ( u, gradTB ) * ( u, gradTB )
		for ( l=0; l<TkSize; l++ ) 
		{
		    Matrix< std::complex<double> > Nj_gradTBl( 1, ctcNumDofs );
			
			for ( j=0; j<mNumNodes; j++ ) 
		    {
				double Nj_gradTB_X = 0.00;
				double Nj_gradTB_Y = 0.00;
				double Nj_gradTB_Z = 0.00;
				
				for ( gp=0; gp<nGaussPoints; gp++ ) 
				{	
		    		Nj_gradTB_X += W[gp] * N[j][gp] * ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] );
					Nj_gradTB_Y += W[gp] * N[j][gp] * ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] );
					Nj_gradTB_Z += W[gp] * N[j][gp] * ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] );
				}

				Nj_gradTBl[0][ j               ] = Nj_gradTB_X;
				Nj_gradTBl[0][ j + mNumNodes   ] = Nj_gradTB_Y;
				Nj_gradTBl[0][ j + mNumNodes*2 ] = Nj_gradTB_Z;
			}
			
			Matrix< std::complex<double> > Nj_gradTB      = Nj_gradTBl;
			Matrix< std::complex<double> > Nj_gradTB_tran = Nj_gradTBl.Transpose();

			Sh_div += ( Nj_gradTB_tran * Nj_gradTB );
		}

		Sh_div /= norm_gradTkBk;

		////////////////////////////////////////////////////////////////////////////////////////

		// Norm curl
		double curlX, curlY, curlZ;
		double norm_curlTkBk = 0.00;
		
		for ( l=0; l<TkSize; l++ ) 
		{
			double norm_curlTkBk_l = 0.00;

			for ( gp=0; gp<nGaussPoints; gp++ ) 
			{	
				curlX = ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] ) - ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] );
				curlY = ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] ) - ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] );
				curlZ = ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] ) - ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] );
				
				norm_curlTkBk_l += W[gp] * ( curlX*curlX + curlY*curlY + curlZ*curlZ );
			}

			norm_curlTkBk += norm_curlTkBk_l;
		}
		
		// ( u, curl TkBk) * ( v, curl TkBk)
		for ( l=0; l<TkSize; l++ ) 
		{
		    Matrix< std::complex<double> > Nj_curlTBl( 1, ctcNumDofs );
			
			for ( j=0; j<mNumNodes; j++ ) 
		    {
				double Nj_curlTBl_X = 0.00;
				double Nj_curlTBl_Y = 0.00;
				double Nj_curlTBl_Z = 0.00;
				
				for ( gp=0; gp<nGaussPoints; gp++ ) 
				{	
		    		Nj_curlTBl_X += W[gp] * N[j][gp] * ( ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] ) - ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] ) );
					Nj_curlTBl_Y += W[gp] * N[j][gp] * ( ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] ) - ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] ) );
					Nj_curlTBl_Z += W[gp] * N[j][gp] * ( ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] ) - ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] ) );
				}

				Nj_curlTBl[ 0 ][ j               ] = Nj_curlTBl_X;
				Nj_curlTBl[ 0 ][ j + mNumNodes   ] = Nj_curlTBl_Y;
				Nj_curlTBl[ 0 ][ j + mNumNodes*2 ] = Nj_curlTBl_Z;
			}

			Matrix< std::complex<double> > Nj_curlTB      = Nj_curlTBl;
			Matrix< std::complex<double> > Nj_curlTB_tran = Nj_curlTBl.Transpose();

			Sh_cur += ( Nj_curlTB_tran * Nj_curlTB );
		}

		Sh_cur /= norm_curlTkBk;
	}

    //**************************************************************************************************
    //* - Adds a stabilization term S_h = 
    //**************************************************************************************************	
    void VolumeElement_1bb::Add_CtCDiv_Stabilization( Matrix< std::complex<double> >& CtC_Sh_div ) 
	{
		CtC_Sh_div.Resize( mNumDofs, mNumDofs, mZero );

		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;

        Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			// 'E' = PEC normal | 'D' = Plasma Sheath | 'H' = PMC | 'T' = TEPMC
			//if ( ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'C' ) && ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'R' ) &&
			//	 ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'H' ) && ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'T' ) ) continue;

			if ( ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'C' ) && ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'R' ) ) continue;

			Matrix< std::complex<double> > Face_CtC_Sh_div; 
			
			Face_CtC_Sh_div.Resize( mNumDofs, mNumDofs, mZero );

	        // Loop counters
		    int i, j, l, gp;
		    
		    // Gaussian derivation weights
		    std::vector<double> cX; 
		    std::vector<double>	cY;
		    std::vector<double> cZ;
		    std::vector<double>  W;
		    
		    int nGaussPoints = GaussPoints3D_Order10( cX, cY, cZ, W );
		    
		    double detJ = 6.00 * mVolume;
		    
		    for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 
		    
            // Derivatives
		    Matrix<double> dNdx, dNdy, dNdz;
		    
		    LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );
		    
		    // Basis functions on Gauss points
		    Matrix<double> N;
		    
		    Lagrange3D_Ni_2bb( N, cX, cY, cZ );
		    
		    // Basis functions T
		    Matrix<double> Tk;
            Matrix<double> dTkdx, dTkdy, dTkdz;
		    
		    //////////////////////////////////////////////////////////////////////////
		    
		    //int TkSize = 4; 
		    //Lagrange3D_Ni_1st        ( Tk,                  cX, cY, cZ );
		    //LagrangeDerivatives3D_1st( dTkdx, dTkdy, dTkdz, cX, cY, cZ );
		    
		    //int TkSize = 10; 
		    //Lagrange3D_Ni_2nd        ( Tk,                  cX, cY, cZ );
		    //LagrangeDerivatives3D_2nd( dTkdx, dTkdy, dTkdz, cX, cY, cZ );
		    
            //int TkSize = 20; 
		    //Lagrange3D_Ni_3th        ( Tk,                  cX, cY, cZ );
		    //LagrangeDerivatives3D_3th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );
		    
		    int TkSize = 35; 
		    Lagrange3D_Ni_4th        ( Tk,                  cX, cY, cZ );
		    LagrangeDerivatives3D_4th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );
		    
		    ////////////////////////////////////////////////////////////////////////////////////////

			// Cambiar:  ahora N[16] sera la bubble asociada a la cara N1*N2*N3

			//Vector<double> Bf(nGaussPoints);

			//for ( gp=0; gp<nGaussPoints; gp++ ) Bf[gp] = N[ nf[0] ][ gp ] * N[ nf[1] ][ gp ] * N[ nf[2] ][ gp ];
            //     
            //Vector<double> dBfdx(nGaussPoints);
			//Vector<double> dBfdy(nGaussPoints);
			//Vector<double> dBfdz(nGaussPoints);

			//for ( gp=0; gp<nGaussPoints; gp++ )
			//{
			//	dBfdx[gp] = dNdx[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
			//		           N[ nf[0] ][ gp ] * dNdx[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
			//		           N[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] * dNdx[ nf[2] ][ gp ];

			//    dBfdy[gp] = dNdy[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
			//		           N[ nf[0] ][ gp ] * dNdy[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
			//		           N[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] * dNdy[ nf[2] ][ gp ];

			//	dBfdz[gp] = dNdz[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
			//		           N[ nf[0] ][ gp ] * dNdz[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
			//		           N[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] * dNdz[ nf[2] ][ gp ];
			//}

			///////////////////////////////////////////////////////////////////////////////////////////
		 //   
		 //   // Norm
		 //   double gradX, gradY, gradZ;
		 //   double norm_gradTkBk = 0.00;
		 //   
		 //   for ( l=0; l<TkSize; l++ ) 
		 //   {
		 //   	double norm_gradTkBk_l = 0.00;
		 //   
		 //   	for ( gp=0; gp<nGaussPoints; gp++ ) 
		 //   	{	
		 //   		gradX = Bf[gp] * dTkdx[l][gp] + Tk[l][gp] * dBfdx[gp];
		 //   		gradY = Bf[gp] * dTkdy[l][gp] + Tk[l][gp] * dBfdy[gp];
		 //   		gradZ = Bf[gp] * dTkdz[l][gp] + Tk[l][gp] * dBfdz[gp];
		 //   		
		 //   		norm_gradTkBk_l += W[gp] * ( gradX*gradX + gradY*gradY + gradZ*gradZ );
		 //   	}
		 //   
		 //   	norm_gradTkBk += norm_gradTkBk_l;
		 //   }
		 //   
		 //   // ( u, gradTB ) * ( u, gradTB )
		 //   for ( l=0; l<TkSize; l++ ) 
		 //   {
		 //       Matrix< std::complex<double> > Nj_gradTBl( 1, ctcNumDofs );
		 //   	
		 //   	for ( j=0; j<mNumNodes; j++ ) 
		 //       {
		 //   		double Nj_gradTB_X = 0.00;
		 //   		double Nj_gradTB_Y = 0.00;
		 //   		double Nj_gradTB_Z = 0.00;
		 //   		
		 //   		for ( gp=0; gp<nGaussPoints; gp++ ) 
		 //   		{	
		 //       		Nj_gradTB_X += W[gp] * N[j][gp] * ( Bf[gp] * dTkdx[l][gp] + Tk[l][gp] * dBfdx[gp] );
		 //   			Nj_gradTB_Y += W[gp] * N[j][gp] * ( Bf[gp] * dTkdy[l][gp] + Tk[l][gp] * dBfdy[gp] );
		 //   			Nj_gradTB_Z += W[gp] * N[j][gp] * ( Bf[gp] * dTkdz[l][gp] + Tk[l][gp] * dBfdz[gp] );
		 //   		}
		 //   
		 //   		Nj_gradTBl[0][ j               ] = Nj_gradTB_X;
		 //   		Nj_gradTBl[0][ j + mNumNodes   ] = Nj_gradTB_Y;
		 //   		Nj_gradTBl[0][ j + mNumNodes*2 ] = Nj_gradTB_Z;
		 //   	}
		 //   	
		 //   	Matrix< std::complex<double> > Nj_gradTB      = Nj_gradTBl;
		 //   	Matrix< std::complex<double> > Nj_gradTB_tran = Nj_gradTBl.Transpose();
		 //   
		 //   	Face_CtC_Sh_div += ( Nj_gradTB_tran * Nj_gradTB );
		 //   }


		    // Norm
		    double gradX, gradY, gradZ;
		    double norm_gradTkBk = 0.00;
		    
		    for ( l=0; l<TkSize; l++ ) 
		    {
				for( int ib=3; ib<6; ib++ )
				{
                    double norm_gradTkBk_l = 0.00;

		    		for ( gp=0; gp<nGaussPoints; gp++ ) 
		    		{	
		    			gradX = N[nf[ib]][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[nf[ib]][gp];
		    			gradY = N[nf[ib]][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[nf[ib]][gp];
		    			gradZ = N[nf[ib]][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[nf[ib]][gp];
		    		
		    			norm_gradTkBk_l += W[gp] * ( gradX*gradX + gradY*gradY + gradZ*gradZ );
		    		}

                    norm_gradTkBk += norm_gradTkBk_l;
			    }
		    }
		    
		    // ( u, gradTB ) * ( u, gradTB )
		    for ( l=0; l<TkSize; l++ ) 
		    {
				for( int ib=3; ib<6; ib++ )
				{
					Matrix< std::complex<double> > Nj_gradTBl( 1, ctcNumDofs );
		    	
		    		for ( j=0; j<mNumNodes; j++ ) 
					{
		    			double Nj_gradTB_X = 0.00;
		    			double Nj_gradTB_Y = 0.00;
		    			double Nj_gradTB_Z = 0.00;
		    		
		    			for ( gp=0; gp<nGaussPoints; gp++ ) 
		    			{	
		        			Nj_gradTB_X += W[gp] * N[j][gp] * ( N[nf[ib]][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[nf[ib]][gp] );
		    				Nj_gradTB_Y += W[gp] * N[j][gp] * ( N[nf[ib]][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[nf[ib]][gp] );
		    				Nj_gradTB_Z += W[gp] * N[j][gp] * ( N[nf[ib]][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[nf[ib]][gp] );
		    			}
		    
		    			Nj_gradTBl[0][ j               ] = Nj_gradTB_X;
		    			Nj_gradTBl[0][ j + mNumNodes   ] = Nj_gradTB_Y;
		    			Nj_gradTBl[0][ j + mNumNodes*2 ] = Nj_gradTB_Z;
					}

		    	    Matrix< std::complex<double> > Nj_gradTB      = Nj_gradTBl;
		    	    Matrix< std::complex<double> > Nj_gradTB_tran = Nj_gradTBl.Transpose();
		    	    
		    	    Face_CtC_Sh_div += ( Nj_gradTB_tran * Nj_gradTB );
			    }
		    }

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// chekear esto : que pasa si hay varias caras del mismo elemento????
		    Face_CtC_Sh_div /= norm_gradTkBk;
		   
			CtC_Sh_div += Face_CtC_Sh_div;
		}
	}


	//**************************************************************************************************
    //* - Adds a stabilization term S_h = 
    //**************************************************************************************************	
    void VolumeElement_1bb::Add_nxnCur_Stabilization( Matrix< std::complex<double> >& nxn_Sh_cur ) 
	{
		nxn_Sh_cur.Resize( mNumDofs, mNumDofs, mZero );

		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;

        Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			// 'E' = PEC normal | 'D' = Plasma Sheath | 'H' = PMC | 'T' = TEPMC
			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) continue;

			Matrix< std::complex<double> > Face_nxn_Sh_cur; 
			
			Face_nxn_Sh_cur.Resize( mNumDofs, mNumDofs, mZero );

	        // Loop counters
		    int i, j, l, gp;
		    
		    // Gaussian derivation weights
		    std::vector<double> cX; 
		    std::vector<double>	cY;
		    std::vector<double> cZ;
		    std::vector<double>  W;
		    
		    int nGaussPoints = GaussPoints3D_Order10( cX, cY, cZ, W );
		    
		    double detJ = 6.00 * mVolume;
		    
		    for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 
		    
            // Derivatives
		    Matrix<double> dNdx, dNdy, dNdz;
		    
		    LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );
		    
		    // Basis functions on Gauss points
		    Matrix<double> N;
		    
		    Lagrange3D_Ni_2bb( N, cX, cY, cZ );
		    
		    // Basis functions T
		    Matrix<double> Tk;
            Matrix<double> dTkdx, dTkdy, dTkdz;
		    
		    //////////////////////////////////////////////////////////////////////////
		    
		    //int TkSize = 4; 
		    //Lagrange3D_Ni_1st        ( Tk,                  cX, cY, cZ );
		    //LagrangeDerivatives3D_1st( dTkdx, dTkdy, dTkdz, cX, cY, cZ );
		    
		    //int TkSize = 10; 
		    //Lagrange3D_Ni_2nd        ( Tk,                  cX, cY, cZ );
		    //LagrangeDerivatives3D_2nd( dTkdx, dTkdy, dTkdz, cX, cY, cZ );
		    
            //int TkSize = 20; 
		    //Lagrange3D_Ni_3th        ( Tk,                  cX, cY, cZ );
		    //LagrangeDerivatives3D_3th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );
		    
		    int TkSize = 35; 
		    Lagrange3D_Ni_4th        ( Tk,                  cX, cY, cZ );
		    LagrangeDerivatives3D_4th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );
		    
		    ////////////////////////////////////////////////////////////////////////////////////////

			// Cambiar:  ahora N[16] sera la bubble asociada a la cara N1*N2*N3

			Vector<double> Bf(nGaussPoints);

			for ( gp=0; gp<nGaussPoints; gp++ ) Bf[gp] = N[ nf[0] ][ gp ] * N[ nf[1] ][ gp ] * N[ nf[2] ][ gp ];
        
            Vector<double> dBfdx(nGaussPoints);
			Vector<double> dBfdy(nGaussPoints);
			Vector<double> dBfdz(nGaussPoints);

			for ( gp=0; gp<nGaussPoints; gp++ )
			{
				dBfdx[gp] = dNdx[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
					           N[ nf[0] ][ gp ] * dNdx[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
					           N[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] * dNdx[ nf[2] ][ gp ];

			    dBfdy[gp] = dNdy[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
					           N[ nf[0] ][ gp ] * dNdy[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
					           N[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] * dNdy[ nf[2] ][ gp ];

				dBfdz[gp] = dNdz[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
					           N[ nf[0] ][ gp ] * dNdz[ nf[1] ][ gp ] *    N[ nf[2] ][ gp ] + 
					           N[ nf[0] ][ gp ] *    N[ nf[1] ][ gp ] * dNdz[ nf[2] ][ gp ];
			}

			/////////////////////////////////////////////////////////////////////////////////////////
		    
		    // Norm curl
		    double curlX, curlY, curlZ;
		    double norm_curlTkBk = 0.00;
		    
		    for ( l=0; l<TkSize; l++ ) 
		    {
		    	double norm_curlTkBk_l = 0.00;
		    
		    	for ( gp=0; gp<nGaussPoints; gp++ ) 
		    	{	
		    		curlX = ( Bf[gp] * dTkdy[l][gp] + Tk[l][gp] * dBfdy[gp] ) - ( Bf[gp] * dTkdz[l][gp] + Tk[l][gp] * dBfdz[gp] );
		    		curlY = ( Bf[gp] * dTkdz[l][gp] + Tk[l][gp] * dBfdz[gp] ) - ( Bf[gp] * dTkdx[l][gp] + Tk[l][gp] * dBfdx[gp] );
		    		curlZ = ( Bf[gp] * dTkdx[l][gp] + Tk[l][gp] * dBfdx[gp] ) - ( Bf[gp] * dTkdy[l][gp] + Tk[l][gp] * dBfdy[gp] );
		    		
		    		norm_curlTkBk_l += W[gp] * ( curlX*curlX + curlY*curlY + curlZ*curlZ );
		    	}
		    
		    	norm_curlTkBk += norm_curlTkBk_l;
		    }
		    
		    // ( u, curl TkBk) * ( v, curl TkBk)
		    for ( l=0; l<TkSize; l++ ) 
		    {
		        Matrix< std::complex<double> > Nj_curlTBl( 1, ctcNumDofs );
		    	
		    	for ( j=0; j<mNumNodes; j++ ) 
		        {
		    		double Nj_curlTBl_X = 0.00;
		    		double Nj_curlTBl_Y = 0.00;
		    		double Nj_curlTBl_Z = 0.00;
		    		
		    		for ( gp=0; gp<nGaussPoints; gp++ ) 
		    		{	
		        		Nj_curlTBl_X += W[gp] * N[j][gp] * ( ( Bf[gp] * dTkdy[l][gp] + Tk[l][gp] * dBfdy[gp] ) - ( Bf[gp] * dTkdz[l][gp] + Tk[l][gp] * dBfdz[gp] ) );
		    			Nj_curlTBl_Y += W[gp] * N[j][gp] * ( ( Bf[gp] * dTkdz[l][gp] + Tk[l][gp] * dBfdz[gp] ) - ( Bf[gp] * dTkdx[l][gp] + Tk[l][gp] * dBfdx[gp] ) );
		    			Nj_curlTBl_Z += W[gp] * N[j][gp] * ( ( Bf[gp] * dTkdx[l][gp] + Tk[l][gp] * dBfdx[gp] ) - ( Bf[gp] * dTkdy[l][gp] + Tk[l][gp] * dBfdy[gp] ) );
		    		}
		    
		    		Nj_curlTBl[ 0 ][ j               ] = Nj_curlTBl_X;
		    		Nj_curlTBl[ 0 ][ j + mNumNodes   ] = Nj_curlTBl_Y;
		    		Nj_curlTBl[ 0 ][ j + mNumNodes*2 ] = Nj_curlTBl_Z;
		    	}
		    
		    	Matrix< std::complex<double> > Nj_curlTB      = Nj_curlTBl;
		    	Matrix< std::complex<double> > Nj_curlTB_tran = Nj_curlTBl.Transpose();
		    
		    	Face_nxn_Sh_cur += ( Nj_curlTB_tran * Nj_curlTB );
		    }
		    
		    Face_nxn_Sh_cur /= norm_curlTkBk;

			nxn_Sh_cur += Face_nxn_Sh_cur;
		}
	}

    //**************************************************************************************************
    //* - Adds a stabilization term S_h = 
    //**************************************************************************************************	
	void VolumeElement_1bb::Add_CurGBk_Stabilization( Matrix< std::complex<double> >& Sh_cur ) 
    { 
		Sh_cur.Resize( mNumDofs, mNumDofs, mZero );

	    // Loop counters
		int i, j, l, gp;

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

		int nGaussPoints = GaussPoints3D_Order10( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Derivatives
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Basis functions T
		Matrix<double> Tk;
        Matrix<double> dTkdx, dTkdy, dTkdz;

		//////////////////////////////////////////////////////////////////////////

		//int TkSize = 4; 
		//Lagrange3D_Ni_1st        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_1st( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		//int TkSize = 10; 
		//Lagrange3D_Ni_2nd        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_2nd( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

        //int TkSize = 20; 
		//Lagrange3D_Ni_3th        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_3th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		int TkSize = 35; 
		Lagrange3D_Ni_4th        ( Tk,                  cX, cY, cZ );
		LagrangeDerivatives3D_4th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		//////////////////////////////////////////////////////////////////////////

		//// Norm grad
		//double gradX, gradY, gradZ;
		//double norm_gradTkBk = 0.00;
		//
		//for ( l=0; l<TkSize; l++ ) 
		//{
		//	double norm_gradTkBk_l = 0.00;

		//	for ( gp=0; gp<nGaussPoints; gp++ ) 
		//	{	
		//		gradX = N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp];
		//		gradY = N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp];
		//		gradZ = N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp];
		//		
		//		norm_gradTkBk_l += W[gp] * ( gradX*gradX + gradY*gradY + gradZ*gradZ );
		//	}

		//	norm_gradTkBk += norm_gradTkBk_l;
		//}
		//
		//// ( curl u, TkBk) * ( curl v, TkBk)
		//for ( l=0; l<TkSize; l++ ) 
		//{
		//    Matrix< std::complex<double> > curNj_TBl( 3, mNumDofs );
		//	
		//	for ( j=0; j<mNumNodes; j++ ) 
		//    {
		//		double dNjdX_TBl = 0.00;
		//		double dNjdY_TBl = 0.00;
		//		double dNjdZ_TBl = 0.00;
		//		
		//		for ( gp=0; gp<nGaussPoints; gp++ ) 
		//		{	
		//    		dNjdX_TBl += W[gp] * dNdx[j][gp] * N[16][gp] * Tk[l][gp] ;
		//			dNjdY_TBl += W[gp] * dNdy[j][gp] * N[16][gp] * Tk[l][gp] ;
		//			dNjdZ_TBl += W[gp] * dNdz[j][gp] * N[16][gp] * Tk[l][gp] ;
		//		}

		//		curNj_TBl[ 0 ][ j               ] =      mZero;
		//		curNj_TBl[ 0 ][ j + mNumNodes   ] = -dNjdZ_TBl;
		//		curNj_TBl[ 0 ][ j + mNumNodes*2 ] = +dNjdY_TBl;
		//				   	 
		//		curNj_TBl[ 1 ][ j               ] = +dNjdZ_TBl;
		//		curNj_TBl[ 1 ][ j + mNumNodes   ] =      mZero;
		//		curNj_TBl[ 1 ][ j + mNumNodes*2 ] = -dNjdX_TBl;
		//				   	 
		//		curNj_TBl[ 2 ][ j               ] = -dNjdY_TBl;
		//		curNj_TBl[ 2 ][ j + mNumNodes   ] = +dNjdX_TBl;
		//		curNj_TBl[ 2 ][ j + mNumNodes*2 ] =      mZero;
		//	}

		//	Matrix< std::complex<double> > curNj_TB      = curNj_TBl;
		//	Matrix< std::complex<double> > curNj_TB_tran = curNj_TBl.Transpose();

		//	Sh_cur += ( curNj_TB_tran * curNj_TB );
		//}

		//Sh_cur /= ( cMu * norm_gradTkBk );

		////////////////////////////////////////////////////////////////////////////////////////

		// Norm curl
		double curlX, curlY, curlZ;
		double norm_curlTkBk = 0.00;
		
		for ( l=0; l<TkSize; l++ ) 
		{
			double norm_curlTkBk_l = 0.00;

			for ( gp=0; gp<nGaussPoints; gp++ ) 
			{	
				curlX = ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] ) - ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] );
				curlY = ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] ) - ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] );
				curlZ = ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] ) - ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] );
				
				norm_curlTkBk_l += W[gp] * ( curlX*curlX + curlY*curlY + curlZ*curlZ );
			}

			norm_curlTkBk += norm_curlTkBk_l;
		}
		
		// ( u, curl TkBk) * ( v, curl TkBk)
		for ( l=0; l<TkSize; l++ ) 
		{
		    Matrix< std::complex<double> > Nj_curlTBl( 1, mNumDofs );
			
			for ( j=0; j<mNumNodes; j++ ) 
		    {
				double Nj_curlTBl_X = 0.00;
				double Nj_curlTBl_Y = 0.00;
				double Nj_curlTBl_Z = 0.00;
				
				for ( gp=0; gp<nGaussPoints; gp++ ) 
				{	
		    		Nj_curlTBl_X += W[gp] * N[j][gp] * ( ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] ) - ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] ) );
					Nj_curlTBl_Y += W[gp] * N[j][gp] * ( ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] ) - ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] ) );
					Nj_curlTBl_Z += W[gp] * N[j][gp] * ( ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] ) - ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] ) );
				}

				Nj_curlTBl[ 0 ][ j               ] = Nj_curlTBl_X;
				Nj_curlTBl[ 0 ][ j + mNumNodes   ] = Nj_curlTBl_Y;
				Nj_curlTBl[ 0 ][ j + mNumNodes*2 ] = Nj_curlTBl_Z;
			}

			Matrix< std::complex<double> > Nj_curlTB      = Nj_curlTBl;
			Matrix< std::complex<double> > Nj_curlTB_tran = Nj_curlTBl.Transpose();

			Sh_cur += ( Nj_curlTB_tran * Nj_curlTB );
		}

		Sh_cur /= ( cMu * norm_curlTkBk );
		Sh_cur *= ( std::exp(- mFreq * mFreq * cEp ) );
	}

	//**************************************************************************************************
    //* - Adds a stabilization term S_h = 
    //**************************************************************************************************	
    void VolumeElement_1bb::Add_DivGBk_Stabilization( Matrix< std::complex<double> >& Sh_div ) 
    {  
		Sh_div.Resize( mNumDofs, mNumDofs, mZero );
		
	    // Loop counters
		int i, j, l, gp;

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

		int nGaussPoints = GaussPoints3D_Order10( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Derivatives
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Basis functions T
		Matrix<double> Tk;
        Matrix<double> dTkdx, dTkdy, dTkdz;

		//////////////////////////////////////////////////////////////////////////

		//int TkSize = 4; 
		//Lagrange3D_Ni_1st        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_1st( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		//int TkSize = 10; 
		//Lagrange3D_Ni_2nd        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_2nd( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

        //int TkSize = 20; 
		//Lagrange3D_Ni_3th        ( Tk,                  cX, cY, cZ );
		//LagrangeDerivatives3D_3th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		int TkSize = 35; 
		Lagrange3D_Ni_4th        ( Tk,                  cX, cY, cZ );
		LagrangeDerivatives3D_4th( dTkdx, dTkdy, dTkdz, cX, cY, cZ );

		/////////////////////////////////////////////////////////////////////////////

		// Norm
		double gradX, gradY, gradZ;
		double norm_gradTkBk = 0.00;
		
		for ( l=0; l<TkSize; l++ ) 
		{
			double norm_gradTkBk_l = 0.00;

			for ( gp=0; gp<nGaussPoints; gp++ ) 
			{	
				gradX = N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp];
				gradY = N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp];
				gradZ = N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp];
				
				norm_gradTkBk_l += W[gp] * ( gradX*gradX + gradY*gradY + gradZ*gradZ );
			}

			norm_gradTkBk += norm_gradTkBk_l;
		}
		
		// ( u, gradTB ) * ( u, gradTB )
		for ( l=0; l<TkSize; l++ ) 
		{
		    Matrix< std::complex<double> > Nj_gradTBl( 1, mNumDofs );
			
			for ( j=0; j<mNumNodes; j++ ) 
		    {
				double Nj_gradTB_X = 0.00;
				double Nj_gradTB_Y = 0.00;
				double Nj_gradTB_Z = 0.00;
				
				for ( gp=0; gp<nGaussPoints; gp++ ) 
				{	
		    		Nj_gradTB_X += W[gp] * N[j][gp] * ( N[16][gp] * dTkdx[l][gp] + Tk[l][gp] * dNdx[16][gp] );
					Nj_gradTB_Y += W[gp] * N[j][gp] * ( N[16][gp] * dTkdy[l][gp] + Tk[l][gp] * dNdy[16][gp] );
					Nj_gradTB_Z += W[gp] * N[j][gp] * ( N[16][gp] * dTkdz[l][gp] + Tk[l][gp] * dNdz[16][gp] );
				}

				Nj_gradTBl[0][ j               ] = Nj_gradTB_X;
				Nj_gradTBl[0][ j + mNumNodes   ] = Nj_gradTB_Y;
				Nj_gradTBl[0][ j + mNumNodes*2 ] = Nj_gradTB_Z;
			}
			
			Matrix< std::complex<double> > Nj_gradTB      = Nj_gradTBl;
			Matrix< std::complex<double> > Nj_gradTB_tran = Nj_gradTBl.Transpose();

			Sh_div += ( Nj_gradTB_tran * Nj_gradTB );
		}

		Sh_div /= ( cMu * norm_gradTkBk );
		Sh_div *= ( std::exp(- mFreq * mFreq * cEp ) );
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
		StiffMatrix.Resize( ctcNumDofs, ctcNumDofs, mZero );
		
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

		//std::complex<double> cEp( mu_real , mu_imag );
        //std::complex<double> cMu( eps_real, eps_imag + ( sigma / mFreq ) );

		// Gaussian derivation weights
		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order8( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Derivatives
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Diameter of the element
		double hk = Calculate_Element_Diameter();
		
		std::complex<double> cteCurlDiv = ( hk * hk ) / cMu;
		//std::complex<double> cteCurlDiv = ( 1.0 / cMu );

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

				StiffMatrix[ i                 ][ j                 ] += cteCurlDiv * intCurlDiv;
				StiffMatrix[ i + ctcNumNodes   ][ j + ctcNumNodes   ] += cteCurlDiv * intCurlDiv;
				StiffMatrix[ i + ctcNumNodes*2 ][ j + ctcNumNodes*2 ] += cteCurlDiv * intCurlDiv;
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
	            
				StiffMatrix[ i ][ j + ctcNumNodes ] += cteCurlDiv * intCurlDiv;
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
	            
				StiffMatrix[ i ][ j + ctcNumNodes*2 ] += cteCurlDiv * intCurlDiv;
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
	            
				StiffMatrix[ i + ctcNumNodes ][ j + ctcNumNodes*2 ] += cteCurlDiv * intCurlDiv;
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
		
		//return;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
			Vector<int> nf = *fs_it;

			//if ( ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) && ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'R' ) ) continue;
			if ( ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) ) continue;

			// Loop indexes
		    int i, j, gp;
		    
            // Gauss points and weights
		    Vector<double> cXs; 
		    Vector<double> cYs;
		    Vector<double>  Ws;
		    
		    int nGaussPoints = GaussPoints2D_Order10( cXs, cYs, Ws );
		    
            // Jacobian
		    double detJs = 2.00 * Calculate_Area( nf );
		    
		    for( gp=0; gp<nGaussPoints; gp++ ) Ws[gp] *= detJs; 
		    
            // Lagrangian base on Gauss points
		    Matrix<double> Ns;
		    
		    Lagrange2D_Ni_2bb( Ns, cXs, cYs );

		    Vector<double> extN; Calculate_ExtNormal( extN, nf );

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

			std::complex<double> cteNxN =  ( hF / cMu );

			for( i=0; i<6; i++ )
			{
				for( j=0; j<6; j++ )
				{
		       		std::complex<double> Sintg_NsiNsj = 0.00;
		           
					for( gp=0; gp<nGaussPoints; gp++ ) 
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

		return;

		//for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		//{
		//	Vector<int> nf = *fs_it;

		//	//if ( ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) && ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'R' ) ) continue;
		//	if ( ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'L' ) ) continue;

		//	// Loop indexes
		//    int i, j, gp;
		//    
  //          // Gauss points and weights
		//    Vector<double> cXs; 
		//    Vector<double> cYs;
		//    Vector<double>  Ws;
		//    
		//    int nGaussPoints = GaussPoints2D_Order10( cXs, cYs, Ws );
		//    
  //          // Jacobian
		//    double detJs = 2.00 * Calculate_Area( nf );
		//    
		//    for( gp=0; gp<nGaussPoints; gp++ ) Ws[gp] *= detJs; 
		//    
  //          // Lagrangian base on Gauss points
		//    Matrix<double> Ns;
		//    
		//    Lagrange2D_Ni_2bb( Ns, cXs, cYs );

		//    Vector<double> extN; Calculate_ExtNormal( extN, nf );

		//	double X1 = mNodes[nf[0]]->X(), Y1 = mNodes[nf[0]]->Y(), Z1 = mNodes[nf[0]]->Z();
		//    double X2 = mNodes[nf[1]]->X(), Y2 = mNodes[nf[1]]->Y(), Z2 = mNodes[nf[1]]->Z();
		//    double X3 = mNodes[nf[2]]->X(), Y3 = mNodes[nf[2]]->Y(), Z3 = mNodes[nf[2]]->Z();
		//    
		//    // Element edge lenghts 
		//    double Lenght_12 = sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) + (Z1-Z2)*(Z1-Z2) );
		//    double Lenght_13 = sqrt( (X1-X3)*(X1-X3) + (Y1-Y3)*(Y1-Y3) + (Z1-Z3)*(Z1-Z3) );
		//    double Lenght_23 = sqrt( (X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3) + (Z2-Z3)*(Z2-Z3) );
		//    
  //          // Maximum lenght 
		//    double hF = Lenght_12;
		//    
		//    if ( Lenght_13 > hF ) hF = Lenght_13;
		//    if ( Lenght_23 > hF ) hF = Lenght_23;

		//	std::complex<double> cteNxN =  ( hF / cMu );

		//	std::complex< double > cEpR( 5.00*mEo, 1.00e6 / mFreq );
		//	std::complex< double > cEpL( 1.00*mEo, 0.00   / mFreq );

		//	std::complex<double> rmat = ( cEpL / cEpR ) - 1.00;
		//	std::complex<double> gamm = ( cEpR / cEpL );
		//	std::complex<double> invL = ( 1.00 / cEpL );

  //          Vector<double> n; Calculate_ExtNormal( n, nf );

		//	Matrix< std::complex<double> > CtCM; CtCM.Resize(3,3);

		//	CtCM[0][0] = (n[0]*n[0]) * rmat + 1.00; CtCM[0][1] = (n[0]*n[1]) * rmat       ; CtCM[0][2] = (n[0]*n[2]) * rmat;
		//	CtCM[1][0] = (n[1]*n[0]) * rmat       ; CtCM[1][1] = (n[1]*n[1]) * rmat + 1.00; CtCM[1][2] = (n[1]*n[2]) * rmat;
		//	CtCM[2][0] = (n[2]*n[0]) * rmat       ; CtCM[2][1] = (n[2]*n[1]) * rmat       ; CtCM[2][2] = (n[2]*n[2]) * rmat + 1.00;

		//	CtCM *= -gamm;

		//	CtCM[0][0] += 1.00;
		//	CtCM[1][1] += 1.00;
		//	CtCM[2][2] += 1.00;

		//	//Matrix< std::complex<double> > n_ctcM; n_ctcM.Resize(1,3);
		//	//
		//	//n_ctcM[0][0] = n[0]*CtCM[0][0] + n[1]*CtCM[1][0] + n[2]*CtCM[2][0];
		//	//n_ctcM[0][1] = n[0]*CtCM[0][1] + n[1]*CtCM[1][1] + n[2]*CtCM[2][1];
		//	//n_ctcM[0][2] = n[0]*CtCM[0][2] + n[1]*CtCM[1][2] + n[2]*CtCM[2][2];

		//	//n_ctcM *= gamm;

		//	Matrix< std::complex<double> > CtCMT; CtCMT.Resize(3,3, mZero);

		//	CtCMT += CtCM.Transpose() * CtCM;

		//	for( i=0; i<6; i++ )
		//	{
		//		for( j=0; j<6; j++ )
		//		{
		//       		std::complex<double> Sintg_NsiNsj = 0.00;
		//           
		//			for( gp=0; gp<nGaussPoints; gp++ ) 
		//			{
		//				Sintg_NsiNsj += Ws[gp] * Ns[i][gp] * Ns[j][gp]; 
		//			}

		//			Sintg_NsiNsj *= cteNxN;

		//			StiffMatrix[ nf[i]      ][ nf[j]      ] += ( CtCMT[0][0] ) * Sintg_NsiNsj; 
		//			StiffMatrix[ nf[i] +  6 ][ nf[j]      ] += ( CtCMT[1][0] ) * Sintg_NsiNsj; 
		//			StiffMatrix[ nf[i] + 12 ][ nf[j]      ] += ( CtCMT[2][0] ) * Sintg_NsiNsj; 

		//			StiffMatrix[ nf[i]      ][ nf[j] +  6 ] += ( CtCMT[0][1] ) * Sintg_NsiNsj; 
		//			StiffMatrix[ nf[i] +  6 ][ nf[j] +  6 ] += ( CtCMT[1][1] ) * Sintg_NsiNsj; 
		//			StiffMatrix[ nf[i] + 12 ][ nf[j] +  6 ] += ( CtCMT[2][1] ) * Sintg_NsiNsj; 

		//			StiffMatrix[ nf[i]      ][ nf[j] + 12 ] += ( CtCMT[0][2] ) * Sintg_NsiNsj; 
		//			StiffMatrix[ nf[i] +  6 ][ nf[j] + 12 ] += ( CtCMT[1][2] ) * Sintg_NsiNsj; 
		//			StiffMatrix[ nf[i] + 12 ][ nf[j] + 12 ] += ( CtCMT[2][2] ) * Sintg_NsiNsj; 
		//		}
		//	}
		//}
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

		int nGaussPoints = GaussPoints3D_Order8( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

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
		//Matrix< std::complex<double> > iAc_x_Ccur = iA_cur * C_cur;
		//Matrix< std::complex<double> > Ccur_tran  = C_cur.Transpose();
		//StiffMatrix += ( Ccur_tran * iAc_x_Ccur );

		Matrix< std::complex<double> > iAd_x_Ddiv = iA_div * D_div;
		Matrix< std::complex<double> > Ddiv_tran  = D_div.Transpose();
		StiffMatrix += ( Ddiv_tran * iAd_x_Ddiv );

		//StiffMatrix -= N_nij;
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
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;
		
		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) continue;

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
		    
		    Lagrange2D_Ni_2bb( N, cX, cY );

		    Vector<double> extN; Calculate_ExtNormal( extN, nf );

			for( i=0; i<3; i++ )
		    {
		        for( j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_NiNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}
   		
                    C_cur[ nf[i]     ][ nf[j] + ctcNumNodes   ] += (extN[2] * Sintg_NiNj); 
                    C_cur[ nf[i]     ][ nf[j] + ctcNumNodes*2 ] -= (extN[1] * Sintg_NiNj); 
		        					    	 					 
		        	C_cur[ nf[i] + 4 ][ nf[j]                 ] -= (extN[2] * Sintg_NiNj); 
                    C_cur[ nf[i] + 4 ][ nf[j] + ctcNumNodes*2 ] += (extN[0] * Sintg_NiNj);
		        					 	 						
		        	C_cur[ nf[i] + 8 ][ nf[j]                 ] += (extN[1] * Sintg_NiNj); 
                    C_cur[ nf[i] + 8 ][ nf[j] + ctcNumNodes   ] -= (extN[0] * Sintg_NiNj); 
		        }
		    }
		}

		/////////////////////////////////////////////////////////////////////////////////////////////
		/*
		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'L' ) continue;

			Vector<int> cf; cf.resize( 6 );
            
			cf[0] = 17; cf[1] = 18; cf[2] = 19; cf[3] = 20; cf[4] = 21; cf[5] = 22;

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
		    
		    Lagrange2D_Ni_2bb( N, cX, cY );

		    Vector<double> extN; Calculate_ExtNormal( extN, nf );

			for( i=0; i<3; i++ )
		    {
		        for( j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_NiNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}
   		
                    C_cur[ nf[i]     ][ nf[j] + ctcNumNodes   ] += (extN[2] * Sintg_NiNj); 
                    C_cur[ nf[i]     ][ nf[j] + ctcNumNodes*2 ] -= (extN[1] * Sintg_NiNj); 
		        					    	 					 
		        	C_cur[ nf[i] + 4 ][ nf[j]                 ] -= (extN[2] * Sintg_NiNj); 
                    C_cur[ nf[i] + 4 ][ nf[j] + ctcNumNodes*2 ] += (extN[0] * Sintg_NiNj);
		        					 	 						
		        	C_cur[ nf[i] + 8 ][ nf[j]                 ] += (extN[1] * Sintg_NiNj); 
                    C_cur[ nf[i] + 8 ][ nf[j] + ctcNumNodes   ] -= (extN[0] * Sintg_NiNj); 


					C_cur[ cf[i]     ][ cf[j] + ctcNumNodes   ] -= (extN[2] * Sintg_NiNj); 
                    C_cur[ cf[i]     ][ cf[j] + ctcNumNodes*2 ] += (extN[1] * Sintg_NiNj); 
		        							 
		        	C_cur[ cf[i] + 4 ][ cf[j]                 ] += (extN[2] * Sintg_NiNj); 
                    C_cur[ cf[i] + 4 ][ cf[j] + ctcNumNodes*2 ] -= (extN[0] * Sintg_NiNj);
		        						
		        	C_cur[ cf[i] + 8 ][ cf[j]                 ] -= (extN[1] * Sintg_NiNj); 
                    C_cur[ cf[i] + 8 ][ cf[j] + ctcNumNodes   ] += (extN[0] * Sintg_NiNj); 
		        }
		    }
		}
		*/
	}

	//******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::Add_C_NxN_Matrix_ctc( Matrix< std::complex<double> >& C_cur ) 
	{
		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;
		
		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) continue;

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
		    
		    Lagrange2D_Ni_2bb( N, cX, cY );


			// Lagrangian base on Gauss points
		    Matrix<double> Xij;
		    
		    Lagrange2D_Ni_1st( Xij, cX, cY );

			int numXSnodes = 3;
			int numXVnodes = 4;


		    Vector<double> extN; Calculate_ExtNormal( extN, nf );

			for( i=0; i<numXSnodes; i++ )
		    {
		        for( j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_XiNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_XiNj += W[gp] * Xij[i][gp] * N[j][gp]; 
					}
   		
                    C_cur[ nf[i]                ][ nf[j] + ctcNumNodes   ] += (extN[2] * Sintg_XiNj); 
                    C_cur[ nf[i]                ][ nf[j] + ctcNumNodes*2 ] -= (extN[1] * Sintg_XiNj); 
		        					    	 					 
		        	C_cur[ nf[i] + numXVnodes   ][ nf[j]                 ] -= (extN[2] * Sintg_XiNj); 
                    C_cur[ nf[i] + numXVnodes   ][ nf[j] + ctcNumNodes*2 ] += (extN[0] * Sintg_XiNj);
		        					 	 						
		        	C_cur[ nf[i] + numXVnodes*2 ][ nf[j]                 ] += (extN[1] * Sintg_XiNj); 
                    C_cur[ nf[i] + numXVnodes*2 ][ nf[j] + ctcNumNodes   ] -= (extN[0] * Sintg_XiNj); 
		        }
		    }
		}
	}

	//******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::Add_D_NxN_Matrix( Matrix< std::complex<double> >& D_div ) 
	{
		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;
		
		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			// 'E' = PEC normal | 'D' = Plasma Sheath | 'H' = PMC | 'T' = TEPMC
			if ( ((*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'H') && ((*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'T') ) continue;
			//if ( ((*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'R') ) continue;

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
		    
		    Lagrange2D_Ni_2bb( N, cX, cY );
		    
		    Vector<double> extN; Calculate_ExtNormal( extN, nf );
		    
		    for( i=0; i<3; i++ )
		    {
		        for( j=0; j<6; j++ )
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
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::Add_CtC_D_NxN_Matrix( Matrix< std::complex<double> >& D_div ) 
	{
		/*
		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;
		
		Vector< Vector<int> >::iterator fs_it;

		for( int iface=0; iface<4; iface++ )
		{
			Vector<int> nf = Faces[iface];

			if ( ((*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'C') ) continue;
			

			Vector<int> cf(6);// =  mCTC_DuplicateFace[iface];

			for (int i=0; i<6; i++) cf[i] = 17+i;

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
		    
		    Lagrange2D_Ni_2bb( N, cX, cY );
		    
		    Vector<double> extN; Calculate_ExtNormal( extN, nf );

			double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
            double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
            double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
			
			std::complex< double > ThisE_cEp( eps_real, eps_imag + ( sigma / mFreq ) );

			std::complex< double > Other_cEp = mCtC_Properties[iface];

			std::complex< double > cEpO_div_cEpT = Other_cEp / ThisE_cEp;

		    for( i=0; i<3; i++ )
		    {
		        for( j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_NiNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}

					D_div[ nf[i] ][ nf[j]                 ] += (extN[0] * Sintg_NiNj);
		    	    D_div[ nf[i] ][ nf[j] + ctcNumNodes   ] += (extN[1] * Sintg_NiNj);
		    	    D_div[ nf[i] ][ nf[j] + ctcNumNodes*2 ] += (extN[2] * Sintg_NiNj);

					D_div[ nf[i] ][ cf[j]                 ] -= (cEpO_div_cEpT * extN[0] * Sintg_NiNj);
		    	    D_div[ nf[i] ][ cf[j] + ctcNumNodes   ] -= (cEpO_div_cEpT * extN[1] * Sintg_NiNj);
		    	    D_div[ nf[i] ][ cf[j] + ctcNumNodes*2 ] -= (cEpO_div_cEpT * extN[2] * Sintg_NiNj);
		        }
		    }
		}
		*/

		//return;

		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;
		
		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			// 'E' = PEC normal | 'D' = Plasma Sheath | 'H' = PMC | 'T' = TEPMC
			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'L' ) continue;

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
		    
		    Lagrange2D_Ni_2bb( N, cX, cY );
		    
		    std::complex< double > cEpR( 5.00*mEo, 1.00e6 / mFreq );
			std::complex< double > cEpL( 1.00*mEo, 0.00   / mFreq );

			std::complex<double> rmat = ( cEpL / cEpR ) - 1.00;
			std::complex<double> gamm = ( cEpR / cEpL );
			std::complex<double> invL = ( 1.00 / cEpL );

            Vector<double> n; Calculate_ExtNormal( n, nf );

			std::vector<double> t(3), b(3); TangencialCoordinates(n,t,b);

			Matrix< std::complex<double> > CtCM; CtCM.Resize(3,3);

			CtCM[0][0] = (n[0]*n[0]) * rmat + 1.00; CtCM[0][1] = (n[0]*n[1]) * rmat       ; CtCM[0][2] = (n[0]*n[2]) * rmat;
			CtCM[1][0] = (n[1]*n[0]) * rmat       ; CtCM[1][1] = (n[1]*n[1]) * rmat + 1.00; CtCM[1][2] = (n[1]*n[2]) * rmat;
			CtCM[2][0] = (n[2]*n[0]) * rmat       ; CtCM[2][1] = (n[2]*n[1]) * rmat       ; CtCM[2][2] = (n[2]*n[2]) * rmat + 1.00;

			//CtCM *= gamm;

			//Matrix< std::complex<double> > sCtCM; sCtCM.Resize(3,3);

			//sCtCM[0][0] = 1.0 - CtCM[0][0]; sCtCM[0][1] =     - CtCM[0][1]; sCtCM[0][2] =     - CtCM[0][2];
			//sCtCM[1][0] =     - CtCM[1][0]; sCtCM[1][1] = 1.0 - CtCM[1][1]; sCtCM[1][2] =     - CtCM[1][2];
			//sCtCM[2][0] =     - CtCM[2][0]; sCtCM[2][1] =     - CtCM[2][1]; sCtCM[2][2] = 1.0 - CtCM[2][2];

			//sCtCM *= gamm;

			Matrix< std::complex<double> > n_ctcM; n_ctcM.Resize(1,3);
			
			n_ctcM[0][0] = n[0]*CtCM[0][0] + n[1]*CtCM[1][0] + n[2]*CtCM[2][0];
			n_ctcM[0][1] = n[0]*CtCM[0][1] + n[1]*CtCM[1][1] + n[2]*CtCM[2][1];
			n_ctcM[0][2] = n[0]*CtCM[0][2] + n[1]*CtCM[1][2] + n[2]*CtCM[2][2];

			n_ctcM *= gamm;

			Matrix< std::complex<double> > t_ctcM; t_ctcM.Resize(1,3);
			
			t_ctcM[0][0] = t[0]*CtCM[0][0] + t[1]*CtCM[1][0] + t[2]*CtCM[2][0];
			t_ctcM[0][1] = t[0]*CtCM[0][1] + t[1]*CtCM[1][1] + t[2]*CtCM[2][1];
			t_ctcM[0][2] = t[0]*CtCM[0][2] + t[1]*CtCM[1][2] + t[2]*CtCM[2][2];

			t_ctcM *= invL;

			Matrix< std::complex<double> > b_ctcM; b_ctcM.Resize(1,3);
			
			b_ctcM[0][0] = b[0]*CtCM[0][0] + b[1]*CtCM[1][0] + b[2]*CtCM[2][0];
			b_ctcM[0][1] = b[0]*CtCM[0][1] + b[1]*CtCM[1][1] + b[2]*CtCM[2][1];
			b_ctcM[0][2] = b[0]*CtCM[0][2] + b[1]*CtCM[1][2] + b[2]*CtCM[2][2];

			b_ctcM *= invL;

		    for( i=0; i<3; i++ )
		    {
		        for( j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_NiNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}

					// - n * UL
					D_div[ nf[i] ][ nf[j]               ] -= ( n[0] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes   ] -= ( n[1] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] -= ( n[2] * Sintg_NiNj );
		            
			        // + n * (ER/EL) * [ ] * UL
			        D_div[ nf[i] ][ nf[j]               ] += ( n_ctcM[0][0] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes   ] += ( n_ctcM[0][1] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] += ( n_ctcM[0][2] * Sintg_NiNj );

					/////////////////////////////////////////////////////////////////////////////

					// - n * UL
					D_div[ nf[i] ][ nf[j]               ] -= ( invL * t[0] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes   ] -= ( invL * t[1] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] -= ( invL * t[2] * Sintg_NiNj );

					// + n * (ER/EL) * [ ] * UL
					D_div[ nf[i] ][ nf[j]               ] += ( t_ctcM[0][0] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes   ] += ( t_ctcM[0][1] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] += ( t_ctcM[0][2] * Sintg_NiNj );


					// - n * UL
					D_div[ nf[i] ][ nf[j]               ] -= ( invL * b[0] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes   ] -= ( invL * b[1] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] -= ( invL * b[2] * Sintg_NiNj );

					// + n * (ER/EL) * [ ] * UL
					D_div[ nf[i] ][ nf[j]               ] += ( b_ctcM[0][0] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes   ] += ( b_ctcM[0][1] * Sintg_NiNj );
		            D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] += ( b_ctcM[0][2] * Sintg_NiNj );

					/////////////////////////////////////////////////////////////////////////////

					//D_div[ nf[i] ][ nf[j]               ] -= ( sCtCM[0][0] * Sintg_NiNj );
		   // 	    D_div[ nf[i] ][ nf[j] + mNumNodes   ] -= ( sCtCM[0][1] * Sintg_NiNj );
		   // 	    D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] -= ( sCtCM[0][2] * Sintg_NiNj );
					//									  
					//D_div[ nf[i] ][ nf[j]               ] -= ( sCtCM[1][0] * Sintg_NiNj );
		   // 	    D_div[ nf[i] ][ nf[j] + mNumNodes   ] -= ( sCtCM[1][1] * Sintg_NiNj );
		   // 	    D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] -= ( sCtCM[1][2] * Sintg_NiNj );
					//									  
					//D_div[ nf[i] ][ nf[j]               ] -= ( sCtCM[2][0] * Sintg_NiNj );
		   // 	    D_div[ nf[i] ][ nf[j] + mNumNodes   ] -= ( sCtCM[2][1] * Sintg_NiNj );
		   // 	    D_div[ nf[i] ][ nf[j] + mNumNodes*2 ] -= ( sCtCM[2][2] * Sintg_NiNj );

		        }
		    }
		}
	}

	//*************************************************************************************************
    //* - Local system of coordinates in PEC surface for n
    //*************************************************************************************************
    void VolumeElement_1bb::TangencialCoordinates( std::vector<double>& n, std::vector<double>& t, std::vector<double>& b )
    {
         double tnorm = sqrt(n[0]*n[0]+n[1]*n[1]);

         if (tnorm > 0.00)
         {
             t[0] = n[1]/tnorm;
             t[1] =-n[0]/tnorm;
             t[2] = 0.00;
         }
         else
         {
             t[0]  = n[2];
             t[1]  = 0.00;
             t[2]  = 0.00;
         }

         b[0] = n[1]*t[2] - n[2]*t[1];
         b[1] = n[2]*t[0] - n[0]*t[2];
         b[2] = n[0]*t[1] - n[1]*t[0];
    }

	//******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::GetStiffMatrix_LL2P_CtC( Matrix< std::complex<double> >& StiffMatrix ) 
	{
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( mNumDofs, mNumDofs, mZero );

		return;

		// Material properties
		double sigma    = (*mProperties)(IHL_ELECTRIC_CONDUCTIVITY);
        double eps_real = (*mProperties)(REAL_ELECTRIC_PERMITTIVITY) * mEo;
        double eps_imag = (*mProperties)(IMAG_ELECTRIC_PERMITTIVITY) * mEo;
        double mu_real  = (*mProperties)(REAL_MAGNETIC_PERMEABILITY) * mMo;
        double mu_imag  = (*mProperties)(IMAG_MAGNETIC_PERMEABILITY) * mMo;

        std::complex<double> cEp( eps_real, eps_imag + ( sigma / mFreq ) );
        std::complex<double> cMu( mu_real , mu_imag );

        // Faces
		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;

		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			// 'E' = PEC normal | 'D' = Plasma Sheath | 'H' = PMC | 'T' = TEPMC
			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'L' ) continue;

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
		    
		    Lagrange2D_Ni_2bb( N, cX, cY );
		    
		    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		    std::complex< double > cEpR( 5.00*mEo, 1.00 / mFreq );
			std::complex< double > cEpL( 1.00*mEo, 0.00 / mFreq );

			std::complex<double> rmat = ( cEpL / cEpR ) - 1.00;
			std::complex<double> gamm = ( cEpR / cEpL );

            Vector<double> n; Calculate_ExtNormal( n, nf );

			Matrix< std::complex<double> > CtCM; CtCM.Resize(3,3);

			CtCM[0][0] = (n[0]*n[0]) * rmat + 1.00; CtCM[0][1] = (n[0]*n[1]) * rmat       ; CtCM[0][2] = (n[0]*n[2]) * rmat       ;
			CtCM[1][0] = (n[1]*n[0]) * rmat       ; CtCM[1][1] = (n[1]*n[1]) * rmat + 1.00; CtCM[1][2] = (n[1]*n[2]) * rmat       ;
			CtCM[2][0] = (n[2]*n[0]) * rmat       ; CtCM[2][1] = (n[2]*n[1]) * rmat       ; CtCM[2][2] = (n[2]*n[2]) * rmat + 1.00;

			//CtCM *= gamm;

			Matrix< std::complex<double> > n_ctcM; n_ctcM.Resize(1,3);
			
			n_ctcM[0][0] = n[0]*CtCM[0][0] + n[1]*CtCM[1][0] + n[2]*CtCM[2][0];
			n_ctcM[0][1] = n[0]*CtCM[0][1] + n[1]*CtCM[1][1] + n[2]*CtCM[2][1];
			n_ctcM[0][2] = n[0]*CtCM[0][2] + n[1]*CtCM[1][2] + n[2]*CtCM[2][2];

			//n_ctcM *= gamm;

            //////////////////////////////////////////////////     Inv(T)*StiffnesMatrix*T     //////////////////////////////////////////////////////////////
			Matrix< std::complex<double> > CtC_surf; CtC_surf.Resize( 3, 18, mZero );

			for( i=0; i<3; i++ )
		    {
		        for( j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_NiNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
					}
   		
                    // - n * UL
					/*CtC_surf[ i ][ j      ] -= ( n[0] * Sintg_NiNj );
		            CtC_surf[ i ][ j + 6  ] -= ( n[1] * Sintg_NiNj );
		            CtC_surf[ i ][ j + 12 ] -= ( n[2] * Sintg_NiNj );*/
		            
			        // + n * (ER/EL) * [ ] * UL
			        CtC_surf[ i ][ j      ] += ( n_ctcM[0][0] * Sintg_NiNj );
		            CtC_surf[ i ][ j + 6  ] += ( n_ctcM[0][1] * Sintg_NiNj );
		            CtC_surf[ i ][ j + 12 ] += ( n_ctcM[0][2] * Sintg_NiNj );
		        }
		    }

			double area =  Calculate_Area( nf );

		    // Projection basis matrices
		    Matrix< std::complex<double> > iA_div; iA_div.Resize( 3, 3, mZero );
		    
		    for( i=0; i<3; i++ )
		    {
                for( j=0; j<3; j++ )
                {
            	     if ( i != j ) iA_div[ i ][ j ] = -3.0 / area;
            		 else          iA_div[ i ][ j ] = +9.0 / area;                          
                }									   
		    }
		    
		    // Scaling the projection matrices
		    iA_div /= cMu ;
		    
		    // Final assembly: [ CT * Ac^(-1) * C ] + [ DT * Ad^(-1) * D ] - [ Ni * Nj ]
		    Matrix< std::complex<double> > iAc_x_CtCs = iA_div * CtC_surf;
		    Matrix< std::complex<double> > Ctcs_tran  = CtC_surf.Transpose();
			Matrix< std::complex<double> > RJump      = Ctcs_tran * iAc_x_CtCs;

			for( i=0; i<6; i++ )
			{
				for( j=0; j<6; j++ )
				{
			        StiffMatrix[ nf[i]      ][ nf[j]      ] += RJump[i    ][j];
			        StiffMatrix[ nf[i] +  6 ][ nf[j]      ] += RJump[i + 6][j];
			        StiffMatrix[ nf[i] + 12 ][ nf[j]      ] += RJump[i +12][j];
			        
			        StiffMatrix[ nf[i]      ][ nf[j] +  6 ] += RJump[i    ][j + 6];
			        StiffMatrix[ nf[i] +  6 ][ nf[j] +  6 ]	+= RJump[i + 6][j + 6];
			        StiffMatrix[ nf[i] + 12 ][ nf[j] +  6 ]	+= RJump[i +12][j + 6];
			        
			        StiffMatrix[ nf[i]      ][ nf[j] + 12 ] += RJump[i     ][j + 12];
			        StiffMatrix[ nf[i] +  6 ][ nf[j] + 12 ]	+= RJump[i +  6][j + 12];
			        StiffMatrix[ nf[i] + 12 ][ nf[j] + 12 ]	+= RJump[i + 12][j + 12];
				}
			}
		}
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

		//std::complex<double> cEp( mu_real , mu_imag );
        //std::complex<double> cMu( eps_real, eps_imag + ( sigma / mFreq ) );
      
		// Gauss points and weights
		std::vector<double> cX, cY, cZ, W;

		int nGaussPoints = GaussPoints3D_Order9( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		/*double cEpr = cEp / mEo;
		
		for ( i=4; i<17; i++ ) 
		{
			i = 16;
			for( gp=0; gp<nGaussPoints; gp++ )
            {
                N   [16][gp] *= (eps_real*mEo);
				dNdx[16][gp] *= (eps_real*mEo);
				dNdy[16][gp] *= (eps_real*mEo);
				dNdz[16][gp] *= (eps_real*mEo);
            } 
		}*/

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
	
		// Surface integral
		Add_C_NxN_Matrix( C_cur ); 
		//Add_D_NxN_Matrix( D_div );
		//Add_CtC_D_NxN_Matrix( D_div ); 

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij; N_nij.Resize( mNumDofs, mNumDofs, mZero );
		
		std::complex<double> cteNiNj = mFreq * mFreq * cEp; 
		//std::complex<double> cteNiNj (0.0, - mFreq * sigma ); 

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
		
		/*
		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_jwe; N_jwe.Resize( 12, mNumDofs, mZero );
		
		for( i=0; i<4; i++ )
		{
            for( j=0; j<mNumNodes; j++ )
            {
                double int_NiNj = 0.00;

                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	int_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
                }
            
            	N_jwe[ i     ][ j               ] = int_NiNj; 
            	N_jwe[ i + 4 ][ j + mNumNodes   ] = int_NiNj; 
            	N_jwe[ i + 8 ][ j + mNumNodes*2 ] = int_NiNj; 
            }
		}
		*/
		
		// Projection basis matrices
		Matrix< std::complex<double> > iA_div; iA_div.Resize(  4,  4, mZero );
		Matrix< std::complex<double> > iA_cur; iA_cur.Resize( 12, 12, mZero );
		//Matrix< std::complex<double> > iA_nij; iA_nij.Resize( 12, 12, mZero );
		
		for( i=0; i<4; i++ )
		{
            for( j=0; j<4; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -4.00 / mVolume;
        		else          iA_div[ i ][ j ] = +16.0 / mVolume;                          
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                iA_cur[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];

				//iA_nij[ i     ][ j     ] = iA_div[ i ][ j ]; 
                //iA_nij[ i + 4 ][ j + 4 ] = iA_div[ i ][ j ];
                //iA_nij[ i + 8 ][ j + 8 ] = iA_div[ i ][ j ];
            }									   
		}
		
		// Scaling the projection matrices
		iA_div /= cMu ;
		iA_cur /= cMu ;
		//iA_nij *= cteNiNj;
		
		// Final assembly: [ CT * Ac^(-1) * C ] + [ DT * Ad^(-1) * D ] - [ Ni * Nj ]
		Matrix< std::complex<double> > iAc_x_Ccur = iA_cur * C_cur;
		Matrix< std::complex<double> > Ccur_tran  = C_cur.Transpose();
		StiffMatrix += ( Ccur_tran * iAc_x_Ccur );

		Matrix< std::complex<double> > iAd_x_Ddiv = iA_div * D_div;
		Matrix< std::complex<double> > Ddiv_tran  = D_div.Transpose();
		StiffMatrix += ( Ddiv_tran * iAd_x_Ddiv );

		//StiffMatrix -= N_jwe.Transpose() * ( iA_nij * N_jwe );
		StiffMatrix -= N_nij;
    }

    //******************************************************************************************************
    //* - Stiffness matrix with local L2 projection method ( div(v), q ) - ( curl(v), w )
    //******************************************************************************************************	
    void VolumeElement_1bb::GetStiffMatrix_LL2P_CurlDivE_CtC( Matrix< std::complex<double> >& StiffMatrix ) 
    { 
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( ctcNumDofs, ctcNumDofs, mZero );

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

		int nGaussPoints = GaussPoints3D_Order9( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< std::complex<double> > D_div; D_div.Resize(  4, ctcNumDofs, mZero );
		Matrix< std::complex<double> > C_cur; C_cur.Resize( 12, ctcNumDofs, mZero );

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
		    	
		    	D_div[ i     ][ j                 ] = +vintg_Ni_dNjdX; 
		    	D_div[ i     ][ j + ctcNumNodes   ] = +vintg_Ni_dNjdY; 
		    	D_div[ i     ][ j + ctcNumNodes*2 ] = +vintg_Ni_dNjdZ; 
		    			    
		    	C_cur[ i     ][ j                 ] =            0.00; 
                C_cur[ i     ][ j + ctcNumNodes   ] = -vintg_Ni_dNjdZ; 
                C_cur[ i     ][ j + ctcNumNodes*2 ] = +vintg_Ni_dNjdY; 
		    								    
		    	C_cur[ i + 4 ][ j                 ] = +vintg_Ni_dNjdZ; 
                C_cur[ i + 4 ][ j + ctcNumNodes   ] =            0.00; 
                C_cur[ i + 4 ][ j + ctcNumNodes*2 ] = -vintg_Ni_dNjdX; 
		    								    
		    	C_cur[ i + 8 ][ j                 ] = -vintg_Ni_dNjdY; 
                C_cur[ i + 8 ][ j + ctcNumNodes   ] = +vintg_Ni_dNjdX; 
                C_cur[ i + 8 ][ j + ctcNumNodes*2 ] =            0.00; 
		    }
		}
	
		// Surface integral
		//Add_C_NxN_Matrix( C_cur ); 
		//Add_D_NxN_Matrix( D_div );
		//Add_CtC_D_NxN_Matrix( D_div ); 

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij; N_nij.Resize( ctcNumDofs, ctcNumDofs, mZero );
		
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
            
            	N_nij[ i                 ][ j                 ] = cteNiNj * int_NiNj; 
            	N_nij[ i + ctcNumNodes   ][ j + ctcNumNodes   ] = cteNiNj * int_NiNj; 
            	N_nij[ i + ctcNumNodes*2 ][ j + ctcNumNodes*2 ] = cteNiNj * int_NiNj; 
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
		Matrix< std::complex<double> > iAc_x_Ccur = iA_cur * C_cur;
		Matrix< std::complex<double> > Ccur_tran  = C_cur.Transpose();
		StiffMatrix += ( Ccur_tran * iAc_x_Ccur );

		Matrix< std::complex<double> > iAd_x_Ddiv = iA_div * D_div;
		Matrix< std::complex<double> > Ddiv_tran  = D_div.Transpose();
		StiffMatrix += ( Ddiv_tran * iAd_x_Ddiv );

		StiffMatrix -= N_nij;
    }

/*
	//******************************************************************************************************
    //* - Stiffness matrix with local L2 projection method ( div(v), q ) - ( curl(v), w )
    //******************************************************************************************************	
    void VolumeElement_1bb::GetStiffMatrix_LL2P_CurlDivE_ctc( Matrix< std::complex<double> >& StiffMatrix ) 
    { 
		// Resizing local stiffness matrix 
        StiffMatrix.Resize( ctcNumDofs, ctcNumDofs, mZero );
		
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

		int nGaussPoints = GaussPoints3D_Order10( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		//----
		Matrix<double> Xij;

		Lagrange3D_Ni_1st( Xij, cX, cY, cZ );

		//int nPoints = cX.size();

		//Xij.Resize( 1, nPoints );
		//
		//for( int p=0; p<nPoints; p++ ) 
  //      {
  //          Xij[0][p] = 1.0;

		//}

		int numXnodes = 4;
		int numXdofs  = numXnodes*3;

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		// Building projection matrices (curl and div contributions)
		Matrix< std::complex<double> > D_div; D_div.Resize(  numXnodes, ctcNumDofs, mZero );
		Matrix< std::complex<double> > C_cur; C_cur.Resize(  numXdofs , ctcNumDofs, mZero );

		for( i=0; i<numXnodes; i++ )
		{
		    for( j=0; j<mNumNodes; j++ )
		    {
		    	double vintg_Xi_dNjdX = 0.00;
		        double vintg_Xi_dNjdY = 0.00;
		        double vintg_Xi_dNjdZ = 0.00;
		    
                for( gp=0; gp<nGaussPoints; gp++ )
                {
                	vintg_Xi_dNjdX += W[gp] * Xij[i][gp] * dNdx[j][gp]; 
		    		vintg_Xi_dNjdY += W[gp] * Xij[i][gp] * dNdy[j][gp];
		    		vintg_Xi_dNjdZ += W[gp] * Xij[i][gp] * dNdz[j][gp];
                }      
		    	
		    	D_div[ i ][ j                 ] = +vintg_Xi_dNjdX; 
		    	D_div[ i ][ j + ctcNumNodes   ] = +vintg_Xi_dNjdY; 
		    	D_div[ i ][ j + ctcNumNodes*2 ] = +vintg_Xi_dNjdZ; 
		    			  
		    	C_cur[ i               ][ j                 ] =           mZero; 
                C_cur[ i               ][ j + ctcNumNodes   ] = -vintg_Xi_dNjdZ; 
                C_cur[ i               ][ j + ctcNumNodes*2 ] = +vintg_Xi_dNjdY; 
		    								    
		    	C_cur[ i + numXnodes   ][ j                 ] = +vintg_Xi_dNjdZ; 
                C_cur[ i + numXnodes   ][ j + ctcNumNodes   ] =           mZero; 
                C_cur[ i + numXnodes   ][ j + ctcNumNodes*2 ] = -vintg_Xi_dNjdX; 
		    								    
		    	C_cur[ i + numXnodes*2 ][ j                 ] = -vintg_Xi_dNjdY; 
                C_cur[ i + numXnodes*2 ][ j + ctcNumNodes   ] = +vintg_Xi_dNjdX; 
                C_cur[ i + numXnodes*2 ][ j + ctcNumNodes*2 ] =           mZero; 
		    }
		}

		// Surface integral
		Add_C_NxN_Matrix_ctc( C_cur ); 
		//Add_D_NxN_Matrix( D_div );
		Add_CtC_D_NxN_Matrix( D_div ); 

		// Volumetric integral of NiNj
		Matrix< std::complex<double> > N_nij; N_nij.Resize( ctcNumDofs, ctcNumDofs, mZero );
		
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
            
            	N_nij[ i                 ][ j                 ] = cteNiNj * int_NiNj; 
            	N_nij[ i + ctcNumNodes   ][ j + ctcNumNodes   ] = cteNiNj * int_NiNj; 
            	N_nij[ i + ctcNumNodes*2 ][ j + ctcNumNodes*2 ] = cteNiNj * int_NiNj; 
            }
		}
				
        Matrix< std::complex<double> > iA_div; iA_div.Resize( numXnodes, numXnodes, mZero );
		Matrix< std::complex<double> > iA_cur; iA_cur.Resize( numXdofs , numXdofs , mZero );

		Matrix<double> XiXj; XiXj.Resize( numXnodes, numXnodes, 0.0 );

		for( int i=0; i<numXnodes; i++ )
		{
            for( int j=0; j<numXnodes; j++ )
            {
				double vintg_Xi_Xj = 0.00;
		    
                for( int gp=0; gp<nGaussPoints; gp++ )
                {
		    		vintg_Xi_Xj += W[gp] * Xij[i][gp] * Xij[j][gp]; 
				}
        		
				XiXj[ i ][ j ] = vintg_Xi_Xj;
            }									   
		}

		Invert_Matrix( XiXj );
		
		for( i=0; i<numXnodes; i++ )
		{
            for( j=0; j<numXnodes; j++ )
            {
        		iA_div[ i ][ j ] = XiXj[i][j];                          
				
				iA_cur[ i               ][ j               ] = iA_div[ i ][ j ]; 
                iA_cur[ i + numXnodes   ][ j + numXnodes   ] = iA_div[ i ][ j ];
                iA_cur[ i + numXnodes*2 ][ j + numXnodes*2 ] = iA_div[ i ][ j ];
            }									   
		}

		// Scaling the projection matrices
		iA_div /= cMu ;
		iA_cur /= cMu ;

		// Final assembly: [ CT * Ac^(-1) * C ] + [ DT * Ad^(-1) * D ] - [ Ni * Nj ]
		StiffMatrix += C_cur.Transpose() * ( iA_cur * C_cur );
		StiffMatrix += D_div.Transpose() * ( iA_div * D_div );
		StiffMatrix -= N_nij;
    }
*/

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

		int nGaussPoints = GaussPoints3D_Order8( cX, cY, cZ, W );

		double detJ = 6.00 * mVolume;

		for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

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
		//StiffMatrix += C_cur.Transpose() * ( iA_cur * C_cur );

		//if ( mPeso != 0.00 ) 
		//{
		//	StiffMatrix += D_div.Transpose() * ( iA_div * D_div );
		//}

		Matrix< std::complex<double> > iAd_x_Ddiv = iA_div * D_div;
		Matrix< std::complex<double> > Ddiv_tran  = D_div.Transpose();
		StiffMatrix += ( Ddiv_tran * iAd_x_Ddiv );

		//StiffMatrix -= N_nij;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

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

		// Surface integral
		//Add_C_NxN_Matrix_edge( C_cur ); 
		//Add_D_NxN_Matrix     ( D_div );

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
	void VolumeElement_1bb::EdgeElementBase_2D( Matrix<double>& Lx, Matrix<double>& Ly, Matrix<double>& Lz, 
		                                        Matrix<double>&  N, Vector<int>& face  ) 
	{
		int nPoints = N[0].size();
		int nEdges  = 3  ;
		
		Lx.Resize( nEdges, nPoints );
		Ly.Resize( nEdges, nPoints );
		Lz.Resize( nEdges, nPoints );

		Matrix<double> DN; 
		
		Calculate_Ni_Derivatives( DN );
		
		double edgeL[3];

		double X0 = mNodes[face[0]]->X(); double Y0 = mNodes[face[0]]->Y(); double Z0 = mNodes[face[0]]->Z();
		double X1 = mNodes[face[1]]->X(); double Y1 = mNodes[face[1]]->Y(); double Z1 = mNodes[face[1]]->Z();
		double X2 = mNodes[face[2]]->X(); double Y2 = mNodes[face[2]]->Y(); double Z2 = mNodes[face[2]]->Z();

		edgeL[0] = sqrt( ( X1 - X0 )*( X1 - X0 ) + ( Y1 - Y0 )*( Y1 - Y0 ) + ( Z1 - Z0 )*( Z1 - Z0 ) );
		edgeL[1] = sqrt( ( X2 - X1 )*( X2 - X1 ) + ( Y2 - Y1 )*( Y2 - Y1 ) + ( Z2 - Z1 )*( Z2 - Z1 ) );
		edgeL[2] = sqrt( ( X2 - X0 )*( X2 - X0 ) + ( Y2 - Y0 )*( Y2 - Y0 ) + ( Z2 - Z0 )*( Z2 - Z0 ) );

        if ( mNodes[face[1]]->Id() > mNodes[face[0]]->Id() ) edgeL[0] *= -1.0;
        if ( mNodes[face[2]]->Id() > mNodes[face[1]]->Id() ) edgeL[1] *= -1.0;
        if ( mNodes[face[2]]->Id() > mNodes[face[0]]->Id() ) edgeL[2] *= -1.0;

		for( int p=0; p<nPoints; p++ ) 
        {
			Lx[ 0 ][ p ] = edgeL[ 0 ] * ( N[ 0 ][ p ] * DN[ 0 ][ face[1] ] - N[ 1 ][ p ] * DN[ 0 ][ face[0] ] );
			Ly[ 0 ][ p ] = edgeL[ 0 ] * ( N[ 0 ][ p ] * DN[ 1 ][ face[1] ] - N[ 1 ][ p ] * DN[ 1 ][ face[0] ] );
			Lz[ 0 ][ p ] = edgeL[ 0 ] * ( N[ 0 ][ p ] * DN[ 2 ][ face[1] ] - N[ 1 ][ p ] * DN[ 2 ][ face[0] ] );
	     	  		  		     	 	 	 	      	   	 
			Lx[ 1 ][ p ] = edgeL[ 1 ] * ( N[ 1 ][ p ] * DN[ 0 ][ face[2] ] - N[ 2 ][ p ] * DN[ 0 ][ face[1] ] );
			Ly[ 1 ][ p ] = edgeL[ 1 ] * ( N[ 1 ][ p ] * DN[ 1 ][ face[2] ] - N[ 2 ][ p ] * DN[ 1 ][ face[1] ] );
			Lz[ 1 ][ p ] = edgeL[ 1 ] * ( N[ 1 ][ p ] * DN[ 2 ][ face[2] ] - N[ 2 ][ p ] * DN[ 2 ][ face[1] ] );
							     	  		  		     	 	 	 	      	   	 
			Lx[ 2 ][ p ] = edgeL[ 2 ] * ( N[ 0 ][ p ] * DN[ 0 ][ face[2] ] - N[ 2 ][ p ] * DN[ 0 ][ face[0] ] );
			Ly[ 2 ][ p ] = edgeL[ 2 ] * ( N[ 0 ][ p ] * DN[ 1 ][ face[2] ] - N[ 2 ][ p ] * DN[ 1 ][ face[0] ] );
			Lz[ 2 ][ p ] = edgeL[ 2 ] * ( N[ 0 ][ p ] * DN[ 2 ][ face[2] ] - N[ 2 ][ p ] * DN[ 2 ][ face[0] ] );
        }	
	}

	//******************************************************************************************************
    //* - 
    //******************************************************************************************************	
	void VolumeElement_1bb::Add_C_NxN_Matrix_edge( Matrix< std::complex<double> >& C_cur ) 
	{
		Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
		Faces[1].resize( 6 ); Faces[1][0] = 1; Faces[1][1] = 3; Faces[1][2] = 2; Faces[1][3] = 10; Faces[1][4] = 11; Faces[1][5] = 12;
		Faces[2].resize( 6 ); Faces[2][0] = 2; Faces[2][1] = 3; Faces[2][2] = 0; Faces[2][3] = 13; Faces[2][4] = 14; Faces[2][5] = 15;
		Faces[3].resize( 6 ); Faces[3][0] = 3; Faces[3][1] = 1; Faces[3][2] = 0; Faces[3][3] =  7; Faces[3][4] =  8; Faces[3][5] =  9;
		
		Vector< Vector<int> >::iterator fs_it;

		for( fs_it = Faces.begin(); fs_it != Faces.end(); fs_it++ )
		{
			Vector<int> nf = *fs_it;

			if ( (*gpNormalType)[ mNodes[ nf[3] ]->Id() ] != 'E' ) continue;

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
		    
		    Lagrange2D_Ni_2bb( N, cX, cY );

            Matrix<double> Lx, Ly, Lz;
		    
		    EdgeElementBase_2D( Lx, Ly, Lz, N, nf );

		    Vector<double> extN; Calculate_ExtNormal( extN, nf );

			for( i=0; i<3; i++ )
		    {
		        for( j=0; j<6; j++ )
		        {
		        	std::complex<double> Sintg_LxiNj(0.00,0.00);
					std::complex<double> Sintg_LyiNj(0.00,0.00);
					std::complex<double> Sintg_LziNj(0.00,0.00);
		            
                    for( gp=0; gp<nGaussPoints; gp++ ) 
					{
						Sintg_LxiNj += W[gp] * Lx[i][gp] * N[j][gp]; 
						Sintg_LyiNj += W[gp] * Ly[i][gp] * N[j][gp]; 
						Sintg_LziNj += W[gp] * Lz[i][gp] * N[j][gp]; 
					}
		        		
					C_cur[ nf[i] ][ nf[j]               ] += ( extN[2] * Sintg_LyiNj - extN[1] * Sintg_LziNj ); 
                    C_cur[ nf[i] ][ nf[j] + mNumNodes   ] += ( extN[0] * Sintg_LziNj - extN[2] * Sintg_LxiNj ); 
                    C_cur[ nf[i] ][ nf[j] + mNumNodes*2 ] += ( extN[1] * Sintg_LxiNj - extN[0] * Sintg_LyiNj ); 
		        }
		    }
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

		Lagrange3D_Ni_2bb( N, cX, cY, cZ );

		// Derivatives on Gauss points
		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_2bb( dNdx, dNdy, dNdz, cX, cY, cZ );

		Matrix<double> Lx, Ly, Lz;

		EdgeElementBase_3D( Lx, Ly, Lz, N, dNdx, dNdy, dNdz );

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
	//=
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
	//=
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

    //===========================================================================
	//=
	//===========================================================================
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

	////////////////////////////////////////////////////////////////////////////////////////

    // prints an arbitrary size matrix to the standard output
    void VolumeElement_1bb::printMatrix(double **a, int rows, int cols) {
    	int i, j;
    
    	for (i = 0; i < rows; i++) {
    		for (j = 0; j < cols; j++) {
    			printf("%+.4lf ", a[i][j]);
    		}
    		printf("\n");
    	}
    	printf("\n");
    }
    
    // prints an arbitrary size vector to the standard output
    void VolumeElement_1bb::printVector(double *v, int size) {
    	int i;
    
    	for (i = 0; i < size; i++) {
    		printf("%+.4lf ", v[i]);
    	}
    	printf("\n\n");
    }
    
    // calculates sqrt( a^2 + b^2 ) with decent precision
    double VolumeElement_1bb::pythag(double a, double b) {
    	double absa, absb;
    
    	absa = fabs(a);
    	absb = fabs(b);
    
    	if (absa > absb)
    		return (absa * sqrt(1.0 + SQR(absb/absa)));
    	else
    		return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
    }
    
    /*
     Modified from Numerical Recipes in C
     Given a matrix a[nRows][nCols], svdcmp() computes its singular value
     decomposition, A = U * W * Vt.  A is replaced by U when svdcmp
     returns.  The diagonal matrix W is output as a vector w[nCols].
     V (not V transpose) is output as the matrix V[nCols][nCols].
     */
    int VolumeElement_1bb::svdcmp(double **a, int nRows, int nCols, double *w, double **v) 
	{
    	int flag, i, its, j, jj, k, l, nm;
    	double anorm, c, f, g, h, s, scale, x, y, z, *rv1;
    
    	rv1 = (double*) malloc(sizeof(double) * nCols);
    	if (rv1 == NULL) {
    		printf("svdcmp(): Unable to allocate vector\n");
    		return (-1);
    	}
    
    	g = scale = anorm = 0.0;
    	for (i = 0; i < nCols; i++) {
    		l = i + 1;
    		rv1[i] = scale * g;
    		g = s = scale = 0.0;
    		if (i < nRows) {
    			for (k = i; k < nRows; k++)
    				scale += fabs(a[k][i]);
    			if (scale) {
    				for (k = i; k < nRows; k++) {
    					a[k][i] /= scale;
    					s += a[k][i] * a[k][i];
    				}
    				f = a[i][i];
    				g = -SIGN(sqrt(s),f);
    				h = f * g - s;
    				a[i][i] = f - g;
    				for (j = l; j < nCols; j++) {
    					for (s = 0.0, k = i; k < nRows; k++)
    						s += a[k][i] * a[k][j];
    					f = s / h;
    					for (k = i; k < nRows; k++)
    						a[k][j] += f * a[k][i];
    				}
    				for (k = i; k < nRows; k++)
    					a[k][i] *= scale;
    			}
    		}
    		w[i] = scale * g;
    		g = s = scale = 0.0;
    		if (i < nRows && i != nCols - 1) {
    			for (k = l; k < nCols; k++)
    				scale += fabs(a[i][k]);
    			if (scale) {
    				for (k = l; k < nCols; k++) {
    					a[i][k] /= scale;
    					s += a[i][k] * a[i][k];
    				}
    				f = a[i][l];
    				g = -SIGN(sqrt(s),f);
    				h = f * g - s;
    				a[i][l] = f - g;
    				for (k = l; k < nCols; k++)
    					rv1[k] = a[i][k] / h;
    				for (j = l; j < nRows; j++) {
    					for (s = 0.0, k = l; k < nCols; k++)
    						s += a[j][k] * a[i][k];
    					for (k = l; k < nCols; k++)
    						a[j][k] += s * rv1[k];
    				}
    				for (k = l; k < nCols; k++)
    					a[i][k] *= scale;
    			}
    		}
    		anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    
    		printf(".");
    		fflush(stdout);
    	}
    
    	for (i = nCols - 1; i >= 0; i--) {
    		if (i < nCols - 1) {
    			if (g) {
    				for (j = l; j < nCols; j++)
    					v[j][i] = (a[i][j] / a[i][l]) / g;
    				for (j = l; j < nCols; j++) {
    					for (s = 0.0, k = l; k < nCols; k++)
    						s += a[i][k] * v[k][j];
    					for (k = l; k < nCols; k++)
    						v[k][j] += s * v[k][i];
    				}
    			}
    			for (j = l; j < nCols; j++)
    				v[i][j] = v[j][i] = 0.0;
    		}
    		v[i][i] = 1.0;
    		g = rv1[i];
    		l = i;
    		printf(".");
    		fflush(stdout);
    	}
    
    	for (i = IMIN(nRows,nCols) - 1; i >= 0; i--) {
    		l = i + 1;
    		g = w[i];
    		for (j = l; j < nCols; j++)
    			a[i][j] = 0.0;
    		if (g) {
    			g = 1.0 / g;
    			for (j = l; j < nCols; j++) {
    				for (s = 0.0, k = l; k < nRows; k++)
    					s += a[k][i] * a[k][j];
    				f = (s / a[i][i]) * g;
    				for (k = i; k < nRows; k++)
    					a[k][j] += f * a[k][i];
    			}
    			for (j = i; j < nRows; j++)
    				a[j][i] *= g;
    		} else
    			for (j = i; j < nRows; j++)
    				a[j][i] = 0.0;
    		++a[i][i];
    		printf(".");
    		fflush(stdout);
    	}
    
    	for (k = nCols - 1; k >= 0; k--) {
    		for (its = 0; its < 30; its++) {
    			flag = 1;
    			for (l = k; l >= 0; l--) {
    				nm = l - 1;
    				if ((fabs(rv1[l]) + anorm) == anorm) {
    					flag = 0;
    					break;
    				}
    				if ((fabs(w[nm]) + anorm) == anorm)
    					break;
    			}
    			if (flag) {
    				c = 0.0;
    				s = 1.0;
    				for (i = l; i <= k; i++) {
    					f = s * rv1[i];
    					rv1[i] = c * rv1[i];
    					if ((fabs(f) + anorm) == anorm)
    						break;
    					g = w[i];
    					h = pythag(f, g);
    					w[i] = h;
    					h = 1.0 / h;
    					c = g * h;
    					s = -f * h;
    					for (j = 0; j < nRows; j++) {
    						y = a[j][nm];
    						z = a[j][i];
    						a[j][nm] = y * c + z * s;
    						a[j][i] = z * c - y * s;
    					}
    				}
    			}
    			z = w[k];
    			if (l == k) {
    				if (z < 0.0) {
    					w[k] = -z;
    					for (j = 0; j < nCols; j++)
    						v[j][k] = -v[j][k];
    				}
    				break;
    			}
    			if (its == 29)
    				printf("no convergence in 30 svdcmp iterations\n");
    			x = w[l];
    			nm = k - 1;
    			y = w[nm];
    			g = rv1[nm];
    			h = rv1[k];
    			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
    			g = pythag(f, 1.0);
    			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,f)))- h)) / x;
    			c = s = 1.0;
    			for (j = l; j <= nm; j++) {
    				i = j + 1;
    				g = rv1[i];
    				y = w[i];
    				h = s * g;
    				g = c * g;
    				z = pythag(f, h);
    				rv1[j] = z;
    				c = f / z;
    				s = h / z;
    				f = x * c + g * s;
    				g = g * c - x * s;
    				h = y * s;
    				y *= c;
    				for (jj = 0; jj < nCols; jj++) {
    					x = v[jj][j];
    					z = v[jj][i];
    					v[jj][j] = x * c + z * s;
    					v[jj][i] = z * c - x * s;
    				}
    				z = pythag(f, h);
    				w[j] = z;
    				if (z) {
    					z = 1.0 / z;
    					c = f * z;
    					s = h * z;
    				}
    				f = c * g + s * y;
    				x = c * y - s * g;
    				for (jj = 0; jj < nRows; jj++) {
    					y = a[jj][j];
    					z = a[jj][i];
    					a[jj][j] = y * c + z * s;
    					a[jj][i] = z * c - y * s;
    				}
    			}
    			rv1[l] = 0.0;
    			rv1[k] = f;
    			w[k] = x;
    		}
    		printf(".");
    		fflush(stdout);
    	}
    	printf("\n");
    
    	free(rv1);
    
    	return (0);
    }

} /* end namespace Kratos */ 