
#include "VolumeElement_4th.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{   
   /**********************************************************************************
    * -Global index of the nodes
    **********************************************************************************/	
    void VolumeElement_4th::GetEquationIdVector(Vector<int>& EquationId)
    {
        const int numNodes = 35;
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
    double VolumeElement_4th::Calculate_Volume()
    {
        double det;

        det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs(det/6.00);
    }

   /**********************************************************
    * - Calculation of dN/dz column vector
    **********************************************************/
    void VolumeElement_4th::Calculate_DN(double DN[3][4])
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

   /**********************************************************************************
    * - Set singularities
    **********************************************************************************/	
    void VolumeElement_4th::SetPeso(double Peso)
    {
        mPeso = Peso;
    }

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Distance to an edge
	////////////////////////////////////////////////////////////////////////////////////////////////////
	double VolumeElement_4th::RtoS(double alpha, double beta, double gamma)
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

   /***********************************************************************************
    * - Calculation of dN/dx, dN/dy, dN/dz on cX,cY,cZ points
    ***********************************************************************************/
    void VolumeElement_4th::LagrangeDerivatives3D_4th(Matrix<double>& dNdx,
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

   /**********************************************************************************
    * - Calculation of the stiffness matrix
    **********************************************************************************/	
    void VolumeElement_4th::GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix) 
    {  
        const int numNodes = 35;
		const int numDofs  = numNodes*3;

		StiffMatrix.Resize(numDofs,numDofs);

        int i,j;
        
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

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int gp;

		int nGaussPoints = GaussPoints3D_Order6(cX, cY, cZ, W);

		Matrix<double> N;

		Lagrange3D_Ni_4th(N, cX, cY, cZ);

		Matrix<double> dNdx, dNdy, dNdz;

		LagrangeDerivatives3D_4th(dNdx, dNdy, dNdz, cX, cY, cZ);

        ///////////////////////////////////////////  Pesos  /////////////////////////////////////////////////
        std::vector<double> RW(nGaussPoints);

		if ( mPeso == 0.00) 
		{
			for ( gp=0; gp<nGaussPoints; gp++ ) RW[gp] = 0.00;
		}
		else
		{
			for ( gp=0; gp<nGaussPoints; gp++ ) RW[gp] = W[gp]; 
		}

        /////////////////////////////////////////////////////////////////////////////////////////////////////
        double jacob = 6.00*Calculate_Volume();

        double intCurl, 
               intDiv, 
               intNN;

        std::complex<double> cteDiv  = jacob / cMu;
        std::complex<double> cteCurl = jacob / cMu;
        std::complex<double> cteNN   = jacob * (freq*freq) * cEps;    

        ////////////////////////////////////// Bloque [Kxx] ////////////////////////////////////////////////
        for (i=0; i<numNodes; i++)
        {
            for (j=i; j<numNodes; j++)
            {
                intCurl = 0.00;
                intDiv  = 0.00; 
                intNN   = 0.00;

                for ( gp=0; gp<nGaussPoints; gp++ )
                {
                    intCurl +=  W[gp] * (dNdy[i][gp] * dNdy[j][gp] + dNdz[i][gp] * dNdz[j][gp]);
                    intDiv  += RW[gp] * (dNdx[i][gp] * dNdx[j][gp]); 
                    intNN   +=  W[gp] * (   N[i][gp] *    N[j][gp]);
                }
                
                StiffMatrix[i][j] = (cteDiv * intDiv) + (cteCurl * intCurl) - (cteNN * intNN);
            }
        }

        ////////////////////////////////////// Bloque [Kyy] ////////////////////////////////////////////////
        for (i=0; i<numNodes; i++)
        {
            for (j=i; j<numNodes; j++)
            {
                intCurl = 0.00;
                intDiv  = 0.00; 
                intNN   = 0.00;

                for ( gp=0; gp<nGaussPoints; gp++ )
                {
                    intCurl +=  W[gp] * (dNdx[i][gp] * dNdx[j][gp] + dNdz[i][gp] * dNdz[j][gp]);
                    intDiv  += RW[gp] * (dNdy[i][gp] * dNdy[j][gp]); 
                    intNN   +=  W[gp] * (   N[i][gp] *    N[j][gp]);
                }
                
                StiffMatrix[i+numNodes][j+numNodes] = (cteDiv * intDiv) + (cteCurl * intCurl) - (cteNN * intNN);
            }
        }

		////////////////////////////////////// Bloque [Kzz] ////////////////////////////////////////////////
        for (i=0; i<numNodes; i++)
        {
            for (j=i; j<numNodes; j++)
            {
                intCurl = 0.00;
                intDiv  = 0.00; 
                intNN   = 0.00;

                for ( gp=0; gp<nGaussPoints; gp++ )
                {
                    intCurl +=  W[gp] * (dNdx[i][gp] * dNdx[j][gp] + dNdy[i][gp] * dNdy[j][gp]);
                    intDiv  += RW[gp] * (dNdz[i][gp] * dNdz[j][gp]); 
                    intNN   +=  W[gp] * (   N[i][gp] *    N[j][gp]);
                }
                
                StiffMatrix[i+(2*numNodes)][j+(2*numNodes)] = (cteDiv * intDiv) + (cteCurl * intCurl) - (cteNN * intNN);
            }
        }

        ////////////////////////////////////// Bloque [Kxy] ////////////////////////////////////////////////
        for (i=0; i<numNodes; i++)
        {
            for (j=0; j<numNodes; j++)
            {
                intCurl = 0.00;
                intDiv  = 0.00; 

                for ( gp=0; gp<nGaussPoints; gp++ )
                {
                    intCurl +=  W[gp] * (dNdy[i][gp] * dNdx[j][gp]);
                    intDiv  += RW[gp] * (dNdx[i][gp] * dNdy[j][gp]); 
                }
                
                StiffMatrix[i][j+numNodes] = (cteDiv * intDiv) - (cteCurl * intCurl);
            }
        }

		////////////////////////////////////// Bloque [Kxz] ////////////////////////////////////////////////
        for (i=0; i<numNodes; i++)
        {
            for (j=0; j<numNodes; j++)
            {
                intCurl = 0.00;
                intDiv  = 0.00; 

                for ( gp=0; gp<nGaussPoints; gp++ )
                {
                    intCurl +=  W[gp] * (dNdz[i][gp] * dNdx[j][gp]);
                    intDiv  += RW[gp] * (dNdx[i][gp] * dNdz[j][gp]); 
                }
                
                StiffMatrix[i][j+(2*numNodes)] = (cteDiv * intDiv) - (cteCurl * intCurl);
            }
        }

		////////////////////////////////////// Bloque [Kyz] ////////////////////////////////////////////////
        for (i=0; i<numNodes; i++)
        {
            for (j=0; j<numNodes; j++)
            {
                intCurl = 0.00;
                intDiv  = 0.00; 

                for ( gp=0; gp<nGaussPoints; gp++ )
                {
                    intCurl +=  W[gp] * (dNdz[i][gp] * dNdy[j][gp]);
                    intDiv  += RW[gp] * (dNdy[i][gp] * dNdz[j][gp]); 
                }
                
                StiffMatrix[i+numNodes][j+(2*numNodes)] = (cteDiv * intDiv) - (cteCurl * intCurl);
            }
        }

        //////////////////////////////////// Bloque Simetrico /////////////////////////////////////////////
        for (i=0; i<numDofs; i++)
        {
            for (j=i+1; j<numDofs; j++)
            {
                StiffMatrix[j][i] =  StiffMatrix[i][j];
            }
        }        
    }

   /*******************************************************************
    * - Calculate derivatives
    ********************************************************************/
    void VolumeElement_4th::DerivativesInNodes(Matrix<double>& ndNdx, 
		                                       Matrix<double>& ndNdy, 
											   Matrix<double>& ndNdz)
    {
		const int numNodes = 35;
		
		std::vector<double> cX(numNodes); 
		std::vector<double>	cY(numNodes);
		std::vector<double> cZ(numNodes);

		//Points 0,1,2,3
		cX[0] = 0.00; cY[0] = 0.00; cZ[0] = 0.00;
		cX[1] = 1.00; cY[1] = 0.00; cZ[1] = 0.00;
		cX[2] = 0.00; cY[2] = 1.00; cZ[2] = 0.00;
		cX[3] = 0.00; cY[3] = 0.00; cZ[3] = 1.00;
		
		//Edge [0 1]
		cX[4]  = 0.25; cY[4]  = 0.00; cZ[4]  = 0.00;
		cX[5]  = 0.50; cY[5]  = 0.00; cZ[5]  = 0.00;
		cX[6]  = 0.75; cY[6]  = 0.00; cZ[6]  = 0.00;
		//Edge [1 2]
		cX[7]  = 0.75; cY[7]  = 0.25; cZ[7]  = 0.00;
		cX[8]  = 0.50; cY[8]  = 0.50; cZ[8]  = 0.00;
		cX[9]  = 0.25; cY[9]  = 0.75; cZ[9]  = 0.00;
		//Edge [2 0] 
		cX[10] = 0.00; cY[10] = 0.75; cZ[10] = 0.00;
		cX[11] = 0.00; cY[11] = 0.50; cZ[11] = 0.00;
		cX[12] = 0.00; cY[12] = 0.25; cZ[12] = 0.00;
		//Edge [0 3]
		cX[13] = 0.00; cY[13] = 0.00; cZ[13] = 0.25;
		cX[14] = 0.00; cY[14] = 0.00; cZ[14] = 0.50;
		cX[15] = 0.00; cY[15] = 0.00; cZ[15] = 0.75;
		//Edge [1 3]
		cX[16] = 0.75; cY[16] = 0.00; cZ[16] = 0.25;
		cX[17] = 0.50; cY[17] = 0.00; cZ[17] = 0.50;
		cX[18] = 0.25; cY[18] = 0.00; cZ[18] = 0.75;
		//Edge [2 3]
		cX[19] = 0.00; cY[19] = 0.75; cZ[19] = 0.25;
		cX[20] = 0.00; cY[20] = 0.50; cZ[20] = 0.50;
		cX[21] = 0.00; cY[21] = 0.25; cZ[21] = 0.75;

		//Face [0 1 2]
		cX[22] = 0.25; cY[22] = 0.25; cZ[22] = 0.00;
		cX[23] = 0.50; cY[23] = 0.25; cZ[23] = 0.00;
		cX[24] = 0.25; cY[24] = 0.50; cZ[24] = 0.00;
		//Face [0 1 3]
		cX[25] = 0.25; cY[25] = 0.00; cZ[25] = 0.25;
		cX[26] = 0.50; cY[26] = 0.00; cZ[26] = 0.25;
		cX[27] = 0.25; cY[27] = 0.00; cZ[27] = 0.50;
		//Face [1 2 3]
		cX[28] = 0.50; cY[28] = 0.25; cZ[28] = 0.25;
		cX[29] = 0.25; cY[29] = 0.50; cZ[29] = 0.25;
		cX[30] = 0.25; cY[30] = 0.25; cZ[30] = 0.50;
		//Face [0 2 3]
		cX[31] = 0.00; cY[31] = 0.25; cZ[31] = 0.25;
		cX[32] = 0.00; cY[32] = 0.50; cZ[32] = 0.25;
		cX[33] = 0.00; cY[33] = 0.25; cZ[33] = 0.50;

		//Volumne [0 1 2 3]
		cX[34] = 0.25; cY[34] = 0.25; cZ[34] = 0.25; 

		LagrangeDerivatives3D_4th(ndNdx, ndNdy, ndNdz, cX, cY, cZ);
	}

   /*******************************************************************
    * - Calculate derivatives
    ********************************************************************/
    void VolumeElement_4th::CalculateDerivatives(Vector<std::complex<double> >& ncHx,
		                                         Vector<std::complex<double> >& ncHy,
												 Vector<std::complex<double> >& ncHz)
    {
        int i,n;

		const int numNodes = 35;
        
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