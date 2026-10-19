
#include "VolumeElement_1st_Electrostatic.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"


namespace Kratos
{   
    //**********************************************************************************************************************
    //* -Global index of the nodes
    //**********************************************************************************************************************	
    void VolumeElement_1st_Electrostatic::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( mNumDofs );

        for( int i=0; i<mNumNodes; i++ ) 
        { 
            EquationId[i] = mNodes[i]->pDofV()->EquationId();
        }
    }

    //**********************************************************************************************************************
    //* - Calculate element volume
    //**********************************************************************************************************************
    double VolumeElement_1st_Electrostatic::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
         
        return fabs( Det / 6.0 );
    }

    //**********************************************************************************************************************
    //* - Calculation of dN/dx, dN/dy, dN/dz column vector
    //**********************************************************************************************************************
    void VolumeElement_1st_Electrostatic::Calculate_DN( double DN[3][4] )
    {
        double IntgCte = 1.0 / ( 6.0 * mVolume ); 

        // dN/dx
        DN[0][0] = IntgCte*( Y(3)*Z(2) + Y(4)*Z(3) + Y(2)*Z(4) - Y(3)*Z(4) - Y(2)*Z(3) - Y(4)*Z(2) );
        DN[0][1] = IntgCte*( Y(3)*Z(4) + Y(4)*Z(1) + Y(1)*Z(3) - Y(3)*Z(1) - Y(1)*Z(4) - Y(4)*Z(3) );
        DN[0][2] = IntgCte*( Y(2)*Z(1) + Y(1)*Z(4) + Y(4)*Z(2) - Y(2)*Z(4) - Y(4)*Z(1) - Y(1)*Z(2) );
        DN[0][3] = IntgCte*( Y(2)*Z(3) + Y(3)*Z(1) + Y(1)*Z(2) - Y(2)*Z(1) - Y(3)*Z(2) - Y(1)*Z(3) );
        
        // dN/dy
        DN[1][0] = IntgCte*( X(3)*Z(4) + X(4)*Z(2) + X(2)*Z(3) - X(3)*Z(2) - X(2)*Z(4) - X(4)*Z(3) ); 
        DN[1][1] = IntgCte*( X(3)*Z(1) + X(1)*Z(4) + X(4)*Z(3) - X(1)*Z(3) - X(3)*Z(4) - X(4)*Z(1) );
        DN[1][2] = IntgCte*( X(2)*Z(4) + X(4)*Z(1) + X(1)*Z(2) - X(2)*Z(1) - X(4)*Z(2) - X(1)*Z(4) );
        DN[1][3] = IntgCte*( X(1)*Z(3) + X(3)*Z(2) + X(2)*Z(1) - X(3)*Z(1) - X(1)*Z(2) - X(2)*Z(3) );
        
        // dN/dz
        DN[2][0] = IntgCte*( X(4)*Y(3) + X(2)*Y(4) + X(3)*Y(2) - X(2)*Y(3) - X(4)*Y(2) - X(3)*Y(4) );
        DN[2][1] = IntgCte*( X(3)*Y(4) + X(4)*Y(1) + X(1)*Y(3) - X(3)*Y(1) - X(4)*Y(3) - X(1)*Y(4) ); 
        DN[2][2] = IntgCte*( X(2)*Y(1) + X(4)*Y(2) + X(1)*Y(4) - X(4)*Y(1) - X(2)*Y(4) - X(1)*Y(2) ); 
        DN[2][3] = IntgCte*( X(2)*Y(3) + X(3)*Y(1) + X(1)*Y(2) - X(2)*Y(1) - X(3)*Y(2) - X(1)*Y(3) );    
    }

    //**********************************************************************************************************************
    //* - Calculation of dN/dx, dN/dy, dN/dz on cX,cY,cZ points
    //**********************************************************************************************************************
    void VolumeElement_1st_Electrostatic::LagrangeDerivatives3D_1st( Matrix<double>& dNdx,
		                                                             Matrix<double>& dNdy,
													                 Matrix<double>& dNdz,
													                 Vector<double>&   cX, 
													                 Vector<double>&   cY, 
													                 Vector<double>&   cZ )
    {
		int nPoints = cX.size();

		dNdx.Resize( mNumNodes, nPoints );
		dNdy.Resize( mNumNodes, nPoints );
		dNdz.Resize( mNumNodes, nPoints );
		
		double DN[3][4];

        Calculate_DN( DN );

		double dL1dx, dL2dx, dL3dx, dL4dx;
        double dL1dy, dL2dy, dL3dy, dL4dy;
		double dL1dz, dL2dz, dL3dz, dL4dz;

        dL1dx = DN[0][0], dL2dx = DN[0][1], dL3dx = DN[0][2], dL4dx = DN[0][3];
        dL1dy = DN[1][0], dL2dy = DN[1][1], dL3dy = DN[1][2], dL4dy = DN[1][3];
		dL1dz = DN[2][0], dL2dz = DN[2][1], dL3dz = DN[2][2], dL4dz = DN[2][3];

		// dNi/dx, dNi/dy, dNi/dz
        for( int p=0; p<nPoints; p++ ) 
        {
			dNdx[0][p] = dL1dx;
			dNdy[0][p] = dL1dy;
			dNdz[0][p] = dL1dz;

			dNdx[1][p] = dL2dx;
			dNdy[1][p] = dL2dy;
			dNdz[1][p] = dL2dz;

			dNdx[2][p] = dL3dx;
			dNdy[2][p] = dL3dy;
			dNdz[2][p] = dL3dz;

			dNdx[3][p] = dL4dx;
			dNdy[3][p] = dL4dy;
			dNdz[3][p] = dL4dz;
        }
	}

    //**********************************************************************************************************************
    //* - Calculation of the stiffness matrix
    //**********************************************************************************************************************	
    void VolumeElement_1st_Electrostatic::GetStiffnessMatrix( Matrix<double>& StiffMatrix ) 
    {  
        StiffMatrix.Resize( mNumDofs, mNumDofs, 0.0 );

        double eo = 8.8541878176e-12;

        double sigma = (*mProperties)( REAL_ELECTRIC_CONDUCTIVITY );
        double epsln = (*mProperties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;

        double MatVol;

        if( abs( sigma ) > 0.0 )
        {
            MatVol = mVolume * sigma;
        }
        else
        {
            MatVol = mVolume * epsln;   
        }

        double DN[ 3 ][ 4 ]; 
        
        Calculate_DN( DN );

        // Diagonal and upper diagonal
		for( int i=0; i<mNumDofs; i++ )
		{
			for( int j=i; j<mNumDofs; j++ ) 
			{
				StiffMatrix[ i ][ j ] = MatVol * ( DN[ 0 ][ i ] * DN[ 0 ][ j ] + DN[ 1 ][ i ] * DN[ 1 ][ j ] + DN[ 2 ][ i ] * DN[ 2 ][ j ] );
			}
		}

        // Lower diagonal
        for( int i=0; i<mNumDofs; i++ )
        {
            for( int j=i+1; j<mNumDofs; j++ )
            {
                StiffMatrix[ j ][ i ] = StiffMatrix[ i ][ j ];
            }
        } 
    }

    //**********************************************************************************************************************
    //* - Dirichlet contribution to the residual vector
    //**********************************************************************************************************************
    void VolumeElement_1st_Electrostatic::GetResidualVector_Dirichlet( std::map<unsigned int, double>& FixedValue, 
                                                                       Matrix<double>& EleStiffMatrix, 
                                                                       Vector<double>& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, 0.0 );

        for( int j=0; j<mNumNodes; j++ ) 
        { 
            if( mNodes[j]->pDofV()->IsFixed() )
            {    
                double DValue = FixedValue[ mNodes[j]->Id() ];

                for( int i=0; i<mNumDofs; i++ ) 
                {
                    ResidualVector[i] -= EleStiffMatrix[i][j] * DValue;
                }
            }
        }
    }

    //**********************************************************************************************************************
    //* - Calculate derivatives
    //**********************************************************************************************************************
    void VolumeElement_1st_Electrostatic::DerivativesInNodes( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz )
    {
		Vector<double> cX( mNumNodes ); 
		Vector<double> cY( mNumNodes );
		Vector<double> cZ( mNumNodes );

		// Points 0, 1, 2, 3
		cX[0] = 0.0; cY[0] = 0.0; cZ[0] = 0.0;
		cX[1] = 1.0; cY[1] = 0.0; cZ[1] = 0.0;
		cX[2] = 0.0; cY[2] = 1.0; cZ[2] = 0.0;
		cX[3] = 0.0; cY[3] = 0.0; cZ[3] = 1.0;

		LagrangeDerivatives3D_1st( ndNdx, ndNdy, ndNdz, cX, cY, cZ );
	}

    //**********************************************************************************************************************
    //* - Calculates gradient of V on nodes
    //**********************************************************************************************************************
    void VolumeElement_1st_Electrostatic::Calculate_Gradient_OnNodes( Vector2Type& grad_OnNodes )
    {
        Matrix<double> ndNdx; 
        Matrix<double> ndNdy; 
        Matrix<double> ndNdz;

        DerivativesInNodes( ndNdx, ndNdy, ndNdz );

        Vector<double> V_OnNodes( mNumNodes );

        for( int n=0; n<mNumNodes; n++ )
        {
            V_OnNodes[n] = (*mProperties)( VOLTAGE, *mNodes[n] );
        }
 
        grad_OnNodes.resize( mNumNodes );

        Vector<double> vZero( 3, 0.0 );

        for( int n=0; n<mNumNodes; n++ )
        {
            grad_OnNodes[n] = vZero;
  
            for( int i=0; i<mNumNodes; i++ )
            {
                grad_OnNodes[n][0] += ndNdx[i][n] * V_OnNodes[i];
                grad_OnNodes[n][1] += ndNdy[i][n] * V_OnNodes[i];
                grad_OnNodes[n][2] += ndNdz[i][n] * V_OnNodes[i];
            }
        }
    }

    //**********************************************************************************************************************
    //* - Calculates gradient of V on Gauss points
    //**********************************************************************************************************************
    void VolumeElement_1st_Electrostatic::Calculate_Gradient_OnGaussPoints( Vector2Type& grad_OnGPs, int numResultsOnGPs )
    {
        // GiD internal Gauss points
		Vector<double> cX, cY, cZ;
	
		int numGPs = GaussPoints3D_InnerGiD_0rd( cX, cY, cZ );

        // Derivatives dN[i]/dx[gp], dN[i]/dy[gp], dN[i]/dz[gp]
        Matrix<double> ndNdx, ndNdy, ndNdz;

        LagrangeDerivatives3D_1st( ndNdx, ndNdy, ndNdz, cX, cY, cZ );

        // Voltage on nodes
        Vector<double> V_OnNodes( mNumNodes );

        for( int n=0; n<mNumNodes; n++ )
        {
            V_OnNodes[n] = (*mProperties)( VOLTAGE, *mNodes[n] );
        }

        // Gradient of V on Gauss points
        grad_OnGPs.resize( numGPs );

        Vector<double> vZero( 3, 0.0 );
      
        for( int gp=0; gp<numGPs; gp++ )
        {
            grad_OnGPs[gp] = vZero;

            for( int i=0; i<mNumNodes; i++ )
            {
                grad_OnGPs[gp][0] += ndNdx[i][gp] * V_OnNodes[i];
                grad_OnGPs[gp][1] += ndNdy[i][gp] * V_OnNodes[i];
                grad_OnGPs[gp][2] += ndNdz[i][gp] * V_OnNodes[i];
            }
        }
    }
} 