
#include "FarField_1eg_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //******************************************************************************************************************************************
    //* - Global DOFs Ids
    //******************************************************************************************************************************************	
    void FarField_1eg_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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

            for( int i=0; i<mNumVerts; i++ ) 
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
    // - Area of the element
    //******************************************************************************************************************************************
    double FarField_1eg_FullWave::Calculate_Area()
    {
		double na[ 3 ], v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ 2 ]->X() - mNodes[ 0 ]->X();
        v2[ 1 ] = mNodes[ 2 ]->Y() - mNodes[ 0 ]->Y();
        v2[ 2 ] = mNodes[ 2 ]->Z() - mNodes[ 0 ]->Z();
		    
        v1[ 0 ] = mNodes[ 1 ]->X() - mNodes[ 0 ]->X();
        v1[ 1 ] = mNodes[ 1 ]->Y() - mNodes[ 0 ]->Y();
        v1[ 2 ] = mNodes[ 1 ]->Z() - mNodes[ 0 ]->Z();

        // Area = 0.5 * ||v2 x v1||
        na[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        na[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        na[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

		return ( 0.5 * sqrt( na[ 0 ]*na[ 0 ] + na[ 1 ]*na[ 1 ] + na[ 2 ]*na[ 2 ] ) );
    }

    //******************************************************************************************************************************************
    //* - Exterior normal 
    //******************************************************************************************************************************************	
    void FarField_1eg_FullWave::Calculate_Exterior_Normal( Vector<double>& Next )
    {
		double v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ 2 ]->X() - mNodes[ 0 ]->X();
        v2[ 1 ] = mNodes[ 2 ]->Y() - mNodes[ 0 ]->Y();
        v2[ 2 ] = mNodes[ 2 ]->Z() - mNodes[ 0 ]->Z();
		    
        v1[ 0 ] = mNodes[ 1 ]->X() - mNodes[ 0 ]->X();
        v1[ 1 ] = mNodes[ 1 ]->Y() - mNodes[ 0 ]->Y();
        v1[ 2 ] = mNodes[ 1 ]->Z() - mNodes[ 0 ]->Z();

        Next.resize( 3 );

        Next[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        Next[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        Next[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

        Next *= ( 1.0 / sqrt( Next[ 0 ]*Next[ 0 ] + Next[ 1 ]*Next[ 1 ] + Next[ 2 ]*Next[ 2 ] ) );
    }

    //******************************************************************************************************************************************
    // - Change coordinate system from ( X, Y, Z ) to ( pX, pY )  
    //******************************************************************************************************************************************
    void FarField_1eg_FullWave::Get_Vertices_Surface_Coordinates( Vector<double>& pX,  Vector<double>& pY )
    {
        // Z axis of the co-planar coordinate system 
        Vector<double> Z_axis; 
        
        Calculate_Exterior_Normal( Z_axis );

        // X axis of the co-planar coordinate system 
        Vector<double> X_axis( 3 );

        X_axis[ 0 ] = X(2) - X(1); 
		X_axis[ 1 ] = Y(2) - Y(1); 
		X_axis[ 2 ] = Z(2) - Z(1);

        X_axis *= ( 1.0 / sqrt( X_axis[ 0 ]*X_axis[ 0 ] + X_axis[ 1 ]*X_axis[ 1 ] + X_axis[ 2 ]*X_axis[ 2 ] ) ); 

        // Y axis of the co-planar coordinate system ( Y_axis = X_axis x Z_axis )
        Vector<double> Y_axis( 3 );

        Y_axis[ 0 ] = X_axis[ 1 ] * Z_axis[ 2 ] - X_axis[ 2 ] * Z_axis[ 1 ];
        Y_axis[ 1 ] = X_axis[ 2 ] * Z_axis[ 0 ] - X_axis[ 0 ] * Z_axis[ 2 ];
        Y_axis[ 2 ] = X_axis[ 0 ] * Z_axis[ 1 ] - X_axis[ 1 ] * Z_axis[ 0 ];

        Y_axis *= ( 1.0 / sqrt( Y_axis[ 0 ]*Y_axis[ 0 ] + Y_axis[ 1 ]*Y_axis[ 1 ] + Y_axis[ 2 ]*Y_axis[ 2 ] ) ); 

        // Position vectors of vertice nodes
        Vector<double> R_2( 3 ); 
        
        R_2[ 0 ] = X(2) - X(1); 
        R_2[ 1 ] = Y(2) - Y(1); 
        R_2[ 2 ] = Z(2) - Z(1);
        
        Vector<double> R_3( 3 ); 

        R_3[ 0 ] = X(3) - X(1); 
        R_3[ 1 ] = Y(3) - Y(1); 
        R_3[ 2 ] = Z(3) - Z(1);

        // Co-planar coordinates of vertice nodes
        pX.resize( mNumVerts ); 
        pY.resize( mNumVerts );

        pX[ 0 ] = 0.0; 
        pY[ 0 ] = 0.0;

        pX[ 1 ] = R_2[ 0 ] * X_axis[ 0 ] + R_2[ 1 ] * X_axis[ 1 ] + R_2[ 2 ] * X_axis[ 2 ]; 
        pY[ 1 ] = R_2[ 0 ] * Y_axis[ 0 ] + R_2[ 1 ] * Y_axis[ 1 ] + R_2[ 2 ] * Y_axis[ 2 ];

        pX[ 2 ] = R_3[ 0 ] * X_axis[ 0 ] + R_3[ 1 ] * X_axis[ 1 ] + R_3[ 2 ] * X_axis[ 2 ]; 
        pY[ 2 ] = R_3[ 0 ] * Y_axis[ 0 ] + R_3[ 1 ] * Y_axis[ 1 ] + R_3[ 2 ] * Y_axis[ 2 ];
    }

    //******************************************************************************************************************************************
    // - Lagrange 1st order base derivatives - DN[ X, Y ][ 0, 1, 2 ]
    //******************************************************************************************************************************************
    void FarField_1eg_FullWave::Calculate_DN( double DN[ 2 ][ 3 ] )
    {
        double Cte = 1.0 / ( 2.0 * mArea ); 
        
        // Plane coordinates
        Vector<double> pX, pY; 
        
        Get_Vertices_Surface_Coordinates( pX, pY );

        // dNi/dx
        DN[ 0 ][ 0 ] = Cte * ( pY[ 1 ] - pY[ 2 ] );
        DN[ 0 ][ 1 ] = Cte * ( pY[ 2 ] - pY[ 0 ] );
        DN[ 0 ][ 2 ] = Cte * ( pY[ 0 ] - pY[ 1 ] );
                               
        // dNi/dy              
        DN[ 1 ][ 0 ] = Cte * ( pX[ 2 ] - pX[ 1 ] ); 
        DN[ 1 ][ 1 ] = Cte * ( pX[ 0 ] - pX[ 2 ] );
        DN[ 1 ][ 2 ] = Cte * ( pX[ 1 ] - pX[ 0 ] );
    }

    //******************************************************************************************************************************************
    // - Edge lengths with sign
    //******************************************************************************************************************************************
    void FarField_1eg_FullWave::Calculate_Edge_Lengths( Vector<double>& eLength )
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

		// Edge lengths
		for( int i=0; i<mNumEdges; i++ )
		{
			eLength[ i ] = sqrt( R[ i ][ 0 ]*R[ i ][ 0 ] + R[ i ][ 1 ]*R[ i ][ 1 ] + R[ i ][ 2 ]*R[ i ][ 2 ] );
		}
		
		// Check directions
		if( mNodes[ 0 ]->Id() > mNodes[ 1 ]->Id() ) eLength[ 0 ] *= -1.0;
		if( mNodes[ 1 ]->Id() > mNodes[ 2 ]->Id() ) eLength[ 1 ] *= -1.0;
		if( mNodes[ 2 ]->Id() > mNodes[ 0 ]->Id() ) eLength[ 2 ] *= -1.0;
	}

    //******************************************************************************************************************************************
    // - 1st order Whitney base for 2D triangles
    //******************************************************************************************************************************************
	void FarField_1eg_FullWave::Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY )
	{
		double L1, L2, L3;
		
		int numPoints = cX.size();

		Nx.Resize( mNumEdges, numPoints );
		Ny.Resize( mNumEdges, numPoints );
		
		double DN[ 2 ][ 3 ]; Calculate_DN( DN );

		double dL1_dX = DN[ 0 ][ 0 ], dL1_dY = DN[ 1 ][ 0 ];
		double dL2_dX = DN[ 0 ][ 1 ], dL2_dY = DN[ 1 ][ 1 ];
		double dL3_dX = DN[ 0 ][ 2 ], dL3_dY = DN[ 1 ][ 2 ];

		Vector<double> eLength( mNumEdges ); Calculate_Edge_Lengths( eLength );

		for( int p=0; p<numPoints; p++ ) 
        {
            L2 = cX[ p ];
            L3 = cY[ p ];
            L1 = 1.0 - L2 - L3;

			// Edge: [0 1]
			Nx[ 0 ][ p ] = ( L1 * dL2_dX - L2 * dL1_dX ) * eLength[ 0 ];
            Ny[ 0 ][ p ] = ( L1 * dL2_dY - L2 * dL1_dY ) * eLength[ 0 ];
		    
			// Edge: [1 2]
			Nx[ 1 ][ p ] = ( L2 * dL3_dX - L3 * dL2_dX ) * eLength[ 1 ];
            Ny[ 1 ][ p ] = ( L2 * dL3_dY - L3 * dL2_dY ) * eLength[ 1 ];
		
			// Edge: [2 0]
			Nx[ 2 ][ p ] = ( L3 * dL1_dX - L1 * dL3_dX ) * eLength[ 2 ];
            Ny[ 2 ][ p ] = ( L3 * dL1_dY - L1 * dL3_dY ) * eLength[ 2 ];
        }
	}

    //******************************************************************************************************************************************
    //* - Stiffness matrix ( Ef and AV formulation )
    //******************************************************************************************************************************************	
    void FarField_1eg_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Initializing matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );
        
        // Wave vector -jK
        double eo = 8.8541878176e-12;
		double pi = 3.14159265358979;
		double mu = pi * 4.0e-7;

        std::complex<double> mjK( 0.0, -mFreq * sqrt( eo / mu ) );		

        // Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.0 * mArea;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

        // 1st order Whitney base functions
        Matrix<double> Nx, Ny;

		Whitney2D_Ni_1st( Nx, Ny, cX, cY );

		// Edge element matrix
		Matrix< std::complex<double> > EdgeBaseMatrix( mNumEdges, mNumEdges );

		// Upper diagonal 
		for( int i=0; i<mNumEdges; i++ )
		{
			for( int j=i; j<mNumEdges; j++ )
			{
				double sIntg_NiNj = 0.0;

				for( int gp=0; gp<numGaussPoints; gp++ )
				{
					sIntg_NiNj += W[ gp ] * ( Nx[ i ][ gp ] * Nx[ j ][ gp ] + Ny[ i ][ gp ] * Ny[ j ][ gp ] );
				}

				EdgeBaseMatrix[ i ][ j ] = mjK * sIntg_NiNj;
			}
		}

		// Lower diagonal 
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
} 