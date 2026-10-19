
#include "GenericRobin_1st_Electrostatic.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //*********************************************************************************************
    //* - Global IDs of the DOFs
    //*********************************************************************************************	
    void GenericRobin_1st_Electrostatic::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( mNumDofs );

        for( int i=0; i<mNumNodes; i++ ) 
        { 
            EquationId[ i ] = mNodes[ i ]->pDofV()->EquationId();
        }
    }

    //*********************************************************************************************
    //* - Area of the element
    //*********************************************************************************************	
    double GenericRobin_1st_Electrostatic::Calculate_Area()
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

    //*********************************************************************************************
    //* - Stiffness matrix
    //*********************************************************************************************
    void GenericRobin_1st_Electrostatic::GetStiffnessMatrix( Matrix<double>& StiffMatrix ) 
    {
		StiffMatrix.Resize( mNumDofs, mNumDofs, 0.0 );

        double GenericRobinCoeff = (*mProperties)( COMPLEX_IBC          )[ 0 ];
        double TypeOfBoundCondtn = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ];

        if( ( GenericRobinCoeff == 0.0 ) || ( TypeOfBoundCondtn != 0.0 ) ) 
        {
            return;
        }

        double GRC_x_Area = GenericRobinCoeff * mArea;

        // Diagonal and upper diagonal
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=i; j<mNumNodes; j++ )
            {
                if( j != i )
                {
                    StiffMatrix[ i ][ j ] = GRC_x_Area / 12.0;
                }
                else
                {
                    StiffMatrix[ i ][ j ] = GRC_x_Area / 6.0;
                }
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

    //*********************************************************************************************
    //* - Residual vector
    //*********************************************************************************************
    void GenericRobin_1st_Electrostatic::GetResidualVector( Vector<double>& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, 0.0 );

        double ESFlux = (*mProperties)( COMPLEX_IBC_2o       )[ 0 ];
        double BCType = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ];

        if( ( ESFlux == 0.0 ) || ( BCType != 0.0 ) ) 
        {
            return;
        }

        double Flux_x_Area = ESFlux * mArea;

        for( int i=0; i<mNumDofs; i++ )
        {
            ResidualVector[ i ] = Flux_x_Area / 3.0;
        } 
    }
} 