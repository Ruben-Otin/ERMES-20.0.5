
#include "FarField_3sb_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //*******************************************************************************************************
    //* - Global DOFs Ids
    //*******************************************************************************************************	
    void FarField_3sb_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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

            for( int i=0; i<3; i++ ) 
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

    //*******************************************************************************************************
    // - Area of the element
    //*******************************************************************************************************
    double FarField_3sb_FullWave::Calculate_Area()
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

    //*******************************************************************************************************
    //* - Element stiffness matrix
    //*******************************************************************************************************	
    void FarField_3sb_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Initializing matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );
        
        // Wave vector K
        double eo = 8.8541878176e-12;
		double pi = 3.14159265358979;
		double mu = pi * 4.0e-7;

        std::complex<double> mjK( 0.0, -mFreq * sqrt( eo / mu ) );		
        
        // Surface integrals
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

        Matrix<double> N; Lagrange2D_Ni_3sb( N, cX, cY );

        Matrix<double> intg_NiNj; intg_NiNj.Resize( mNumNodes, mNumNodes, 0.0 );

        double Jacob = 2.0 * mArea;

        for( int i=0; i<mNumNodes; i++ )
        {
        	for( int j=0; j<mNumNodes; j++ )
        	{
                for( int gp=0; gp<numGaussPoints; gp++ ) 
                {
                    intg_NiNj[ i ][ j ] += Jacob * W[ gp ] * N[ i ][ gp ] * N[ j ][ gp ];
                }
        	}
        }

        // Stiffness matrix
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                StiffMatrix[ i               ][ j               ] = mjK * intg_NiNj[ i ][ j ];
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = mjK * intg_NiNj[ i ][ j ];
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = mjK * intg_NiNj[ i ][ j ];
            }
        }          
    } 
} 