
#if !defined(KRATOS_FAR_FIELD_2EG_COLDPLASMA)
#define KRATOS_FAR_FIELD_2EG_COLDPLASMA

#include "../includes/element.h"

namespace Kratos
{
    class FarField_2eg_ColdPlasma : public Element
    {
      public:

        FarField_2eg_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumFaces = 1;  
            mNumEdges = 3;
            mNumVerts = 3;
            
            mNumBases = 2*mNumEdges + 2*mNumFaces;
            
            mNumNodes = 10;
            mNumDofs  = 30; 

            mPotentialsOn = false;

            mArea = Calculate_Area();
        }
        
        void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ) 
        { 
            mpColdPlasma = pColdPlasma; 
        }
        
        void SetFrequency( double NewFrequency )             
        { 
            mFreq = NewFrequency; 
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 33;
                mPotentialsOn = true; 
            }
        }

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ); 

      private:

        int mNumFaces;
        int mNumEdges;
        int mNumVerts;
        int mNumBases;
          
        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;

        double mFreq;
        double mArea;

        ColdPlasma::Pointer mpColdPlasma;
        
        String Get_Wave_Type();
        
        double Calculate_Area(); 

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void Calculate_DN( double DN[ 2 ][ 3 ] );

        void Calculate_Edge_Lengths( Vector<double>& eLength );

        void Calculate_Exterior_Normal( Vector<double>& Next );

        void Get_Vertices_Surface_Coordinates( Vector<double>& pX,  Vector<double>& pY );

        void jwK_OnGaussPoints( Vector< std::complex<double> >& jwKgp, Matrix<double>& N );

        void Order_Face_Bases( int n[ 2 ], int e[ 2 ] );

        void Face_Base( int edgeBase[ 2 ], int nodeBase[ 2 ], int nodeID[ 3 ], int edgeID[ 3 ], int oppnID[ 3 ] );

        void Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY );

        void Whitney2D_Ni_2nd( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY );
    }; 
} 

#endif 