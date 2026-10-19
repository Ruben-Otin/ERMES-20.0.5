
#if !defined(KRATOS_FAR_FIELD_1EG_COLDPLASMA)
#define KRATOS_FAR_FIELD_1EG_COLDPLASMA

#include "../includes/element.h"

namespace Kratos
{
    class FarField_1eg_ColdPlasma : public Element
    {
      public:

        FarField_1eg_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumEdges = 3;
            mNumVerts = 3;
            
            mNumNodes = 6;
            mNumDofs  = 18;

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
                mNumDofs = 21;
                mPotentialsOn = true; 
            }
        }

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ); 

      private:

        int mNumEdges;
        int mNumVerts;
          
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

        void Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY );
    }; 
} 

#endif 