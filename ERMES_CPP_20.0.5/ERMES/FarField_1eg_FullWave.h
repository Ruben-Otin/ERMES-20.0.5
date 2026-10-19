
#if !defined(KRATOS_FAR_FIELD_1EG_FULLWAVE)
#define KRATOS_FAR_FIELD_1EG_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class FarField_1eg_FullWave : public Element
    {
      public:

        FarField_1eg_FullWave( NodesArrayType& rNodes ) : Element( rNodes )
        {
            mNumEdges = 3;
            mNumVerts = 3;
            
            mNumNodes = 6;
            mNumDofs  = 18; 

            mPotentialsOn = false;
			
            mArea = Calculate_Area();
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

        double Calculate_Area(); 

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void Calculate_DN( double DN[ 2 ][ 3 ] );

        void Calculate_Edge_Lengths( Vector<double>& eLength );

        void Calculate_Exterior_Normal( Vector<double>& Next );

        void Get_Vertices_Surface_Coordinates( Vector<double>& pX,  Vector<double>& pY );

        void Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY );
    }; 
} 

#endif 