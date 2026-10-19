
#if !defined(KRATOS_RWPORT_TE10_2EG_FULLWAVE)
#define KRATOS_RWPORT_TE10_2EG_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class RWPortTE10_2eg_FullWave : public Element
    {
      public:

        RWPortTE10_2eg_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
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

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 33;
                mPotentialsOn = true; 
            }
        }

        double GetPortNumber();

        std::complex<double> GetProjection();

        std::complex<double> GetNormalization();

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetStiffnessMatrix( Matrix<std::complex<double> >& StiffMatrix ); 

        void GetResidualVector( Vector<std::complex<double> >& ResidualVector );

      private:

        int mNumFaces;
        int mNumEdges;
        int mNumVerts;
        int mNumBases;

        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;

        double mArea;

        double Calculate_Area();
        double Calculate_Width();
        double Calculate_Height();

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void Calculate_r0( double* r0 );
        void Calculate_b ( double* b  ); 
        void Calculate_t ( double* t  ); 
        
        std::complex<double> Calculate_N10();
        std::complex<double> Calculate_BetaMuCoef(); 

        void Calculate_DN( double DN[ 2 ][ 3 ] );

        void Calculate_Edge_Lengths( Vector<double>& eLength );

        void Calculate_Exterior_Normal( Vector<double>& Next );

        void Get_CoPlanar_XY_Axis( Vector<double>& X_axis,  Vector<double>& Y_axis );

        void Get_Vertices_Surface_Coordinates( Vector<double>& pX,  Vector<double>& pY );

        void Order_Face_Bases( int n[ 2 ], int e[ 2 ] );

        void Face_Base( int edgeBase[ 2 ], int nodeBase[ 2 ], int nodeID[ 3 ], int edgeID[ 3 ], int oppnID[ 3 ] );

        void Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY );

        void Whitney2D_Ni_2nd( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY );

        void Input_Field_Vector( Vector<double>& Fx, Vector<double>& Fy, Vector<double>& cX, Vector<double>& cY );
    }; 
} 

#endif 