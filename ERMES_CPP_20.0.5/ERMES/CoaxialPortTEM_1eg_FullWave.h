
#if !defined(KRATOS_COAXIAL_PORT_TEM_1EG_FULLWAVE)
#define KRATOS_COAXIAL_PORT_TEM_1EG_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class CoaxialPortTEM_1eg_FullWave : public Element
    {
      public:

        CoaxialPortTEM_1eg_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumEdges = 3;
            mNumVerts = 3;
            
            mNumNodes = 6;
            mNumDofs  = 18; 

            mPotentialsOn = false;

            mArea = Calculate_Area();
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 21;
                mPotentialsOn = true; 
            }
        }

        double GetPortNumber();

        std::complex<double> GetProjection();

        std::complex<double> GetNormalization();

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void GetResidualVector( Vector< std::complex<double> >& ResidualVector );
          
      private:

        int mNumEdges;
        int mNumVerts;

        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;

        double mArea;

	    double Calculate_Area();

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

	    std::complex<double> Calculate_KMuCoeff(); 
	    std::complex<double> Calculate_NrmCoeff();

	    void Calculate_r0( double* r0 ); 
	    
        void Calculate_DN( double DN[ 2 ][ 3 ] );

        void Calculate_Edge_Lengths( Vector<double>& eLength );

        void Calculate_Exterior_Normal( Vector<double>& Next );

        void Get_CoPlanar_XY_Axis( Vector<double>& X_axis,  Vector<double>& Y_axis );

        void Get_Vertices_Surface_Coordinates( Vector<double>& pX,  Vector<double>& pY );

        void Whitney2D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Vector<double>& cX, Vector<double>& cY );

        void Input_Field_Vector( Vector<double>& Fx, Vector<double>& Fy, Vector<double>& cX, Vector<double>& cY );
    }; 
} 

#endif 