
#if !defined(KRATOS_VOLUME_ELEMENT_1ST_ELECTROSTATIC)
#define KRATOS_VOLUME_ELEMENT_1ST_ELECTROSTATIC

#include "../includes/element.h"

namespace Kratos
{
    class VolumeElement_1st_Electrostatic : public Element
    {
      public:

        VolumeElement_1st_Electrostatic( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 4;
            mNumDofs  = 4;

            mVolume = Calculate_Volume();
        }
        
        void GetStiffnessMatrix( Matrix<double>& StiffMatrix );

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetResidualVector_Dirichlet( std::map<unsigned int, double>& FixedValue, Matrix<double>& EleStiffMatrix, Vector<double>& ResidualVector ); 
        
        void Calculate_Gradient_OnGaussPoints( Vector2Type& grad_OnGPs, int numResultsOnGPs );

        void Calculate_Gradient_OnNodes( Vector2Type& grad_OnNodes );
       
      private:

        int mNumNodes;
        int mNumDofs ;

        double mVolume;

        double Calculate_Volume(); 

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

		void Calculate_DN( double DN[3][4] );

		void LagrangeDerivatives3D_1st( Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

		void DerivativesInNodes( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz );
    }; 
} 

#endif