
#if !defined(KRATOS_COAXIAL_PORT_TEM_1ST_FULLWAVE)
#define KRATOS_COAXIAL_PORT_TEM_1ST_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class CoaxialPortTEM_1st_FullWave : public Element
    {
      public:
	  
        CoaxialPortTEM_1st_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 3;
            mNumDofs  = 9; 

            mPotentialsOn = false;

            mArea = Calculate_Area();
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 12;
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

        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;

        double mArea;

	    double Calculate_Area();

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

	    void Calculate_r0( double* r0 ); 
	    
	    std::complex<double> Calculate_KMuCoeff(); 
	    std::complex<double> Calculate_NrmCoeff();
    }; 
} 

#endif 