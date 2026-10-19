
#if !defined(KRATOS_RWPORT_TE10_3SB_FULLWAVE)
#define KRATOS_RWPORT_TE10_3SB_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class RWPortTE10_3sb_FullWave : public Element
    {
      public:

        RWPortTE10_3sb_FullWave ( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
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

        void GetStiffnessMatrix( Matrix<std::complex<double> >& StiffMatrix ); 

        void GetResidualVector( Vector<std::complex<double> >& ResidualVector );
          
      private:

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
    }; 
} 

#endif 