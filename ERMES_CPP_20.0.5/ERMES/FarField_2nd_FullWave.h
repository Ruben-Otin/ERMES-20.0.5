
#if !defined(KRATOS_FAR_FIELD_2ND_FULLWAVE)
#define KRATOS_FAR_FIELD_2ND_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class FarField_2nd_FullWave : public Element
    {
      public:

        FarField_2nd_FullWave( NodesArrayType& rNodes ) : Element( rNodes )
        {
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

        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;

        double mFreq;
        double mArea;

        double Calculate_Area(); 
    }; 
} 

#endif 