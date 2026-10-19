
#if !defined(KRATOS_FAR_FIELD_1ST_FULLWAVE)
#define KRATOS_FAR_FIELD_1ST_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class FarField_1st_FullWave : public Element
    {
      public:

        FarField_1st_FullWave( NodesArrayType& rNodes ) : Element( rNodes )
        {
            mNumNodes = 3;
            mNumDofs  = 9; 

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
                mNumDofs = 12;
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