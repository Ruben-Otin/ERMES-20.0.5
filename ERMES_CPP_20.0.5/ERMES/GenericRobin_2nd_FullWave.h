
#if !defined(KRATOS_GENERIC_ROBIN_2ND_FULLWAVE)
#define KRATOS_GENERIC_ROBIN_2ND_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class GenericRobin_2nd_FullWave : public Element
    {
      public:

        GenericRobin_2nd_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
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

        void GetResidualVector_Flux( Vector< Vector<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector );

        void GetResidualVector_Flux( std::map< unsigned int, std::complex<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector );

      private:

        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;

        double mFreq;
        double mArea;

        double Calculate_Area(); 

		std::complex<double> Calculate_jK();

        std::complex<double> Calculate_jK_div_Mu();

		void Calculate_Exterior_Normal( Vector<double>& Next );

		void Calculate_NiNj_Matrix( Matrix<double>& NiNj_Matrix ); 

        void GetResidualVector_Gauss_Beams( Vector< std::complex<double> >& ResidualVector );

		void Add_PlaneWave_Source( Vector< std::complex<double> >& cF, Vector<double>& Kdir, Vector< std::complex<double> >& ResidualVector );

		void GetResidualVector_Single_Wave( Vector< std::complex<double> >& ResidualVector );

		void GetResidualVector_Multpl_Waves( Vector< Vector<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector );
    }; 
} 

#endif 