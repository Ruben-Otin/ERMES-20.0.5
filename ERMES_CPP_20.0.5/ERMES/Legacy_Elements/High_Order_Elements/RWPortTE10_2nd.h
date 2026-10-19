
#if !defined(KRATOS_RWPORT_TE10_2ND)
#define  KRATOS_RWPORT_TE10_2ND

#include "../includes/element.h"

namespace Kratos
{
    class RWPortTE10_2nd : public Element
    {
        public:

            RWPortTE10_2nd( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
            {
			}

            void GetStiffMatrix     ( Matrix<std::complex<double> >& StiffMatrix    ); 
      		  
			void GetResidualVector  ( Vector<std::complex<double> >& ResidualVector );
									  
            void GetEquationIdVector( Vector<int>& EquationId );

			std::complex<double> getProjection   ();

	        std::complex<double> getNormalization();

			double getPortNumber();

        private:

			void Calculate_r0( double* r0 );
			void Calculate_b ( double*  b ); 
			void Calculate_t ( double*  t ); 

			double CalculateArea  ();
			double CalculateWidth ();
			double CalculateHeight();
			
			std::complex<double> CalculateBetaMuCoef(); 
			std::complex<double> CalculateN10       ();
    }; 
} 

#endif 