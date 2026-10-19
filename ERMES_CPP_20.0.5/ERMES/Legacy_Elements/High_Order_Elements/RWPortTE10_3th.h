
#if !defined(KRATOS_RWPORT_TE10_3TH)
#define  KRATOS_RWPORT_TE10_3TH

#include "../includes/element.h"

namespace Kratos
{

    class RWPortTE10_3th : public Element
    {
        public:

            RWPortTE10_3th(NodesArrayType& rNodes, 
                           Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}

            void GetStiffMatrix     (Matrix<std::complex<double> >& StiffMatrix); 
              
			void GetResidualVector  (Vector<std::complex<double> >& ResidualVector);

            void GetEquationIdVector(Vector<int>& EquationId);

			double getPortNumber();

			std::complex<double> getProjection();

	        std::complex<double> getNormalization();

        private:

			void Calculate_r0(double* r0);
			void Calculate_b (double* b); 
			void Calculate_t (double* t); 

			double CalculateArea();
			double CalculateWidth();
			double CalculateHeight();
			
			std::complex<double> CalculateBetaMuCoef(); 
			std::complex<double> CalculateN10();
    }; 

} 

#endif 