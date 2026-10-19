
#if !defined(KRATOS_COAXIAL_PORT_TEM_4TH)
#define  KRATOS_COAXIAL_PORT_TEM_4TH

#include "../includes/element.h"

namespace Kratos
{
    class CoaxialPortTEM_4th : public Element
    {
        public:

            CoaxialPortTEM_4th(NodesArrayType& rNodes, 
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
			
			double CalculateArea();

			std::complex<double> CalculateKMuCoef(); 
			
			std::complex<double> CalculateNormalizationCoef(); 
    }; 
} 

#endif 