
#if !defined(KRATOS_COAXIAL_PORT_TEM_2NDQ)
#define  KRATOS_COAXIAL_PORT_TEM_2NDQ

#include "../includes/element.h"

namespace Kratos
{
    class CoaxialPortTEM_2ndQ : public Element
    {
        public:

            CoaxialPortTEM_2ndQ(NodesArrayType& rNodes, 
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
			
			std::complex<double> CalculateKMuCoef(); 
			
			std::complex<double> CalculateNormalizationCoef(); 

			void Calculate_detJ(std::vector<double>& detJ,
		                        std::vector<double>& cX, 
								std::vector<double>& cY);
    }; 
} 

#endif 