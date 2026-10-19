

#if !defined(KRATOS_RWPORT_TE10_ERMES)
#define  KRATOS_RWPORT_TE10_ERMES


#include "../includes/element.h"

namespace Kratos
{

    class RWPortTE10_ERMES : public Element
    {
        public:

            RWPortTE10_ERMES(NodesArrayType& rNodes, 
                             Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}

            void GetStiffMatrix9x9(std::complex<double> StiffMatrix_9x9[9][9]);

            void GetResidualVector(std::complex<double>* ResidualVector_9); 

            void GetEquationIdVector(int* EqIdVector_9);

            void GetMagnitudeVector(std::complex<double>* MagnitudVector); 


        private:

			void SetUpVariables();

            double NiNj_matrix[3][3];

            double n[3];

            std::complex<double> cbetaMu;
    }; 

} 


#endif 