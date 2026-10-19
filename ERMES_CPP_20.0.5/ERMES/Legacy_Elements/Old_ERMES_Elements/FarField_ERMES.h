
#if !defined(KRATOS_FAR_FIELD_ERMES)
#define  KRATOS_FAR_FIELD_ERMES


#include "../includes/element.h"

namespace Kratos
{

    class FarField_ERMES : public Element
    {
        public:

            FarField_ERMES(NodesArrayType& rNodes): Element(rNodes)
            {}

            void GetStiffMatrix9x9(std::complex<double> StiffMatrix_9x9[9][9]);

            void GetResidualVector(std::complex<double>* ResidualVector_9); 

            void GetEquationIdVector(int* EqIdVector_9);

            void SetFrequency(double NewFrequency);

        private:

            double freq;
    }; 

} 


#endif 