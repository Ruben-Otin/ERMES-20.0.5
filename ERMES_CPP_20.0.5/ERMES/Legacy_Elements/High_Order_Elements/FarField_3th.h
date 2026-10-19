
#if !defined(KRATOS_FAR_FIELD_3TH)
#define  KRATOS_FAR_FIELD_3TH

#include "../includes/element.h"

namespace Kratos
{
    class FarField_3th: public Element
    {
        public:

            FarField_3th(NodesArrayType& rNodes):Element(rNodes)
            {}

            void GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix); 
              
            void GetEquationIdVector(Vector<int>& EquationId);

            void SetFrequency(double NewFrequency);

        private:

            double freq;
    }; 
} 

#endif 