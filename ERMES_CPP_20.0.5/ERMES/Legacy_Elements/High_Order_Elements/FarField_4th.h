
#if !defined(KRATOS_FAR_FIELD_4TH)
#define  KRATOS_FAR_FIELD_4TH

#include "../includes/element.h"

namespace Kratos
{
    class FarField_4th: public Element
    {
        public:

            FarField_4th(NodesArrayType& rNodes):Element(rNodes)
            {}

            void GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix); 
              
            void GetEquationIdVector(Vector<int>& EquationId);

            void SetFrequency(double NewFrequency);

        private:

            double freq;
    }; 
} 

#endif 