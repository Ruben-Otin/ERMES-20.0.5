
#if !defined(KRATOS_FAR_FIELD_2ND)
#define  KRATOS_FAR_FIELD_2ND

#include "../includes/element.h"

namespace Kratos
{
    class FarField_2nd: public Element
    {
        public:

            FarField_2nd(NodesArrayType& rNodes):Element(rNodes)
            {}

            void GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix); 
              
            void GetEquationIdVector(Vector<int>& EquationId);

            void SetFrequency(double NewFrequency);


        private:

            double freq;

    }; 
} 

#endif 