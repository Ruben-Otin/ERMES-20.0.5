
#if !defined(KRATOS_FAR_FIELD_2NDQ)
#define  KRATOS_FAR_FIELD_2NDQ

#include "../includes/element.h"

namespace Kratos
{
    class FarField_2ndQ: public Element
    {
        public:

            FarField_2ndQ(NodesArrayType& rNodes):Element(rNodes)
            {}

            void GetStiffMatrix(Matrix<std::complex<double> >& StiffMatrix); 
              
            void GetEquationIdVector(Vector<int>& EquationId);

            void SetFrequency(double NewFrequency);


        private:

            double freq;

			void Calculate_detJ(std::vector<double>& detJ,
		                        std::vector<double>& cX, 
								std::vector<double>& cY);

    }; 
} 

#endif 