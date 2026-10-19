
#if !defined(KRATOS_GENERIC_ROBIN_2NDQ)
#define  KRATOS_GENERIC_ROBIN_2NDQ

#include "../includes/element.h"

namespace Kratos
{
    class GenericRobin_2ndQ : public Element
    {
        public:

            GenericRobin_2ndQ(NodesArrayType& rNodes, 
                              Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}

            void GetStiffMatrix     (Matrix<std::complex<double> >& StiffMatrix); 

            void GetEquationIdVector(Vector<int>& EquationId);
			

		private:

			void Calculate_detJ(std::vector<double>& detJ,
		                        std::vector<double>& cX, 
								std::vector<double>& cY);
    }; 
} 

#endif 