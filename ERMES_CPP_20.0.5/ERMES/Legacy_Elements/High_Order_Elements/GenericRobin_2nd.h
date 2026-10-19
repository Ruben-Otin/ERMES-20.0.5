
#if !defined(KRATOS_GENERIC_ROBIN_2ND)
#define  KRATOS_GENERIC_ROBIN_2ND

#include "../includes/element.h"

namespace Kratos
{
    class GenericRobin_2nd : public Element
    {
        public:

            GenericRobin_2nd(NodesArrayType& rNodes, 
                             Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}

            void GetStiffMatrix     (Matrix<std::complex<double> >& StiffMatrix); 

            void GetEquationIdVector(Vector<int>& EquationId);		
    }; 
} 

#endif 