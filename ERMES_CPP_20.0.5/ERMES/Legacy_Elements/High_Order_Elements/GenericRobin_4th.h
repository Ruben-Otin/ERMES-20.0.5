
#if !defined(KRATOS_GENERIC_ROBIN_4TH)
#define  KRATOS_GENERIC_ROBIN_4TH

#include "../includes/element.h"

namespace Kratos
{
    class GenericRobin_4th : public Element
    {
        public:

            GenericRobin_4th(NodesArrayType& rNodes, 
                             Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}

            void GetStiffMatrix     (Matrix<std::complex<double> >& StiffMatrix); 

            void GetEquationIdVector(Vector<int>& EquationId);		
    }; 
} 

#endif 