
#if !defined(KRATOS_GENERIC_ROBIN_3TH)
#define  KRATOS_GENERIC_ROBIN_3TH

#include "../includes/element.h"

namespace Kratos
{
    class GenericRobin_3th : public Element
    {
        public:

            GenericRobin_3th(NodesArrayType& rNodes, 
                             Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}

            void GetStiffMatrix     (Matrix<std::complex<double> >& StiffMatrix); 

            void GetEquationIdVector(Vector<int>& EquationId);		
    }; 
} 

#endif 