
#if !defined(DEFINE_KRATOS_COMP_PDOF)
#define DEFINE_KRATOS_COMP_PDOF

#include "../includes/dof.h"

namespace Kratos
{
    struct ComparePDof
    {
        // const: required by std::set/std::map comparators since C++17 (GCC >= 11 rejects it otherwise)
        bool operator()( const Dof::Pointer& a1, const Dof::Pointer& a2 ) const
        {
            if ( a1->NodeIndex() ==  a2->NodeIndex() ) 
            {
                return ( a1->GetVariableKey() < a2->GetVariableKey() ); 
            }
            else 
            {
                return ( a1->NodeIndex() <  a2->NodeIndex() );
            }
        }
    };
}  

#endif 

