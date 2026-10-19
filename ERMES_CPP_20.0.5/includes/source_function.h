
#if !defined(KRATOS_SOURCE_FUNCTION)
#define KRATOS_SOURCE_FUNCTION

#include "../includes/geometry_function.h"
#include "../includes/results.h"

namespace Kratos
{
    template<class TDataType>
    class SourceFunction : public GeometryFunction<TDataType, TDataType>
    {
      public:
       
        typedef std::shared_ptr<SourceFunction<TDataType> > Pointer;
    
        SourceFunction()
        {
        }
        
        virtual ~SourceFunction()
        {
        }
        
        virtual void Update( Results& rResults )
        {
        }
    }; 
}  

#endif 

