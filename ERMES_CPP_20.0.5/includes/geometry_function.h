
#if !defined(KRATOS_GEOMETRY_FUNCTION)
#define KRATOS_GEOMETRY_FUNCTION

#include "../includes/point.h"

namespace Kratos
{
    // Base class for all the geometry functions.
    // Geometry functions are group of functions which get a point as input 
    // argument. This function could be set by their constructor and () operator
    // will be override with a point as argument. This implementation is not
    // good in performance aspects and must change in the case of high
    // performance usage. This class is abstract.
    template<class TArgumentType, class TResultType>
    class GeometryFunction
    {
      public:
      
		typedef std::shared_ptr<GeometryFunction<TArgumentType, TResultType> > Pointer;
          
        GeometryFunction()
        {
        }

        virtual ~GeometryFunction()
        {
        }
              
		virtual TResultType operator()( const Point<double>& rThisPoint )
        {
            return TResultType();
        }
            
        virtual TResultType getValue()
        {
            return TResultType();
        }
    }; 
} 

#endif 

