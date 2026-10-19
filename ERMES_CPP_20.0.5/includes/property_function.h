
#if !defined(KRATOS_PROPERTY_FUNCTION)
#define KRATOS_PROPERTY_FUNCTION

#include "../includes/geometry_function.h"

namespace Kratos
{
    // Base class for all property function classes.
    // This function add an Update method to geometry function. This Update
    // method will use Result database to update the function respect to results
    // obtained from previous iteration. This ability make this function useful
    // in the case of nonlinear or transient or other iterative solutions. 
    template<class TDataType>
    class PropertyFunction : public GeometryFunction<TDataType, TDataType>
    {
      public:

        typedef std::shared_ptr<PropertyFunction<TDataType> > Pointer;
      
        PropertyFunction()
        {
        }

        virtual ~PropertyFunction()
        {
        }

	    virtual TDataType operator()( const Point<double>& rThisPoint, double Time = double() )
        {
            return TDataType();
        }

        virtual TDataType getValue()
        {
            return TDataType();
        }
    }; 
}  

#endif 

