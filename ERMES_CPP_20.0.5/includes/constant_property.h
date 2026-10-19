
#if !defined(KRATOS_CONSTANT_PROPERTY)
#define KRATOS_CONSTANT_PROPERTY

#include <memory>
#include "../includes/property_function.h"

namespace Kratos
{
    template<class TDataType>
    class ConstantProperty : public PropertyFunction<TDataType>
    {
      public:
      
	    typedef std::shared_ptr<ConstantProperty<TDataType> > Pointer;
        
        ConstantProperty( const TDataType& rData = TDataType() ) : mData(rData) 
        {
        }
        
        virtual ~ConstantProperty()
        {
        }
        
        TDataType operator()( const Point<double>& rThisPoint, double Time = double() )
        {
            return mData;
        }
        
        TDataType getValue()
        {
            return mData;
        }
        
      private:
     
        TDataType mData;
    }; 
} 

#endif 

