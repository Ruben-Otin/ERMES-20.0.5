
#if !defined(KRATOS_CONSTANT_SOURCE)
#define KRATOS_CONSTANT_SOURCE

#include "../includes/source_function.h"

namespace Kratos
{
    template<class TDataType>
    class ConstantSource : public SourceFunction<TDataType>
    {
      public:
      
        typedef std::shared_ptr<ConstantSource<TDataType> > Pointer;
        
        ConstantSource( const TDataType& rData = TDataType() ) : mData(rData) 
        {
        }
        
        virtual ~ConstantSource()
        {
        }
        
        TDataType operator()( const Point<double>& rThisPoint )
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

