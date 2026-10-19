
#if !defined(KRATOS_SOURCES_BASE)
#define KRATOS_SOURCES_BASE

#include <map>

#include "../includes/constant_source.h"
#include "../includes/variable.h"

namespace Kratos
{
    template<class TDataType>
    class SourcesBase
    {
      public:
             
        SourcesBase() : mpDefaultSource( new ConstantSource<TDataType>( TDataType() ) )
        {
        }
    
        virtual ~SourcesBase()
        {
        }
          
        typename SourceFunction<TDataType>::Pointer GetSource( const Variable<TDataType>& rVariable )
        {
            if( HasSourceFor( rVariable ) )
            {
                return mData[rVariable.getKey()];
            }
	 	    else   
            {
                return mpDefaultSource;
            }
        }
    
        void AddSource( const Variable<TDataType>& rVariable, typename SourceFunction<TDataType>::Pointer pNewFunction )
        {
            mData[rVariable.getKey()] = pNewFunction;
        }
       
	    bool HasSourceFor( const Variable<TDataType>& rVariable )
        {
            return ( mData.find( rVariable.getKey() ) != mData.end() );
        }
           
      private:
        
	    typename SourceFunction<TDataType>::Pointer mpDefaultSource;
         
        std::map<int, typename SourceFunction<TDataType>::Pointer> mData;
    }; 
} 

#endif 

