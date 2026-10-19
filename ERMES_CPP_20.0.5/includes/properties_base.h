
#if !defined(KRATOS_PROPERTIES_BASE)
#define KRATOS_PROPERTIES_BASE

#include <map>

#include "../includes/constant_property.h"
#include "../includes/variable.h"
#include "../includes/properties_default_function.h"

namespace Kratos
{
    template<class TDataType>
    class PropertiesBase
    {
      public:

        typedef PropertyFunction<TDataType>                   FunctionType       ;
        typedef PropertiesDefaultFunction<TDataType>          DefaultFunctionType;
        typedef std::map<int, typename FunctionType::Pointer> ContainerType      ;
        typedef Variable<TDataType>                           VariableType       ;
 
        PropertiesBase() : mpDefaultFunction( new DefaultFunctionType() )
        {
        }
        
        PropertiesBase( typename PropertiesDefaultFunction<TDataType>::Pointer pNewFunction ) : mpDefaultFunction( pNewFunction )
        {
        }
        
        virtual ~PropertiesBase()
        {
        }
        
        // Find the Function respect to the given variable
        void SetProperty( const VariableType& rVariable, typename FunctionType::Pointer pNewFunction )
        {
            mData[rVariable.getKey()] = pNewFunction;
        }

        TDataType GetProperty( const VariableType& rVariable, const Node& rThisNode, double Time )
        {
            typename ContainerType::iterator result;
        
            if( ( result = mData.find( rVariable.getKey() ) ) != mData.end() )
            {
                return ( *(result->second) )( rThisNode );
            }
            else
            {
                return ( *mpDefaultFunction )( rVariable, rThisNode, Time );
            }
        }
        
        TDataType GetProperty( const VariableType& rVariable, const Node& rThisNode )
        {
            typename ContainerType::iterator result;
            
            if( ( result = mData.find( rVariable.getKey() ) ) != mData.end() )
            {
                return ( *(result->second) )( rThisNode );
            }
            else
            {
                return ( *mpDefaultFunction )( rVariable, rThisNode );
            }
        }
 
        TDataType GetProperty( const VariableType& rVariable )
        {
            typename ContainerType::iterator result;
            
            if( ( result = mData.find( rVariable.getKey() ) ) != mData.end() )
            {
                return ( result->second )->getValue();
            }
            else
            {
                return ( *mpDefaultFunction )( rVariable );
            }
        }

      private:

        // Data stored in this memeber variable
        ContainerType mData;

        typename DefaultFunctionType::Pointer mpDefaultFunction;
    }; 
}  

#endif 

