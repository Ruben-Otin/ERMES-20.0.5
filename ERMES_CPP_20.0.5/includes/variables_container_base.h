
#if !defined(KRATOS_VARIABLES_CONTAINER_BASE)
#define KRATOS_VARIABLES_CONTAINER_BASE

#include "../includes/variable.h"

namespace Kratos
{
    template<class TDataType>
    class VariablesContainerBase
    {
      public:

        VariablesContainerBase()
        {
        }

	    VariablesContainerBase( const VariablesContainerBase& rOther ) : mData( rOther.mData )
        {
        }

        virtual ~VariablesContainerBase()
        {
        }

        /** Find data respect to the given variable. */
        TDataType& GetValue( const Variable<TDataType>& rVariable )
        {
            return mData[rVariable.getKey()];
        }

        /** Find data respect to the given variable. */
        const TDataType& GetValue( const Variable<TDataType>& rVariable ) const
        {
            typename std::map<int, TDataType>::const_iterator i = mData.find( rVariable.getKey() );

            if( i != mData.end() ) 
            {
                return i->second;
            }
            else 
            {
                return TDataType();
            }
        }

      protected:

        std::map<int, TDataType> mData;
    }; 
}  

#endif 

