
#if !defined(KRATOS_VARIABLE)
#define KRATOS_VARIABLE

#include "../includes/kratos_string.h"
#include <typeinfo>

namespace Kratos
{
    template<class TDataType>
    class Variable
    {
      public:

        Variable( const String& Name, unsigned short int Key ) : mName( Name ), mKey( Key )
        {
        }
        
        Variable( const String& Name, unsigned short int Key, const Variable<TDataType>& rTimeDerivative ) : mName( Name ), mKey( Key )
        {
        }
        
        Variable( const Variable<TDataType>& rOtherVariable ) : mName( rOtherVariable.mName ), mKey( rOtherVariable.mKey )
        {
        }
        
        virtual ~Variable()
        {
        }
        
        unsigned short int getKey() const
        {
            return mKey;
        }
        
        Variable<TDataType>& operator=( const Variable<TDataType>& rOtherVariable )
        {
            mName = rOtherVariable.mName;
            mKey  = rOtherVariable.mKey;
            return *this;
        }
        
        friend bool operator==( const Variable<TDataType>& rFirstVariable, const Variable<TDataType>& rSecondVariable )
        {
            return ( rFirstVariable.mName == rSecondVariable.mName ) && ( rFirstVariable.mKey  == rSecondVariable.mKey );
        }
        
	    static const Variable<TDataType> NONE()
        {
            return Variable<TDataType>( "NONE", 0 );
        }
        
        const String& Name() const
        {
            return mName;
        }
        
        const std::type_info& getType()
        {
            return typeid( TDataType );
        }

      private:

        String mName;
        
        unsigned short int mKey;
        
        Variable();
    };
}  

#endif 

