
#if !defined(KRATOS_STRING)
#define KRATOS_STRING

#include <string>
#include <sstream>

#include "../includes/define.h"  

namespace Kratos
{
    // String class used in Kratos.
	// This class derived from STL string for portability porpose and add some
	// functionality to it for assigning numeric values to it.
	// It has constructor and operator necessary for convert to or from C
	// null terminated character array.*/
    class String : public std::string
    {
      public:
    
        String( void )
        {
        }
        
        // This constructor initialize string with a C const char null terminated array
        String( PointerToConstCharType pUserString ) : std::basic_string<char>( pUserString )
        {
        }
        
        virtual ~String()
        {
        }
        
        // This operator do casting to the const char*
        operator PointerToConstCharType() const 
        {
            return c_str();
        }
        
        // This assignment operator put given Value into the string.
        // Previous data in string will be erased. 
        template<class TDataType>
        String& operator = ( const TDataType& Value )
        {
            std::stringstream Buffer;
            Buffer << Value;
            assign( Buffer.str() );
            return *this;
        }
    
	    template<class TDataType>
	    String& operator << ( const TDataType& Value )
        {
	        std::stringstream buffer;
	        buffer << *this;
	        buffer << Value;
	        assign( buffer.str() );
	        return *this;
        }
    }; 
}  

#endif 