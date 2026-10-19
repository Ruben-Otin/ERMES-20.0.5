
#if !defined(KRATOS_EXCEPTION)
#define KRATOS_EXCEPTION

#include <stdexcept>

#include "../includes/kratos_string.h"

namespace Kratos
{
    template<class TExceptionType> TExceptionType Exception( const String& rWhere, const TExceptionType& rException )
    {
        return rException;
    }
    
    template<class TExceptionType, class TInfoType> TExceptionType Exception( const String& rWhere, const TExceptionType& rException, const TInfoType& rInfo )
    {
        String buffer;
        buffer << "in " << rInfo << " " << rWhere << "\nwith subject    :  " << rException.what();
        return TExceptionType(buffer);
    }
    
    template<class TExceptionType, class TInfoType> TExceptionType AddMoreExceptionInfo( const String& rWhere, const TExceptionType& rException, const TInfoType& rInfo )
    {
        String buffer;
        buffer << rException.what() << "\nwhile executing : " << rInfo << " " << rWhere;
        return TExceptionType(buffer);
    }
}  

#endif 

