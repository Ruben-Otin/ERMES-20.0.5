
#if !defined(KRATOS_DEFINE)
#define KRATOS_DEFINE

namespace Kratos
{
typedef const char* PointerToConstCharType;
  
#if defined(_MSC_VER)
#pragma warning (disable: 4355)
#pragma warning (disable: 4503)
#pragma warning (disable: 4786)
#endif

//Exception handling  
#define KRATOS_TRY try {
#define KRATOS_ERROR(ExceptionType, ErrorPlace, ErrorMessage, MoreInfo) \
      throw Exception(ErrorPlace, ExceptionType(ErrorMessage), MoreInfo);

#define KRATOS_CATCH_WITH_BLOCK(CatchPlace, MoreInfo,Block) \
} \
catch(std::overflow_error& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::underflow_error& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::range_error& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::out_of_range& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::length_error& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::invalid_argument& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::domain_error& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::logic_error& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::runtime_error& e) { Block throw AddMoreExceptionInfo(CatchPlace, e, MoreInfo); } \
catch(std::exception& e) { Block KRATOS_ERROR(std::runtime_error, CatchPlace, e.what(), MoreInfo) } \
catch(...) { Block KRATOS_ERROR(std::runtime_error, CatchPlace, "Unknown error", MoreInfo) } 

#define KRATOS_CATCH_BLOCK_BEGIN class ExceptionBlock{public: void operator()(void){
#define KRATOS_CATCH_BLOCK_END }} exception_block; exception_block(); 
  
#define KRATOS_CATCH(CatchPlace, MoreInfo) \
KRATOS_CATCH_WITH_BLOCK(CatchPlace, MoreInfo,{})
  
#ifdef KRATOS_PRINT_TRACE
  
#define KRATOS_TRACE(A,B) gTrace.Inform(A,B)

#else

#define KRATOS_TRACE(A,B)
#endif

#define KRATOS_MAJOR_VERSION  0

#define KRATOS_MINOR_VERSION  1

#define KRATOS_CORRECTION     0
}  
 
#endif 

