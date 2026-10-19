
#if !defined(KRATOS_INPUT_C_SCANNER)
#define KRATOS_INPUT_C_SCANNER

#include <malloc.h>
#include <FlexLexer.h>

namespace Kratos
{
    class InputCScanner : public yyFlexLexer
    {
      public:
        
        InputCScanner( std::istream* NewInput = 0, std::ostream* NewOutput = 0 ) 
        : yyFlexLexer( NewInput, NewOutput ), mNumberOfLines( 1 ), mNumberOfNodes( 1 ), mNumberOfElements( 1 ){}
      
        virtual ~InputCScanner(){}
        
        int yylex( yy_InputCParser_stype &yylval );
        
        int& rNumberOfNodes   (){ return mNumberOfNodes;   }
        int& rNumberOfLines   (){ return mNumberOfLines;   }
        int& rNumberOfElements(){ return mNumberOfElements;}
 
        char* GetYYText(){ return yytext; }
           
      private:
          
        int mNumberOfLines;
        int mNumberOfNodes;
        int mNumberOfElements;        
    }; 
}  

#endif 

