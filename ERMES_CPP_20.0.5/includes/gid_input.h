
#if !defined(KRATOS_GID_INPUT)
#define KRATOS_GID_INPUT

#include "../includes/input_c_parser.h"
#include "../includes/input.h"

namespace Kratos
{
    class GidInput : public Input
    {
      public:
    
        typedef std::shared_ptr<GidInput> Pointer;
        
        GidInput( const String& InputFileName ) : Input( InputFileName )
        {
        }
        
        virtual ~GidInput()
        {
        }
        
        void Read( Kernel* pKernel )
        {
            InputCParser parser( pKernel, mpInputFile, 0 );
            parser.yyparse();
        } 
    }; 
} 

#endif 

