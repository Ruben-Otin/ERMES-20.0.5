
#if !defined(KRATOS_INPUT)
#define KRATOS_INPUT

#include <memory>

namespace Kratos
{
    class Kernel;
    
    class Input
    {
      public:
        
        typedef std::shared_ptr<Input> Pointer;
    
        Input( const String& InputFileName )
        {
            mpInputFile = new std::ifstream( InputFileName );
          
            if( mpInputFile->fail() )
            {
                delete mpInputFile;
                String buffer;
                buffer << "error opening data file: " <<  InputFileName ;
                throw Exception( "Input::SetInput", std::invalid_argument( buffer ) );
            }        
        }
    
        void SetNewInput( const String& InputFileName )
        {
            delete mpInputFile;
            mpInputFile = new std::ifstream( InputFileName );
        }
    
        virtual ~Input()
        {
            delete mpInputFile;
        }
             
        virtual void Read( Kernel* pKernel ){}
        
      protected:
         
        std::ifstream* mpInputFile;        
    }; 
}  

#endif 

