
#include "../includes/define.h"
#include "../includes/process.h"

namespace Kratos
{
    Process::Process()
    {
        KRATOS_TRACE( "Process::Process()", "Process Created" );  
    }
    
    Process::~Process()
    {
        KRATOS_TRACE( "Process::~Process()", "Process destroyed" );
    }
    
    void  Process::Execute()
    {
        KRATOS_TRACE( "Process::Execute()", "Process Executed" );  
    }
    
    void Process::AddProcess( Process::Pointer  pNewProcess )
    {
        KRATOS_TRACE( "Process::AddProcess( Process::Pointer  pNewProcess )", "Process added" );

        mProcessArray.push_back( pNewProcess );
    }
} 
