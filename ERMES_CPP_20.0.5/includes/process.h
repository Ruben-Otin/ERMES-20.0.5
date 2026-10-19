
#if !defined(KRATOS_PROCESS)
#define KRATOS_PROCESS

#include <vector>

#include <memory>

namespace Kratos
{
    // Process base class.
    // Process is the entity that controls the flow in Kratos, 
    // giving to it a great flexibility. 
    // Every general action on the flow can be derived as a process. 
    // Actions like solve a equation system, write results, read some data, 
    // output some other can be defined as processes. 
    // Acording to program  structure Kernel is the entity who manages this actions, 
    // just generating the model and executing the processes stored in model, 
    // or maybe can decide by some "inteligent" way the processes are going to be executed.
    class Process
    {
      public:
    
        // Counted pointer of Process 
        typedef std::shared_ptr<Process> Pointer;
        
	    // Vector of Process::Pointer 
	    typedef std::vector<Pointer> ProcessArrayType;
    
        // Constructor
        Process();
    
        virtual ~Process();
             
        // Execute the process. 
        // This is the virtual funcion of every process.
        virtual void Execute();
             
        // Adds a process in a loop process. 
        // This function only makes sense in the case of LoopProcess, 
        // a Process that is made of other processes
        virtual  void AddProcess( Process::Pointer pNewProcess );
        
	  protected:
	   
	    // Array of processes to be executed 
	    Process::ProcessArrayType  mProcessArray;
    };  
}  

#endif 

