
#include <iostream>
#include <complex>
#include <exception>
#include <cstdlib>
#include <string>

#include "../includes/kernel.h"
#include "../includes/variables.h"


// ERMES version ( printed at start-up and by "ERMES --version" )
#define ERMES_VERSION "20.0.5"

int main( int argc, char *argv[ ], char *envp[ ] )
{
    if( ( argc >= 2 ) && ( std::string( argv[ 1 ] ) == "--version" ) )
    {
        std::cout << "ERMES " << ERMES_VERSION << std::endl;

        return EXIT_SUCCESS;
    }

    if( argc < 2 ) 
    {
        std::cerr << "ERMES " << ERMES_VERSION << std::endl;
        std::cerr << "Usage: ERMES <problem_name>   ( reads <problem_name>.dat and the related GiD files )" << std::endl;

        return EXIT_FAILURE;
    }

    // Any exception escaping the analysis ends ERMES with a readable message and a failure exit status
    // ( previously an uncaught exception aborted the program, e.g. when the .dat file did not exist ).
    try
    {
        Kratos::Kernel* p_kernel = new Kratos::Kernel;	      

        std::cout << "ERMES " << ERMES_VERSION << " analysis started..." << std::endl; 
        std::cout << std::endl;

        p_kernel->Execute( argv[ 1 ] );

        std::cout << std::endl;
	    std::cout << "ERMES analysis finished." << std::endl;
        std::cout << std::endl;

        delete p_kernel;
    }
    catch( std::exception& e )
    {
        std::cout << std::endl << "!!!!!!!!!! ERMES error: " << e.what() << " !!!!!!!!!!" << std::endl;
        std::cerr << "ERMES error: " << e.what() << std::endl;

        return EXIT_FAILURE;
    }
    catch( ... )
    {
        std::cout << std::endl << "!!!!!!!!!! ERMES error: unknown exception !!!!!!!!!!" << std::endl;
        std::cerr << "ERMES error: unknown exception" << std::endl;

        return EXIT_FAILURE;
    }

	return EXIT_SUCCESS;
}
