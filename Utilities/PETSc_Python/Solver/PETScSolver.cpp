static char help[] = "Solves a complex sparse linear system in parallel with KSP.\n\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// - PETSc solver C++ source code. 
//
// - Check the "README.txt" file provided in "/External_Solvers/PETSc" folder for instructions about how to install PETSc.
//
// - Once PETSc has been installed, compile this file.
//
// - Type "make PETScSolver" on this same folder path to compile "PETScSolver.cpp".
//
// - The executable can be called from the Python scripts located on the "/External_Solvers/PETSc/Scripts" folder.
//
// - The executable can also be called directly from a Linux bash script (see an example in file "/Documents/Batch_Scripts/PETSc_batch.sh").
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <petscksp.h>
#include <iostream>
#include <fstream>
#include <complex>

int main( int argc, char **args )
{
	Vec            x, b;       // Approx solution, RHS 
	Mat            A;          // Linear system matrix 
	KSP            ksp;        // Linear solver context 
	PetscMPIInt    rank, size; // PETSc mpi rank and size object 
	PetscViewer    viewer;     // PETSc object for visualization and reading 
	PetscErrorCode ierr;       // Error code 

	// Initialise PETSc, includes MPI initialisation
	ierr = PetscInitialize( &argc, &args, (char*)0, help ); if ( ierr ) return ierr;

	MPI_Comm_rank( PETSC_COMM_WORLD, &rank );
	MPI_Comm_size( PETSC_COMM_WORLD, &size );

	// Problem folder full path ( ending with "/" )
	char Problem_Full_Path[ PETSC_MAX_PATH_LEN ]; PetscBool path_given;

	ierr = PetscOptionsGetString( NULL, NULL, "-my_folder", Problem_Full_Path, sizeof( Problem_Full_Path ), &path_given ); CHKERRQ( ierr );

	// If no path is given, give it here
	if ( !path_given ) strcpy( Problem_Full_Path, "" );

	// File names 
	const char Matrix_A_File_Name[] = "Matrix_A_PETSc.dat";
	const char Vector_B_File_Name[] = "Vector_B_PETSc.dat";
	const char Vector_X_File_Name[] = "Vector_X_PETSc.dat";
	
	// Matrices and vector full paths
	char Matrix_A_Full_Path[ strlen( Problem_Full_Path ) + strlen ( Matrix_A_File_Name ) + 1 ];
	char Vector_B_Full_Path[ strlen( Problem_Full_Path ) + strlen ( Vector_B_File_Name ) + 1 ];
	char Vector_X_Full_Path[ strlen( Problem_Full_Path ) + strlen ( Vector_X_File_Name ) + 1 ];
	
	strcpy( Matrix_A_Full_Path, Problem_Full_Path ); strcat( Matrix_A_Full_Path, Matrix_A_File_Name );
	strcpy( Vector_B_Full_Path, Problem_Full_Path ); strcat( Vector_B_Full_Path, Vector_B_File_Name );
	strcpy( Vector_X_Full_Path, Problem_Full_Path ); strcat( Vector_X_Full_Path, Vector_X_File_Name );
	
	// Reading the matrix that defines the linear system Ax = b
	ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD, Matrix_A_Full_Path, FILE_MODE_READ, &viewer ); CHKERRQ( ierr );
	ierr = MatCreate            ( PETSC_COMM_WORLD, &A                                          ); CHKERRQ( ierr );
	ierr = MatLoad              ( A, viewer                                                     ); CHKERRQ( ierr );
	ierr = PetscViewerDestroy   ( &viewer                                                       ); CHKERRQ( ierr );

	// Reading the right-hand-side vector that defines the linear system Ax = b
	ierr = PetscViewerBinaryOpen( PETSC_COMM_WORLD, Vector_B_Full_Path, FILE_MODE_READ, &viewer ); CHKERRQ( ierr );
	ierr = VecCreate            ( PETSC_COMM_WORLD, &b                                          ); CHKERRQ( ierr );
	ierr = VecLoad              ( b, viewer                                                     ); CHKERRQ( ierr );
	ierr = PetscViewerDestroy   ( &viewer                                                       ); CHKERRQ( ierr );
	
	// Defining X vector solution ( duplicating B vector size )
	ierr = VecDuplicate( b, &x ); CHKERRQ( ierr );
	
	// Create linear solver and set options
	ierr = KSPCreate        ( PETSC_COMM_WORLD, &ksp ); CHKERRQ( ierr );
	ierr = KSPSetOperators  ( ksp, A, A              ); CHKERRQ( ierr );
	ierr = KSPSetFromOptions( ksp                    ); CHKERRQ( ierr );
	
	// Solve linear system
	ierr = PetscPrintf( PETSC_COMM_WORLD, "Solving linear system A*X = B...\n\n" ); CHKERRQ( ierr );
	ierr = KSPSolve   ( ksp, b, x                                                ); CHKERRQ( ierr );
	ierr = PetscPrintf( PETSC_COMM_WORLD, "\nLinear system solved.\n"            ); CHKERRQ( ierr );

	// Save solution in file
	ierr = PetscViewerBinaryOpen         ( PETSC_COMM_WORLD, Vector_X_Full_Path, FILE_MODE_WRITE, &viewer ); CHKERRQ( ierr );
	ierr = PetscViewerBinarySetSkipHeader( viewer, PETSC_TRUE                                             ); CHKERRQ( ierr );
	ierr = VecView                       ( x, viewer                                                      ); CHKERRQ( ierr );
	ierr = PetscViewerDestroy            ( &viewer                                                        ); CHKERRQ( ierr );
	
	// Clean PETSc objects 
	ierr = KSPDestroy( &ksp ); CHKERRQ( ierr );
	ierr = VecDestroy( &x   ); CHKERRQ( ierr );
	ierr = VecDestroy( &b   ); CHKERRQ( ierr ); 
	ierr = MatDestroy( &A   ); CHKERRQ( ierr );
	
	// Finalize PETSc
	ierr = PetscFinalize();

	return ierr;
}