
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include "complex_solver.h"

namespace Kratos
{
    //*********************************************************************
    //* - Initializing private properties
    //*********************************************************************    
    ComplexSolver::ComplexSolver( unsigned int NumThreads, unsigned int MaxIterations, unsigned int StepIterations, double Tolerance, String& Precond )
    {
        mNumThreads      = NumThreads;
		mMaxIterations   = MaxIterations;
		mStepIterations  = StepIterations;
        mTolerance       = Tolerance; 
        mResidualNorm    = 0.0;
		mResEstimation   = 0.0;
        msPreconditioner = Precond;

		mWriteSolutionEveryStep = false;
    }

    //*********************************************************************
    //* - Show info
    //*********************************************************************    
    void ComplexSolver::PrintInfo()
    {
        // mResEstimation : relative residual estimated by the iterative recurrence ( used as stopping criterion )
        // mResidualNorm  : true residual ||b - A*x|| recomputed after the iterations
        const double TrueRelResidual = ( mBNorm != 0.00 ) ? ( mResidualNorm / mBNorm ) : 0.00;

        const bool Breakdown = !std::isfinite( mResEstimation ) || !std::isfinite( mResidualNorm );

        if (mBNorm == 0.00)
        {
            std::cout << "!!!!!!!!!!!!!!! ||b|| = 0 !!!!!!!!!!!!!!!" << std::endl;
        }
        else if ( Breakdown )
        {
            // NaN / Inf: previous versions reported this case as "CONVERGENCE ACHIEVED" because NaN > tol is false
            std::cout << "!!!!!!!!!!! SOLVER BREAKDOWN ( NaN / Inf ) !!!!!!!!!!!" << std::endl;
			std::cout << std::endl;
            std::cout << "||b||             = " << mBNorm          << std::endl;
            std::cout << "||Ax-b||          = " << mResidualNorm   << std::endl;
			std::cout << "Solver estimation = " << mResEstimation  << std::endl;
        }
        else if (mResEstimation > mTolerance)
        {
            std::cout << "!!!!!!!!!!! NO CONVERGENCE !!!!!!!!!!!" << std::endl;
			std::cout << std::endl;
            std::cout << "||b||             = " << mBNorm          << std::endl;
            std::cout << "||Ax-b||          = " << mResidualNorm   << std::endl;
            std::cout << "||Ax-b||/||b||    = " << TrueRelResidual << std::endl;
			std::cout << "Solver estimation = " << mResEstimation  << std::endl;
        }
        else
        {
            std::cout << "--------> CONVERGENCE ACHIEVED <-------- " << std::endl;  
            std::cout << std::endl;
            std::cout << "||b||             = " << mBNorm          << std::endl;
            std::cout << "||Ax-b||          = " << mResidualNorm   << std::endl;
            std::cout << "||Ax-b||/||b||    = " << TrueRelResidual << std::endl;
			std::cout << "Solver estimation = " << mResEstimation  << std::endl;

            // Residual gap: the recursively updated residual can drift away from the true one ( rounding in BiCG/QMR ).
            if( TrueRelResidual > 10.0 * mTolerance )
            {
                std::cout << std::endl;
                std::cout << "Warning: true residual ||Ax-b||/||b|| is more than 10x the tolerance ( residual gap )." << std::endl;
                std::cout << "         The solution may be less accurate than the solver estimation suggests."      << std::endl;
            }
        }

        std::cout << std::endl;
        std::cout << "Max. iterations: " << mMaxIterations << std::endl;
        std::cout << "Tolerance      : " << mTolerance     << std::endl;
        std::cout << "Problem size   : " << mProblemSize   << std::endl;
        std::cout << std::endl;
        std::cout << "Iterations     : " << mIterationNumber                      << std::endl;
        std::cout << "Time spent     : " << difftime(mcEnd,mcStart) << " seconds" << std::endl;
        std::cout << std::endl; 
    }

    //***********************************************************
    //* - Initialization of iterative solver
    //***********************************************************
    bool ComplexSolver::IniSolver( const char* SolverName, ComplexMatrixType& A, ComplexVectorType& b )
    {
        // The direct solver ( "Direct sparse LU..." ) uses no preconditioner
        const bool Direct = ( std::string( SolverName ).find( "Direct" ) != std::string::npos );

        std::cout << ( Direct ? "Direct solver   : " : "Iterative solver: " ) << SolverName << std::endl;

        if( !Direct ) std::cout << "Preconditioner  : " << msPreconditioner << std::endl;
        std::cout << std::endl;

        mProblemSize     = b.size();
        mIterationNumber = 0;
        mInfoCount       = 0;
 
        time(&mcStart);

        mBNorm = sqrt(cTwoNorm(b));

        if (mBNorm == 0.00) 
        {
            time(&mcEnd);
            return false;
        }

        return true;
    }

	//***********************************************************
    //* - Write the x_vector every step iteration
    //***********************************************************
    void ComplexSolver::SetWriteStepResult( bool ToWriteStepResult )
    {
		mWriteSolutionEveryStep = ToWriteStepResult;
    }

    //***********************************************************
    //* - Finishing iterative solver
    //***********************************************************
    void ComplexSolver::FinishSolver()
    {
        time(&mcEnd);
        std::cout << std::endl;
        PrintInfo();

        // A NaN/Inf solution is useless: stop with a failure status instead of writing garbage results
        if( ( mBNorm != 0.00 ) && ( !std::isfinite( mResEstimation ) || !std::isfinite( mResidualNorm ) ) )
        {
            std::cout << "ERMES analysis finished." << std::endl;
            std::cerr << "ERMES error: iterative solver breakdown ( NaN / Inf )" << std::endl;
            exit( EXIT_FAILURE );
        }
    }

}