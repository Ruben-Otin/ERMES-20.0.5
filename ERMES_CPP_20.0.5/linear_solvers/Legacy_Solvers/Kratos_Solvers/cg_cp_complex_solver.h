
#if !defined(CG_CP_COMPLEX_SOLVER_H_INCLUDED)
#define CG_CP_COMPLEX_SOLVER_H_INCLUDED

#include "complex_solver.h"

namespace Kratos
{
    //**********************************************************************************
    //* - Conjugate Gradient for cold plasma
    //**********************************************************************************    
    class CG_CP_ComplexSolver : public ComplexSolver
    {
        public:

            CG_CP_ComplexSolver( int NumThreads    , 
                                 int MaxIterations , 
                                 int StepIterations, 
                                 double Tolerance  , 
                                 String& Precond   ) : ComplexSolver( NumThreads    , 
                                                                      MaxIterations , 
                                                                      StepIterations, 
                                                                      Tolerance     , 
                                                                      Precond       )
            {};
    
            void Solve( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x );

	    private:

            std::vector<ComplexType> mDiagonal;

            void InitializePrecond    ( ComplexMatrixType& rA );

            void ApplyPrecondToVector ( ComplexVectorType& rV, ComplexVectorType& rY );

		    void Solve_PrecondCGS     ( ComplexMatrixType& A , ComplexVectorType& b, ComplexVectorType& x );
            void Solve_PrecondBiCGSTAB( ComplexMatrixType& A , ComplexVectorType& b, ComplexVectorType& x );
    };
}

#endif 
