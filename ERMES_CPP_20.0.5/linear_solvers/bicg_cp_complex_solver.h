
#if !defined(BiCG_CP_COMPLEX_SOLVER_H_INCLUDED)
#define BiCG_CP_COMPLEX_SOLVER_H_INCLUDED

#include "complex_solver.h"

namespace Kratos
{
    //*******************************************************************************
    //* - Bi-Conjugate Gradient  for cold plasma
    //*******************************************************************************    
    class BiCG_CP_ComplexSolver : public ComplexSolver
    {
        public:

            BiCG_CP_ComplexSolver( int NumThreads    ,
			                       int MaxIterations , 
			                       int StepIterations,
                                   double Tolerance  ,
                                   String& Precond   ): ComplexSolver( NumThreads    ,
						                                               MaxIterations ,
						                                               StepIterations,
                                                                       Tolerance     ,
                                                                       Precond       )
            {}
    
            void Solve         ( ComplexMatrixType& A,                           ComplexVectorType& b, ComplexVectorType& x );
            void Solve_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x );
            void Solve_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x );

	    private:

		    std::vector<ComplexType> mDiagonal;

            void InitializePrecond( ComplexMatrixType& rA                            );
            void InitializePrecond( ComplexMatrixType& rA, ComplexMatrixType& rA_aux );

            void ApplyPrecondToVector( ComplexVectorType& rV, ComplexVectorType& rY );

		    void Solve_PrecondBiCG         ( ComplexMatrixType& A,                           ComplexVectorType& b, ComplexVectorType& x );
            void Solve_PrecondBiCG_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x );
            void Solve_PrecondBiCG_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x );
    };
}

#endif 
