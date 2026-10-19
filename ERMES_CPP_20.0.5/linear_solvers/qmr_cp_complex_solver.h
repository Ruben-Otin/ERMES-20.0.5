
#if !defined(QMR_CP_COMPLEX_SOLVER_H_INCLUDED)
#define QMR_CP_COMPLEX_SOLVER_H_INCLUDED

#include "complex_solver.h"

namespace Kratos
{
	//*********************************************************************************
	//* - Quasi-Minimal Residual for cold plasma
	//*********************************************************************************    
	class QMR_CP_ComplexSolver : public ComplexSolver
	{
	    public:

		    QMR_CP_ComplexSolver( int NumThreads    , 
			                      int MaxIterations , 
			                      int StepIterations, 
			                      double Tolerance  , 
			                      String& Precond   ) : ComplexSolver( NumThreads    , 
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

            void ApplyPrecond_M1( ComplexVectorType& rV, ComplexVectorType& rY );
            void ApplyPrecond_M2( ComplexVectorType& rV, ComplexVectorType& rY );
		    
		    void Solve_PrecondQMR         ( ComplexMatrixType& A,                           ComplexVectorType& b, ComplexVectorType& x );
            void Solve_PrecondQMR_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x );
            void Solve_PrecondQMR_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x );
	};
}

#endif 
