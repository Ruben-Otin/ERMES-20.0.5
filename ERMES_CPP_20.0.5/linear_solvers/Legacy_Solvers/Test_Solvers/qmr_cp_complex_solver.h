
#if !defined(QMR_CP_COMPLEX_SOLVER_H_INCLUDED)
#define QMR_CP_COMPLEX_SOLVER_H_INCLUDED

#include "complex_solver.h"

namespace Kratos
{
	//*********************************************************************
	//* - Quasi-Minimal Residual for cold plasma
	//*********************************************************************    
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
            {};

		    void Solve( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x );

	    private:

		    ComplexVectorType        mTemp;
		    std::vector<ComplexType> mDiagonal;

		    void InitializePrecond(ComplexMatrixType& rA);
		    void PrecondMultFull  (ComplexMatrixType& rA, ComplexVectorType& rV, ComplexVectorType& rY);

		    void ApplyInversePrecond(ComplexVectorType& rV);
		    void ApplyPrecond       (ComplexVectorType& rV);

		    void Solve_QMR       (ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x);
		    void Solve_PrecondQMR(ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x);
	};
}

#endif 
