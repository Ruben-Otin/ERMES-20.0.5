
#if !defined(CG_COMPLEX_SOLVER_H_INCLUDED)
#define CG_COMPLEX_SOLVER_H_INCLUDED


#include "complex_solver.h"


namespace Kratos
{

    //*********************************************************************
    //* - Conjugate Gradient 
    //*********************************************************************    
    class CGComplexSolver : public ComplexSolver
    {
     public:

         CGComplexSolver(int NumThreads,
			             int MaxIterations, 
			             int StepIterations,
                         double Tolerance,
                         String& Precond): ComplexSolver(NumThreads,
						                                 MaxIterations,
						                                 StepIterations,
                                                         Tolerance,
                                                         Precond){};
    
         void Solve(ComplexMatrixType& A, 
                    ComplexVectorType& b, 
                    ComplexVectorType& x);

	 private:

		 ComplexVectorType        mTemp;
		 std::vector<ComplexType> mDiagonal;

		 void InitializePrecond  (ComplexMatrixType& rA);
		 void PrecondMultSym     (ComplexMatrixType& rA, ComplexVectorType& rV, ComplexVectorType& rY);
		 void PrecondConjgMultSym(ComplexMatrixType& rA, ComplexVectorType& rX, ComplexVectorType& rY);

		 void ApplyInversePrecond(ComplexVectorType& rV);
		 void ApplyPrecond       (ComplexVectorType& rV);

		 void Solve_CG       (ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x);
		 void Solve_PrecondCG(ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x);
    };

}

#endif 
