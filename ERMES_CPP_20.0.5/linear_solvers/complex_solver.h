

#if !defined(COMPLEX_SOLVER_H_INCLUDED)
#define  COMPLEX_SOLVER_H_INCLUDED

#include <complex>
#include <time.h>

#include <memory>

#include "../includes/vector.h"
#include "../includes/csr_matrix.h"
#include "../includes/kratos_string.h"
#include "complex_matrix_operations.h"

#include "../ERMES/Utils/MatrixFileWritting.h"

namespace Kratos
{
    //*********************************************************************
    //* - Base class for complex solvers. 
    //*********************************************************************    
    class ComplexSolver 
    {
        public:

            typedef std::shared_ptr<ComplexSolver> Pointer;

            typedef std::complex<double> ComplexType;

            typedef CSRMatrix<ComplexType> ComplexMatrixType;

            typedef Vector<ComplexType> ComplexVectorType;

            ComplexSolver(unsigned int NumThreads, unsigned int MaxIterations, unsigned int StepIterations, double Tolerance, String& Precond);
    
            virtual void Solve         ( ComplexMatrixType& A,                           ComplexVectorType& b, ComplexVectorType& x ){};
            virtual void Solve_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x ){};
            virtual void Solve_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x ){};

            virtual ~ComplexSolver(){}

            virtual void PrintInfo();

            bool IniSolver( const char* SolverName, ComplexMatrixType& A, ComplexVectorType& b );

            void SetWriteStepResult( bool ToWriteStepResult );

            void FinishSolver();

        protected:

            double mTolerance;
            double mResidualNorm;
            double mBNorm;
            double mResEstimation;

            unsigned int mNumThreads;
            unsigned int mMaxIterations;
            unsigned int mStepIterations;
            unsigned int mIterationNumber;
            unsigned int mProblemSize;

            unsigned short int mInfoCount;

            bool mWriteSolutionEveryStep;

            time_t mcStart;
            time_t mcEnd;

            String msPreconditioner;
    };
}

#endif 