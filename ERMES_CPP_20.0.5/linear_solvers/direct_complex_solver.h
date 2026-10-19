
#if !defined(DIRECT_COMPLEX_SOLVER_H_INCLUDED)
#define DIRECT_COMPLEX_SOLVER_H_INCLUDED

#include <vector>
#include <string>

#include "complex_solver.h"

namespace Kratos
{
    // Approximate minimum degree ordering of the pattern of A + A^T ( A in CSC ). Q[k] = k-th pivot.
    void ApproximateMinimumDegree( int N, const std::vector<int>& Ap, const std::vector<int>& Ai, std::vector<int>& Q );

    //*****************************************************************************************************
    //* - Reference implementation, not used by DirectComplexSolver since the multifrontal solver: kept for
    //*   the unit tests ( singular matrix detection ) and as a simple, serial, robust LU for small matrices.
    //*
    //* - Sparse LU factorization of a general complex matrix stored in compressed sparse columns ( CSC ).
    //*
    //*   Ordering : approximate minimum degree on the pattern of A + A^T ( quotient graph with element
    //*              absorption ), applied symmetrically to columns and, preferentially, to rows.
    //*   Numeric  : left-looking LU ( Gilbert-Peierls ) with threshold partial pivoting. The diagonal of the
    //*              permuted matrix is kept as pivot whenever |a_kk| >= PivotTol * max_i |a_ik|, which keeps
    //*              the fill of the minimum degree ordering; otherwise the largest entry is used.
    //*
    //*   P * A * Q = L * U ,  L unit lower triangular, U upper triangular.
    //*****************************************************************************************************
    class ComplexSparseLU
    {
        public:

            typedef std::complex<double> ComplexType;

            ComplexSparseLU() : mN( 0 ) {}

            // Takes ownership of the CSC arrays ( they are swapped, the inputs are left empty )
            void SetMatrix( int n, std::vector<int>& Ap, std::vector<int>& Ai, std::vector<ComplexType>& Ax );

            // Ordering + numeric factorization. Returns false ( and an explanation ) if A is singular.
            bool Factorize( double PivotTol, std::string& ErrorMsg );

            // b <- inv(A) * b
            void Solve( std::vector<ComplexType>& b );

            // r <- b - A*x
            void Residual( const std::vector<ComplexType>& b, const std::vector<ComplexType>& x, std::vector<ComplexType>& r );

            long long NnzA() const { return mAp.empty() ? 0 : mAp[ mN ]; }
            long long NnzL() const { return mLp.empty() ? 0 : mLp[ mN ]; }
            long long NnzU() const { return mUp.empty() ? 0 : mUp[ mN ]; }

            int NumOffDiagonalPivots() const { return mOffDiagPivots; }

            void Clear();

            // Physical memory of the machine in bytes ( 0 if unknown )
            static double PhysicalMemoryBytes();

        private:

            int mN;

            int mOffDiagPivots;

            // A in CSC
            std::vector<int>         mAp, mAi;
            std::vector<ComplexType> mAx;

            // Column permutation ( fill-reducing ) and inverse row permutation ( pivoting )
            std::vector<int> mQ, mPinv;

            // L ( unit diagonal stored first in each column ) and U ( diagonal stored last ) in CSC
            std::vector<int>         mLp, mLi, mUp, mUi;
            std::vector<ComplexType> mLx, mUx;

            void MinimumDegreeOrdering();

            long long SymbolicAnalysis( std::vector<double>& ColWork );

            int Reach( int col, std::vector<int>& xi, std::vector<int>& mark, std::vector<int>& pstack, int stamp );
    };

    //*****************************************************************************************************
    //* - Internal direct solver ( sparse LU ) for ERMES complex systems.
    //*
    //*   Works with every matrix storage used by ERMES:
    //*     Solve          : A complex symmetric ( upper triangle stored ), or full when FullStorage = true
    //*                      ( cold plasma 'Full_matrix' mode )
    //*     Solve_HermSymm : A Hermitian ( upper ) + A_aux complex symmetric ( upper )
    //*     Solve_HermFull : A Hermitian ( upper ) + A_aux full
    //*
    //*   Steps ( DirectComplexSolver::SolveCombined ):
    //*     1) conversion of the ERMES row storage to a full CSC matrix;
    //*     2) symmetric equilibration ( Ruiz scaling );
    //*     3) supernodal multifrontal factorization ( ComplexMultifrontal ), LDL^T for complex symmetric
    //*        matrices, LU otherwise, OpenMP parallel;
    //*     4) GMRES preconditioned by the factorization ( 1-2 iterations when it is exact );
    //*     5) if LDL^T fails or GMRES does not converge, refactorization as a general matrix ( LU ).
    //*   Memory grows quickly with the problem size; the predicted memory is printed and checked before
    //*   factorizing. For very large problems use the iterative solvers or the External_solver ( PETSc/MUMPS ).
    //*****************************************************************************************************
    class DirectComplexSolver : public ComplexSolver
    {
        public:

            DirectComplexSolver( int NumThreads    ,
                                 int MaxIterations ,
                                 int StepIterations,
                                 double Tolerance  ,
                                 String& Precond   ,
                                 bool FullStorage  ) : ComplexSolver( NumThreads, MaxIterations, StepIterations, Tolerance, Precond ),
                                                       mFullStorage( FullStorage )
            {}

            void Solve         ( ComplexMatrixType& A,                           ComplexVectorType& b, ComplexVectorType& x );
            void Solve_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x );
            void Solve_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x );

            // Storage of an ERMES row matrix when converted to CSC
            enum StorageType { FULL = 0, SYMMETRIC_UPPER = 1, HERMITIAN_UPPER = 2 };

        private:

            bool mFullStorage;

            void SolveCombined( ComplexMatrixType& A, StorageType sA, ComplexMatrixType* pA_aux, StorageType sAux,
                                ComplexVectorType& b, ComplexVectorType& x, const char* SolverName );
    };

    //*****************************************************************************************************
    //* - Converts one or two ERMES row matrices ( summed ) to a full CSC matrix ( exposed for unit tests ).
    //*****************************************************************************************************
    void BuildFullCSC( CSRMatrix< std::complex<double> >& A, DirectComplexSolver::StorageType sA,
                       CSRMatrix< std::complex<double> >* pA_aux, DirectComplexSolver::StorageType sAux,
                       int n, std::vector<int>& Ap, std::vector<int>& Ai, std::vector< std::complex<double> >& Ax );
}

#endif
