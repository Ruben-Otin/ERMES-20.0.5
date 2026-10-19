
#include <iostream>
#include <math.h>
#include <cmath>
#include "bicg_complex_solver.h"

namespace Kratos
{
    //*********************************************************************
    //* - Iterative Bi-Conjugate Gradient solver 
    //*********************************************************************    
    void BiCGComplexSolver::Solve(ComplexMatrixType& A, 
                                  ComplexVectorType& b, 
                                  ComplexVectorType& x)
    {
        if (!IniSolver("Bi-Conjugate Gradient", A, b))
        {
            FinishSolver();
            return;
        }

		ComplexVectorType v_aux(mProblemSize);

		cMultSymm(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

		mResidualNorm  = sqrt( cTwoNorm(v_aux) );
		mResEstimation = mResidualNorm / mBNorm;

		if ( mResEstimation > mTolerance ) 
		{
			if (msPreconditioner == "Diagonal")
			{
				Solve_PrecondBiCG(A,b,x);
			}
			else 
			{
				Solve_BiCG(A,b,x);
			}
		}
		else 
		{
			FinishSolver();
		}
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Biconjugate Gradient solver. Unpreconditioned 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//*********************************************************************
    //* - Iterative Bi-Conjugate Gradient solver. Unpreconditioned 
    //*********************************************************************    
    void BiCGComplexSolver::Solve_BiCG(ComplexMatrixType& A, 
                                       ComplexVectorType& b, 
                                       ComplexVectorType& x)
    {
        ComplexVectorType v_aux(mProblemSize);

        cMultSymm(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

        ComplexVectorType r(v_aux);
        ComplexVectorType p(r);

        ComplexType rnorm = cDot(r,r);
        ComplexType alpha;
        ComplexType gamma;

		double tol_X_bnorm = mTolerance*mBNorm;

        do
        {
            cMultSymm(A, p, v_aux, mNumThreads);

            alpha = rnorm/cDot(v_aux,p);

            // x = (alpha * p) + x 
            cScaleFirstAndAdd(alpha, p, x);

            // r = (-alpha * v_aux) + r;
            cScaleFirstAndAdd(-alpha, v_aux, r);

            gamma = cDot(r,r)/rnorm;

            // p = r + (gamma * p); 
            cScaleSecondAndAdd(r, gamma, p);

            rnorm = gamma * rnorm;

            mResidualNorm = sqrt( cTwoNorm(r) );

            mIterationNumber++;
            mInfoCount++;

            if( mInfoCount == mStepIterations )
            {
                 std::cout << "it = "<< mIterationNumber << "   res = " << ( mResidualNorm / mBNorm );
				 mInfoCount = 0;

				 if( mWriteSolutionEveryStep )
				 {
					 WriteStepResultsInFile( x );
					 std::cout << "     <Write result on file done>";
				 }

				 std::cout << std::endl;
            }

        }while( ( mIterationNumber < mMaxIterations ) && 
                ( mResidualNorm    > tol_X_bnorm    ) && 
                std::isfinite( mResidualNorm )            );

		// Residual estimation
		mResEstimation = mResidualNorm / mBNorm;

		// Calculation of the real residual norm
		cMultSymm(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

		mResidualNorm = sqrt( cTwoNorm(v_aux) );

        FinishSolver();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Bi-Conjugate Gradient solver. Preconditioned 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//*********************************************************************
    //* - Iterative Bi-Conjugate Gradient solver working
    //*********************************************************************    
    void BiCGComplexSolver::Solve_PrecondBiCG(ComplexMatrixType& A, 
                                              ComplexVectorType& b, 
                                              ComplexVectorType& x)
    {
        // Initializing preconditioner. 
		// We will solve ( (P*A*P)*(inv(P)*x) = P*b ) instead of ( A*x = b )     
		InitializePrecond(A);

		// inv(P)*x
		ApplyInversePrecond(x);

		// P*b
		ApplyPrecond(b);

		// Now we solve with (P*b)
		mBNorm = sqrt( cTwoNorm(b) );

        ComplexVectorType v_aux(mProblemSize);
        
		// v_aux = (P*A*P)*x donde x = inv(P)*x
		PrecondMultSym(A, x, v_aux);

        cDeduct(b, v_aux);

        ComplexVectorType r(v_aux);
        ComplexVectorType p(r);

        ComplexType rnorm = cDot(r,r);
        ComplexType alpha;
        ComplexType gamma;

		double tol_X_bnorm = mTolerance*mBNorm;

        do
        {
			// v_aux = (P*A*P)*p
			PrecondMultSym(A, p, v_aux);

            alpha = rnorm/cDot(v_aux,p);

            // x = (alpha * p) + x 
            cScaleFirstAndAdd(alpha, p, x);

            // r = (-alpha * v_aux) + r;
            cScaleFirstAndAdd(-alpha, v_aux, r);

            gamma = cDot(r,r)/rnorm;

            // p = r + (gamma * p); 
            cScaleSecondAndAdd(r, gamma, p);

            rnorm = gamma * rnorm;

            mResidualNorm = sqrt( cTwoNorm(r) );

            mIterationNumber++;
            mInfoCount++;

            if( mInfoCount == mStepIterations )
            {
                 std::cout << "it = " << mIterationNumber << "   res = " << ( mResidualNorm / mBNorm );
				 mInfoCount = 0;

				 if( mWriteSolutionEveryStep )
				 {
					 ApplyPrecond( x );
					 WriteStepResultsInFile( x );
					 ApplyInversePrecond( x );
					 std::cout << "     <Write result on file done>";
				 }

				 std::cout << std::endl;
            }

        }while( ( mIterationNumber < mMaxIterations ) && 
                ( mResidualNorm    > tol_X_bnorm    ) && 
                std::isfinite( mResidualNorm )            );

		// Residual estimation
		mResEstimation = mResidualNorm / mBNorm;
		
		// x = P*y. 
		// We have find y = inv(P)*x
		ApplyPrecond(x);

		// Calculation of the real residual norm
		// inv(P)*b
		ApplyInversePrecond(b);

		cMultSymm(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

		mResidualNorm = sqrt( cTwoNorm(v_aux) );
		mBNorm        = sqrt( cTwoNorm(b)     );

        FinishSolver();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Diagonal Preconditioner
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//**********************************************************************
    //* - Initializing diagonal precondiotioner P = ( 1.00/sqrt( A(i,i) ) )
    //**********************************************************************    
	void BiCGComplexSolver::InitializePrecond(ComplexMatrixType& rA)
    {
        mDiagonal.resize(mProblemSize);
        mTemp.resize(mProblemSize);

		for(int i = 0 ; i < mProblemSize ; ++i) 
		{
			std::complex<double> sqrtDiag = std::sqrt( rA(i,i) );

			if ( std::abs(sqrtDiag) != 0.00 ) mDiagonal[i] = 1.00 / sqrtDiag;
			else                              mDiagonal[i] = 1.00;
		}
    }

	//*********************************************************************
    //* - inv(P)*X
    //*********************************************************************    
	void BiCGComplexSolver::ApplyInversePrecond(ComplexVectorType& rV)
    {    
	     for(int i = 0 ; i < mProblemSize ; ++i) rV[i] /= mDiagonal[i];
    }

	//*********************************************************************
    //* - P*X
    //*********************************************************************    
	void BiCGComplexSolver::ApplyPrecond(ComplexVectorType& rV)
    {
		for(int i = 0 ; i < mProblemSize ; ++i) rV[i] *= mDiagonal[i]; 
    }

	//*********************************************************************
    //* - Y = (P*A*P)*X , donde X = inv(P)*x
    //*********************************************************************  
	void BiCGComplexSolver::PrecondMultSym(ComplexMatrixType& rA, ComplexVectorType& rX, ComplexVectorType& rY)
    {
		 for(int i = 0 ; i < mProblemSize ; ++i) mTemp[i] = rX[i] * mDiagonal[i];

		 cMultSymm(rA,mTemp,rY, mNumThreads);

         ApplyPrecond(rY);
    }
}


