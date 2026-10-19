
#include <iostream>
#include <math.h>
#include "cg_complex_solver.h"


namespace Kratos
{
    //*********************************************************************
    //* - Iterative Conjugate Gradient solver 
    //*********************************************************************    
    void CGComplexSolver::Solve(ComplexMatrixType& A, 
                                ComplexVectorType& b, 
                                ComplexVectorType& x)
    {
        if (!IniSolver("Conjugate Gradient", A, b))
        {
            FinishSolver();
            return;
        }

		ComplexVectorType v_aux(mProblemSize);
		cMultSym(A, x, v_aux, mNumThreads);
        cDeduct(b, v_aux);

		mResidualNorm  = sqrt( cTwoNorm(v_aux) );
		mResEstimation = mResidualNorm / mBNorm;

		if ( mResEstimation > mTolerance ) 
		{
			if (msPreconditioner == "Diagonal")
			{
				Solve_PrecondCG(A,b,x);
			}
			else 
			{
				Solve_CG(A,b,x);
			}
		}
		else 
		{
			FinishSolver();
		}
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Conjugate Gradient solver. Unpreconditioned 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//*********************************************************************
    //* - Iterative Conjugate Gradient solver. Unpreconditioned 
    //*********************************************************************    
    void CGComplexSolver::Solve_CG(ComplexMatrixType& A, 
                                   ComplexVectorType& b, 
                                   ComplexVectorType& x)
    {
        ComplexVectorType r(b);

        ComplexVectorType v_aux(mProblemSize);

        cConjgMultSym(A, r, v_aux, mNumThreads);
		
        double beta = 1.00 / cTwoNorm(v_aux);

        ComplexVectorType p(mProblemSize);

        cScale(beta, v_aux, p);

        double alpha = 0.00;

		double tol_X_bnorm = mTolerance*mBNorm;

        do
        {
            cMultSym(A, p, v_aux, mNumThreads);
           
            alpha = 1.00 / cTwoNorm(v_aux);

            cScaleFirstAndAdd( alpha, p    , x);

            cScaleFirstAndAdd(-alpha, v_aux, r);

            cConjgMultSym(A, r, v_aux, mNumThreads);
		
            beta = 1.00 / cTwoNorm(v_aux);

            cScaleFirstAndAdd(beta, v_aux, p);

            mResidualNorm = sqrt( cTwoNorm(r) );

            mIterationNumber++;
            mInfoCount++;

            if ( mInfoCount == mStepIterations )
            {
                 std::cout<<"it = "<<mIterationNumber<<"   res = "<< (mResidualNorm / mBNorm);
				 mInfoCount = 0;

				 if (mWriteResultEveryStep)
				 {
					 WriteStepResultsInFile(x);
					 std::cout<<"     <Write result in file done>";
				 }

				 std::cout << std::endl;
            }

        }while( ( mIterationNumber < mMaxIterations ) && 
                ( mResidualNorm    > tol_X_bnorm    )    );

		// Residual estimation
		mResEstimation = mResidualNorm / mBNorm;

		// Calculation of the real residual norm
		cMultSym(A, x, v_aux, mNumThreads);
        cDeduct(b, v_aux);
		mResidualNorm = sqrt( cTwoNorm(v_aux) );

        FinishSolver();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Conjugate Gradient solver. Preconditioned 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//*********************************************************************
    //* - Iterative Conjugate Gradient solver working
    //*********************************************************************    
    void CGComplexSolver::Solve_PrecondCG(ComplexMatrixType& A, 
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

        ComplexVectorType r(b);

        ComplexVectorType v_aux(mProblemSize);

        PrecondConjgMultSym(A, r, v_aux);
		
        double beta = 1.00 / cTwoNorm(v_aux);

        ComplexVectorType p(mProblemSize);

        cScale(beta, v_aux, p);

        double alpha = 0.00;

		double tol_X_bnorm = mTolerance*mBNorm;

        do
        {
		    // v_aux = (P*A*P)*p
			PrecondMultSym(A, p, v_aux);
           
            alpha = 1.00 / cTwoNorm(v_aux);

            cScaleFirstAndAdd( alpha, p    , x);

            cScaleFirstAndAdd(-alpha, v_aux, r);

            PrecondConjgMultSym(A, r, v_aux);
		
            beta = 1.00 / cTwoNorm(v_aux);

            cScaleFirstAndAdd(beta, v_aux, p);

            mResidualNorm = sqrt( cTwoNorm(r) );

            mIterationNumber++;
            mInfoCount++;

            if ( mInfoCount == mStepIterations )
            {
                 std::cout<<"it = "<<mIterationNumber<<"   res = "<< (mResidualNorm / mBNorm);
				 mInfoCount = 0;

				 if (mWriteResultEveryStep)
				 {
					 ApplyPrecond(x);
					 WriteStepResultsInFile(x);
					 ApplyInversePrecond(x);
					 std::cout<<"     <Write result in file done>";
				 }

				 std::cout << std::endl;
            }

        }while( ( mIterationNumber < mMaxIterations ) && 
                ( mResidualNorm    > tol_X_bnorm    )    );

		// Residual estimation
		mResEstimation = mResidualNorm / mBNorm;
		
		// x = P*y. 
		// We have find y = inv(P)*x
		ApplyPrecond(x);

		// Calculation of the real residual norm
		// inv(P)*b
		ApplyInversePrecond(b);
		cMultSym(A, x, v_aux, mNumThreads);
        cDeduct(b, v_aux);
		mResidualNorm = sqrt( cTwoNorm(v_aux) );
		mBNorm        = sqrt( cTwoNorm(b) );

        FinishSolver();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Diagonal Preconditioner
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//**********************************************************************
    //* - Initializing diagonal precondiotioner P = ( 1.00/sqrt( A(i,i) ) )
    //**********************************************************************    
	void CGComplexSolver::InitializePrecond(ComplexMatrixType& rA)
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
	void CGComplexSolver::ApplyInversePrecond(ComplexVectorType& rV)
    {    
	     for(int i = 0 ; i < mProblemSize ; ++i) rV[i] /= mDiagonal[i];
    }

	//*********************************************************************
    //* - P*X
    //*********************************************************************    
	void CGComplexSolver::ApplyPrecond(ComplexVectorType& rV)
    {
		for(int i = 0 ; i < mProblemSize ; ++i) rV[i] *= mDiagonal[i]; 
    }

	//*********************************************************************
    //* - Y = (P*A*P)*X , donde X = inv(P)*x
    //*********************************************************************  
	void CGComplexSolver::PrecondMultSym(ComplexMatrixType& rA, ComplexVectorType& rX, ComplexVectorType& rY)
    {
		 for(int i = 0 ; i < mProblemSize ; ++i) mTemp[i] = rX[i] * mDiagonal[i];

		 cMultSym(rA,mTemp,rY, mNumThreads);

         ApplyPrecond(rY);
    }

	//*********************************************************************
    //* - Y = conj(P*A*P)*X, donde X = inv(P)*x, conj(A*B)=conj(A)*conj(B)
    //*********************************************************************  
	void CGComplexSolver::PrecondConjgMultSym(ComplexMatrixType& rA, ComplexVectorType& rX, ComplexVectorType& rY)
    {
		 int i; 
		
		 for(i = 0 ; i < mProblemSize ; ++i) mTemp[i] = rX[i] * std::conj(mDiagonal[i]);

		 cConjgMultSym(rA, mTemp, rY, mNumThreads);

         for(i = 0 ; i < mProblemSize ; ++i) rY[i] *= std::conj(mDiagonal[i]); 
    }
}

