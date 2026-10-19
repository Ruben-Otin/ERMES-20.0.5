
#include <iostream>
#include <math.h>
#include "bicg_cp_complex_solver.h"

namespace Kratos
{
    //*********************************************************************
    //* - Iterative Bi-Conjugate Gradient solver for cold plasma  
    //*********************************************************************    
    void BiCG_CP_ComplexSolver::Solve(ComplexMatrixType& A, 
                                      ComplexVectorType& b, 
                                      ComplexVectorType& x)
    {
        if (!IniSolver("Bi-Conjugate Gradient (Cold Plasma)", A, b))
        {
            FinishSolver();
            return;
        }

		ComplexVectorType v_aux(mProblemSize);

        cMultFull(A, x, v_aux, mNumThreads);
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
    // - Bi-conjugate Gradient solver. Unpreconditioned 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//*********************************************************************
    //* - Iterative Bi-Conjugate Gradient solver. Unpreconditioned 
    //*********************************************************************    
    void BiCG_CP_ComplexSolver::Solve_BiCG(ComplexMatrixType& A, 
                                           ComplexVectorType& b, 
                                           ComplexVectorType& x)
    {
        ComplexVectorType v_aux(mProblemSize);

        cMultFull(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

        ComplexVectorType r(v_aux);
        ComplexVectorType p(r);

        ComplexType rnorm = cDot(r,r);
        ComplexType alpha;
        ComplexType gamma;

		double tol_X_bnorm = mTolerance*mBNorm;

        do
        {
            cMultFull(A, p, v_aux, mNumThreads);

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
		cMultFull(A, x, v_aux, mNumThreads);
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
    void BiCG_CP_ComplexSolver::Solve_PrecondBiCG(ComplexMatrixType& A, 
                                              ComplexVectorType& b, 
                                              ComplexVectorType& x)
    {
        // Initializing preconditioner. 
		// We will solve ( (P*A*P)*(inv(P)*x) = P*b ) instead of ( A*x = b )     
		InitializePrecond(A);

        // Saving original b
        ComplexVectorType b0(b);

        // b = At * b0
        cTransposeMult(A, b0, b);

		// inv(P)*x
		ApplyInversePrecond(x);

		// P*b
		ApplyPrecond(b);

		// Now we solve with (P*b)
		mBNorm = sqrt( cTwoNorm(b) );

        ComplexVectorType v_aux(mProblemSize);
        
		// v_aux = (P*A*P)*x donde x = inv(P)*x
		PrecondMultFull(A, x, v_aux);

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
			PrecondMultFull(A, p, v_aux);

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
		//ApplyInversePrecond(b);

        cMultFull( A , x, v_aux, mNumThreads );
        cDeduct  ( b0, v_aux );

        mResidualNorm = sqrt( cTwoNorm(v_aux) );
        mBNorm        = sqrt( cTwoNorm(b0)    );

        FinishSolver();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Diagonal Preconditioner
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//**********************************************************************
    //* - Initializing diagonal precondiotioner P = ( 1.00/sqrt( A(i,i) ) )
    //**********************************************************************    
	void BiCG_CP_ComplexSolver::InitializePrecond(ComplexMatrixType& rA)
    {
        //mDiagonal.resize(mProblemSize);
        //mTemp.resize    (mProblemSize);

        //for(int i = 0 ; i < mProblemSize ; ++i) 
        //{
        //    std::complex<double> sqrtDiag = std::sqrt( rA(i,i) );

        //    if ( std::abs(sqrtDiag) != 0.00 ) mDiagonal[i] = 1.00 / sqrtDiag;
        //    else                              mDiagonal[i] = 1.00;
        //}

        mDiagonal.resize(mProblemSize,0.00);
        mTemp.resize    (mProblemSize,0.00);

        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        ComplexType aij;
        ComplexType czero(0.0, 0.0);

        for(int i = 0 ; i < mProblemSize ; ++i)  
        {
            mDiagonal[i] = czero;
            mTemp    [i] = czero;
        }

        int i_col;

        for (int i = 0; i < rA.RowsNumber(); i++)
        {
            data_iterator     = rA[i].begin();
            end_data_iterator = rA[i].end();

            while (data_iterator != end_data_iterator)
            {
                i_col = data_iterator->first;
                aij   = data_iterator->second;

                mTemp[i_col] += aij*aij;

                data_iterator++;
            }
        }

        for(int i = 0 ; i < mProblemSize ; ++i) 
        {
            std::complex<double> sqrtDiag = std::sqrt( mTemp[i] );

            if ( std::abs(sqrtDiag) != 0.00 ) mDiagonal[i] = 1.00 / sqrtDiag;
            else                              mDiagonal[i] = 1.00;
        }
    }

	//*********************************************************************
    //* - inv(P)*X
    //*********************************************************************    
	void BiCG_CP_ComplexSolver::ApplyInversePrecond(ComplexVectorType& rV)
    {    
	     for(int i = 0 ; i < mProblemSize ; ++i) rV[i] /= mDiagonal[i];
    }

	//*********************************************************************
    //* - P*X
    //*********************************************************************    
	void BiCG_CP_ComplexSolver::ApplyPrecond(ComplexVectorType& rV)
    {
		for(int i = 0 ; i < mProblemSize ; ++i) rV[i] *= mDiagonal[i]; 
    }

	//*********************************************************************
    //* - Y = (P*A*P)*X , donde X = inv(P)*x
    //*********************************************************************  
	void BiCG_CP_ComplexSolver::PrecondMultFull(ComplexMatrixType& rA, ComplexVectorType& rX, ComplexVectorType& rY)
    {
        //for(int i = 0 ; i < mProblemSize ; ++i) mTemp[i] = rX[i] * mDiagonal[i];

        //cMultFull(rA, mTemp, rY, mNumThreads);

        //ApplyPrecond(rY);

        for(int i = 0 ; i < mProblemSize ; ++i) mTemp[i] = rX[i] * mDiagonal[i];

        ComplexVectorType vTemp(mProblemSize);

        cMultFull(rA, mTemp, vTemp, mNumThreads);

        cTransposeMult(rA, vTemp, rY);

        ApplyPrecond(rY);
    }
}


