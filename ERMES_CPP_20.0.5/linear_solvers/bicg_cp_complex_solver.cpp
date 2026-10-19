
#include <iostream>
#include <math.h>
#include <cmath>
#include "bicg_cp_complex_solver.h"

namespace Kratos
{
    //****************************************************************************************************
    //* - Iterative Bi-Conjugate Gradient solver for cold plasma (Full matrix)
    //****************************************************************************************************  
    void BiCG_CP_ComplexSolver::Solve( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x )
    {
        msPreconditioner = "Diagonal";
        
        if ( !IniSolver( "Bi-Conjugate Gradient (Full matrix)", A, b ) )
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
            Solve_PrecondBiCG(A,b,x);
		}
		else 
		{
			FinishSolver();
		}
    }

    //**************************************************************************************************************************************
    //* - Iterative Bi-Conjugate Gradient solver for cold plasma (A hermitic, A_aux symmetric)
    //**************************************************************************************************************************************
    void BiCG_CP_ComplexSolver::Solve_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        msPreconditioner = "Diagonal";

        if ( !IniSolver( "Bi-Conjugate Gradient (Herm-Symm matrix)", A, b ) )
        {
            FinishSolver();
            return;
        }

        ComplexVectorType v_her(mProblemSize);
        ComplexVectorType v_aux(mProblemSize);

        cMultHerm(A    , x, v_her, mNumThreads);
        cMultSymm(A_aux, x, v_aux, mNumThreads);
        cAdd     (v_her, v_aux); 
        cDeduct  (b    , v_aux);

        mResidualNorm  = sqrt( cTwoNorm(v_aux) );
        mResEstimation = mResidualNorm / mBNorm;

        if ( mResEstimation > mTolerance ) 
        {
            Solve_PrecondBiCG_HermSymm( A, A_aux, b, x );
        }
        else 
        {
            FinishSolver();
        }
    }

    //**************************************************************************************************************************************
    //* - Iterative Bi-Conjugate Gradient solver for cold plasma (A hermitic, A_aux full) 
    //**************************************************************************************************************************************
    void BiCG_CP_ComplexSolver::Solve_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        msPreconditioner = "Diagonal";

        if ( !IniSolver( "Bi-Conjugate Gradient (Herm-Full matrix)", A, b ) )
        {
            FinishSolver();
            return;
        }

        ComplexVectorType v_her(mProblemSize);
        ComplexVectorType v_aux(mProblemSize);

        cMultHerm(A    , x, v_her, mNumThreads);
        cMultFull(A_aux, x, v_aux, mNumThreads);
        cAdd     (v_her, v_aux); 
        cDeduct  (b    , v_aux);

        mResidualNorm  = sqrt( cTwoNorm(v_aux) );
        mResEstimation = mResidualNorm / mBNorm;

        if ( mResEstimation > mTolerance ) 
        {
            Solve_PrecondBiCG_HermFull( A, A_aux, b, x );
        }
        else 
        {
            FinishSolver();
        }
    }

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Preconditioned Bi-Conjugate Gradient solver
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//****************************************************************************************************************
    //* - Iterative Bi-Conjugate Gradient solver (Full matrix)
    //****************************************************************************************************************    
    void BiCG_CP_ComplexSolver::Solve_PrecondBiCG( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x )
    {
        InitializePrecond( A );

        double resid = 1.0;   // relative residual reported if no iteration is run
        double normb = sqrt( cTwoNorm(b) );

        ComplexType rho_1, rho_2;
        ComplexType alpha, beta;

        ComplexVectorType z(mProblemSize), ztilde(mProblemSize); 
        ComplexVectorType p(mProblemSize), ptilde(mProblemSize); 
        ComplexVectorType q(mProblemSize), qtilde(mProblemSize);

        ComplexVectorType v_aux(mProblemSize);

        cMultFull(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

        ComplexVectorType r     (v_aux);
        ComplexVectorType rtilde(r    );

        for ( int i = 1; i <= mMaxIterations; i++ ) 
        {
            ApplyPrecondToVector(r     , z     );
            ApplyPrecondToVector(rtilde, ztilde);

            rho_1 = cDot(z, rtilde);

            if (i == 1) 
            {
                p      = z     ;
                ptilde = ztilde;
            } 
            else 
            {
                beta = rho_1 / rho_2;

                //p = z + beta * p;
                cScaleSecondAndAdd(z, beta, p);

                //ptilde = ztilde + beta * ptilde;
                cScaleSecondAndAdd(ztilde, beta, ptilde);
            }

            // q = A * p
            cMultFull(A, p, q, mNumThreads);
            
            // qtilde = Atrans * ptilde
            cTransMultFull(A, ptilde, qtilde, mNumThreads);

            alpha = rho_1 / cDot(ptilde, q);

            // x = x + alpha * p;
            cScaleFirstAndAdd(alpha, p, x);

            // r = r - alpha * q;
            cScaleFirstAndAdd(-alpha, q, r);

            // rtilde = rtilde - alpha * qtilde;
            cScaleFirstAndAdd(-alpha, qtilde, rtilde);

            rho_2 = rho_1;

            resid = sqrt( cTwoNorm(r) ) / normb;

            mIterationNumber++; 
            
            mInfoCount++;

            if( mInfoCount == mStepIterations )
            {
                std::cout << "it = " << mIterationNumber << "   res = " << resid;
               
                mInfoCount = 0;

                if( mWriteSolutionEveryStep )
                {
                    WriteStepResultsInFile( x );
                    std::cout << "     <Write result on file done>";
                }

                std::cout << std::endl;
            }

            if( ( resid < mTolerance ) || !std::isfinite( resid ) ) 
            {
                break;
            }
        }

        // Residual estimation
        mResEstimation = resid;

        // Real residual
        cMultFull(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

        mResidualNorm = sqrt( cTwoNorm(v_aux) );
        mBNorm        = sqrt( cTwoNorm(b)     );

        FinishSolver();
    }

    //***************************************************************************************************************************************************
    //* - Iterative Bi-Conjugate Gradient solver (A hermitic, A_aux symmetric)
    //***************************************************************************************************************************************************    
    void BiCG_CP_ComplexSolver::Solve_PrecondBiCG_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        InitializePrecond( A, A_aux );

        double resid = 1.0;   // relative residual reported if no iteration is run
        double normb = sqrt( cTwoNorm(b) );

        ComplexType rho_1, rho_2;
        ComplexType alpha, beta;

        ComplexVectorType z(mProblemSize), ztilde(mProblemSize); 
        ComplexVectorType p(mProblemSize), ptilde(mProblemSize); 
        ComplexVectorType q(mProblemSize), qtilde(mProblemSize);

        ComplexVectorType v_her(mProblemSize);
        ComplexVectorType v_aux(mProblemSize);

        cMultHerm(A    , x, v_her, mNumThreads);
        cMultSymm(A_aux, x, v_aux, mNumThreads);
        cAdd     (v_her, v_aux); 
        cDeduct  (b    , v_aux);

        ComplexVectorType r     (v_aux);
        ComplexVectorType rtilde(r    );

        for ( int i = 1; i <= mMaxIterations; i++ ) 
        {
            ApplyPrecondToVector(r     , z     );
            ApplyPrecondToVector(rtilde, ztilde);

            rho_1 = cDot(z, rtilde);

            if (i == 1) 
            {
                p      = z     ;
                ptilde = ztilde;
            } 
            else 
            {
                beta = rho_1 / rho_2;

                //p = z + beta * p;
                cScaleSecondAndAdd(z, beta, p);

                //ptilde = ztilde + beta * ptilde;
                cScaleSecondAndAdd(ztilde, beta, ptilde);
            }

            // q = A * p
            cMultHerm(A    , p, v_her, mNumThreads);
            cMultSymm(A_aux, p, v_aux, mNumThreads);
            cAdd     (v_her, v_aux, q); 

            // qtilde = Atrans * ptilde
            cTransMultHerm(A    , ptilde, v_her, mNumThreads);
            cTransMultSymm(A_aux, ptilde, v_aux, mNumThreads);
            cAdd          (v_her, v_aux, qtilde); 

            alpha = rho_1 / cDot(ptilde, q);

            // x = x + alpha * p;
            cScaleFirstAndAdd(alpha, p, x);

            // r = r - alpha * q;
            cScaleFirstAndAdd(-alpha, q, r);

            // rtilde = rtilde - alpha * qtilde;
            cScaleFirstAndAdd(-alpha, qtilde, rtilde);

            rho_2 = rho_1;

            resid = sqrt( cTwoNorm(r) ) / normb;

            mIterationNumber++; 

            mInfoCount++;

            if( mInfoCount == mStepIterations )
            {
                std::cout << "it = " << mIterationNumber << "   res = " << resid;

                mInfoCount = 0;

                if( mWriteSolutionEveryStep )
                {
                    WriteStepResultsInFile( x );
                    std::cout << "     <Write result on file done>";
                }

                std::cout << std::endl;
            }

            if( ( resid < mTolerance ) || !std::isfinite( resid ) ) 
            {
                break;
            }
        }

        // Residual estimation
        mResEstimation = resid;

        // Real residual
        cMultHerm(A    , x, v_her, mNumThreads);
        cMultSymm(A_aux, x, v_aux, mNumThreads);
        cAdd     (v_her, v_aux); 
        cDeduct  (b    , v_aux);

        mResidualNorm = sqrt( cTwoNorm(v_aux) );
        mBNorm        = sqrt( cTwoNorm(b)     );

        FinishSolver();
    }

    //***************************************************************************************************************************************************
    //* - Iterative Bi-Conjugate Gradient solver (A hermitic, A_aux full)
    //***************************************************************************************************************************************************    
    void BiCG_CP_ComplexSolver::Solve_PrecondBiCG_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        InitializePrecond( A, A_aux );

        double resid = 1.0;   // relative residual reported if no iteration is run
        double normb = sqrt( cTwoNorm(b) );

        ComplexType rho_1, rho_2;
        ComplexType alpha, beta;

        ComplexVectorType z(mProblemSize), ztilde(mProblemSize); 
        ComplexVectorType p(mProblemSize), ptilde(mProblemSize); 
        ComplexVectorType q(mProblemSize), qtilde(mProblemSize);

        ComplexVectorType v_her(mProblemSize);
        ComplexVectorType v_aux(mProblemSize);

        cMultHerm(A    , x, v_her, mNumThreads);
        cMultFull(A_aux, x, v_aux, mNumThreads);
        cAdd     (v_her, v_aux); 
        cDeduct  (b    , v_aux);

        ComplexVectorType r     (v_aux);
        ComplexVectorType rtilde(r    );

        for ( int i = 1; i <= mMaxIterations; i++ ) 
        {
            ApplyPrecondToVector(r     , z     );
            ApplyPrecondToVector(rtilde, ztilde);

            rho_1 = cDot(z, rtilde);

            if (i == 1) 
            {
                p      = z     ;
                ptilde = ztilde;
            } 
            else 
            {
                beta = rho_1 / rho_2;

                //p = z + beta * p;
                cScaleSecondAndAdd(z, beta, p);

                //ptilde = ztilde + beta * ptilde;
                cScaleSecondAndAdd(ztilde, beta, ptilde);
            }

            // q = A * p
            cMultHerm(A    , p, v_her, mNumThreads);
            cMultFull(A_aux, p, v_aux, mNumThreads);
            cAdd     (v_her, v_aux, q); 

            // qtilde = Atrans * ptilde
            cTransMultHerm(A    , ptilde, v_her, mNumThreads);
            cTransMultFull(A_aux, ptilde, v_aux, mNumThreads);
            cAdd          (v_her, v_aux, qtilde); 

            alpha = rho_1 / cDot(ptilde, q);

            // x = x + alpha * p;
            cScaleFirstAndAdd(alpha, p, x);

            // r = r - alpha * q;
            cScaleFirstAndAdd(-alpha, q, r);

            // rtilde = rtilde - alpha * qtilde;
            cScaleFirstAndAdd(-alpha, qtilde, rtilde);

            rho_2 = rho_1;

            resid = sqrt( cTwoNorm(r) ) / normb;

            mIterationNumber++; 

            mInfoCount++;

            if( mInfoCount == mStepIterations )
            {
                std::cout << "it = " << mIterationNumber << "   res = " << resid;

                mInfoCount = 0;

                if( mWriteSolutionEveryStep )
                {
                    WriteStepResultsInFile( x );
                    std::cout << "     <Write result on file done>";
                }

                std::cout << std::endl;
            }

            if( ( resid < mTolerance ) || !std::isfinite( resid ) ) 
            {
                break;
            }
        }

        // Residual estimation
        mResEstimation = resid;

        // Real residual
        cMultHerm(A    , x, v_her, mNumThreads);
        cMultFull(A_aux, x, v_aux, mNumThreads);
        cAdd     (v_her, v_aux); 
        cDeduct  (b    , v_aux);

        mResidualNorm = sqrt( cTwoNorm(v_aux) );
        mBNorm        = sqrt( cTwoNorm(b)     );

        FinishSolver();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////
    // - Diagonal Preconditioner
	///////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************************
    //* - Initializing diagonal precondiotioner P = [ 1.00 / A(i,i) ]
    //*************************************************************************************************************
    void BiCG_CP_ComplexSolver::InitializePrecond( ComplexMatrixType& rA )
    {
        mDiagonal.resize(mProblemSize);

        for(int i = 0 ; i < mProblemSize ; ++i) 
        {
            std::complex<double> diag_ii = rA(i,i);

            if ( std::abs(diag_ii) != 0.00 ) 
            {
                mDiagonal[i] = 1.00 / diag_ii;
            }
            else
            {
                mDiagonal[i] = 1.00;
            }
        }
    }

    //*************************************************************************************************************
    //* - Initializing diagonal precondiotioner P = [ 1.00 / ( A(i,i) + A_aux(i,i) ) ]
    //*************************************************************************************************************
    void BiCG_CP_ComplexSolver::InitializePrecond( ComplexMatrixType& rA, ComplexMatrixType& rA_aux )
    {
        mDiagonal.resize(mProblemSize);

        for(int i = 0 ; i < mProblemSize ; ++i) 
        {
            std::complex<double> diag_ii = rA(i,i) + rA_aux(i,i);

            if ( std::abs(diag_ii) != 0.00 ) 
            {
                mDiagonal[i] = 1.00 / diag_ii;
            }
            else
            {
                mDiagonal[i] = 1.00;
            }
        }
    }

    //********************************************************************************************
    //* - Y = P * X
    //********************************************************************************************    
    void BiCG_CP_ComplexSolver::ApplyPrecondToVector(ComplexVectorType& rV, ComplexVectorType& rY)
    {
        for( int i = 0 ; i < mProblemSize ; ++i ) rY[i] = mDiagonal[i] * rV[i]; 
    }
}


