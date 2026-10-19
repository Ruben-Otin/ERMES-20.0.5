
#include <iostream>
#include <math.h>
#include "cg_cp_complex_solver.h"

namespace Kratos
{
    //*********************************************************************
    //* - Iterative Conjugate Gradient solver for cold plasma
    //*********************************************************************    
    void CG_CP_ComplexSolver::Solve( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x )
    {
        msPreconditioner = "Diagonal";

        if ( !IniSolver( "Conjugate Gradient (Cold Plasma)", A, b ) )
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
			//Solve_PrecondCGS(A,b,x);
            Solve_PrecondBiCGSTAB(A,b,x);
		}
		else 
		{
			FinishSolver();
		}
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Preconditioned Conjugate Gradient iterative solver type
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //**************************************************************************************************************
    //* - Preconditioned Conjugate Gradient Squared iterative solver 
    //**************************************************************************************************************    
    void CG_CP_ComplexSolver::Solve_PrecondCGS( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x )
    {
        InitializePrecond(A);

        ComplexType rho_1, rho_2;
        ComplexType alpha, beta;

        ComplexVectorType p(mProblemSize), phat(mProblemSize); 
        ComplexVectorType q(mProblemSize), qhat(mProblemSize); 
        ComplexVectorType u(mProblemSize), uhat(mProblemSize); 
        ComplexVectorType vhat(mProblemSize);

        ComplexVectorType v_aux(mProblemSize);

        cMultFull(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

        ComplexVectorType r     (v_aux);
        ComplexVectorType rtilde(r    );

        double normb = sqrt( cTwoNorm(b) );
        double resid = sqrt( cTwoNorm(r) ) / normb;

        for ( int i = 1; i <= mMaxIterations; i++ ) 
        {
            rho_1 = cDot( rtilde, r );

            if ( std::abs(rho_1) == 0 ) 
            {
                std::cout << "!!!!!  Solver error: |rho_1| = 0  !!!!!" << std::endl;
                break;
            }

            if ( i == 1 ) 
            {
                u = r;
                p = u;
            } 
            else 
            {
                beta = rho_1 / rho_2;

                // u = r + beta * q;
                cScaleFirstAndAdd( beta, q, r, u );

                // p = u + beta * (q + beta * p);
                cScaleFirstAndAdd( beta, p    , q, v_aux );
                cScaleFirstAndAdd( beta, v_aux, u, p     );
            }

            // phat = M * p
            ApplyPrecondToVector( p, phat );

            // vhat = A * phat
            cMultFull(A, phat, vhat, mNumThreads);

            alpha = rho_1 / cDot( rtilde, vhat );

            // q = u - alpha * vhat
            cScaleFirstAndAdd( -alpha, vhat, u, q );

            // uhat = M.solve(u + q)
            cAdd                ( u, q, v_aux );
            ApplyPrecondToVector( v_aux, uhat );

            // x = x + alpha * uhat
            cScaleFirstAndAdd( alpha, uhat, x );

            // qhat = A * uhat
            cMultFull(A, uhat, qhat, mNumThreads);

            // r = r - alpha * qhat
            cScaleFirstAndAdd( -alpha, qhat, r );
            
            rho_2 = rho_1;

            resid = sqrt( cTwoNorm(r) ) / normb;

            mIterationNumber++; 

            mInfoCount++;

            if ( mInfoCount == mStepIterations )
            {
                std::cout<<"it = "<< mIterationNumber <<"   res = "<< resid;

                mInfoCount = 0;

                if (mWriteResultEveryStep)
                {
                    WriteStepResultsInFile(x);
                    std::cout<<"     <Write result in file done>";
                }

                std::cout << std::endl;
            }

            if ( resid < mTolerance ) 
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

    //**************************************************************************************************************
    //* - Preconditioned Bi-Conjugate Gradient Stabilized iterative solver 
    //**************************************************************************************************************    
    void CG_CP_ComplexSolver::Solve_PrecondBiCGSTAB( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x )
    {
        InitializePrecond(A);

        ComplexType rho_1, rho_2;
        ComplexType alpha, beta, omega;

        ComplexVectorType p(mProblemSize), phat(mProblemSize); 
        ComplexVectorType s(mProblemSize), shat(mProblemSize); 
        ComplexVectorType t(mProblemSize);
        ComplexVectorType v(mProblemSize);

        ComplexVectorType v_aux(mProblemSize);

        cMultFull(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

        ComplexVectorType r     (v_aux);
        ComplexVectorType rtilde(r    );

        double normb = sqrt( cTwoNorm(b) );
        double resid = sqrt( cTwoNorm(r) ) / normb;

        for ( int i = 1; i <= mMaxIterations; i++ ) 
        {
            rho_1 = cDot(rtilde, r);

            if ( std::abs(rho_1) == 0 ) 
            {
                std::cout << "!!!!!  Solver error: |rho_1| = 0  !!!!!" << std::endl;
                break;
            }

            if (i == 1)
            {
                p = r;
            }
            else 
            {
                beta = ( rho_1 / rho_2 ) * ( alpha / omega );

                // v_aux = p - omega * v;
                cScaleFirstAndAdd( -omega, v, p, v_aux );

                // p = r + beta * (p - omega * v);
                cScaleFirstAndAdd( beta, v_aux, r, p );
            }

            // phat = M * p
            ApplyPrecondToVector( p, phat );

            // v = A * phat;
            cMultFull(A, phat, v, mNumThreads);

            alpha = rho_1 / cDot( rtilde, v );

            // s = r - alpha * v;
            cScaleFirstAndAdd( -alpha, v, r, s );

            resid = sqrt( cTwoNorm(s) ) / normb;

            if ( resid < mTolerance ) 
            {
                // x = x + alpha * phat;
                cScaleFirstAndAdd( alpha, phat, x);
                break;
            }

            // shat = M * s
            ApplyPrecondToVector( s, shat );

            // t = A * shat
            cMultFull(A, shat, t, mNumThreads);

            omega = cDot(t,s) / cDot(t,t);

            // x = x + ( alpha * phat ) + ( omega * shat )
            cScaleAndAdd( alpha, phat, omega, shat, v_aux );
            cAdd        ( v_aux, x );

            // r = s - omega * t;
            cScaleFirstAndAdd( -omega, t, s, r );

            rho_2 = rho_1;

            resid = sqrt( cTwoNorm(r) ) / normb;

            mIterationNumber++; 

            mInfoCount++;

            if ( mInfoCount == mStepIterations )
            {
                std::cout<<"it = "<< mIterationNumber <<"   res = "<< resid;

                mInfoCount = 0;

                if (mWriteResultEveryStep)
                {
                    WriteStepResultsInFile(x);
                    std::cout<<"     <Write result in file done>";
                }

                std::cout << std::endl;
            }

            if ( resid < mTolerance ) 
            {
                break;
            }
            
            if ( std::abs(omega) == 0 ) 
            {
                std::cout << "!!!!!  Solver error: |omega| = 0  !!!!!" << std::endl;
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

	///////////////////////////////////////////////////////////////////////////////////////////////
    // - Diagonal Preconditioner
	///////////////////////////////////////////////////////////////////////////////////////////////

	//**************************************************************************
    //* - Initializing diagonal precondiotioner P = [ 1.00 / A(i,i) ]
    //**************************************************************************   
	void CG_CP_ComplexSolver::InitializePrecond(ComplexMatrixType& rA)
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

    //********************************************************************************************
    //* - Y = P * X
    //********************************************************************************************    
    void CG_CP_ComplexSolver::ApplyPrecondToVector(ComplexVectorType& rV, ComplexVectorType& rY)
    {
        for( int i = 0 ; i < mProblemSize ; ++i ) rY[i] = mDiagonal[i] * rV[i]; 
    }
}

