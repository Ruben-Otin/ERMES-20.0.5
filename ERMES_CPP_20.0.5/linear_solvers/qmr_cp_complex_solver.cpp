
#include <iostream>
#include <math.h>
#include <cmath>
#include "qmr_cp_complex_solver.h"

namespace Kratos
{
    //****************************************************************************************************
    //* - Iterative Quasi-Minimal Residual solver for cold plasma (Full matrix)
    //****************************************************************************************************
    void QMR_CP_ComplexSolver::Solve( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x )
    {
        msPreconditioner = "Diagonal";
        
        if ( !IniSolver( "Quasi-Minimal Residual (Full matrix)", A, b ) )
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
            Solve_PrecondQMR(A,b,x);
        }
        else 
        {
            FinishSolver();
        }
    }

    //**************************************************************************************************************************************
    //* - Iterative Quasi-Minimal Residual solver for cold plasma (A hermitic, A_aux symmetric)
    //**************************************************************************************************************************************
    void QMR_CP_ComplexSolver::Solve_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        msPreconditioner = "Diagonal";

        if ( !IniSolver( "Quasi-Minimal Residual (Herm-Symm matrix)", A, b ) )
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
            Solve_PrecondQMR_HermSymm( A, A_aux, b, x );
        }
        else 
        {
            FinishSolver();
        }
    }

    //**************************************************************************************************************************************
    //* - Iterative Quasi-Minimal Residual solver for cold plasma (A hermitic, A_aux full) 
    //**************************************************************************************************************************************
    void QMR_CP_ComplexSolver::Solve_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        msPreconditioner = "Diagonal";

        if ( !IniSolver( "Quasi-Minimal Residual (Herm-Full matrix)", A, b ) )
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
            Solve_PrecondQMR_HermFull( A, A_aux, b, x );
        }
        else 
        {
            FinishSolver();
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Preconditioned Quasi-Minimal Residual solver
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************************
    //* - Iterative Quasi-Minimal Residual solver (Full matrix)
    //*************************************************************************************************************    
    void QMR_CP_ComplexSolver::Solve_PrecondQMR( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x )
    {
        InitializePrecond( A );

        ComplexType rho  , rho_1  ; 
        ComplexType gamma, gamma_1;
        ComplexType theta, theta_1;
        ComplexType xi   , eta    ;
        ComplexType ep   , beta   ;
        ComplexType delta         ;

        ComplexVectorType v(mProblemSize), v_tld(mProblemSize); 
        ComplexVectorType w(mProblemSize), w_tld(mProblemSize); 
        ComplexVectorType y(mProblemSize), y_tld(mProblemSize); 
        ComplexVectorType z(mProblemSize), z_tld(mProblemSize); 
        ComplexVectorType p(mProblemSize), p_tld(mProblemSize); 
        ComplexVectorType q(mProblemSize), d    (mProblemSize); 
        ComplexVectorType r(mProblemSize), s    (mProblemSize); 

        ComplexVectorType v_aux(mProblemSize);

        cMultFull(A, x, v_aux, mNumThreads);
        cDeduct  (b, v_aux);

        r = v_aux;

        double normb = sqrt( cTwoNorm(b) );
        double resid = sqrt( cTwoNorm(r) ) / normb;

        v_tld = r;

        // y = M1.solve(v_tld);
        ApplyPrecond_M1(v_tld, y);

        rho   = sqrt( cTwoNorm(y) );
        w_tld = r;

        // z = M2.trans_solve(w_tld);
        ApplyPrecond_M2(w_tld, z);

        xi    = sqrt( cTwoNorm(z) );
        gamma = 1.0;
        eta   =-1.0;
        theta = 0.0;

        for ( int i = 1; i <= mMaxIterations; i++ ) 
        {
            if ( std::abs(rho) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |rho| = 0  !!!!!" << std::endl;
                break;
            }

            if ( std::abs(xi) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |xi| = 0  !!!!!" << std::endl;
                break;
            }
                
            // v = (1.0/rho) * v_tld;
            cScale( (1.0/rho), v_tld, v);

            // y = (1.0/rho) * y;
            cScale( (1.0/rho), y );

            // w = (1.0/xi) * w_tld;
            cScale( (1.0/xi), w_tld, w);

            // z = (1.0/xi) * z;
            cScale( (1.0/xi), z );

            delta = cDot( z, y );

            if ( std::abs(delta) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |delta| = 0  !!!!!" << std::endl;
                break;
            }

            // y_tld = M2.solve(y);               
            ApplyPrecond_M2(y, y_tld);

            // z_tld = M1.trans_solve(z);
            ApplyPrecond_M1(z, z_tld);

            if (i > 1) 
            {
                // p = y_tld - (xi * delta / ep) * p;
                cScaleSecondAndAdd( y_tld, -(xi*delta/ep), p );

                // q = z_tld - (rho * delta / ep) * q;
                cScaleSecondAndAdd( z_tld, -(rho*delta/ep), q );
            } 
            else 
            {
                p = y_tld;
                q = z_tld;
            }

            // p_tld = A * p;
            cMultFull( A, p, p_tld, mNumThreads );

            ep = cDot( q, p_tld );

            if ( std::abs(ep) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |ep| = 0  !!!!!" << std::endl;
                break;
            }

            beta = ep / delta;

            if ( std::abs(beta) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |beta| = 0  !!!!!" << std::endl;
                break;
            }

            // v_tld = p_tld - beta * v;
            cScaleFirstAndAdd( -beta, v, p_tld, v_tld );

            // y = M1.solve(v_tld);
            ApplyPrecond_M1(v_tld, y);

            rho_1 = rho;
            rho   = sqrt( cTwoNorm(y) );

            // w_tld = A.trans_mult(q) - beta(0) * w;
            cTransMultFull   (A, q, v_aux, mNumThreads);
            cScaleFirstAndAdd(-beta, w, v_aux, w_tld);

            // z = M2.trans_solve(w_tld);
            ApplyPrecond_M2(w_tld, z);

            xi = sqrt( cTwoNorm(z) );

            gamma_1 = gamma;
            theta_1 = theta;

            theta = rho / (gamma_1 * beta);
            gamma = 1.0 / sqrt(1.0 + theta * theta);

            if ( std::abs(gamma) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |gamma| = 0  !!!!!" << std::endl;
                break;
            }

            eta = ( -eta * rho_1 * gamma * gamma ) / ( beta * gamma_1 * gamma_1 );

            if (i > 1) 
            {
                // d = eta * p + (theta_1 * theta_1 * gamma * gamma) * d;
                cScale(eta, p, v_aux);
                cScaleSecondAndAdd( v_aux, (theta_1*theta_1*gamma*gamma), d );

                // s = eta * p_tld + (theta_1*theta_1*gamma*gamma) * s;
                cScale(eta, p_tld, v_aux);
                cScaleSecondAndAdd( v_aux, (theta_1*theta_1*gamma*gamma), s );
            } 
            else 
            {
                // d = eta * p;
                cScale(eta, p, d);

                // s = eta * p_tld;
                cScale(eta, p_tld, s);
            }

            // x = x + d;  
            cAdd(d, x);
            
            // r = r - s;  
            cScaleFirstAndAdd(-1.00, s, r);

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

    //*************************************************************************************************************************************************
    //* - Iterative Quasi-Minimal Residual solver (A hermitic, A_aux symmetric)
    //*************************************************************************************************************************************************    
    void QMR_CP_ComplexSolver::Solve_PrecondQMR_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        InitializePrecond( A, A_aux );

        ComplexType rho  , rho_1  ; 
        ComplexType gamma, gamma_1;
        ComplexType theta, theta_1;
        ComplexType xi   , eta    ;
        ComplexType ep   , beta   ;
        ComplexType delta         ;

        ComplexVectorType v(mProblemSize), v_tld(mProblemSize); 
        ComplexVectorType w(mProblemSize), w_tld(mProblemSize); 
        ComplexVectorType y(mProblemSize), y_tld(mProblemSize); 
        ComplexVectorType z(mProblemSize), z_tld(mProblemSize); 
        ComplexVectorType p(mProblemSize), p_tld(mProblemSize); 
        ComplexVectorType q(mProblemSize), d    (mProblemSize); 
        ComplexVectorType r(mProblemSize), s    (mProblemSize); 

        ComplexVectorType v_her(mProblemSize);
        ComplexVectorType v_aux(mProblemSize);

        cMultHerm(A    , x, v_her, mNumThreads);
        cMultSymm(A_aux, x, v_aux, mNumThreads);
        cAdd     (v_her, v_aux); 
        cDeduct  (b    , v_aux);

        r = v_aux;

        double normb = sqrt( cTwoNorm(b) );
        double resid = sqrt( cTwoNorm(r) ) / normb;

        v_tld = r;

        // y = M1.solve(v_tld);
        ApplyPrecond_M1(v_tld, y);

        rho   = sqrt( cTwoNorm(y) );
        w_tld = r;

        // z = M2.trans_solve(w_tld);
        ApplyPrecond_M2(w_tld, z);

        xi    = sqrt( cTwoNorm(z) );
        gamma = 1.0;
        eta   =-1.0;
        theta = 0.0;

        for ( int i = 1; i <= mMaxIterations; i++ ) 
        {
            if ( std::abs(rho) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |rho| = 0  !!!!!" << std::endl;
                break;
            }

            if ( std::abs(xi) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |xi| = 0  !!!!!" << std::endl;
                break;
            }

            // v = (1.0/rho) * v_tld;
            cScale( (1.0/rho), v_tld, v);

            // y = (1.0/rho) * y;
            cScale( (1.0/rho), y );

            // w = (1.0/xi) * w_tld;
            cScale( (1.0/xi), w_tld, w);

            // z = (1.0/xi) * z;
            cScale( (1.0/xi), z );

            delta = cDot( z, y );

            if ( std::abs(delta) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |delta| = 0  !!!!!" << std::endl;
                break;
            }

            // y_tld = M2.solve(y);               
            ApplyPrecond_M2(y, y_tld);

            // z_tld = M1.trans_solve(z);
            ApplyPrecond_M1(z, z_tld);

            if (i > 1) 
            {
                // p = y_tld - (xi * delta / ep) * p;
                cScaleSecondAndAdd( y_tld, -(xi*delta/ep), p );

                // q = z_tld - (rho * delta / ep) * q;
                cScaleSecondAndAdd( z_tld, -(rho*delta/ep), q );
            } 
            else 
            {
                p = y_tld;
                q = z_tld;
            }

            // p_tld = A * p;
            cMultHerm(A    , p, v_her, mNumThreads);
            cMultSymm(A_aux, p, v_aux, mNumThreads);
            cAdd     (v_her, v_aux, p_tld); 
            
            ep = cDot( q, p_tld );

            if ( std::abs(ep) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |ep| = 0  !!!!!" << std::endl;
                break;
            }

            beta = ep / delta;

            if ( std::abs(beta) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |beta| = 0  !!!!!" << std::endl;
                break;
            }

            // v_tld = p_tld - beta * v;
            cScaleFirstAndAdd( -beta, v, p_tld, v_tld );

            // y = M1.solve(v_tld);
            ApplyPrecond_M1(v_tld, y);

            rho_1 = rho;
            rho   = sqrt( cTwoNorm(y) );

            // w_tld = A.trans_mult(q) - beta(0) * w;
            cTransMultHerm(A    , q, v_her, mNumThreads);
            cTransMultSymm(A_aux, q, v_aux, mNumThreads);
            cAdd          (v_her, v_aux); 
   
            cScaleFirstAndAdd(-beta, w, v_aux, w_tld);

            // z = M2.trans_solve(w_tld);
            ApplyPrecond_M2(w_tld, z);

            xi = sqrt( cTwoNorm(z) );

            gamma_1 = gamma;
            theta_1 = theta;

            theta = rho / (gamma_1 * beta);
            gamma = 1.0 / sqrt(1.0 + theta * theta);

            if ( std::abs(gamma) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |gamma| = 0  !!!!!" << std::endl;
                break;
            }

            eta = ( -eta * rho_1 * gamma * gamma ) / ( beta * gamma_1 * gamma_1 );

            if (i > 1) 
            {
                // d = eta * p + (theta_1 * theta_1 * gamma * gamma) * d;
                cScale(eta, p, v_aux);
                cScaleSecondAndAdd( v_aux, (theta_1*theta_1*gamma*gamma), d );

                // s = eta * p_tld + (theta_1*theta_1*gamma*gamma) * s;
                cScale(eta, p_tld, v_aux);
                cScaleSecondAndAdd( v_aux, (theta_1*theta_1*gamma*gamma), s );
            } 
            else 
            {
                // d = eta * p;
                cScale(eta, p, d);

                // s = eta * p_tld;
                cScale(eta, p_tld, s);
            }

            // x = x + d;  
            cAdd(d, x);

            // r = r - s;  
            cScaleFirstAndAdd(-1.0, s, r);

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

    //*************************************************************************************************************
    //* - Iterative Quasi-Minimal Residual solver (A hermitic, A_aux full)
    //*************************************************************************************************************    
    void QMR_CP_ComplexSolver::Solve_PrecondQMR_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        InitializePrecond( A, A_aux );

        ComplexType rho  , rho_1  ; 
        ComplexType gamma, gamma_1;
        ComplexType theta, theta_1;
        ComplexType xi   , eta    ;
        ComplexType ep   , beta   ;
        ComplexType delta         ;

        ComplexVectorType v(mProblemSize), v_tld(mProblemSize); 
        ComplexVectorType w(mProblemSize), w_tld(mProblemSize); 
        ComplexVectorType y(mProblemSize), y_tld(mProblemSize); 
        ComplexVectorType z(mProblemSize), z_tld(mProblemSize); 
        ComplexVectorType p(mProblemSize), p_tld(mProblemSize); 
        ComplexVectorType q(mProblemSize), d    (mProblemSize); 
        ComplexVectorType r(mProblemSize), s    (mProblemSize); 

        ComplexVectorType v_her(mProblemSize);
        ComplexVectorType v_aux(mProblemSize);

        cMultHerm(A    , x, v_her, mNumThreads);
        cMultFull(A_aux, x, v_aux, mNumThreads);
        cAdd     (v_her, v_aux); 
        cDeduct  (b    , v_aux);

        r = v_aux;

        double normb = sqrt( cTwoNorm(b) );
        double resid = sqrt( cTwoNorm(r) ) / normb;

        v_tld = r;

        // y = M1.solve(v_tld);
        ApplyPrecond_M1(v_tld, y);

        rho   = sqrt( cTwoNorm(y) );
        w_tld = r;

        // z = M2.trans_solve(w_tld);
        ApplyPrecond_M2(w_tld, z);

        xi    = sqrt( cTwoNorm(z) );
        gamma = 1.0;
        eta   =-1.0;
        theta = 0.0;

        for ( int i = 1; i <= mMaxIterations; i++ ) 
        {
            if ( std::abs(rho) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |rho| = 0  !!!!!" << std::endl;
                break;
            }

            if ( std::abs(xi) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |xi| = 0  !!!!!" << std::endl;
                break;
            }

            // v = (1.0/rho) * v_tld;
            cScale( (1.0/rho), v_tld, v);

            // y = (1.0/rho) * y;
            cScale( (1.0/rho), y );

            // w = (1.0/xi) * w_tld;
            cScale( (1.0/xi), w_tld, w);

            // z = (1.0/xi) * z;
            cScale( (1.0/xi), z );

            delta = cDot( z, y );

            if ( std::abs(delta) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |delta| = 0  !!!!!" << std::endl;
                break;
            }

            // y_tld = M2.solve(y);               
            ApplyPrecond_M2(y, y_tld);

            // z_tld = M1.trans_solve(z);
            ApplyPrecond_M1(z, z_tld);

            if (i > 1) 
            {
                // p = y_tld - (xi * delta / ep) * p;
                cScaleSecondAndAdd( y_tld, -(xi*delta/ep), p );

                // q = z_tld - (rho * delta / ep) * q;
                cScaleSecondAndAdd( z_tld, -(rho*delta/ep), q );
            } 
            else 
            {
                p = y_tld;
                q = z_tld;
            }

            // p_tld = A * p;
            cMultHerm(A    , p, v_her, mNumThreads);
            cMultFull(A_aux, p, v_aux, mNumThreads);
            cAdd     (v_her, v_aux, p_tld); 

            ep = cDot( q, p_tld );

            if ( std::abs(ep) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |ep| = 0  !!!!!" << std::endl;
                break;
            }

            beta = ep / delta;

            if ( std::abs(beta) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |beta| = 0  !!!!!" << std::endl;
                break;
            }

            // v_tld = p_tld - beta * v;
            cScaleFirstAndAdd( -beta, v, p_tld, v_tld );

            // y = M1.solve(v_tld);
            ApplyPrecond_M1(v_tld, y);

            rho_1 = rho;
            rho   = sqrt( cTwoNorm(y) );

            // w_tld = A.trans_mult(q) - beta(0) * w;
            cTransMultHerm(A    , q, v_her, mNumThreads);
            cTransMultFull(A_aux, q, v_aux, mNumThreads);
            cAdd          (v_her, v_aux); 

            cScaleFirstAndAdd(-beta, w, v_aux, w_tld);

            // z = M2.trans_solve(w_tld);
            ApplyPrecond_M2(w_tld, z);

            xi = sqrt( cTwoNorm(z) );

            gamma_1 = gamma;
            theta_1 = theta;

            theta = rho / (gamma_1 * beta);
            gamma = 1.0 / sqrt(1.0 + theta * theta);

            if ( std::abs(gamma) == 0.0 )
            {
                std::cout << "!!!!!  Solver error: |gamma| = 0  !!!!!" << std::endl;
                break;
            }

            eta = ( -eta * rho_1 * gamma * gamma ) / ( beta * gamma_1 * gamma_1 );

            if (i > 1) 
            {
                // d = eta * p + (theta_1 * theta_1 * gamma * gamma) * d;
                cScale(eta, p, v_aux);
                cScaleSecondAndAdd( v_aux, (theta_1*theta_1*gamma*gamma), d );

                // s = eta * p_tld + (theta_1*theta_1*gamma*gamma) * s;
                cScale(eta, p_tld, v_aux);
                cScaleSecondAndAdd( v_aux, (theta_1*theta_1*gamma*gamma), s );
            } 
            else 
            {
                // d = eta * p;
                cScale(eta, p, d);

                // s = eta * p_tld;
                cScale(eta, p_tld, s);
            }

            // x = x + d;  
            cAdd(d, x);

            // r = r - s;  
            cScaleFirstAndAdd(-1.00, s, r);

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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Diagonal Preconditioner
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************************
    //* - Initializing diagonal precondiotioner P = [ 1.00 / A(i,i) ]
    //*************************************************************************************************************
    void QMR_CP_ComplexSolver::InitializePrecond( ComplexMatrixType& rA )
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
    void QMR_CP_ComplexSolver::InitializePrecond( ComplexMatrixType& rA, ComplexMatrixType& rA_aux )
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

    //*************************************************************************************************************
    //* - Y = M1 * X
    //*************************************************************************************************************    
    void QMR_CP_ComplexSolver::ApplyPrecond_M1( ComplexVectorType& rV, ComplexVectorType& rY )
    {
        for( int i = 0 ; i < mProblemSize ; ++i ) rY[i] = mDiagonal[i] * rV[i]; 
    }

    //*************************************************************************************************************
    //* - Y = M2 * X
    //*************************************************************************************************************    
    void QMR_CP_ComplexSolver::ApplyPrecond_M2( ComplexVectorType& rV, ComplexVectorType& rY )
    {
        for( int i = 0 ; i < mProblemSize ; ++i ) rY[i] = rV[i]; 
    }
}