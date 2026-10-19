
#include <iostream>
#include <math.h>
#include "qmr_cp_complex_solver.h"

namespace Kratos
{
    //*********************************************************************
    //* - Iterative Quasi-Minimal Residual solver for cold plasma  
    //*********************************************************************    
    void QMR_CP_ComplexSolver::Solve(ComplexMatrixType& A, 
                                     ComplexVectorType& b, 
                                     ComplexVectorType& x)
    {
        if (!IniSolver("Quasi-Minimal Residual (Cold Plasma)", A, b))
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
                Solve_PrecondQMR(A,b,x);
            }
            else 
            {
                Solve_QMR(A,b,x);
            }
        }
        else 
        {
            FinishSolver();
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // - Quasi-Minimal Residual. Unpreconditioned 
    ////////////////////////////////////////////////////////////////////////

    //*********************************************************************
    //* - Iterative Quasi-Minimal Residual solver. Unpreconditioned 
    //*********************************************************************    
    void QMR_CP_ComplexSolver::Solve_QMR(ComplexMatrixType& A, 
                                         ComplexVectorType& b, 
                                         ComplexVectorType& x)
    {
        ComplexVectorType v_aux(b);

        cMultFull(A, x, v_aux, mNumThreads);

        // v_aux = b - v_aux
        cDeduct(b, v_aux);

        double resd = sqrt( cTwoNorm(v_aux) );
        double rho  = resd;

        // v = v_aux;
        ComplexVectorType v(v_aux);

        // v = (1/rho) * v
        cScale((1.00/rho), v); 

        // p = 0.00
        ComplexVectorType p(mProblemSize);

        cScale(0.00, p); 

        // d = p
        ComplexVectorType d(p);

        double c     = 1.00;
        double theta = 1.00;

        double c_m1;
        double rho_1;
        double theta_m1;
        double abs_beta;

        ComplexType ep    ( 1.00, 0.00);
        ComplexType s     ( 1.00, 0.00);
        ComplexType nu    (-1.00, 0.00);
        ComplexType delta ( 0.00, 0.00);
        ComplexType beta  ( 0.00, 0.00);

        mResidualNorm = resd;

        double tol_X_bnorm = mTolerance*mBNorm;

        do
        {
            // delta = v * v
            delta = cDot(v,v); 

            // p = v - (rho*delta/ep) * p
            cScaleSecondAndAdd(v, -(rho*delta/ep), p);

            // v_aux = A*p
            cMultFull(A, p, v_aux, mNumThreads);

            // ep = p * ( A*p )
            ep = cDot(p, v_aux);

            // beta = ep / delta
            beta     = ep / delta;
            abs_beta = std::abs(beta);

            // v = v_aux - beta * v
            cScaleSecondAndAdd(v_aux, -beta, v);

            theta_m1 = theta;

            // rho_n+1 = sqrt( conj(v) * v )
            rho_1 = sqrt( cTwoNorm(v) );

            // theta = rho/(c*|beta|)
            theta = rho_1 / ( c * abs_beta );

            c_m1 = c;

            // c_n = 1.00 / sqrt( 1.00 + theta^2 )
            c = 1.00 / sqrt( 1.00 + (theta*theta) );

            s = ( beta / abs_beta ) * ( (rho_1 * c) / (c_m1 * abs_beta) );

            nu = - nu * ( (rho * c * c) / (beta * c_m1 * c_m1) );

            // d = nu * p + (theta*c_aux)^2 * d
            cScaleAndAdd(nu, p, (theta_m1*c)*(theta_m1*c), d);  

            // x = d + x
            cAdd(d, x);

            // v = (1/rho_n+1) * v
            cScale((1.00/rho_1), v); 

            rho = rho_1;

            mIterationNumber++;
            mInfoCount++;

            mResidualNorm = mResidualNorm * sqrt( 1.00 + (1.00/mIterationNumber) ) * std::abs(s);

            if ( mInfoCount == mStepIterations ) 
            {
                cMultFull(A, x, v_aux, mNumThreads);
                cDeduct  (b, v_aux);

                mResidualNorm = sqrt( cTwoNorm(v_aux) );

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

    ///////////////////////////////////////////////////////////////////////
    // - Quasi-Minimal Residual solver. Preconditioned 
    ///////////////////////////////////////////////////////////////////////

    //*********************************************************************
    //* - Iterative Quasi-Minimal Residual solver working
    //*********************************************************************    
    void QMR_CP_ComplexSolver::Solve_PrecondQMR(ComplexMatrixType& A, 
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
        mBNorm = sqrt(cTwoNorm(b));

        ComplexVectorType v_aux(b);

        // v_aux = (P*A*P)*x, donde x = inv(P)*x
        PrecondMultFull(A, x, v_aux);

        // v_aux = b - v_aux
        cDeduct(b, v_aux);

        double resd = sqrt( cTwoNorm(v_aux) );
        double rho  = resd;

        // v = v_aux;
        ComplexVectorType v(v_aux);

        // v = (1/rho) * v
        cScale((1.00/rho), v); 

        // p = 0.00
        ComplexVectorType p(mProblemSize);

        cScale(0.00, p); 

        // d = p
        ComplexVectorType d(p);

        double c     = 1.00;
        double theta = 1.00;

        double c_m1;
        double rho_1;
        double theta_m1;
        double abs_beta;

        ComplexType ep    ( 1.00, 0.00);
        ComplexType s     ( 1.00, 0.00);
        ComplexType nu    (-1.00, 0.00);
        ComplexType delta ( 0.00, 0.00);
        ComplexType beta  ( 0.00, 0.00);

        mResidualNorm = resd;

        double tol_X_bnorm = mTolerance*mBNorm;

        do
        {
            // delta = v * v
            delta = cDot(v,v); 

            // p = v - (rho*delta/ep) * p
            cScaleSecondAndAdd(v, -(rho*delta/ep), p);

            // v_aux = (P*A*P)*p
            PrecondMultFull(A, p, v_aux);

            // ep = p * (P*A*P)*p
            ep = cDot(p, v_aux);

            // beta = ep / delta
            beta     = ep / delta;
            abs_beta = std::abs(beta);

            // v = v_aux - beta * v
            cScaleSecondAndAdd(v_aux, -beta, v);

            theta_m1 = theta;

            // rho_n+1 = sqrt( conj(v) * v )
            rho_1 = sqrt(cTwoNorm(v));

            // theta = rho/(c*|beta|)
            theta = rho_1 / ( c * abs_beta );

            c_m1 = c;

            // c_n = 1.00 / sqrt( 1.00 + theta^2 )
            c  = 1.00 / sqrt( 1.00 + (theta*theta) );

            s  = ( beta / abs_beta ) * ( (rho_1 * c) / (c_m1 * abs_beta) );

            nu = - nu * ( (rho * c * c) / (beta * c_m1 * c_m1) );

            // d = nu * p + (theta*c_aux)^2 * d
            cScaleAndAdd(nu, p, (theta_m1*c)*(theta_m1*c), d);  

            // x = d + x
            cAdd(d, x);

            // v = (1/rho_n+1) * v
            cScale((1.00/rho_1), v); 

            rho = rho_1;

            mIterationNumber++;
            mInfoCount++;

            mResidualNorm = mResidualNorm * sqrt( 1.00 + (1.00/mIterationNumber) ) * std::abs(s);

            if ( mInfoCount == mStepIterations ) 
            {
                PrecondMultFull(A, x, v_aux);
                cDeduct        (b, v_aux);

                mResidualNorm = sqrt( cTwoNorm(v_aux) );

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

    /////////////////////////////////////////////////////////////////////////
    // - Diagonal Preconditioner
    /////////////////////////////////////////////////////////////////////////

    //**********************************************************************
    //* - Initializing diagonal precondiotioner P = ( 1.00/sqrt( A(i,i) ) )
    //**********************************************************************    
    void QMR_CP_ComplexSolver::InitializePrecond(ComplexMatrixType& rA)
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
    void QMR_CP_ComplexSolver::ApplyInversePrecond(ComplexVectorType& rV)
    {    
        for(int i = 0 ; i < mProblemSize ; ++i) rV[i] /= mDiagonal[i];
    }

    //*********************************************************************
    //* - P*X
    //*********************************************************************    
    void QMR_CP_ComplexSolver::ApplyPrecond(ComplexVectorType& rV)
    {
        for(int i = 0 ; i < mProblemSize ; ++i) rV[i] *= mDiagonal[i]; 
    }

    //*********************************************************************
    //* - Y = (P*A*P)*X , donde X = inv(P)*x
    //*********************************************************************  
    void QMR_CP_ComplexSolver::PrecondMultFull(ComplexMatrixType& rA, ComplexVectorType& rX, ComplexVectorType& rY)
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