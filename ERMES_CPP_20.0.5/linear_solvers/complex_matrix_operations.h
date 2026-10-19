
#if !defined(COMPLEX_MATRIX_OPERATIONS_H_INCLUDED)
#define COMPLEX_MATRIX_OPERATIONS_H_INCLUDED

#include <complex>
#include "../includes/vector.h"
#include "../includes/csr_matrix.h"

#include <math.h>
#include <numeric>

// Vectors shorter than this are processed serially ( OpenMP start-up cost dominates for small sizes )
#if !defined(ERMES_OMP_MIN_SIZE)
#define ERMES_OMP_MIN_SIZE 20000
#endif

namespace Kratos
{
    typedef std::complex<double> ComplexType;

    typedef CSRMatrix<ComplexType> ComplexMatrixType;

    typedef std::vector<ComplexType> ComplexVectorType;

    //---------------------------------------------------------------------------------------------
    // Sparse matrix-vector products ( see complex_matrix_operations.cpp for storage conventions )
    //---------------------------------------------------------------------------------------------
    void cMultSymm     (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cConjgMultSymm(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);

    void cMultHerm     (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cConjgMultHerm(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);

    void cMultFull     (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cConjgMultFull(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cTransMultFull(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cTransMultHerm(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
    void cTransMultSymm(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads);
  
    void cMult         (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y);
    void cConjgMult    (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y);

    void cTransposeMult     (ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y);
    void cTransposeConjgMult(ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y);

    // Frees the per-thread buffers used by the symmetric / Hermitian / transpose products
    void cReleaseWorkspace();

    //---------------------------------------------------------------------------------------------
    // Vector operations ( OpenMP parallel above ERMES_OMP_MIN_SIZE entries )
    //---------------------------------------------------------------------------------------------

    // NOTE: returns the SQUARED 2-norm, sum |x_i|^2 ( callers apply sqrt )
    double cTwoNorm(ComplexVectorType& X);

    ComplexType cDot     (ComplexVectorType& X, ComplexVectorType& Y);  // sum x_i*y_i       ( bilinear )
    ComplexType cConjgDot(ComplexVectorType& X, ComplexVectorType& Y);  // sum x_i*conj(y_i) ( sesquilinear )

    void cAdd(ComplexVectorType& X, ComplexVectorType& Y);                         // Y = X + Y
    void cAdd(ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z);   // Z = X + Y

    void cDeduct(ComplexVectorType& X, ComplexVectorType& Y);                       // Y = X - Y
    void cDeduct(ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z); // Z = X - Y

    //*****************************************
    //* Y = (c * Y) 
    //*****************************************
    template<class TDataType>
    void cScale(TDataType c, ComplexVectorType& Y) 
    {
        const long long n = Y.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Y[ i ] = c * Y[ i ];
    } 

    //*****************************************
    //*  Y = (a * X)
    //*****************************************
    template<class TDataType>
    void cScale(TDataType a, ComplexVectorType& X, ComplexVectorType& Y) 
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Y[ i ] = a * X[ i ];
    } 

    //*****************************************
    //*  Z = (a * X) + (b * Y)
    //*****************************************
    template<class TDataType, class TDataType2>
    void cScaleAndAdd(TDataType a, ComplexVectorType& X, TDataType2 b, ComplexVectorType& Y, ComplexVectorType& Z)  
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Z[ i ] = ( a * X[ i ] ) + ( b * Y[ i ] );
    } 

    //***********************************************
    //*  Y = (a * X) + (b * Y) 
    //***********************************************
    template<class TDataType, class TDataType2>
    void cScaleAndAdd(TDataType a, ComplexVectorType& X, TDataType2 b, ComplexVectorType& Y) 
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Y[ i ] = ( a * X[ i ] ) + ( b * Y[ i ] );
    } 

    //***********************************************
    //*  Y = (a * X) +  Y 
    //***********************************************
    template<class TDataType>
    void cScaleFirstAndAdd(TDataType a, ComplexVectorType& X, ComplexVectorType& Y) 
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Y[ i ] = ( a * X[ i ] ) + Y[ i ];
    } 

	//*****************************************
    //*  Z = (a * X) +  Y
    //*****************************************
    template<class TDataType>
    void cScaleFirstAndAdd(TDataType a, ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z)  
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Z[ i ] = ( a * X[ i ] ) + Y[ i ];
    } 
 
    //***********************************************
    //*  Y =  X + (b * Y) 
    //***********************************************
    template<class TDataType>
    void cScaleSecondAndAdd(ComplexVectorType& X, TDataType b, ComplexVectorType& Y) 
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Y[ i ] = X[ i ] + ( b * Y[ i ] );
    } 
}

#endif
