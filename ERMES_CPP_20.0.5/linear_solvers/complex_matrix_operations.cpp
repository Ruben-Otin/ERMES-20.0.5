#include <math.h>
#include <numeric>
#include <complex>
#include <vector>
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <algorithm>

#include "complex_matrix_operations.h"

//*****************************************************************************************************************
//* Sparse matrix-vector products for the ERMES iterative solvers.
//*
//* Storage conventions ( see Modeler::Assemble_ElementMatrix_to_Global_* ):
//*  - Full : every row i stores all its non-zero entries a_ij.
//*  - Symm / Herm : only the upper triangle ( j >= i ) is stored. The lower triangle is a_ji = a_ij ( complex
//*    symmetric ) or a_ji = conj(a_ij) ( Hermitian ).
//*
//* Products with an implicit lower triangle need a "scatter" to other rows ( Y[j] += ... ), which is a data race
//* when rows are processed in parallel. Each thread therefore accumulates its scatter contributions in a private
//* buffer, and the buffers are reduced in parallel at the end.
//*
//* Performance changes with respect to ERMES 20.0.4:
//*  - Scatter buffers are allocated once and reused across calls ( previously every thread copied a full
//*    vector on every product: one heap allocation of N complex values per thread per iteration ).
//*  - Rows are split in contiguous static blocks. With upper-triangular storage, the scatter of the thread owning
//*    rows [r0, r1) only touches columns >= r0, so its buffer only spans [r0, N) ( ~half the memory on average ).
//*  - The final reduction is parallel over the vector index ( previously serial inside an omp critical section ).
//*  - Products without scatter ( cMultFull, cConjgMultFull ) write Y[i] directly, without buffers.
//*  - Dot products, norms and vector updates are parallel ( see also complex_matrix_operations.h ).
//*
//* The results are mathematically identical to the previous implementation; only the floating point summation
//* order changes.
//*****************************************************************************************************************

namespace Kratos
{
    namespace
    {
        //*********************************************************************************************************
        //* Entry transformations used to build the different product variants
        //*********************************************************************************************************
        struct OpId    { inline ComplexType operator()( const ComplexType& a ) const { return a;            } };
        struct OpConj  { inline ComplexType operator()( const ComplexType& a ) const { return std::conj( a ); } };

        //*********************************************************************************************************
        //* Reusable per-thread scatter buffers ( the solvers call the products sequentially, never concurrently )
        //*********************************************************************************************************
        std::vector< ComplexVectorType > gScatterBuffers;

        //*********************************************************************************************************
        //* Contiguous block [r0, r1) of rows owned by thread tid
        //*********************************************************************************************************
        inline void RowBlock( long long n, int nthreads, int tid, long long& r0, long long& r1 )
        {
            long long chunk = n / nthreads;
            long long rest  = n % nthreads;

            r0 = tid * chunk + ( tid < rest ? tid : rest );
            r1 = r0 + chunk + ( tid < rest ? 1 : 0 );
        }

        inline int ValidThreads( int NumThreads )
        {
            return ( NumThreads > 0 ) ? NumThreads : 1;
        }

        //*********************************************************************************************************
        //* Generic product with gather ( row i ) and scatter ( implicit lower triangle ):
        //*
        //*   Y[i]  = fd( a_ii ) * X[i] + sum_{j != i} fo( a_ij ) * X[j]      ( stored entries of row i )
        //*   Y[j] +=                     g( a_ij )  * X[i]   for j != i       ( implicit entries of column i )
        //*
        //* If UpperOnly is true the matrix stores only j >= i, so thread scatter ranges start at its first row.
        //* If DoGather is false only the scatter part is computed ( transpose of a full matrix ).
        //*********************************************************************************************************
        template< class TFd, class TFo, class TG, bool DoGather, bool UpperOnly >
        void GatherScatterProduct( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
        {
            const long long n = A.RowsNumber();

            const int nthreads = ValidThreads( NumThreads );

            omp_set_num_threads( nthreads );

            TFd fd; TFo fo; TG g;

            if( (int)gScatterBuffers.size() < nthreads ) gScatterBuffers.resize( nthreads );

            // Real number of threads may be lower than requested ( nested / limited OpenMP ): the block partition
            // uses the number of threads actually running the region.
            #pragma omp parallel
            {
                const int tid = omp_get_thread_num();
                const int nth = omp_get_num_threads();

                long long r0, r1;
                RowBlock( n, nth, tid, r0, r1 );

                // Scatter range of this thread: [ s0, n )
                const long long s0 = UpperOnly ? r0 : 0;

                ComplexVectorType& buf = gScatterBuffers[ tid ];

                if( (long long)buf.size() != n - s0 ) buf.resize( n - s0 );

                std::fill( buf.begin(), buf.end(), ComplexType( 0.0, 0.0 ) );

                for( long long i = r0; i < r1; i++ )
                {
                    const ComplexType xi = X[ i ];

                    ComplexType result( 0.0, 0.0 );

                    std::vector< std::pair< unsigned int, ComplexType > >& row = A[ i ];

                    const std::size_t nnz = row.size();

                    for( std::size_t k = 0; k < nnz; k++ )
                    {
                        const long long   j   = row[ k ].first;
                        const ComplexType aij = row[ k ].second;

                        // Symmetric / Hermitian storage must only hold the upper triangle ( Modeler assembly ensures it )
                        if( UpperOnly && ( j < s0 ) )
                        {
                            std::cerr << "ERMES error: lower-triangular entry in a symmetric/Hermitian matrix ( row " << i << ", column " << j << " )" << std::endl;
                            std::abort();
                        }

                        if( j == i )
                        {
                            if( DoGather ) result   += fd( aij ) * xi;
                            else           buf[ i - s0 ] += g( aij ) * xi;
                        }
                        else
                        {
                            if( DoGather ) result += fo( aij ) * X[ j ];

                            buf[ j - s0 ] += g( aij ) * xi;
                        }
                    }

                    if( DoGather ) Y[ i ] = result;
                }

                #pragma omp barrier

                // Y[j] += sum over threads of their scatter buffers. Each thread reduces its own block of rows,
                // reading the buffers of the threads whose scatter range covers it.
                for( int t = 0; t < nth; t++ )
                {
                    long long t0, t1;
                    RowBlock( n, nth, t, t0, t1 );

                    const long long ts0 = UpperOnly ? t0 : 0;

                    ComplexVectorType& tbuf = gScatterBuffers[ t ];

                    const long long j0 = ( r0 > ts0 ) ? r0 : ts0;

                    for( long long j = j0; j < r1; j++ )
                    {
                        if( DoGather || t > 0 ) Y[ j ] += tbuf[ j - ts0 ];
                        else                    Y[ j ]  = tbuf[ j - ts0 ];
                    }
                }
            }
        }

        //*********************************************************************************************************
        //* Row-wise product without scatter ( full storage ): Y[i] = sum_j f( a_ij ) * X[j]
        //*********************************************************************************************************
        template< class TF >
        void RowProduct( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
        {
            const long long n = A.RowsNumber();

            TF f;

            omp_set_num_threads( ValidThreads( NumThreads ) );

            #pragma omp parallel for schedule(static)
            for( long long i = 0; i < n; i++ )
            {
                ComplexType result( 0.0, 0.0 );

                std::vector< std::pair< unsigned int, ComplexType > >& row = A[ i ];

                const std::size_t nnz = row.size();

                for( std::size_t k = 0; k < nnz; k++ )
                {
                    result += f( row[ k ].second ) * X[ row[ k ].first ];
                }

                Y[ i ] = result;
            }
        }
    }

    //*************************************************************************************************
    //* Frees the scatter buffers ( call after the solver to return the memory before post-processing )
    //*************************************************************************************************
    void cReleaseWorkspace()
    {
        std::vector< ComplexVectorType >().swap( gScatterBuffers );
    }

    //*************************************************************************************************
    //*  Y = A*X with A complex symmetric ( upper triangle stored )
    //*************************************************************************************************
    void cMultSymm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        GatherScatterProduct< OpId, OpId, OpId, true, true >( A, X, Y, NumThreads );
    }

	//*************************************************************************************************
	//*  Y = A*X with A Hermitian ( aij = conj(aji), upper triangle stored )
	//*************************************************************************************************
	void cMultHerm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
	{
        GatherScatterProduct< OpId, OpId, OpConj, true, true >( A, X, Y, NumThreads );
	}

    //*************************************************************************************************
    //*  Y = A*X with A full non-symmetric
    //*************************************************************************************************
    void cMultFull( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        RowProduct< OpId >( A, X, Y, NumThreads );
    }

    //*************************************************************************************************
    //*  Y = Atrans*X with A full non-symmetric ( scatter only )
    //*************************************************************************************************
    void cTransMultFull( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        GatherScatterProduct< OpId, OpId, OpId, false, false >( A, X, Y, NumThreads );
    }

    //*************************************************************************************************
    //*  Y = Atrans*X with A Hermitian ( upper triangle stored ). Atrans = conj(A):
    //*  row i: a_ii*X[i] + sum_{j>i} conj(a_ij)*X[j] ;  column part: Y[j] += a_ij*X[i]
    //*  ( the diagonal is used as stored, exactly as in ERMES 20.0.4 )
    //*************************************************************************************************
    void cTransMultHerm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        GatherScatterProduct< OpId, OpConj, OpId, true, true >( A, X, Y, NumThreads );
    }

    //*************************************************************************************************
    //*  Y = Atrans*X with A complex symmetric ( same as cMultSymm )
    //*************************************************************************************************
    void cTransMultSymm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        GatherScatterProduct< OpId, OpId, OpId, true, true >( A, X, Y, NumThreads );
    }

    //*************************************************************************************************
    //*  Y = conj(A)*X with A full non-symmetric
    //*************************************************************************************************
    void cConjgMultFull( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        RowProduct< OpConj >( A, X, Y, NumThreads );
    }

    //*************************************************************************************************
    //*  Y = conj(A)*X with A complex symmetric ( upper triangle stored )
    //*************************************************************************************************
    void cConjgMultSymm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
    {
        GatherScatterProduct< OpConj, OpConj, OpConj, true, true >( A, X, Y, NumThreads );
	}

	//*************************************************************************************************
	//*  Y = conj(A)*X with A Hermitian ( upper triangle stored ): conj(A)_ji = a_ij
	//*************************************************************************************************
	void cConjgMultHerm( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y, int NumThreads )
	{
        GatherScatterProduct< OpConj, OpConj, OpId, true, true >( A, X, Y, NumThreads );
	}

    //*****************************************************************************
    //*  Y = A*X ( serial )
    //*****************************************************************************
    void cMult( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    {
        RowProduct< OpId >( A, X, Y, 1 );
    }

    //*****************************************************************************
    //*  Y = conj(A)*X ( serial )
    //*****************************************************************************
    void cConjgMult( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    {
        RowProduct< OpConj >( A, X, Y, 1 );
    }

    //*****************************************************************************
    //*  Y = Atrans*X ( serial )
    //*****************************************************************************
    void cTransposeMult( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    {
        std::fill( Y.begin(), Y.end(), ComplexType() );

        const long long n = A.RowsNumber();

        for( long long i = 0; i < n; i++ )
        {
            std::vector< std::pair< unsigned int, ComplexType > >& row = A[ i ];

            for( std::size_t k = 0; k < row.size(); k++ ) Y[ row[ k ].first ] += row[ k ].second * X[ i ];
        }
    }

    //*****************************************************************************
    //*  Y = Aherm*X ( serial )
    //*****************************************************************************
    void cTransposeConjgMult( ComplexMatrixType& A, ComplexVectorType& X, ComplexVectorType& Y )
    {
        std::fill( Y.begin(), Y.end(), ComplexType() );

        const long long n = A.RowsNumber();

        for( long long i = 0; i < n; i++ )
        {
            std::vector< std::pair< unsigned int, ComplexType > >& row = A[ i ];

            for( std::size_t k = 0; k < row.size(); k++ ) Y[ row[ k ].first ] += std::conj( row[ k ].second ) * X[ i ];
        }
    }

    //*****************************************************************************
    //*  Squared 2-norm: sum |x_i|^2  ( callers take the sqrt )
    //*****************************************************************************
    double cTwoNorm( ComplexVectorType& X )
    {
        const long long n = X.size();

        double result = 0.00;

        #pragma omp parallel for reduction(+:result) schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ )
        {
            result += std::norm( X[ i ] );
        }

        return result;
    } 

    //*****************************************************************************
    //*  Bilinear ( unconjugated ) dot product: sum x_i * y_i  ( COCG / complex symmetric QMR )
    //*****************************************************************************
    ComplexType cDot( ComplexVectorType& X, ComplexVectorType& Y ) 
    {
        const long long n = X.size();

        double re = 0.00, im = 0.00;

        #pragma omp parallel for reduction(+:re,im) schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ )
        {
            const ComplexType p = X[ i ] * Y[ i ];

            re += p.real();
            im += p.imag();
        }

        return ComplexType( re, im );
    } 

    //*****************************************************************************
    //*  Sesquilinear dot product: sum x_i * conj(y_i)
    //*****************************************************************************
    ComplexType cConjgDot( ComplexVectorType& X, ComplexVectorType& Y ) 
    {
        const long long n = X.size();

        double re = 0.00, im = 0.00;

        #pragma omp parallel for reduction(+:re,im) schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ )
        {
            const ComplexType p = X[ i ] * std::conj( Y[ i ] );

            re += p.real();
            im += p.imag();
        }

        return ComplexType( re, im );
    } 

    //*******************************************************************************
    //*  Y = X + Y 
    //*******************************************************************************
    void cAdd( ComplexVectorType& X, ComplexVectorType& Y ) 
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Y[ i ] = X[ i ] + Y[ i ];
    } 

    //*******************************************************************************
    //*  Z = X + Y 
    //*******************************************************************************
    void cAdd( ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z ) 
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Z[ i ] = X[ i ] + Y[ i ];
    } 

    //*******************************************************************************
    //*  Y = X - Y 
    //*******************************************************************************
    void cDeduct( ComplexVectorType& X, ComplexVectorType& Y ) 
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Y[ i ] = X[ i ] - Y[ i ];
    } 

    //*******************************************************************************
    //*  Z = X - Y 
    //*******************************************************************************
    void cDeduct( ComplexVectorType& X, ComplexVectorType& Y, ComplexVectorType& Z ) 
    {
        const long long n = X.size();

        #pragma omp parallel for schedule(static) if( n > ERMES_OMP_MIN_SIZE )
        for( long long i = 0; i < n; i++ ) Z[ i ] = X[ i ] - Y[ i ];
    } 
}
