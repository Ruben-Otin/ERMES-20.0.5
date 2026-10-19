
#include <iostream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <new>
#include <omp.h>

#include "multifrontal_complex_solver.h"
#include "direct_complex_solver.h"

// Runtime CPU dispatch for the dense kernel ( GCC on x86-64 ): a generic version plus versions compiled
// for AVX2/FMA ( Haswell and newer ) and AVX-512 ( Skylake-SP and newer ). The best one is selected when ERMES
// starts, so the executable stays portable to any x86-64 CPU and still uses the vector units of modern ones.
#if defined(__GNUC__) && !defined(__clang__) && defined(__x86_64__) && !defined(ERMES_NO_CLONES)
#define ERMES_KERNEL_CLONES __attribute__(( target_clones( "default", "arch=haswell", "arch=skylake-avx512" ) ))
#else
#define ERMES_KERNEL_CLONES
#endif

//*****************************************************************************************************************
//* Optional external libraries ( see linear_solvers/README.txt and unix/Makefile ):
//*
//*  - ERMES_USE_BLAS  : dense updates with cblas_zgemm from an optimized BLAS ( OpenBLAS, BLIS, ... ).
//*                      The BLAS must run single-threaded inside ERMES ( e.g. OPENBLAS_NUM_THREADS=1 ):
//*                      ERMES already distributes the work over its OpenMP threads.
//*  - ERMES_USE_METIS : nested dissection ordering ( METIS_NodeND, METIS 5 ) for large matrices, which
//*                      reduces fill-in and operations for 3D meshes compared with minimum degree.
//*
//* Without them ERMES uses its own kernel and its approximate minimum degree ordering.
//*****************************************************************************************************************
#if defined(ERMES_USE_BLAS)
extern "C" void cblas_zgemm( int Order, int TransA, int TransB, int M, int N, int K, const void* alpha,
                             const void* A, int lda, const void* B, int ldb, const void* beta, void* C, int ldc );
#endif

#if defined(ERMES_USE_METIS)
#include <metis.h>
#endif

// The BLAS is called only outside ERMES's OpenMP regions ( top fronts ), with all threads; OpenBLAS and BLIS
// are told the number of threads directly, so no environment variable is needed.
#if defined(ERMES_USE_OPENBLAS)
extern "C" void openblas_set_num_threads( int num_threads );
extern "C" int  openblas_get_parallel();                       // 0: sequential build, 1: pthreads, 2: OpenMP
#endif

#if defined(ERMES_USE_BLIS)
#include <stdint.h>
extern "C" void    bli_thread_set_num_threads( int64_t n_threads );
extern "C" int64_t bli_info_get_enable_threading();
#endif

namespace
{
    // External BLAS used for the top fronts ( set by ComplexMultifrontal::SelectDenseKernel )
    bool gUseBlas = false;
}

#include <cstdlib>

namespace Kratos
{
    namespace
    {
        typedef std::complex<double> CT;

        //*************************************************************************************************
        //* Dense kernels on interleaved complex arrays ( column major ). Plain double arithmetic is used:
        //* std::complex multiplication calls a slow NaN-checking routine unless -ffast-math is used.
        //*************************************************************************************************

        // y[0..len) -= a * x[0..len)
        inline void Axpy( double* y, const double* x, double ar, double ai, long long len )
        {
            for( long long i=0; i<len; i++ )
            {
                const double xr = x[ 2*i ], xi = x[ 2*i + 1 ];

                y[ 2*i     ] -= xr*ar - xi*ai;
                y[ 2*i + 1 ] -= xr*ai + xi*ar;
            }
        }

        //*************************************************************************************************
        //* Dense update kernel ( the core of the factorization, "GEMM" ).
        //*************************************************************************************************
        ERMES_KERNEL_CLONES
        void SchurUpdate( double* C, long long ldc, const double* A, long long lda, const double* B, long long ldb,
                          long long mr, long long nc, long long k, bool Lower, int NumThreads )
        {
            // C(i,j) -= sum_p A(i,p) * B(p,j),  i in [0,mr), j in [0,nc), p in [0,k)  ( i >= j only if Lower:
            // entries above the diagonal of the current 4-column group may also be updated; they are unused )
            //
            // Blocking: 32 columns per task ( OpenMP ), 128 rows x 128 inner products per cache block, 4-column
            // register groups. A and C blocks are packed with real and imaginary parts in separate arrays, so the
            // inner loops are plain real multiply-adds that the compiler vectorizes ( SSE2 / AVX2 / AVX-512 ).
#if defined(ERMES_USE_BLAS)
            // Optimized BLAS: C = C - A * B with cblas_zgemm ( column major, no transpose ). It is only called
            // from outside OpenMP parallel regions ( the top fronts ), one call at a time, and the BLAS uses
            // its own threads there. Inside parallel regions ( independent subtrees, one per thread ) ERMES's own
            // kernel is used instead: single-threaded BLAS builds are often not thread-safe when called
            // concurrently, which can corrupt results ( NaN ).
            if( gUseBlas && !omp_in_parallel() )
            {
                const int CblasColMajor = 102, CblasNoTrans = 111;

                const std::complex<double> MinusOne( -1.0, 0.0 ), One( 1.0, 0.0 );

                // Lower: column blocks, each updating only the rows from its first column down
                const long long JBB = Lower ? 256 : nc;

                for( long long j0=0; j0<nc; j0+=JBB )
                {
                    const long long nb = std::min( JBB, nc - j0 );
                    const long long i0 = Lower ? j0 : 0;

                    if( i0 >= mr ) continue;

                    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, (int)( mr - i0 ), (int)nb, (int)k, &MinusOne,
                                 A + 2 * i0, (int)lda, B + 2 * ( j0 * ldb ), (int)ldb, &One, C + 2 * ( i0 + j0 * ldc ), (int)ldc );
                }

                return;
            }
#endif

            const long long JB = 32, PB = 128, IB = 128;

            const long long nJB = ( nc + JB - 1 ) / JB;

            const bool Parallel = ( NumThreads > 1 ) && ( (double)mr * (double)nc * (double)k > 2.0e6 );

            #pragma omp parallel for schedule(dynamic,1) num_threads( NumThreads ) if( Parallel )
            for( long long jb=0; jb<nJB; jb++ )
            {
                std::vector<double> Ar( IB * PB ), Ai( IB * PB ), Cr( 4 * IB ), Ci( 4 * IB );

                const long long j0 = jb * JB;
                const long long j1 = std::min( nc, j0 + JB );

                const long long iStart = Lower ? j0 : 0;

                for( long long i0=iStart; i0<mr; i0+=IB )
                {
                    const long long i1 = std::min( mr, i0 + IB );
                    const long long li = i1 - i0;

                    for( long long p0=0; p0<k; p0+=PB )
                    {
                        const long long p1 = std::min( k, p0 + PB );

                        // Pack A( i0:i1, p0:p1 )
                        for( long long p=p0; p<p1; p++ )
                        {
                            const double* a = A + 2 * ( i0 + p * lda );

                            double* ar = &Ar[ ( p - p0 ) * IB ];
                            double* ai = &Ai[ ( p - p0 ) * IB ];

                            for( long long i=0; i<li; i++ ) { ar[ i ] = a[ 2*i ]; ai[ i ] = a[ 2*i + 1 ]; }
                        }

                        for( long long j=j0; j<j1; j+=4 )
                        {
                            const long long nj = std::min( (long long)4, j1 - j );

                            const long long ib = Lower ? std::max( i0, j ) : i0;

                            if( ib >= i1 ) continue;

                            const long long off = ib - i0, len = i1 - ib;

                            // Unpack C( ib:i1, j:j+nj )
                            for( long long q=0; q<nj; q++ )
                            {
                                const double* c = C + 2 * ( ib + ( j + q ) * ldc );

                                for( long long i=0; i<len; i++ ) { Cr[ q * IB + i ] = c[ 2*i ]; Ci[ q * IB + i ] = c[ 2*i + 1 ]; }
                            }

                            for( long long q=nj; q<4; q++ )
                                for( long long i=0; i<len; i++ ) { Cr[ q * IB + i ] = 0.0; Ci[ q * IB + i ] = 0.0; }

                            double* cr0 = &Cr[ 0 ];        double* ci0 = &Ci[ 0 ];
                            double* cr1 = &Cr[ IB ];       double* ci1 = &Ci[ IB ];
                            double* cr2 = &Cr[ 2 * IB ];   double* ci2 = &Ci[ 2 * IB ];
                            double* cr3 = &Cr[ 3 * IB ];   double* ci3 = &Ci[ 3 * IB ];

                            for( long long p=p0; p<p1; p++ )
                            {
                                const double* bp = B + 2 * ( p + j * ldb );

                                double br[ 4 ] = { 0.0, 0.0, 0.0, 0.0 }, bi[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };

                                for( long long q=0; q<nj; q++ ) { br[ q ] = bp[ 2*q*ldb ]; bi[ q ] = bp[ 2*q*ldb + 1 ]; }

                                const double b0r = br[ 0 ], b0i = bi[ 0 ], b1r = br[ 1 ], b1i = bi[ 1 ];
                                const double b2r = br[ 2 ], b2i = bi[ 2 ], b3r = br[ 3 ], b3i = bi[ 3 ];

                                const double* ar = &Ar[ ( p - p0 ) * IB + off ];
                                const double* ai = &Ai[ ( p - p0 ) * IB + off ];

                                for( long long i=0; i<len; i++ )
                                {
                                    const double xr = ar[ i ], xi = ai[ i ];

                                    cr0[ i ] -= xr*b0r - xi*b0i;  ci0[ i ] -= xr*b0i + xi*b0r;
                                    cr1[ i ] -= xr*b1r - xi*b1i;  ci1[ i ] -= xr*b1i + xi*b1r;
                                    cr2[ i ] -= xr*b2r - xi*b2i;  ci2[ i ] -= xr*b2i + xi*b2r;
                                    cr3[ i ] -= xr*b3r - xi*b3i;  ci3[ i ] -= xr*b3i + xi*b3r;
                                }
                            }

                            // Pack C back
                            for( long long q=0; q<nj; q++ )
                            {
                                double* c = C + 2 * ( ib + ( j + q ) * ldc );

                                for( long long i=0; i<len; i++ ) { c[ 2*i ] = Cr[ q * IB + i ]; c[ 2*i + 1 ] = Ci[ q * IB + i ]; }
                            }
                        }
                    }
                }
            }
        }



        // Elimination tree ( Liu, path compression ) of the lower pattern given by rows ( rp, ri: j < k in row k )
        void EliminationTree( int n, const std::vector<int>& rp, const std::vector<int>& ri, std::vector<int>& parent )
        {
            std::vector<int> ancestor( n, -1 );

            parent.assign( n, -1 );

            for( int k=0; k<n; k++ )
            {
                for( int p=rp[ k ]; p<rp[ k + 1 ]; p++ )
                {
                    int r = ri[ p ];

                    while( ( ancestor[ r ] != -1 ) && ( ancestor[ r ] != k ) )
                    {
                        int next = ancestor[ r ];
                        ancestor[ r ] = k;
                        r = next;
                    }

                    if( ancestor[ r ] == -1 ) { ancestor[ r ] = k; parent[ r ] = k; }
                }
            }
        }

        // Lower pattern of the permuted A + A^T by rows: row k holds the columns j < k ( duplicates allowed )
        void LowerRows( int n, const std::vector<int>& Ap, const std::vector<int>& Ai, const std::vector<int>& iperm,
                        std::vector<int>& rp, std::vector<int>& ri )
        {
            rp.assign( n + 1, 0 );

            for( int j=0; j<n; j++ )
                for( int p=Ap[ j ]; p<Ap[ j + 1 ]; p++ )
                {
                    int a = iperm[ Ai[ p ] ], b = iperm[ j ];
                    if( a != b ) rp[ std::max( a, b ) + 1 ]++;
                }

            for( int k=0; k<n; k++ ) rp[ k + 1 ] += rp[ k ];

            ri.resize( rp[ n ] );

            std::vector<int> pos( rp.begin(), rp.end() - 1 );

            for( int j=0; j<n; j++ )
                for( int p=Ap[ j ]; p<Ap[ j + 1 ]; p++ )
                {
                    int a = iperm[ Ai[ p ] ], b = iperm[ j ];
                    if( a != b ) ri[ pos[ std::max( a, b ) ]++ ] = std::min( a, b );
                }
        }
    }

    //*****************************************************************************************************
    //* - Nested dissection ordering with METIS ( Q[k] = k-th pivot ). Returns false if not available.
    //*****************************************************************************************************
    static bool MetisOrdering( int n, const std::vector<int>& Ap, const std::vector<int>& Ai, std::vector<int>& Q )
    {
#if defined(ERMES_USE_METIS)
        // Graph of A + A^T without diagonal ( CSR, METIS index type )
        std::vector< std::vector<idx_t> > adj( n );

        for( int j=0; j<n; j++ )
            for( int p=Ap[ j ]; p<Ap[ j + 1 ]; p++ )
            {
                int i = Ai[ p ];
                if( i == j ) continue;
                adj[ i ].push_back( j );
                adj[ j ].push_back( i );
            }

        std::vector<idx_t> xadj( n + 1, 0 ), adjncy;

        for( int i=0; i<n; i++ )
        {
            std::sort( adj[ i ].begin(), adj[ i ].end() );
            adj[ i ].erase( std::unique( adj[ i ].begin(), adj[ i ].end() ), adj[ i ].end() );
            xadj[ i + 1 ] = xadj[ i ] + adj[ i ].size();
        }

        adjncy.reserve( xadj[ n ] );

        for( int i=0; i<n; i++ ) { adjncy.insert( adjncy.end(), adj[ i ].begin(), adj[ i ].end() ); std::vector<idx_t>().swap( adj[ i ] ); }

        if( adjncy.empty() ) return false;   // no edges ( diagonal matrix )

        idx_t nv = n;
        idx_t options[ METIS_NOPTIONS ];

        METIS_SetDefaultOptions( options );
        options[ METIS_OPTION_NUMBERING ] = 0;

        std::vector<idx_t> perm( n ), iperm( n );

        if( METIS_NodeND( &nv, xadj.data(), adjncy.data(), NULL, options, perm.data(), iperm.data() ) != METIS_OK ) return false;

        // METIS: row k of the permuted matrix is row perm[k] of A ( same convention as Q )
        Q.assign( perm.begin(), perm.end() );

        return true;
#else
        (void)n; (void)Ap; (void)Ai; (void)Q;
        return false;
#endif
    }

    //*****************************************************************************************************
    //* - Chooses the dense kernel for the top fronts and returns a description for the .info file.
    //*   The external BLAS is used when it can run with all the threads: a multithreaded build, or one thread.
    //*   A sequential BLAS build with several ERMES threads would run the largest fronts on a single core,
    //*   so ERMES's own parallel kernel is faster there. ERMES_BLAS=OFF forces ERMES's kernel.
    //*****************************************************************************************************
    std::string ComplexMultifrontal::SelectDenseKernel( int NumThreads )
    {
        gUseBlas = false;

#if defined(ERMES_USE_BLAS)
        const char* pEnv = std::getenv( "ERMES_BLAS" );

        if( ( pEnv != NULL ) && ( ( std::string( pEnv ) == "OFF" ) || ( std::string( pEnv ) == "off" ) ) )
        {
            return "ERMES kernel ( external BLAS disabled by ERMES_BLAS=OFF )";
        }

        bool Threaded = false;

#if defined(ERMES_USE_OPENBLAS)
        Threaded = ( openblas_get_parallel() > 0 );
#elif defined(ERMES_USE_BLIS)
        Threaded = ( bli_info_get_enable_threading() > 0 );
#endif

        if( Threaded || ( NumThreads <= 1 ) )
        {
            gUseBlas = true;

            return Threaded ? "external BLAS ( multithreaded )" : "external BLAS ( 1 thread )";
        }

#if defined(ERMES_USE_OPENBLAS) || defined(ERMES_USE_BLIS)
        return "ERMES parallel kernel ( the external BLAS is a sequential build )";
#else
        return "ERMES parallel kernel ( unknown BLAS threading: external BLAS used only with 1 thread )";
#endif
#else
        (void)NumThreads;

        return "ERMES kernel";
#endif
    }

    void ComplexMultifrontal::Clear()
    {
        mN = 0; mNs = 0; mPerturbed = 0;

        std::vector<int>().swap( mPerm ); std::vector<int>().swap( mIperm );
        std::vector<int>().swap( mSfirst ); std::vector<int>().swap( mSparent );
        std::vector<long long>().swap( mStructPtr ); std::vector<int>().swap( mStructIdx );
        std::vector<long long>().swap( mLptr ); std::vector<long long>().swap( mUptr );
        std::vector<CT>().swap( mLx ); std::vector<CT>().swap( mUx );
        std::vector<int>().swap( mPAp ); std::vector<int>().swap( mPAi ); std::vector<int>().swap( mPAmap );
        std::vector<int>().swap( mPiv );
    }

    //*****************************************************************************************************
    //* - Symbolic analysis
    //*****************************************************************************************************
    void ComplexMultifrontal::Analyze( int n, const std::vector<int>& Ap, const std::vector<int>& Ai, bool Symmetric )
    {
        Clear();

        mN   = n;
        mSym = Symmetric;

        // 1) Fill-reducing ordering: METIS nested dissection for large matrices when available ( environment
        //    variable ERMES_ORDERING=AMD forces minimum degree ), otherwise approximate minimum degree
        std::vector<int> Q;

        const char* pOrd = std::getenv( "ERMES_ORDERING" );

        const bool ForceAMD = ( pOrd != NULL ) && ( std::string( pOrd ) == "AMD" );

        mOrdering = "approximate minimum degree";

        if( ForceAMD || ( n < 2000 ) || !MetisOrdering( n, Ap, Ai, Q ) )
        {
            ApproximateMinimumDegree( n, Ap, Ai, Q );
        }
        else
        {
            mOrdering = "METIS nested dissection";
        }

        std::vector<int> iq( n );

        for( int k=0; k<n; k++ ) iq[ Q[ k ] ] = k;

        // 2) Elimination tree and postorder ( supernodes need children numbered just before parents )
        std::vector<int> rp, ri, parent;

        LowerRows( n, Ap, Ai, iq, rp, ri );

        EliminationTree( n, rp, ri, parent );

        std::vector<int> head( n, -1 ), next( n, -1 ), post;

        post.reserve( n );

        for( int j=n-1; j>=0; j-- )
        {
            if( parent[ j ] >= 0 ) { next[ j ] = head[ parent[ j ] ]; head[ parent[ j ] ] = j; }
        }

        std::vector<int> stack;

        for( int r=0; r<n; r++ )
        {
            if( parent[ r ] != -1 ) continue;

            stack.push_back( r );

            while( !stack.empty() )
            {
                int top = stack.back();
                int child = head[ top ];

                if( child == -1 )
                {
                    stack.pop_back();
                    post.push_back( top );
                }
                else
                {
                    head[ top ] = next[ child ];
                    stack.push_back( child );
                }
            }
        }

        mPerm.resize( n );
        mIperm.resize( n );

        for( int k=0; k<n; k++ ) mPerm[ k ] = Q[ post[ k ] ];
        for( int k=0; k<n; k++ ) mIperm[ mPerm[ k ] ] = k;

        // 3) Tree and column counts in the final order
        LowerRows( n, Ap, Ai, mIperm, rp, ri );

        EliminationTree( n, rp, ri, parent );

        std::vector<long long> cc( n, 1 );
        std::vector<int>       mark( n, -1 );

        for( int k=0; k<n; k++ )
        {
            mark[ k ] = k;

            for( int p=rp[ k ]; p<rp[ k + 1 ]; p++ )
            {
                for( int r = ri[ p ]; ( r != -1 ) && ( mark[ r ] != k ); r = parent[ r ] )
                {
                    cc[ r ]++;
                    mark[ r ] = k;
                }
            }
        }

        std::vector<int> nchild( n, 0 );

        for( int j=0; j<n; j++ ) if( parent[ j ] >= 0 ) nchild[ parent[ j ] ]++;

        // 4) Fundamental supernodes + relaxed amalgamation
        struct SN { int f; int nc; long long m; double z; };

        std::vector<SN> sn;

        int j = 0;

        while( j < n )
        {
            SN s; s.f = j; s.nc = 1; s.z = 0.0;

            while( ( j + 1 < n ) && ( parent[ j ] == j + 1 ) && ( cc[ j ] == cc[ j + 1 ] + 1 ) && ( nchild[ j + 1 ] == 1 ) )
            {
                j++;
                s.nc++;
            }

            s.m = cc[ s.f ] - s.nc;

            // Merge with the previous supernodes while they are children of this one and the extra zeros are few
            while( !sn.empty() )
            {
                SN& c = sn.back();

                const int cl = c.f + c.nc - 1;

                if( ( parent[ cl ] < s.f ) || ( parent[ cl ] > s.f + s.nc - 1 ) ) break;

                const double ncNew = c.nc + s.nc;
                const double zNew  = c.z + s.z + (double)c.nc * (double)( s.nc + s.m - c.m );
                const double total = ncNew * ( ncNew + 1.0 ) / 2.0 + ncNew * (double)s.m;
                const double frac  = zNew / total;

                bool merge = ( ncNew <= 4 ) || ( ( ncNew <= 16 ) && ( frac < 0.8 ) ) ||
                             ( ( ncNew <= 48 ) && ( frac < 0.1 ) ) || ( frac < 0.05 );

                if( !merge ) break;

                s.f  = c.f;
                s.nc = (int)ncNew;
                s.z  = zNew;

                sn.pop_back();
            }

            sn.push_back( s );

            j++;
        }

        mNs = sn.size();

        mSfirst.resize( mNs + 1 );

        std::vector<int> snode( n );

        for( int s=0; s<mNs; s++ )
        {
            mSfirst[ s ] = sn[ s ].f;

            for( int c=sn[ s ].f; c<sn[ s ].f + sn[ s ].nc; c++ ) snode[ c ] = s;
        }

        mSfirst[ mNs ] = n;

        mSparent.assign( mNs, -1 );

        for( int s=0; s<mNs; s++ )
        {
            int l = mSfirst[ s + 1 ] - 1;

            if( parent[ l ] >= 0 ) mSparent[ s ] = snode[ parent[ l ] ];
        }

        // 5) Row structure of each supernode ( below its diagonal block ), bottom-up
        std::vector<int> lcp( n + 1, 0 ), lci( rp[ n ] );

        for( int k=0; k<n; k++ ) for( int p=rp[ k ]; p<rp[ k + 1 ]; p++ ) lcp[ ri[ p ] + 1 ]++;
        for( int c=0; c<n; c++ ) lcp[ c + 1 ] += lcp[ c ];
        {
            std::vector<int> pos( lcp.begin(), lcp.end() - 1 );
            for( int k=0; k<n; k++ ) for( int p=rp[ k ]; p<rp[ k + 1 ]; p++ ) lci[ pos[ ri[ p ] ]++ ] = k;
        }

        std::vector<int> shead( mNs, -1 ), snext( mNs, -1 );

        for( int s=mNs-1; s>=0; s-- )
        {
            if( mSparent[ s ] >= 0 ) { snext[ s ] = shead[ mSparent[ s ] ]; shead[ mSparent[ s ] ] = s; }
        }

        mStructPtr.assign( mNs + 1, 0 );
        mStructIdx.clear();

        std::fill( mark.begin(), mark.end(), -1 );

        std::vector<int> rows;

        mPredEntries = 0.0; mPredFront = 0.0; mPredOps = 0.0; mMaxFront = 0;

        for( int s=0; s<mNs; s++ )
        {
            const int f = mSfirst[ s ], l = mSfirst[ s + 1 ] - 1;

            rows.clear();

            for( int c=f; c<=l; c++ )
                for( int p=lcp[ c ]; p<lcp[ c + 1 ]; p++ )
                {
                    int i = lci[ p ];
                    if( ( i > l ) && ( mark[ i ] != s ) ) { mark[ i ] = s; rows.push_back( i ); }
                }

            for( int c=shead[ s ]; c!=-1; c=snext[ c ] )
                for( long long p=mStructPtr[ c ]; p<mStructPtr[ c + 1 ]; p++ )
                {
                    int i = mStructIdx[ p ];
                    if( ( i > l ) && ( mark[ i ] != s ) ) { mark[ i ] = s; rows.push_back( i ); }
                }

            std::sort( rows.begin(), rows.end() );

            mStructIdx.insert( mStructIdx.end(), rows.begin(), rows.end() );

            mStructPtr[ s + 1 ] = mStructIdx.size();

            const double ns = l - f + 1, m = rows.size(), nf = ns + m;

            mPredEntries += ns * nf + ( mSym ? 0.0 : ns * m );

            mPredFront = std::max( mPredFront, nf * nf );

            mMaxFront = std::max( mMaxFront, (int)nf );

            mPredOps += mSym ? ( ns * ns * nf / 2.0 + ns * m * m / 2.0 ) : ( ns * ns * nf + ns * m * m );
        }

        // 6) Permuted matrix for assembly ( CSC, values taken through mPAmap )
        mPAp.assign( n + 1, 0 );

        for( int c=0; c<n; c++ ) mPAp[ mIperm[ c ] + 1 ] += Ap[ c + 1 ] - Ap[ c ];
        for( int c=0; c<n; c++ ) mPAp[ c + 1 ] += mPAp[ c ];

        mPAi.resize( mPAp[ n ] );
        mPAmap.resize( mPAp[ n ] );

        for( int c=0; c<n; c++ )
        {
            int q = mPAp[ mIperm[ c ] ];

            for( int p=Ap[ c ]; p<Ap[ c + 1 ]; p++, q++ )
            {
                mPAi  [ q ] = mIperm[ Ai[ p ] ];
                mPAmap[ q ] = p;
            }
        }
    }

    //*****************************************************************************************************
    //* - Numeric factorization
    //*****************************************************************************************************
    bool ComplexMultifrontal::Factorize( const std::vector<int>& Ap, const std::vector<int>& Ai, const std::vector<CT>& Ax,
                                         int NumThreads, bool Verbose, std::string& ErrorMsg )
    {
        const int n = mN;

        if( NumThreads < 1 ) NumThreads = 1;

// The BLAS is only called for the top fronts, outside OpenMP regions: it may use all threads there
#if defined(ERMES_USE_OPENBLAS)
        openblas_set_num_threads( NumThreads );
#endif
#if defined(ERMES_USE_BLIS)
        bli_thread_set_num_threads( NumThreads );
#endif

        mPerturbed = 0;

        // Static pivot threshold
        double Amax = 0.0;

        for( std::size_t p=0; p<Ax.size(); p++ ) Amax = std::max( Amax, std::abs( Ax[ p ] ) );

        // Static pivot threshold ( relative to the largest entry; DirectComplexSolver equilibrates the matrix
        // so that the largest entry of every row is ~1 ). sqrt( machine epsilon ) as in PARDISO: a larger
        // perturbation limits the growth of the factors, the accuracy is recovered by GMRES.
        const double tau = 1.0e-8 * ( ( Amax > 0.0 ) ? Amax : 1.0 );

        // Upper entries by rows ( LU only ): transpose of the permuted CSC
        std::vector<int> PRp, PRj, PRmap;

        if( !mSym )
        {
            PRp.assign( n + 1, 0 );

            for( std::size_t p=0; p<mPAi.size(); p++ ) PRp[ mPAi[ p ] + 1 ]++;
            for( int i=0; i<n; i++ ) PRp[ i + 1 ] += PRp[ i ];

            PRj.resize( mPAi.size() ); PRmap.resize( mPAi.size() );

            std::vector<int> pos( PRp.begin(), PRp.end() - 1 );

            for( int c=0; c<n; c++ )
                for( int p=mPAp[ c ]; p<mPAp[ c + 1 ]; p++ )
                {
                    int q = pos[ mPAi[ p ] ]++;
                    PRj[ q ] = c; PRmap[ q ] = mPAmap[ p ];
                }
        }

        // Supernodal tree: children lists, subtree ranges and costs
        std::vector<int> shead( mNs, -1 ), snext( mNs, -1 ), sfirst( mNs );
        std::vector<double> nodeOps( mNs ), subOps( mNs );

        for( int s=mNs-1; s>=0; s-- )
            if( mSparent[ s ] >= 0 ) { snext[ s ] = shead[ mSparent[ s ] ]; shead[ mSparent[ s ] ] = s; }

        for( int s=0; s<mNs; s++ )
        {
            const double ns = mSfirst[ s + 1 ] - mSfirst[ s ];
            const double m  = mStructPtr[ s + 1 ] - mStructPtr[ s ];
            const double nf = ns + m;

            nodeOps[ s ] = mSym ? ( ns * ns * nf / 2.0 + ns * m * m / 2.0 ) : ( ns * ns * nf + ns * m * m );
            subOps [ s ] = nodeOps[ s ];
            sfirst [ s ] = s;
        }

        for( int s=0; s<mNs; s++ )   // children have lower indices ( postorder )
        {
            int p = mSparent[ s ];

            if( p >= 0 ) { subOps[ p ] += subOps[ s ]; sfirst[ p ] = std::min( sfirst[ p ], sfirst[ s ] ); }
        }

        try
        {
            mLptr.assign( mNs + 1, 0 );
            mUptr.assign( mNs + 1, 0 );

            for( int s=0; s<mNs; s++ )
            {
                const long long ns = mSfirst[ s + 1 ] - mSfirst[ s ];
                const long long m  = mStructPtr[ s + 1 ] - mStructPtr[ s ];

                mLptr[ s + 1 ] = mLptr[ s ] + ns * ( ns + m );
                mUptr[ s + 1 ] = mUptr[ s ] + ( mSym ? 0 : ns * m );
            }

            mLx.assign( mLptr[ mNs ], CT( 0.0, 0.0 ) );
            mUx.assign( mUptr[ mNs ], CT( 0.0, 0.0 ) );
        }
        catch( std::bad_alloc& )
        {
            Clear();
            ErrorMsg = "not enough memory for the factors";
            return false;
        }

        std::vector< std::vector<CT> > cb( mNs );

        double DoneOps = 0.0;
        int    NextReport = 10;
        bool   OutOfMemory = false;
        bool   Breakdown   = false;

        mPiv.assign( n, 0 );

        time_t t0; time( &t0 );

        //*************************************************************************************************
        //* Factorization of one supernode ( KThreads: threads for its dense kernels )
        //*************************************************************************************************
        auto ProcessNode = [&]( int s, std::vector<CT>& F, std::vector<CT>& W, std::vector<int>& relpos, int KThreads )
        {
            const int f  = mSfirst[ s ];
            const int l  = mSfirst[ s + 1 ] - 1;
            const long long ns = l - f + 1;
            const long long m  = mStructPtr[ s + 1 ] - mStructPtr[ s ];
            const long long nf = ns + m;

            const int* sidx = mStructIdx.data() + mStructPtr[ s ];

            if( (long long)F.size() < nf * nf ) F.resize( nf * nf );

            for( long long i=0; i<ns; i++ ) relpos[ f + i ] = i;
            for( long long i=0; i<m ; i++ ) relpos[ sidx[ i ] ] = ns + i;

            std::fill( F.begin(), F.begin() + nf * nf, CT( 0.0, 0.0 ) );

            // Assembly: original entries
            for( int c=f; c<=l; c++ )
                for( int p=mPAp[ c ]; p<mPAp[ c + 1 ]; p++ )
                {
                    int i = mPAi[ p ];
                    if( i >= c ) F[ relpos[ i ] + relpos[ c ] * nf ] += Ax[ mPAmap[ p ] ];
                }

            if( !mSym )
            {
                for( int i=f; i<=l; i++ )
                    for( int p=PRp[ i ]; p<PRp[ i + 1 ]; p++ )
                    {
                        int c = PRj[ p ];
                        if( c > i ) F[ relpos[ i ] + relpos[ c ] * nf ] += Ax[ PRmap[ p ] ];
                    }
            }

            // Assembly: contribution blocks of the children ( extend-add )
            for( int c=shead[ s ]; c!=-1; c=snext[ c ] )
            {
                const long long mc = mStructPtr[ c + 1 ] - mStructPtr[ c ];
                const int* cidx = mStructIdx.data() + mStructPtr[ c ];
                const CT* B = cb[ c ].data();

                for( long long jj=0; jj<mc; jj++ )
                {
                    const long long cj = relpos[ cidx[ jj ] ] * nf;

                    for( long long ii=( mSym ? jj : 0 ); ii<mc; ii++ ) F[ relpos[ cidx[ ii ] ] + cj ] += B[ ii + jj * mc ];
                }

                std::vector<CT>().swap( cb[ c ] );
            }

            double* Fd = reinterpret_cast<double*>( F.data() );

            int LocalPerturbed = 0;

            // Local pivot order of the ns fully summed variables ( original local index at each position )
            std::vector<int> LocalPiv( ns );

            for( long long i=0; i<ns; i++ ) LocalPiv[ i ] = (int)i;

            // Blocked right-looking factorization of the ns pivots. Each block of NB columns is factorized
            // column by column ( updates restricted to the block ), then the whole trailing part of the front,
            // including the Schur complement F22, is updated with one call to the dense kernel.
            const long long NB = 64;

            for( long long k0=0; k0<ns; k0+=NB )
            {
                const long long k1 = std::min( ns, k0 + NB );

                for( long long k=k0; k<k1; k++ )
                {
                    // Pivoting restricted to the fully summed variables of the front ( no delayed pivots ):
                    //  LU   : partial pivoting, largest entry of column k among the fully summed rows [k,ns),
                    //         with a row interchange over the whole front ( as LAPACK getrf ).
                    //  LDL^T: symmetric 1x1 pivoting, largest diagonal entry among the block columns [k,k1),
                    //         with a symmetric interchange ( lower triangle storage ).
                    long long p = k;

                    if( !mSym )
                    {
                        double best = std::abs( F[ k + k * nf ] );

                        for( long long i=k+1; i<ns; i++ ) { double a = std::abs( F[ i + k * nf ] ); if( a > best ) { best = a; p = i; } }

                        if( p != k )
                        {
                            for( long long c=0; c<nf; c++ ) std::swap( F[ k + c * nf ], F[ p + c * nf ] );

                            std::swap( LocalPiv[ k ], LocalPiv[ p ] );
                        }
                    }
                    else
                    {
                        double best = std::abs( F[ k + k * nf ] );

                        for( long long i=k+1; i<k1; i++ ) { double a = std::abs( F[ i + i * nf ] ); if( a > best ) { best = a; p = i; } }

                        if( p != k )
                        {
                            for( long long j=0; j<k; j++ ) std::swap( F[ k + j * nf ], F[ p + j * nf ] );   // rows of L

                            std::swap( F[ k + k * nf ], F[ p + p * nf ] );                                    // diagonal

                            for( long long i=k+1; i<p; i++ ) std::swap( F[ i + k * nf ], F[ p + i * nf ] );  // (i,k) <-> (p,i)

                            for( long long i=p+1; i<nf; i++ ) std::swap( F[ i + k * nf ], F[ i + p * nf ] ); // (i,k) <-> (i,p)

                            std::swap( LocalPiv[ k ], LocalPiv[ p ] );
                        }
                    }

                    CT d = F[ k + k * nf ];

                    if( !std::isfinite( d.real() ) || !std::isfinite( d.imag() ) ) { Breakdown = true; }

                    if( std::abs( d ) < tau )
                    {
                        d = ( std::abs( d ) > 0.0 ) ? d / std::abs( d ) * tau : CT( tau, 0.0 );
                        F[ k + k * nf ] = d;
                        LocalPerturbed++;
                    }

                    const CT dinv = 1.0 / d;

                    if( mSym )
                    {
                        // F(i,j) -= F(i,k) * F(j,k) / d   for j in (k,k1), i >= j
                        for( long long jj=k+1; jj<k1; jj++ )
                        {
                            const CT fj = F[ jj + k * nf ] * dinv;

                            Axpy( Fd + 2 * ( jj + jj * nf ), Fd + 2 * ( jj + k * nf ), fj.real(), fj.imag(), nf - jj );
                        }

                        for( long long i=k+1; i<nf; i++ ) F[ i + k * nf ] *= dinv;
                    }
                    else
                    {
                        for( long long i=k+1; i<nf; i++ ) F[ i + k * nf ] *= dinv;

                        // F(i,j) -= L(i,k) * U(k,j)   for j in (k,k1), i > k
                        for( long long jj=k+1; jj<k1; jj++ )
                        {
                            const CT u = F[ k + jj * nf ];

                            Axpy( Fd + 2 * ( k + 1 + jj * nf ), Fd + 2 * ( k + 1 + k * nf ), u.real(), u.imag(), nf - k - 1 );
                        }
                    }
                }

                const long long nt = nf - k1;   // trailing size

                if( nt == 0 ) continue;

                const long long kb = k1 - k0;

                if( mSym )
                {
                    // W(p,j) = d_p * L(k1+j, k0+p)
                    if( (long long)W.size() < kb * nt ) W.resize( kb * nt );

                    for( long long jj=0; jj<nt; jj++ )
                        for( long long p=0; p<kb; p++ ) W[ p + jj * kb ] = F[ ( k0 + p ) + ( k0 + p ) * nf ] * F[ k1 + jj + ( k0 + p ) * nf ];

                    SchurUpdate( Fd + 2 * ( k1 + k1 * nf ), nf, Fd + 2 * ( k1 + k0 * nf ), nf,
                                 reinterpret_cast<double*>( W.data() ), kb, nt, nt, kb, true, KThreads );
                }
                else
                {
                    // U rows of the block for the trailing columns: unit lower triangular solve with L11 of the block
                    const bool Par = ( KThreads > 1 ) && ( nt * kb * kb > 200000 );

                    #pragma omp parallel for schedule(dynamic,16) num_threads( KThreads ) if( Par )
                    for( long long jj=k1; jj<nf; jj++ )
                    {
                        for( long long k=k0; k<k1; k++ )
                        {
                            const CT u = F[ k + jj * nf ];

                            if( u == 0.0 ) continue;

                            Axpy( Fd + 2 * ( k + 1 + jj * nf ), Fd + 2 * ( k + 1 + k * nf ), u.real(), u.imag(), k1 - k - 1 );
                        }
                    }

                    SchurUpdate( Fd + 2 * ( k1 + k1 * nf ), nf, Fd + 2 * ( k1 + k0 * nf ), nf,
                                 Fd + 2 * ( k0 + k1 * nf ), nf, nt, nt, kb, false, KThreads );
                }
            }

            if( LocalPerturbed > 0 )
            {
                #pragma omp atomic
                mPerturbed += LocalPerturbed;
            }

            for( long long i=0; i<ns; i++ ) mPiv[ f + i ] = LocalPiv[ i ];

            // Store the factors
            CT* Lp = mLx.data() + mLptr[ s ];

            for( long long c=0; c<ns; c++ ) std::copy( F.begin() + c * nf, F.begin() + c * nf + nf, Lp + c * nf );

            if( !mSym )
            {
                CT* Up = mUx.data() + mUptr[ s ];

                for( long long jj=0; jj<m; jj++ ) std::copy( F.begin() + ( ns + jj ) * nf, F.begin() + ( ns + jj ) * nf + ns, Up + jj * ns );
            }

            // Contribution block for the parent
            if( ( m > 0 ) && ( mSparent[ s ] >= 0 ) )
            {
                std::vector<CT> B( m * m );

                for( long long jj=0; jj<m; jj++ )
                    std::copy( F.begin() + ns + ( ns + jj ) * nf, F.begin() + ns + ( ns + jj ) * nf + m, B.begin() + jj * m );

                cb[ s ].swap( B );
            }

            for( long long i=0; i<ns; i++ ) relpos[ f + i ] = -1;
            for( long long i=0; i<m ; i++ ) relpos[ sidx[ i ] ] = -1;

            // Progress
            #pragma omp atomic
            DoneOps += nodeOps[ s ];

            if( Verbose && ( omp_get_thread_num() == 0 ) && ( mPredOps > 0.0 ) )
            {
                double done;

                // "atomic read" needs OpenMP 3.1; Visual C++ ( OpenMP 2.0 ) reads directly ( progress display only )
#if defined(_OPENMP) && ( _OPENMP >= 201107 )
                #pragma omp atomic read
                done = DoneOps;
#else
                done = DoneOps;
#endif

                const double fr = done / mPredOps;

                if( 100.0 * fr >= NextReport )
                {
                    time_t tn; time( &tn );

                    const double el = difftime( tn, t0 );

                    std::cout << "  " << NextReport << "%   elapsed " << el << " s";

                    if( ( fr > 0.05 ) && ( fr < 1.0 ) ) std::cout << ",  estimated remaining " << el * ( 1.0 - fr ) / fr << " s";

                    std::cout << std::endl;

                    while( NextReport <= 100.0 * fr ) NextReport += 10;
                }
            }
        };

        //*************************************************************************************************
        //* Phase 1 ( several threads ): independent subtrees in parallel, one thread per subtree.
        //* Phase 2: the top of the tree ( large fronts ) in order, with parallel dense kernels.
        //*************************************************************************************************
        std::vector<char> inSubtree( mNs, 0 );
        std::vector<int>  roots;

        if( NumThreads > 1 )
        {
            const double Threshold = mPredOps / ( 8.0 * NumThreads );

            std::vector<int> work;

            for( int s=0; s<mNs; s++ ) if( mSparent[ s ] < 0 ) work.push_back( s );

            while( !work.empty() )
            {
                int s = work.back(); work.pop_back();

                if( subOps[ s ] <= Threshold )
                {
                    roots.push_back( s );

                    for( int q=sfirst[ s ]; q<=s; q++ ) inSubtree[ q ] = 1;
                }
                else
                {
                    for( int c=shead[ s ]; c!=-1; c=snext[ c ] ) work.push_back( c );
                }
            }

            // Largest subtrees first ( better load balance )
            std::sort( roots.begin(), roots.end(), [&]( int a, int b ) { return subOps[ a ] > subOps[ b ]; } );

            #pragma omp parallel num_threads( NumThreads )
            {
                std::vector<CT>  F, W;
                std::vector<int> relpos;

                try { relpos.assign( n, -1 ); } catch( std::bad_alloc& ) { OutOfMemory = true; }

                #pragma omp for schedule(dynamic,1)
                for( long long r=0; r<(long long)roots.size(); r++ )
                {
                    if( OutOfMemory ) continue;

                    try
                    {
                        const int root = roots[ r ];

                        for( int q=sfirst[ root ]; q<=root; q++ ) ProcessNode( q, F, W, relpos, 1 );
                    }
                    catch( std::bad_alloc& ) { OutOfMemory = true; }
                }
            }
        }

        if( !OutOfMemory )
        {
            try
            {
                std::vector<CT>  F, W;
                std::vector<int> relpos( n, -1 );

                F.reserve( (std::size_t)mMaxFront * (std::size_t)mMaxFront );

                for( int s=0; s<mNs; s++ )
                {
                    if( !inSubtree[ s ] ) ProcessNode( s, F, W, relpos, NumThreads );
                }
            }
            catch( std::bad_alloc& ) { OutOfMemory = true; }
        }

        if( OutOfMemory )
        {
            Clear();
            ErrorMsg = "not enough memory for the factorization";
            return false;
        }

        if( Breakdown )
        {
            ErrorMsg = "numerical breakdown ( NaN / Inf pivot )";
            return false;
        }

        // The assembly arrays are not needed by Solve()
        std::vector<int>().swap( mPAp ); std::vector<int>().swap( mPAi ); std::vector<int>().swap( mPAmap );

        return true;
    }


    //*****************************************************************************************************
    //* - b <- inv(A) * b
    //*****************************************************************************************************
    void ComplexMultifrontal::Solve( std::vector<CT>& b ) const
    {
        const int n = mN;

        std::vector<CT> y( n ), t( mMaxFront );

        for( int k=0; k<n; k++ ) y[ k ] = b[ mPerm[ k ] ];

        // Forward: L ( unit diagonal ), in the pivot order of each supernode
        for( int s=0; s<mNs; s++ )
        {
            const int f = mSfirst[ s ];
            const long long ns = mSfirst[ s + 1 ] - f;
            const long long m  = mStructPtr[ s + 1 ] - mStructPtr[ s ];
            const long long nf = ns + m;
            const int* sidx = mStructIdx.data() + mStructPtr[ s ];
            const int* piv  = mPiv.data() + f;
            const CT* L = mLx.data() + mLptr[ s ];

            for( long long k=0; k<ns; k++ ) t[ k ] = y[ f + piv[ k ] ];

            for( long long k=0; k<ns; k++ )
            {
                const CT tk = t[ k ];

                if( tk == 0.0 ) continue;

                for( long long i=k+1; i<ns; i++ ) t[ i ] -= L[ i + k * nf ] * tk;
                for( long long i=0; i<m; i++ )    y[ sidx[ i ] ] -= L[ ns + i + k * nf ] * tk;
            }

            // LU: rows were interchanged, the values now belong to the pivot ( = column ) order.
            // LDL^T: rows and columns were interchanged symmetrically, keep the original positions.
            if( mSym ) for( long long k=0; k<ns; k++ ) y[ f + piv[ k ] ] = t[ k ];
            else       for( long long k=0; k<ns; k++ ) y[ f + k ]        = t[ k ];
        }

        // Backward
        for( int s=mNs-1; s>=0; s-- )
        {
            const int f = mSfirst[ s ];
            const long long ns = mSfirst[ s + 1 ] - f;
            const long long m  = mStructPtr[ s + 1 ] - mStructPtr[ s ];
            const long long nf = ns + m;
            const int* sidx = mStructIdx.data() + mStructPtr[ s ];
            const int* piv  = mPiv.data() + f;
            const CT* L = mLx.data() + mLptr[ s ];

            if( mSym )
            {
                // D^-1 then L^T, in the pivot order
                for( long long k=0; k<ns; k++ ) t[ k ] = y[ f + piv[ k ] ] / L[ k + k * nf ];

                for( long long k=ns-1; k>=0; k-- )
                {
                    CT sum( 0.0, 0.0 );

                    for( long long i=k+1; i<ns; i++ ) sum += L[ i + k * nf ] * t[ i ];
                    for( long long i=0; i<m; i++ )    sum += L[ ns + i + k * nf ] * y[ sidx[ i ] ];

                    t[ k ] -= sum;
                }

                for( long long k=0; k<ns; k++ ) y[ f + piv[ k ] ] = t[ k ];
            }
            else
            {
                // U ( columns in natural order )
                const CT* U = mUx.data() + mUptr[ s ];

                for( long long k=ns-1; k>=0; k-- )
                {
                    CT sum( 0.0, 0.0 );

                    for( long long jj=k+1; jj<ns; jj++ ) sum += L[ k + jj * nf ] * y[ f + jj ];
                    for( long long jj=0; jj<m; jj++ )    sum += U[ k + jj * ns ] * y[ sidx[ jj ] ];

                    y[ f + k ] = ( y[ f + k ] - sum ) / L[ k + k * nf ];
                }
            }
        }

        for( int k=0; k<n; k++ ) b[ mPerm[ k ] ] = y[ k ];
    }
}
