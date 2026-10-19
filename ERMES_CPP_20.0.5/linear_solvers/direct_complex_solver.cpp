
#include <iostream>
#include <algorithm>
#include <set>
#include <cmath>
#include <ctime>
#include <limits>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "direct_complex_solver.h"
#include "multifrontal_complex_solver.h"

namespace Kratos
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ComplexSparseLU
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ComplexSparseLU::SetMatrix( int n, std::vector<int>& Ap, std::vector<int>& Ai, std::vector<ComplexType>& Ax )
    {
        Clear();

        mN = n;

        mAp.swap( Ap );
        mAi.swap( Ai );
        mAx.swap( Ax );
    }

    void ComplexSparseLU::Clear()
    {
        mN = 0;
        mOffDiagPivots = 0;

        std::vector<int>().swap( mAp ); std::vector<int>().swap( mAi ); std::vector<ComplexType>().swap( mAx );
        std::vector<int>().swap( mQ  ); std::vector<int>().swap( mPinv );
        std::vector<int>().swap( mLp ); std::vector<int>().swap( mLi ); std::vector<ComplexType>().swap( mLx );
        std::vector<int>().swap( mUp ); std::vector<int>().swap( mUi ); std::vector<ComplexType>().swap( mUx );
    }

    namespace
    {
        //*************************************************************************************************
        //* Weighted approximate minimum degree on a graph given by sorted adjacency lists ( no self loops ).
        //* w[i] = number of matrix variables represented by node i ( see ApproximateMinimumDegree ).
        //*
        //* Quotient graph: every eliminated node becomes an "element" ( the clique of its neighbours ) stored
        //* as a list of nodes; elements adjacent to the pivot are absorbed, so memory stays O( edges ).
        //* Invariant: node v is in element e  <=>  e is in the element list of v.
        //* Degree of i ( AMD approximate external degree, weighted ):
        //*     w( adjacent nodes ) + ( w(Lp) - w_i ) + sum over other adjacent elements e of w( Le \ Lp )
        //* Elements contained in the new element are absorbed ( aggressive absorption ).
        //* Degree lists are buckets ( doubly linked lists ), so updates cost O(1).
        //*************************************************************************************************
        void WeightedAMD( int n, std::vector< std::vector<int> >& vars, const std::vector<int>& w, std::vector<int>& order )
        {
            std::vector< std::vector<int> > elems( n ), evars( n );

            std::vector<long long> ew( n, 0 );        // element weights
            std::vector<char>      state( n, 0 );     // 0 node, 1 live element, 2 absorbed element
            std::vector<long long> deg( n, 0 );
            std::vector<int>       mark( n, -1 ), wtag( n, -1 );
            std::vector<long long> wext( n, 0 );

            long long TotalW = 0;

            for( int i=0; i<n; i++ ) TotalW += w[ i ];

            // Degree buckets ( degree clamped to n_buckets - 1 )
            const long long nb = TotalW + 1;
            std::vector<int> head( nb, -1 ), next( n, -1 ), prev( n, -1 );
            long long mindeg = 0;

            auto Insert = [&]( int i )
            {
                long long d = deg[ i ];
                next[ i ] = head[ d ]; prev[ i ] = -1;
                if( head[ d ] >= 0 ) prev[ head[ d ] ] = i;
                head[ d ] = i;
                if( d < mindeg ) mindeg = d;
            };

            auto Remove = [&]( int i )
            {
                if( prev[ i ] >= 0 ) next[ prev[ i ] ] = next[ i ]; else head[ deg[ i ] ] = next[ i ];
                if( next[ i ] >= 0 ) prev[ next[ i ] ] = prev[ i ];
            };

            for( int i=0; i<n; i++ )
            {
                long long d = 0;
                for( std::size_t a=0; a<vars[ i ].size(); a++ ) d += w[ vars[ i ][ a ] ];
                deg[ i ] = d;
                Insert( i );
            }

            order.clear();
            order.reserve( n );

            std::vector<int> Lp;

            long long RemainingW = TotalW;

            for( int k=0; k<n; k++ )
            {
                while( head[ mindeg ] < 0 ) mindeg++;

                const int p = head[ mindeg ];

                Remove( p );

                order.push_back( p );

                RemainingW -= w[ p ];

                mark[ p ] = k;

                // New element Lp = nodes adjacent to p directly or through its elements
                Lp.clear();

                long long LpW = 0;

                for( std::size_t a=0; a<vars[ p ].size(); a++ )
                {
                    int v = vars[ p ][ a ];
                    if( ( state[ v ] == 0 ) && ( mark[ v ] != k ) ) { mark[ v ] = k; Lp.push_back( v ); LpW += w[ v ]; }
                }

                for( std::size_t a=0; a<elems[ p ].size(); a++ )
                {
                    int e = elems[ p ][ a ];

                    if( state[ e ] != 1 ) continue;

                    for( std::size_t b=0; b<evars[ e ].size(); b++ )
                    {
                        int v = evars[ e ][ b ];
                        if( ( state[ v ] == 0 ) && ( mark[ v ] != k ) ) { mark[ v ] = k; Lp.push_back( v ); LpW += w[ v ]; }
                    }

                    state[ e ] = 2;
                    std::vector<int>().swap( evars[ e ] );
                }

                state[ p ] = 1;
                evars[ p ] = Lp;
                ew[ p ] = LpW;

                std::vector<int>().swap( vars [ p ] );
                std::vector<int>().swap( elems[ p ] );

                // Pass 1: remove covered node edges and absorbed elements; w( Le \ Lp ) for adjacent elements
                for( std::size_t a=0; a<Lp.size(); a++ )
                {
                    const int i = Lp[ a ];

                    std::vector<int>& vi = vars[ i ];
                    std::size_t q = 0;
                    for( std::size_t b=0; b<vi.size(); b++ ) { int v = vi[ b ]; if( ( state[ v ] == 0 ) && ( mark[ v ] != k ) ) vi[ q++ ] = v; }
                    vi.resize( q );

                    std::vector<int>& ei = elems[ i ];
                    q = 0;
                    for( std::size_t b=0; b<ei.size(); b++ )
                    {
                        const int e = ei[ b ];
                        if( state[ e ] != 1 ) continue;
                        ei[ q++ ] = e;
                        if( wtag[ e ] != k ) { wtag[ e ] = k; wext[ e ] = ew[ e ]; }
                        wext[ e ] -= w[ i ];
                    }
                    ei.resize( q );
                }

                // Pass 2: degrees ( elements inside Lp are absorbed )
                for( std::size_t a=0; a<Lp.size(); a++ )
                {
                    const int i = Lp[ a ];

                    long long d = LpW - w[ i ];

                    for( std::size_t b=0; b<vars[ i ].size(); b++ ) d += w[ vars[ i ][ b ] ];

                    std::vector<int>& ei = elems[ i ];
                    std::size_t q = 0;

                    for( std::size_t b=0; b<ei.size(); b++ )
                    {
                        const int e = ei[ b ];
                        if( state[ e ] != 1 ) continue;
                        if( wext[ e ] <= 0 ) { state[ e ] = 2; std::vector<int>().swap( evars[ e ] ); continue; }
                        d += wext[ e ];
                        ei[ q++ ] = e;
                    }

                    ei.resize( q );
                    ei.push_back( p );

                    const long long cap = RemainingW - w[ i ];
                    if( d > cap ) d = cap;
                    if( d < 0 ) d = 0;

                    Remove( i );
                    deg[ i ] = d;
                    Insert( i );
                }
            }
        }
    }

    //*****************************************************************************************************
    //* - Approximate minimum degree ordering on the pattern of A + A^T ( Q[k] = k-th pivot ).
    //*
    //*   Graph compression: variables with identical closed neighbourhoods ( e.g. the 3 or 4 unknowns of
    //*   the same FEM node: Ax, Ay, Az, V ) are merged into one weighted node before ordering, and expanded
    //*   afterwards. This makes the ordering many times faster on FEM matrices ( the work grows with the
    //*   square of the unknowns per node ) and gives the same kind of ordering.
    //*****************************************************************************************************
    void ApproximateMinimumDegree( int N, const std::vector<int>& Ap, const std::vector<int>& Ai, std::vector<int>& Q )
    {
        const int n = N;

        // Symmetric adjacency of A + A^T without diagonal ( sorted, unique )
        std::vector< std::vector<int> > adj( n );

        for( int j=0; j<n; j++ )
            for( int p=Ap[ j ]; p<Ap[ j + 1 ]; p++ )
            {
                int i = Ai[ p ];
                if( i == j ) continue;
                adj[ i ].push_back( j );
                adj[ j ].push_back( i );
            }

        for( int i=0; i<n; i++ )
        {
            std::sort( adj[ i ].begin(), adj[ i ].end() );
            adj[ i ].erase( std::unique( adj[ i ].begin(), adj[ i ].end() ), adj[ i ].end() );
        }

        // Indistinguishable variables: equal closed neighbourhood adj(i) U {i}. Hash, then verify.
        std::vector< std::pair<unsigned long long, int> > keys( n );

        for( int i=0; i<n; i++ )
        {
            unsigned long long h = 1469598103934665603ULL ^ (unsigned long long)( adj[ i ].size() + 1 );
            unsigned long long s = (unsigned long long)i;
            for( std::size_t a=0; a<adj[ i ].size(); a++ ) s += (unsigned long long)adj[ i ][ a ];
            h = ( h ^ s ) * 1099511628211ULL;
            keys[ i ] = std::make_pair( h, i );
        }

        std::sort( keys.begin(), keys.end() );

        std::vector<int> rep( n, -1 );   // representative of each variable

        auto SameClosed = [&]( int a, int b )
        {
            const std::vector<int>& A = adj[ a ];
            const std::vector<int>& B = adj[ b ];
            if( A.size() != B.size() ) return false;
            // closed neighbourhoods: A U {a} == B U {b}  <=>  A \ {b} == B \ {a} and a in B, b in A
            if( !std::binary_search( A.begin(), A.end(), b ) || !std::binary_search( B.begin(), B.end(), a ) ) return false;
            std::size_t x = 0, y = 0;
            while( ( x < A.size() ) && ( y < B.size() ) )
            {
                if( A[ x ] == b ) { x++; continue; }
                if( B[ y ] == a ) { y++; continue; }
                if( A[ x ] != B[ y ] ) return false;
                x++; y++;
            }
            return true;
        };

        for( int s=0; s<n; )
        {
            int e = s;
            while( ( e < n ) && ( keys[ e ].first == keys[ s ].first ) ) e++;

            for( int a=s; a<e; a++ )
            {
                const int va = keys[ a ].second;
                if( rep[ va ] >= 0 ) continue;
                rep[ va ] = va;
                for( int b=a+1; b<e; b++ )
                {
                    const int vb = keys[ b ].second;
                    if( ( rep[ vb ] < 0 ) && SameClosed( va, vb ) ) rep[ vb ] = va;
                }
            }

            s = e;
        }

        // Compressed graph
        std::vector<int> cid( n, -1 ), crep;

        for( int i=0; i<n; i++ ) if( rep[ i ] == i ) { cid[ i ] = crep.size(); crep.push_back( i ); }

        for( int i=0; i<n; i++ ) cid[ i ] = cid[ rep[ i ] ];

        const int nc = crep.size();

        std::vector<int> w( nc, 0 );

        for( int i=0; i<n; i++ ) w[ cid[ i ] ]++;

        std::vector< std::vector<int> > cadj( nc );

        for( int c=0; c<nc; c++ )
        {
            const std::vector<int>& A = adj[ crep[ c ] ];
            for( std::size_t a=0; a<A.size(); a++ ) { int d = cid[ A[ a ] ]; if( d != c ) cadj[ c ].push_back( d ); }
            std::sort( cadj[ c ].begin(), cadj[ c ].end() );
            cadj[ c ].erase( std::unique( cadj[ c ].begin(), cadj[ c ].end() ), cadj[ c ].end() );
        }

        std::vector< std::vector<int> >().swap( adj );

        std::vector<int> corder;

        WeightedAMD( nc, cadj, w, corder );

        // Expand: members of each compressed node are consecutive in the order
        std::vector< std::vector<int> > members( nc );

        for( int i=0; i<n; i++ ) members[ cid[ i ] ].push_back( i );

        Q.clear();
        Q.reserve( n );

        for( int k=0; k<nc; k++ )
        {
            const std::vector<int>& M = members[ corder[ k ] ];
            Q.insert( Q.end(), M.begin(), M.end() );
        }
    }


    void ComplexSparseLU::MinimumDegreeOrdering()
    {
        ApproximateMinimumDegree( mN, mAp, mAi, mQ );
    }

    //*****************************************************************************************************
    //* - Non-zero pattern of x = L \ A(:,col) in topological order ( xi[top..n-1] ). Iterative DFS on the
    //*   graph of the L columns computed so far ( row i -> column pinv[i] ).
    //*****************************************************************************************************
    int ComplexSparseLU::Reach( int col, std::vector<int>& xi, std::vector<int>& mark, std::vector<int>& pstack, int stamp )
    {
        const int n = mN;

        int top = n;

        // xi[0..n-1] is used as DFS stack ( bottom part ) and output ( top part ); they never overlap
        for( int p=mAp[ col ]; p<mAp[ col + 1 ]; p++ )
        {
            int start = mAi[ p ];

            if( mark[ start ] == stamp ) continue;

            int head = 0;

            xi[ 0 ] = start;

            while( head >= 0 )
            {
                int i = xi[ head ];
                int j = mPinv[ i ];

                if( mark[ i ] != stamp )
                {
                    mark  [ i ] = stamp;
                    pstack[ head ] = ( j < 0 ) ? 0 : mLp[ j ] + 1;   // skip unit diagonal
                }

                bool done = true;

                if( j >= 0 )
                {
                    int pend = mLp[ j + 1 ];

                    for( int q=pstack[ head ]; q<pend; q++ )
                    {
                        int r = mLi[ q ];

                        if( mark[ r ] == stamp ) continue;

                        pstack[ head ] = q + 1;

                        xi[ ++head ] = r;

                        done = false;

                        break;
                    }
                }

                if( done )
                {
                    head--;

                    xi[ --top ] = i;
                }
            }
        }

        return top;
    }

    //*****************************************************************************************************
    //* - Symbolic Cholesky analysis of the pattern of A + A^T in the minimum degree order ( Q ). Gives the
    //*   predicted number of entries of L ( U has the same pattern when the diagonal pivots are kept ) and the
    //*   work per column ( ~ squared column count ), used for the memory check and the progress report.
    //*   Elimination tree with path compression ( Liu ) and column counts from the row subtrees, O( nnz(L) )
    //*   integer operations.
    //*****************************************************************************************************
    long long ComplexSparseLU::SymbolicAnalysis( std::vector<double>& ColWork )
    {
        const int n = mN;

        std::vector<int> pinvQ( n );

        for( int k=0; k<n; k++ ) pinvQ[ mQ[ k ] ] = k;

        // Lower pattern of the permuted A + A^T, row by row ( CSR )
        std::vector<int> rp( n + 1, 0 );

        for( int j=0; j<n; j++ )
        {
            for( int p=mAp[ j ]; p<mAp[ j + 1 ]; p++ )
            {
                int a = pinvQ[ mAi[ p ] ], b = pinvQ[ j ];

                if( a != b ) rp[ std::max( a, b ) + 1 ]++;
            }
        }

        for( int k=0; k<n; k++ ) rp[ k + 1 ] += rp[ k ];

        std::vector<int> ri( rp[ n ] ), pos( rp.begin(), rp.end() - 1 );

        for( int j=0; j<n; j++ )
        {
            for( int p=mAp[ j ]; p<mAp[ j + 1 ]; p++ )
            {
                int a = pinvQ[ mAi[ p ] ], b = pinvQ[ j ];

                if( a != b ) ri[ pos[ std::max( a, b ) ]++ ] = std::min( a, b );
            }
        }

        // Elimination tree
        std::vector<int> parent( n, -1 ), ancestor( n, -1 );

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

        // Column counts ( below the diagonal ) from the row subtrees
        std::vector<long long> count( n, 0 );
        std::vector<int>       mark ( n, -1 );

        for( int k=0; k<n; k++ )
        {
            mark[ k ] = k;

            for( int p=rp[ k ]; p<rp[ k + 1 ]; p++ )
            {
                for( int r = ri[ p ]; ( r != -1 ) && ( mark[ r ] != k ); r = parent[ r ] )
                {
                    count[ r ]++;
                    mark [ r ] = k;
                }
            }
        }

        ColWork.assign( n, 0.0 );

        long long NnzL = 0;

        for( int k=0; k<n; k++ )
        {
            NnzL += 1 + count[ k ];

            ColWork[ k ] = 1.0 + (double)count[ k ] * (double)count[ k ];
        }

        return NnzL;
    }

    //*****************************************************************************************************
    //* - Physical memory of the machine in bytes ( 0 if unknown )
    //*****************************************************************************************************
    double ComplexSparseLU::PhysicalMemoryBytes()
    {
#if defined(_WIN32)
        MEMORYSTATUSEX status;
        status.dwLength = sizeof( status );
        if( GlobalMemoryStatusEx( &status ) ) return (double)status.ullTotalPhys;
        return 0.0;
#else
        long pages = sysconf( _SC_PHYS_PAGES );
        long psize = sysconf( _SC_PAGE_SIZE  );
        if( ( pages > 0 ) && ( psize > 0 ) ) return (double)pages * (double)psize;
        return 0.0;
#endif
    }

    //*****************************************************************************************************
    //* - P*A*Q = L*U
    //*****************************************************************************************************
    bool ComplexSparseLU::Factorize( double PivotTol, std::string& ErrorMsg )
    {
        const int n = mN;

        time_t tStart; time( &tStart );

        MinimumDegreeOrdering();

        // Symbolic analysis: predicted size and cost of the factors ( assuming diagonal pivots )
        std::vector<double> ColWork;

        long long PredNnzL = SymbolicAnalysis( ColWork );

        double TotalWork = 0.0;

        for( int k=0; k<n; k++ ) TotalWork += ColWork[ k ];

        const double PredBytes = 2.0 * (double)PredNnzL * ( sizeof( ComplexType ) + sizeof( int ) );
        const double RAMBytes  = PhysicalMemoryBytes();

        time_t tSym; time( &tSym );

        std::cout << "Ordering + symbolic analysis: " << difftime( tSym, tStart ) << " s" << std::endl;
        std::cout << "Predicted nnz(L+U) = " << 2 * PredNnzL - n << "  ( " << PredBytes / 1073741824.0 << " GB, fill ratio "
                  << ( 2.0 * PredNnzL - n ) / double( NnzA() ) << " ),  operations ~ " << TotalWork << std::endl;

        // Rough single-core estimate ( ~2e8 complex operations per second measured on a 3D FEM-like pattern )
        std::cout << "Estimated factorization time ~ " << TotalWork / 2.0e8 << " s ( rough, single core )" << std::endl;

        if( RAMBytes > 0.0 )
        {
            std::cout << "Physical memory    = " << RAMBytes / 1073741824.0 << " GB" << std::endl;

            if( PredBytes > 0.8 * RAMBytes )
            {
                ErrorMsg = "the LU factors would need more memory than available ( see the prediction above ). "
                           "The problem is too large for the internal direct solver: use Quasi_Minimal_Residual, "
                           "Bi_Conjugate_Gradient or the External_solver ( e.g. PETSc + MUMPS )";

                return false;
            }
        }

        std::cout << std::endl << "Factorizing ( progress in % of the predicted operations )..." << std::endl;

        mPinv.assign( n, -1 );

        mLp.assign( n + 1, 0 );
        mUp.assign( n + 1, 0 );

        // Reserve the predicted size: avoids repeated reallocations, which temporarily double the memory
        const std::size_t guess = (std::size_t)( 1.05 * PredNnzL ) + n;

        mLi.clear(); mLx.clear(); mLi.reserve( guess ); mLx.reserve( guess );
        mUi.clear(); mUx.clear(); mUi.reserve( guess ); mUx.reserve( guess );

        double DoneWork = 0.0;
        int    NextReport = 5;

        time_t tFact0; time( &tFact0 );

        mOffDiagPivots = 0;

        std::vector<ComplexType> x( n, ComplexType( 0.0, 0.0 ) );

        std::vector<int> xi( n ), mark( n, -1 ), pstack( n );

        for( int k=0; k<n; k++ )
        {
            mLp[ k ] = mLi.size();
            mUp[ k ] = mUi.size();

            const int col = mQ[ k ];

            // Progress report every 5% of the predicted work ( the last columns are the most expensive )
            DoneWork += ColWork[ k ];

            if( ( TotalWork > 0.0 ) && ( 100.0 * DoneWork / TotalWork >= NextReport ) )
            {
                time_t tNow; time( &tNow );

                const double Elapsed  = difftime( tNow, tFact0 );
                const double Fraction = DoneWork / TotalWork;

                std::cout << "  " << NextReport << "%   elapsed " << Elapsed << " s";

                if( ( Fraction > 0.02 ) && ( Fraction < 1.0 ) ) std::cout << ",  estimated remaining " << Elapsed * ( 1.0 - Fraction ) / Fraction << " s";

                std::cout << std::endl;

                while( NextReport <= 100.0 * Fraction ) NextReport += 5;
            }

            // Pattern and values of x = L \ A(:,col)
            const int top = Reach( col, xi, mark, pstack, k );

            for( int p=mAp[ col ]; p<mAp[ col + 1 ]; p++ ) x[ mAi[ p ] ] = mAx[ p ];

            for( int t=top; t<n; t++ )
            {
                const int i = xi[ t ];
                const int j = mPinv[ i ];

                if( j < 0 ) continue;

                const ComplexType xj = x[ i ];

                for( int p=mLp[ j ] + 1; p<mLp[ j + 1 ]; p++ ) x[ mLi[ p ] ] -= mLx[ p ] * xj;
            }

            // U part ( already pivoted rows ) and pivot search
            int    ipiv = -1;
            double amax = -1.0;

            for( int t=top; t<n; t++ )
            {
                const int i = xi[ t ];

                if( mPinv[ i ] < 0 )
                {
                    double a = std::abs( x[ i ] );

                    if( a > amax ) { amax = a; ipiv = i; }
                }
                else
                {
                    mUi.push_back( mPinv[ i ] );
                    mUx.push_back( x[ i ] );
                }
            }

            if( ( ipiv < 0 ) || !( amax > 0.0 ) )
            {
                ErrorMsg = "matrix is structurally or numerically singular ( check boundary conditions and unconnected regions )";

                return false;
            }

            // Prefer the diagonal ( keeps the fill-reducing ordering )
            if( ( mPinv[ col ] < 0 ) && ( std::abs( x[ col ] ) >= PivotTol * amax ) )
            {
                ipiv = col;
            }
            else
            {
                mOffDiagPivots++;
            }

            const ComplexType pivot = x[ ipiv ];

            mUi.push_back( k );
            mUx.push_back( pivot );

            mPinv[ ipiv ] = k;

            mLi.push_back( ipiv );
            mLx.push_back( ComplexType( 1.0, 0.0 ) );

            for( int t=top; t<n; t++ )
            {
                const int i = xi[ t ];

                if( mPinv[ i ] < 0 )
                {
                    mLi.push_back( i );
                    mLx.push_back( x[ i ] / pivot );
                }

                x[ i ] = ComplexType( 0.0, 0.0 );
            }
        }

        mLp[ n ] = mLi.size();
        mUp[ n ] = mUi.size();

        // L rows to pivot order
        for( std::size_t p=0; p<mLi.size(); p++ ) mLi[ p ] = mPinv[ mLi[ p ] ];

        return true;
    }

    //*****************************************************************************************************
    //* - b <- inv(A)*b   ( x = Q * inv(U) * inv(L) * P * b )
    //*****************************************************************************************************
    void ComplexSparseLU::Solve( std::vector<ComplexType>& b )
    {
        const int n = mN;

        std::vector<ComplexType> y( n );

        for( int i=0; i<n; i++ ) y[ mPinv[ i ] ] = b[ i ];

        // L y = P b
        for( int j=0; j<n; j++ )
        {
            const ComplexType yj = y[ j ];

            if( yj == 0.0 ) continue;

            for( int p=mLp[ j ] + 1; p<mLp[ j + 1 ]; p++ ) y[ mLi[ p ] ] -= mLx[ p ] * yj;
        }

        // U z = y
        for( int j=n-1; j>=0; j-- )
        {
            y[ j ] /= mUx[ mUp[ j + 1 ] - 1 ];

            const ComplexType yj = y[ j ];

            if( yj == 0.0 ) continue;

            for( int p=mUp[ j ]; p<mUp[ j + 1 ] - 1; p++ ) y[ mUi[ p ] ] -= mUx[ p ] * yj;
        }

        for( int k=0; k<n; k++ ) b[ mQ[ k ] ] = y[ k ];
    }

    //*****************************************************************************************************
    //* - r <- b - A*x
    //*****************************************************************************************************
    void ComplexSparseLU::Residual( const std::vector<ComplexType>& b, const std::vector<ComplexType>& x, std::vector<ComplexType>& r )
    {
        r = b;

        for( int j=0; j<mN; j++ )
        {
            const ComplexType xj = x[ j ];

            for( int p=mAp[ j ]; p<mAp[ j + 1 ]; p++ ) r[ mAi[ p ] ] -= mAx[ p ] * xj;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ERMES matrices -> CSC
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        typedef std::pair< int, std::complex<double> > RowValue;

        void AddRowMatrix( CSRMatrix< std::complex<double> >& A, DirectComplexSolver::StorageType s, int n,
                           std::vector< std::vector<RowValue> >& cols )
        {
            for( int i=0; i<n; i++ )
            {
                std::vector< std::pair< unsigned int, std::complex<double> > >& row = A[ i ];

                for( std::size_t k=0; k<row.size(); k++ )
                {
                    const int j = row[ k ].first;

                    const std::complex<double> a = row[ k ].second;

                    cols[ j ].push_back( RowValue( i, a ) );

                    if( ( j != i ) && ( s == DirectComplexSolver::SYMMETRIC_UPPER ) ) cols[ i ].push_back( RowValue( j, a ) );
                    if( ( j != i ) && ( s == DirectComplexSolver::HERMITIAN_UPPER ) ) cols[ i ].push_back( RowValue( j, std::conj( a ) ) );
                }
            }
        }

        bool RowLess( const RowValue& a, const RowValue& b ) { return a.first < b.first; }
    }

    void BuildFullCSC( CSRMatrix< std::complex<double> >& A, DirectComplexSolver::StorageType sA,
                       CSRMatrix< std::complex<double> >* pA_aux, DirectComplexSolver::StorageType sAux,
                       int n, std::vector<int>& Ap, std::vector<int>& Ai, std::vector< std::complex<double> >& Ax )
    {
        std::vector< std::vector<RowValue> > cols( n );

        AddRowMatrix( A, sA, n, cols );

        if( pA_aux != NULL ) AddRowMatrix( *pA_aux, sAux, n, cols );

        Ap.assign( n + 1, 0 );
        Ai.clear();
        Ax.clear();

        std::size_t nnz = 0;

        for( int j=0; j<n; j++ ) nnz += cols[ j ].size();

        Ai.reserve( nnz );
        Ax.reserve( nnz );

        for( int j=0; j<n; j++ )
        {
            std::vector<RowValue>& c = cols[ j ];

            std::sort( c.begin(), c.end(), RowLess );

            // Sum duplicates ( entries present in both A and A_aux )
            for( std::size_t k=0; k<c.size(); k++ )
            {
                if( ( k > 0 ) && ( c[ k ].first == c[ k - 1 ].first ) ) Ax.back() += c[ k ].second;
                else                                                   { Ai.push_back( c[ k ].first ); Ax.push_back( c[ k ].second ); }
            }

            std::vector<RowValue>().swap( c );

            Ap[ j + 1 ] = Ai.size();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // DirectComplexSolver
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    void DirectComplexSolver::Solve( ComplexMatrixType& A, ComplexVectorType& b, ComplexVectorType& x )
    {
        if( mFullStorage ) SolveCombined( A, FULL           , NULL, FULL, b, x, "Direct sparse LU (Full matrix)" );
        else               SolveCombined( A, SYMMETRIC_UPPER, NULL, FULL, b, x, "Direct sparse LU" );
    }

    void DirectComplexSolver::Solve_HermSymm( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        SolveCombined( A, HERMITIAN_UPPER, &A_aux, SYMMETRIC_UPPER, b, x, "Direct sparse LU (Herm-Symm matrix)" );
    }

    void DirectComplexSolver::Solve_HermFull( ComplexMatrixType& A, ComplexMatrixType& A_aux, ComplexVectorType& b, ComplexVectorType& x )
    {
        SolveCombined( A, HERMITIAN_UPPER, &A_aux, FULL, b, x, "Direct sparse LU (Herm-Full matrix)" );
    }

    namespace
    {
        // r = b - A*x  ( A in CSC )
        void ResidualCSC( int n, const std::vector<int>& Ap, const std::vector<int>& Ai, const std::vector< std::complex<double> >& Ax,
                          const std::vector< std::complex<double> >& b, const std::vector< std::complex<double> >& x,
                          std::vector< std::complex<double> >& r )
        {
            r = b;

            for( int j=0; j<n; j++ )
            {
                const std::complex<double> xj = x[ j ];

                if( xj == 0.0 ) continue;

                for( int p=Ap[ j ]; p<Ap[ j + 1 ]; p++ ) r[ Ai[ p ] ] -= Ax[ p ] * xj;
            }
        }

        double Norm2( const std::vector< std::complex<double> >& v )
        {
            double s = 0.0;

            for( std::size_t i=0; i<v.size(); i++ ) s += std::norm( v[ i ] );

            return std::sqrt( s );
        }

        //*************************************************************************************************
        //* Right-preconditioned restarted GMRES for A x = b ( A in CSC ), preconditioner Prec( v ) ~ inv(A) v.
        //* With an exact factorization it converges in one iteration; with perturbed pivots it recovers the
        //* accuracy where plain iterative refinement may diverge. Returns the true residual norm.
        //*************************************************************************************************
        template< class TPrec >
        double PreconditionedGMRES( int n, const std::vector<int>& Ap, const std::vector<int>& Ai,
                                    const std::vector< std::complex<double> >& Ax, const std::vector< std::complex<double> >& b,
                                    std::vector< std::complex<double> >& x, TPrec& Prec, double RelTol, int MaxIter, int Restart,
                                    unsigned int& Iterations )
        {
            typedef std::complex<double> C;

            const double bnorm = Norm2( b );

            Iterations = 0;

            std::vector<C> r( n ), w( n );

            ResidualCSC( n, Ap, Ai, Ax, b, x, r );

            double beta = Norm2( r );

            if( !( beta > RelTol * bnorm ) || ( bnorm == 0.0 ) ) return beta;

            std::vector< std::vector<C> > V( Restart + 1, std::vector<C>( n ) ), Z( Restart, std::vector<C>( n ) );
            std::vector< std::vector<C> > H( Restart + 1, std::vector<C>( Restart, C( 0.0, 0.0 ) ) );
            std::vector<double> cs( Restart );
            std::vector<C>      sn( Restart ), g( Restart + 1 ), yv( Restart );

            while( (int)Iterations < MaxIter )
            {
                for( int i=0; i<n; i++ ) V[ 0 ][ i ] = r[ i ] / beta;

                std::fill( g.begin(), g.end(), C( 0.0, 0.0 ) );
                g[ 0 ] = beta;

                int jEnd = 0;

                for( int j=0; ( j<Restart ) && ( (int)Iterations < MaxIter ); j++ )
                {
                    Z[ j ] = V[ j ];
                    Prec( Z[ j ] );

                    // w = A * z_j
                    std::fill( w.begin(), w.end(), C( 0.0, 0.0 ) );

                    for( int c=0; c<n; c++ )
                    {
                        const C zc = Z[ j ][ c ];
                        if( zc == 0.0 ) continue;
                        for( int p=Ap[ c ]; p<Ap[ c + 1 ]; p++ ) w[ Ai[ p ] ] += Ax[ p ] * zc;
                    }

                    // Modified Gram-Schmidt ( Hermitian inner product )
                    for( int i=0; i<=j; i++ )
                    {
                        C h( 0.0, 0.0 );
                        for( int q=0; q<n; q++ ) h += std::conj( V[ i ][ q ] ) * w[ q ];
                        H[ i ][ j ] = h;
                        for( int q=0; q<n; q++ ) w[ q ] -= h * V[ i ][ q ];
                    }

                    const double hn = Norm2( w );

                    H[ j + 1 ][ j ] = hn;

                    if( hn > 0.0 ) for( int q=0; q<n; q++ ) V[ j + 1 ][ q ] = w[ q ] / hn;

                    // Previous Givens rotations on column j
                    for( int i=0; i<j; i++ )
                    {
                        const C a = H[ i ][ j ], bb = H[ i + 1 ][ j ];
                        H[ i     ][ j ] =  cs[ i ] * a + sn[ i ] * bb;
                        H[ i + 1 ][ j ] = -std::conj( sn[ i ] ) * a + cs[ i ] * bb;
                    }

                    // New rotation zeroing H(j+1,j)
                    {
                        const C a = H[ j ][ j ], bb = H[ j + 1 ][ j ];
                        const double rr = std::sqrt( std::norm( a ) + std::norm( bb ) );

                        if( rr == 0.0 ) { cs[ j ] = 1.0; sn[ j ] = 0.0; }
                        else if( std::abs( a ) == 0.0 ) { cs[ j ] = 0.0; sn[ j ] = std::conj( bb ) / std::abs( bb ); }
                        else
                        {
                            cs[ j ] = std::abs( a ) / rr;
                            sn[ j ] = ( a / std::abs( a ) ) * std::conj( bb ) / rr;
                        }

                        H[ j     ][ j ] = cs[ j ] * a + sn[ j ] * bb;
                        H[ j + 1 ][ j ] = 0.0;

                        g[ j + 1 ] = -std::conj( sn[ j ] ) * g[ j ];
                        g[ j     ] =  cs[ j ] * g[ j ];
                    }

                    Iterations++;

                    jEnd = j + 1;

                    if( !( std::abs( g[ j + 1 ] ) > RelTol * bnorm ) || ( hn == 0.0 ) ) break;
                }

                // x += Z * inv(H) * g
                for( int i=jEnd-1; i>=0; i-- )
                {
                    C s = g[ i ];
                    for( int k=i+1; k<jEnd; k++ ) s -= H[ i ][ k ] * yv[ k ];
                    yv[ i ] = s / H[ i ][ i ];
                }

                for( int i=0; i<jEnd; i++ ) for( int q=0; q<n; q++ ) x[ q ] += yv[ i ] * Z[ i ][ q ];

                ResidualCSC( n, Ap, Ai, Ax, b, x, r );

                const double NewBeta = Norm2( r );

                if( !std::isfinite( NewBeta ) ) return NewBeta;

                if( !( NewBeta > RelTol * bnorm ) ) return NewBeta;

                if( !( NewBeta < 0.999 * beta ) && ( jEnd < Restart ) ) return NewBeta;   // stagnation

                beta = NewBeta;
            }

            return beta;
        }
    }

    //*****************************************************************************************************
    //* - Equilibrate, factorize ( multifrontal LDL^T or LU ), solve with GMRES preconditioned by the
    //*   factorization; refactorize as LU if LDL^T fails.
    //*   mResEstimation = final relative residual ||Ax-b||/||b||, mIterationNumber = GMRES iterations.
    //*****************************************************************************************************
    void DirectComplexSolver::SolveCombined( ComplexMatrixType& A, StorageType sA, ComplexMatrixType* pA_aux, StorageType sAux,
                                             ComplexVectorType& b, ComplexVectorType& x, const char* SolverName )
    {
        if( !IniSolver( SolverName, A, b ) )
        {
            std::fill( x.begin(), x.end(), ComplexType( 0.0, 0.0 ) );
            FinishSolver();
            return;
        }

        const int n = mProblemSize;

        const int NumThreads = ( mNumThreads > 0 ) ? (int)mNumThreads : 1;

        std::vector<int> Ap, Ai; std::vector<ComplexType> Ax;

        BuildFullCSC( A, sA, pA_aux, sAux, n, Ap, Ai, Ax );

        //*************************************************************************************************
        //* Symmetric equilibration ( Ruiz ): As = D A D with the largest entry of every row and column ~1.
        //* Makes the static pivot threshold meaningful for matrices with very different scales ( materials,
        //* element sizes, boundary terms ). The factorization is computed for As, inv(A) = D inv(As) D.
        //*************************************************************************************************
        std::vector<double> D( n, 1.0 ), rmax( n );
        std::vector<ComplexType> As( Ax );

        for( int it=0; it<5; it++ )
        {
            std::fill( rmax.begin(), rmax.end(), 0.0 );

            for( int j=0; j<n; j++ )
                for( int p=Ap[ j ]; p<Ap[ j + 1 ]; p++ )
                {
                    const double a = std::abs( Ax[ p ] ) * D[ Ai[ p ] ] * D[ j ];
                    rmax[ Ai[ p ] ] = std::max( rmax[ Ai[ p ] ], a );
                    rmax[ j ]       = std::max( rmax[ j ], a );
                }

            for( int i=0; i<n; i++ ) if( rmax[ i ] > 0.0 ) D[ i ] /= std::sqrt( rmax[ i ] );
        }

        for( int j=0; j<n; j++ ) for( int p=Ap[ j ]; p<Ap[ j + 1 ]; p++ ) As[ p ] = Ax[ p ] * ( D[ Ai[ p ] ] * D[ j ] );

        // Complex symmetric only for the full-wave / electrostatic storage ( one symmetric matrix )
        bool Symmetric = ( sA == SYMMETRIC_UPPER ) && ( pA_aux == NULL );

        const double GB  = 1073741824.0;
        const double RAM = ComplexSparseLU::PhysicalMemoryBytes() / GB;

        ComplexMultifrontal MF;

        std::string ErrorMsg;

        // Solution: GMRES preconditioned with the factorization ( one iteration if the factorization is exact )
        struct Preconditioner
        {
            ComplexMultifrontal& rMF; const std::vector<double>& rD;

            void operator()( std::vector<ComplexType>& v )
            {
                for( std::size_t i=0; i<v.size(); i++ ) v[ i ] *= rD[ i ];
                rMF.Solve( v );
                for( std::size_t i=0; i<v.size(); i++ ) v[ i ] *= rD[ i ];
            }
        } Prec = { MF, D };

        std::vector<ComplexType> rhs( b.begin(), b.end() ), sol( n, ComplexType( 0.0, 0.0 ) );

        const double RelTol = std::min( 1.0e-10, mTolerance );

        bool Solved = false;

        mResidualNorm = std::numeric_limits<double>::quiet_NaN();

        // Attempt 1: LDL^T ( symmetric matrices ) or LU. Attempt 2 ( symmetric matrices only ): LU with partial
        // pivoting, if LDL^T broke down or GMRES did not reach the tolerance with it ( strongly indefinite
        // matrices, where symmetric 1x1 pivots are not enough ).
        for( int attempt=0; ( attempt<2 ) && !Solved; attempt++ )
        {
            if( attempt == 1 )
            {
                if( !Symmetric ) break;

                Symmetric = false;

                std::cout << std::endl << "Refactorizing as a general matrix ( LU with partial pivoting )..." << std::endl;
            }

            std::cout << "Direct solver: n = " << n << ", nnz(A) = " << Ap[ n ] << ", "
                      << ( Symmetric ? "complex symmetric ( LDL^T )" : "general ( LU )" ) << ", threads = " << NumThreads << std::endl;
            std::cout << "Ordering and symbolic analysis..." << std::endl;

            time_t t0; time( &t0 );

            MF.Analyze( n, Ap, Ai, Symmetric );

            time_t t1; time( &t1 );

            const double FactorGB = MF.FactorEntries() * sizeof( ComplexType ) / GB;
            const double WorkGB   = 2.0 * MF.PeakFrontEntries() * sizeof( ComplexType ) / GB;

            std::cout << "Ordering: " << MF.Ordering() << ",  dense kernels: " << ComplexMultifrontal::SelectDenseKernel( NumThreads ) << std::endl;
            std::cout << "Analysis time = " << difftime( t1, t0 ) << " s,  supernodes = " << MF.NumSupernodes()
                      << ",  largest front = " << MF.MaxFrontSize() << std::endl;
            std::cout << "Predicted memory: factors " << FactorGB << " GB + work " << WorkGB << " GB";
            if( RAM > 0.0 ) std::cout << "   ( physical memory " << RAM << " GB )";
            std::cout << std::endl;
            std::cout << "Operations ~ " << MF.Operations() << ",  estimated factorization time ~ "
                      << MF.Operations() / ( 1.5e9 * NumThreads ) << " s ( rough )" << std::endl;

            if( ( RAM > 0.0 ) && ( FactorGB + WorkGB > 0.8 * RAM ) )
            {
                std::cout << "Direct solver: the factorization needs more memory than available. Use Quasi_Minimal_Residual, "
                             "Bi_Conjugate_Gradient or the External_solver ( e.g. PETSc + MUMPS )." << std::endl;
                break;
            }

            std::cout << std::endl << "Factorizing..." << std::endl;

            bool ok = MF.Factorize( Ap, Ai, As, NumThreads, true, ErrorMsg );

            time_t t2; time( &t2 );

            if( !ok )
            {
                std::cout << "Direct solver: " << ErrorMsg << std::endl;
                MF.Clear();
                continue;
            }

            std::cout << "Factorization time = " << difftime( t2, t1 ) << " s,  perturbed pivots = " << MF.NumPerturbedPivots() << std::endl;

            std::fill( sol.begin(), sol.end(), ComplexType( 0.0, 0.0 ) );

            mResidualNorm = PreconditionedGMRES( n, Ap, Ai, Ax, rhs, sol, Prec, RelTol, 500, 30, mIterationNumber );

            time_t t3; time( &t3 );

            std::cout << "GMRES iterations = " << mIterationNumber << ",  ||Ax-b||/||b|| = " << mResidualNorm / mBNorm
                      << ",  solution time = " << difftime( t3, t2 ) << " s" << std::endl;

            MF.Clear();

            Solved = std::isfinite( mResidualNorm ) && ( mResidualNorm <= mTolerance * mBNorm );
        }

        for( int i=0; i<n; i++ ) x[ i ] = sol[ i ];

        mResEstimation = mResidualNorm / mBNorm;

        FinishSolver();
    }
}
