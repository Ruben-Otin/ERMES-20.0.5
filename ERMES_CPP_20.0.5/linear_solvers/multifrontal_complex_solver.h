
#if !defined(MULTIFRONTAL_COMPLEX_SOLVER_H_INCLUDED)
#define MULTIFRONTAL_COMPLEX_SOLVER_H_INCLUDED

#include <vector>
#include <complex>
#include <string>

namespace Kratos
{
    //*****************************************************************************************************
    //* - Supernodal multifrontal factorization of a sparse complex matrix ( OpenMP parallel ).
    //*
    //*   Symmetric = true : complex symmetric A ( not Hermitian ), A = L D L^T  ( half memory and work )
    //*   Symmetric = false: general A with symmetric non-zero pattern ( FEM matrices ), A = L U
    //*
    //*   Analysis: approximate minimum degree ordering, elimination tree postorder, fundamental supernodes
    //*   and relaxed amalgamation ( larger dense blocks, as in CHOLMOD / MUMPS ).
    //*   Numeric : multifrontal method. Each supernode assembles a dense frontal matrix from the original
    //*   entries and the contribution blocks of its children, factorizes its pivot columns and passes the
    //*   Schur complement to its parent. The dense updates use cache-blocked kernels, parallel with OpenMP.
    //*   Each front is factorized by blocks of 64 columns, so that almost all the work is done by the dense
    //*   kernel ( ERMES's vectorized kernel, or an external BLAS for the top fronts if it is multithreaded ).
    //*   Parallelism: independent subtrees in parallel ( one thread each ), then the top fronts with
    //*   parallel dense kernels.
    //*   Pivoting: inside each front, among its fully summed variables ( no delayed pivots ): partial
    //*   pivoting for LU, symmetric 1x1 pivoting for LDL^T. Pivots below 1e-8 x max|A| are perturbed
    //*   ( static pivoting ); DirectComplexSolver equilibrates A before and recovers the accuracy with
    //*   GMRES preconditioned by the factorization.
    //*   Ordering: METIS nested dissection when ERMES is built with METIS ( n >= 2000 ), otherwise
    //*   approximate minimum degree.
    //*****************************************************************************************************
    class ComplexMultifrontal
    {
        public:

            typedef std::complex<double> ComplexType;

            ComplexMultifrontal() : mN( 0 ), mSym( true ), mPerturbed( 0 ) {}

            // A in CSC with both triangles ( pattern symmetrized internally )
            void Analyze( int n, const std::vector<int>& Ap, const std::vector<int>& Ai, bool Symmetric );

            bool Factorize( const std::vector<int>& Ap, const std::vector<int>& Ai, const std::vector<ComplexType>& Ax,
                            int NumThreads, bool Verbose, std::string& ErrorMsg );

            // b <- inv(A) * b
            void Solve( std::vector<ComplexType>& b ) const;

            // Predictions from the analysis
            double FactorEntries() const { return mPredEntries; }
            double PeakFrontEntries() const { return mPredFront; }
            double Operations() const { return mPredOps; }
            int    NumSupernodes() const { return mNs; }
            int    MaxFrontSize() const { return mMaxFront; }

            const std::string& Ordering() const { return mOrdering; }

            int NumPerturbedPivots() const { return mPerturbed; }

            void Clear();

            // Chooses the dense kernel ( external BLAS or ERMES's ) and returns its description
            static std::string SelectDenseKernel( int NumThreads );

        private:

            int  mN;
            bool mSym;
            int  mPerturbed;

            // Permutation: new index k <-> original mPerm[ k ]
            std::vector<int> mPerm, mIperm;

            // Supernodes: columns [ mSfirst[s], mSfirst[s+1] ), parent supernode, row structure below the block
            int mNs;
            std::vector<int>       mSfirst, mSparent;
            std::vector<long long> mStructPtr;
            std::vector<int>       mStructIdx;

            // Factors. Panel of supernode s ( nf x ns, column major ): D or U11 diagonal, L11 / L21 below.
            // LU only: U12 block ( ns x m, column major ) and U11 in the upper part of the panel.
            std::vector<long long>   mLptr, mUptr;
            std::vector<ComplexType> mLx, mUx;

            // Analysis predictions
            double mPredEntries, mPredFront, mPredOps;
            int    mMaxFront;

            std::string mOrdering;

            // Permuted full matrix ( CSC ) for assembly
            std::vector<int> mPAp, mPAi, mPAmap;

            // Pivot order inside each supernode ( mPiv[ f + k ] = original local index of position k )
            std::vector<int> mPiv;
    };
}

#endif
