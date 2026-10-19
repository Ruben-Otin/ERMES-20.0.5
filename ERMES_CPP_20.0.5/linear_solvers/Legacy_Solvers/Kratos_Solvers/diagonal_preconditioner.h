//***************************************************************   
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: pooyan $
//*   Date:                $Date: 2003/10/06 08:45:05 $
//*   Revision:            $Revision: 1.6 $
//****************************************************************



#if !defined(KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED



// System includes 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"


namespace Kratos
{  
    //DiagonalPreconditioner class. 
    /** DiagonalPreconditioner for linesr system solvers*/
    template<class TSparseSpaceType, class TDenseSpaceType>
    class DiagonalPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
    {
    public:
      
        //Counted pointer of DiagonalPreconditioner
        typedef boost::shared_ptr<DiagonalPreconditioner> Pointer;

        typedef  Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

        typedef typename TSparseSpaceType::DataType DataType;

        typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

        typedef typename TSparseSpaceType::VectorType VectorType;

        typedef typename TDenseSpaceType::MatrixType DenseMatrixType;


        //Default constructor.
        DiagonalPreconditioner(){}

        /// Destructor.
        virtual ~DiagonalPreconditioner(){}


        /** DiagonalPreconditioner Initialize
        Initialize preconditioner for linear system rA*rX=rB
        @param rA  system matrix.
        @param rX Unknows vector
        @param rB Right side linear system of equations.
        */
        void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) 
        {
            mDiagonal.resize(TSparseSpaceType::Size(rX));
            mTemp.resize(TSparseSpaceType::Size(rX));

            int i;
            int j;

            const DataType zero = DataType();
            const DataType one = DataType(1.00);

            for(i = 0 ; i < rA.RowsNumber() ; ++i)
                if(rA(i,i) != zero)
                    mDiagonal[i] = 1.00 / sqrt(fabs(rA(i,i)));
                else
                    mDiagonal[i] = one; 
        }

        void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) 
        {
            BaseType::Initialize(rA, rX, rB);
        }

        void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
        {
            for(int i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
            mTemp[i] = rX[i] * mDiagonal[i];
            TSparseSpaceType::Mult(rA,mTemp, rY);
            ApplyLeft(rY);
        }

        void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
        {
            for(int i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
            mTemp[i] = rX[i] * mDiagonal[i];
            TSparseSpaceType::TransposeMult(rA,mTemp, rY);
            ApplyRight(rY);
        }

        VectorType& ApplyLeft(VectorType& rX)
        {
            for(int i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
            rX[i] *= mDiagonal[i];

            return rX;
        }

        VectorType& ApplyRight(VectorType& rX)
        {
            for(int i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
            rX[i] *= mDiagonal[i];

            return rX;
        }

        /** DiagonalPreconditioner transpose solver.
        Solving tranpose preconditioner system M^T*x=y, where m^T means transpose.
        @param rMatrix   DiagonalPreconditioner system matrix.
        @param rXVector  Unknows of preconditioner suystem
        @param rYVector  Right side of preconditioner system.
        */    
        VectorType& ApplyTransposeLeft(VectorType& rX)
        {
            for(int i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
            rX[i] *= mDiagonal[i];

            return rX;
        }

        VectorType& ApplyTransposeRight(VectorType& rX)
        {
            for(int i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
            rX[i] *= mDiagonal[i];

            return rX;
        }

        VectorType& ApplyInverseRight(VectorType& rX)
        {
            const DataType zero = DataType();

            for(int i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
            if(mDiagonal[i] != zero)
            rX[i] /= mDiagonal[i];

            return rX;
        }

        VectorType& Finalize(VectorType& rX)
        {
            for(int i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
            rX[i] *= mDiagonal[i];

            return rX;
        }


        /// Return information about this object.
        String Info() const
        {
            return "Diagonal preconditioner";
        }

        /// Print information about this object.
        void  PrintInfo(std::ostream& OStream) const
        {
             OStream << "Diagonal preconditioner";
        }



      
   
      
    private:
      
        VectorType mDiagonal;

        VectorType mTemp;

        /// Assignment operator.
        DiagonalPreconditioner& operator=(const DiagonalPreconditioner& Other);

        /// Copy constructor.
        DiagonalPreconditioner(const DiagonalPreconditioner& Other);
      
    }; // Class DiagonalPreconditioner 
  
 
    //input stream function
    template<class TSparseSpaceType, class TDenseSpaceType>
    inline std::istream& operator >> (std::istream& IStream, 
 				                      DiagonalPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
    {
    }


    /// output stream function
    template<class TSparseSpaceType, class TDenseSpaceType>
    inline std::ostream& operator << (std::ostream& OStream, 
 				                      const DiagonalPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
    {
        rThis.PrintInfo(OStream);
        OStream << std::endl;
        rThis.PrintData(OStream);
        return OStream;
    } 
  
}  // namespace Kratos.

#endif // KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED  defined 

