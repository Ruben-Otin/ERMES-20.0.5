//**************************************************************   
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: pooyan $
//*   Date:                $Date: 2003/06/04 11:38:33 $
//*   Revision:            $Revision: 1.4 $
//**************************************************************


#if !defined(KRATOS_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_PRECONDITIONER_H_INCLUDED



// System includes 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"


namespace Kratos
{
 
  /// Preconditioner class. 
  /* 
      Base class for preconditioners for linesr system solvers which defining
      standard interface for all the preconditioners derived from it.

      Considering a linear solver FooSolver with a method FooSolver::Solve.
      A typical code using this type of Preconditioners would be:
      
      \begin{verbatim}
       FooSolver::Solve(A,b,x,preconditioner) 
       {
          preconditioner.Initialize(A,x,b);
	  ...
	  ...
	  while(...) // Start iteration.
	  {
              preconditioner.ApplyLeft(x);
	      mult(a,x)
	      preconditioner.ApplyRight(x)
	  } // End iteration

	  preconditioner.Finalize(A,x,b);
       }
       \end{verbatim}
  */

    template<class TSparseSpaceType, class TDenseSpaceType>
    class Preconditioner
    {
    public:

        //Counted pointer of Preconditioner
        typedef boost::shared_ptr<Preconditioner> Pointer;

        typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

        typedef typename TSparseSpaceType::VectorType VectorType;

        typedef typename TDenseSpaceType::MatrixType DenseMatrixType;


        //Default constructor.
        Preconditioner(){}

        //Destructor.
        virtual ~Preconditioner(){}


        /** Preconditioner Initialize
        Initialize preconditioner for linear system rA*rX=rB
        @param rA  system matrix.
        @param rX Unknows vector
        @param rB Right side linear system of equations.
        */
        virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) {}

        virtual void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) 
        {
            VectorType x(TDenseSpaceType::Size1(rX));
            VectorType b(TDenseSpaceType::Size1(rB));

            TDenseSpaceType::GetColumn(0,rX, x);
            TDenseSpaceType::GetColumn(0,rB, b);

            Initialize(rA, x, b);
        }

        virtual void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
        {
            VectorType z = rX;
            ApplyRight(z);
            TSparseSpaceType::Mult(rA,z, rY);
            ApplyLeft(rY);
        }

        virtual void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
        {
            VectorType z = rX;
            ApplyTransposeLeft(z);
            TSparseSpaceType::TransposeMult(rA,z, rY);
            ApplyTransposeRight(rY);
        }

        virtual VectorType& ApplyLeft(VectorType& rX){return rX;}

        virtual VectorType& ApplyRight(VectorType& rX){return rX;}

        /** Preconditioner transpose solver.
        Solving tranpose preconditioner system M^T*x=y, where m^T means transpose.
        @param rX  Unknows of preconditioner suystem
        */    
        virtual VectorType& ApplyTransposeLeft(VectorType& rX)
        {
           return rX;
        }

        virtual VectorType& ApplyTransposeRight(VectorType& rX)
        {     
           return rX;
        }

        virtual VectorType& ApplyInverseRight(VectorType& rX){return rX;}

        /* The method Finalize is used to recover the value of rX.
        In principle, it is enough to multiply by the right preconditioner.
        See the diagoinal preconditioner for a nontrivial example. */
        virtual VectorType& Finalize(VectorType& rX)
        {
            return ApplyRight(rX);
        }


   

        /// Return information about this object.
        virtual String Info() const
        {
            return "Preconditioner";
        }

        /// Print information about this object.
        virtual void  PrintInfo(std::ostream& OStream) const
        {
            OStream << "Preconditioner";
        }

        virtual void PrintData(std::ostream& OStream) const
        {
        }




    private:


        //Assignment operator.
        Preconditioner& operator=(const Preconditioner& Other);

        //Copy constructor.
        Preconditioner(const Preconditioner& Other);

         
        
    }; // Class Preconditioner 
  
 
 
    //input stream function
    template<class TSparseSpaceType, class TDenseSpaceType>
    inline std::istream& operator >> (std::istream& IStream,
                                      Preconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
    {
    }

    //output stream function
    template<class TSparseSpaceType, class TDenseSpaceType>
    inline std::ostream& operator << (std::ostream& OStream, 
 				                      const Preconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
    {
        rThis.PrintInfo(OStream);
        OStream << std::endl;
        rThis.PrintData(OStream);
        return OStream;
    }
   
}  // namespace Kratos.

#endif // KRATOS_PRECONDITIONER_H_INCLUDED  defined 

