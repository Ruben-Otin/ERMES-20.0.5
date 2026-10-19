//*******************************************************************   
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: pooyan $
//*   Date:                $Date: 2003/06/04 11:38:33 $
//*   Revision:            $Revision: 1.4 $
//*******************************************************************



#if !defined(KRATOS_CG_SOLVER_H_INCLUDED )
#define  KRATOS_CG_SOLVER_H_INCLUDED



// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "iterative_solver.h"


namespace Kratos
{
 
    template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType, class TReordererType>
    class CGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
    public:
     
      typedef boost::shared_ptr<CGSolver> Pointer;

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
     
      
      /// Default constructor.
      CGSolver(){}

      CGSolver(double NewMaxTolerance) : BaseType(NewMaxTolerance){}

      CGSolver(double NewMaxTolerance, 
               unsigned int NewMaxIterationsNumber) : BaseType(NewMaxTolerance, NewMaxIterationsNumber){}

      CGSolver(double NewMaxTolerance, 
               unsigned int NewMaxIterationsNumber, 
               typename TPreconditionerType::Pointer pNewPreconditioner) : 
	           BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Destructor.
      virtual ~CGSolver(){}
      
      
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	      if(IsNotConsistent(rA, rX, rB))
	        return false;

  	      GetPreconditioner()->Initialize(rA,rX,rB);
 	      GetPreconditioner()->ApplyInverseRight(rX);
	      GetPreconditioner()->ApplyLeft(rB);

	      bool is_solved = IterativeSolve(rA,rX,rB);

 	      GetPreconditioner()->Finalize(rX);

	      return is_solved;
      }
      
      
      /** Multi solve method for solving a set of linear systems with same coefficient matrix.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
      {
  	      GetPreconditioner()->Initialize(rA,rX,rB);

 	      bool is_solved = true;
	      VectorType x(TDenseSpaceType::Size1(rX));
	      VectorType b(TDenseSpaceType::Size1(rB));
	      for(int i = 0 ; i < TDenseSpaceType::Size2(rX) ; i++)
	      {
	          TDenseSpaceType::GetColumn(i,rX, x);
	          TDenseSpaceType::GetColumn(i,rB, b);
	          
	          GetPreconditioner()->ApplyInverseRight(x);
	          GetPreconditioner()->ApplyLeft(b);

	          is_solved &= IterativeSolve(rA,x,b);

	          GetPreconditioner()->Finalize(x);
	      }

	      return is_solved;
      }
     

      /// Return information about this object.
      String Info() const
      {
	      String buffer;
	      buffer << "Conjugate gradient linear solver with " << GetPreconditioner()->Info();
	      return  buffer;
      }
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
      {
	      OStream << "Conjugate gradient linear solver with ";
	      GetPreconditioner()->PrintInfo(OStream);
      }

      /// Print object's data. Using PrintData method of the base class.
      void  PrintData(std::ostream& OStream) const 
      {
	      BaseType::PrintData(OStream);
      }
      
            
     
   private:
      
        

        bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
        {
            const int size = TSparseSpaceType::Size(rX);

            int info_count    = 0;
            mIterationsNumber = 0;

            VectorType r(size);

            PreconditionedMult(rA,rX,r);
           // TSparseSpaceType::Mult(rA, rX, r);

            TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

            mBNorm = TSparseSpaceType::TwoNorm(rB);

            VectorType p(r);
            VectorType q(size);

            double roh0 = TSparseSpaceType::Dot(r, r);
            double roh1 = roh0;
            double beta = 0;

            if(roh0 == 0.00)
            return false;

            do
            {
                PreconditionedMult(rA,p,q);
                //TSparseSpaceType::Mult(rA, p, q);

                double pq = TSparseSpaceType::Dot(p,q);

                if(pq == 0.00)
                break;

                double alpha = roh0 / pq;

                TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
                TSparseSpaceType::ScaleAndAdd(-alpha, q, 1.00, r);

                roh1 = TSparseSpaceType::Dot(r,r);

                beta = (roh1 / roh0);
                TSparseSpaceType::ScaleAndAdd(1.00, r, beta, p);

                roh0 = roh1;

                mResidualNorm = sqrt(roh1);
                mIterationsNumber++;

                info_count++;
                if (info_count==100)
                {
                   std::cout<<"it = "<<mIterationsNumber<<"   res = "<<(mResidualNorm/mBNorm)<<std::endl;
                   info_count=0;
                }

            } while(IterationNeeded() && (roh0 != 0.00));

            std::cout<<std::endl;

            return IsConverged();
        }
        
     
        /// Assignment operator.
        CGSolver& operator=(const CGSolver& Other);

        /// Copy constructor.
        CGSolver(const CGSolver& Other);
   
    }; // Class CGSolver 

  
 
    /// input stream function
    template<class TSparseSpaceType, 
             class TDenseSpaceType, 
             class TPreconditionerType, 
             class TReordererType>
    inline std::istream& operator >> (std::istream& IStream, 
				      CGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
    }

     /// output stream function
    template<class TSparseSpaceType, 
             class TDenseSpaceType, 
             class TPreconditionerType, 
             class TReordererType>
    inline std::ostream& operator << (std::ostream& OStream, 
				    const CGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
        rThis.PrintInfo(OStream);
        OStream << std::endl;
        rThis.PrintData(OStream);
        return OStream;
    }

  
  
}  // namespace Kratos.

#endif // KRATOS_CG_SOLVER_H_INCLUDED  defined 


