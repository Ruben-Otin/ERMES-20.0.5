//*****************************************************************   
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: santi $
//*   Date:                $Date: 2003/09/25 14:38:51 $
//*   Revision:            $Revision: 1.5 $
//*****************************************************************


#if !defined(KRATOS_BICGSTAB_SOLVER_H_INCLUDED )
#define  KRATOS_BICGSTAB_SOLVER_H_INCLUDED



// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "iterative_solver.h"


namespace Kratos
{


  template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType,class TReordererType>
  class BICGSTABSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
  {
    public:
     
      
      /// Counted pointer of BICGSTABSolver
      typedef boost::shared_ptr<BICGSTABSolver> Pointer;

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
   
      /// Default constructor.
      BICGSTABSolver(){}

      BICGSTABSolver(double NewTolerance) : BaseType(NewTolerance){}

      BICGSTABSolver(double NewTolerance, 
                     unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber){}

      BICGSTABSolver(double NewMaxTolerance, 
                     unsigned int NewMaxIterationsNumber, 
                     typename TPreconditionerType::Pointer pNewPreconditioner) : 
	               BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Destructor.
      virtual ~BICGSTABSolver(){}
      

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

     
      // Return information about this object.
      String Info() const
      {
	      String buffer;
	      buffer << "Biconjugate gradient stabilized linear solver with " << GetPreconditioner()->Info();
	      return  buffer;
      }
      
      // Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
      {
	      OStream << "Biconjugate gradient stabilized linear solver with ";
	      GetPreconditioner()->PrintInfo(OStream);
      }

      // Print object's data.
      void  PrintData(std::ostream& OStream) const 
      {
	      BaseType::PrintData(OStream);
      }
      
                
    private:

      
      bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
        const int size = TSparseSpaceType::Size(rX);

        int info_count = 0;
        mIterationsNumber = 0;

        VectorType r(size);

        PreconditionedMult(rA,rX,r);
        //TSparseSpaceType::Mult(rA,rX,r);

        TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

        mBNorm = TSparseSpaceType::TwoNorm(rB);

        VectorType p(r);
        VectorType s(size);
        VectorType q(size);

        VectorType rs(r); 
        VectorType qs(size); 

        double roh0 = TSparseSpaceType::Dot(r, rs);
        double roh1 = roh0;
        double alpha = 0.00;
        double beta = 0.00;
        double omega = 0.00;


        do
        {
	        PreconditionedMult(rA,p,q);
            //TSparseSpaceType::Mult(rA,p,q);

	        alpha = roh0 / TSparseSpaceType::Dot(rs,q);

	        TSparseSpaceType::ScaleAndAdd(1.00, r, -alpha, q, s);

	        PreconditionedMult(rA,s,qs);
            //TSparseSpaceType::Mult(rA,s,qs);

	        omega = TSparseSpaceType::Dot(qs,qs);

	        if(omega == 0.00)
	          break;

	        omega = TSparseSpaceType::Dot(qs,s) / omega;

	        TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
	        TSparseSpaceType::ScaleAndAdd(omega, s, 1.00, rX);
	        TSparseSpaceType::ScaleAndAdd(-omega, qs, 1.00, s, r);

	        roh1 = TSparseSpaceType::Dot(r,rs);

	        if((roh0 == 0.00) || (omega == 0.00))
	          break;
	        
	        beta = (roh1 * alpha) / (roh0 * omega);
	        
	        TSparseSpaceType::ScaleAndAdd(1.00, p, -omega, q);
	        TSparseSpaceType::ScaleAndAdd(1.00, r, beta, q, p);
	          
	        roh0 = roh1;

	        mResidualNorm =TSparseSpaceType::TwoNorm(r);
	        mIterationsNumber++;

            //Print iteration info
            info_count++;
            if (info_count==100)
            {
               std::cout<<"it = "<<mIterationsNumber<<"   res = "<<(mResidualNorm/mBNorm)<<std::endl;
               info_count=0;
            }

        }while(IterationNeeded());

        std::cout<<std::endl;

	    return IsConverged();
      }  
        
      
      /// Assignment operator.
      BICGSTABSolver& operator=(const BICGSTABSolver& Other);

      /// Copy constructor.
      BICGSTABSolver(const BICGSTABSolver& Other);
       
    }; // Class BICGSTABSolver 

 
        
 
    /// input stream function
    template<class TSparseSpaceType, 
             class TDenseSpaceType, 
             class TPreconditionerType, 
             class TReordererType>
    inline std::istream& operator >> (std::istream& IStream, 
				      BICGSTABSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
    }

    /// output stream function
    template<class TSparseSpaceType, 
             class TDenseSpaceType, 
             class TPreconditionerType, 
             class TReordererType>
    inline std::ostream& operator << (std::ostream& OStream, 
				    const BICGSTABSolver<TSparseSpaceType, TDenseSpaceType, 
				    TPreconditionerType, TReordererType>& rThis)
    {
        rThis.PrintInfo(OStream);
        OStream << std::endl;
        rThis.PrintData(OStream);
        return OStream;
    }

  
  
}  // namespace Kratos.

#endif // KRATOS_BICGSTAB_SOLVER_H_INCLUDED  defined 


