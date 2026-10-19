//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: santi $
//   Date:                $Date: 2003/09/25 14:38:51 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_ITERATIVE_SOLVER_H_INCLUDED )
#define  KRATOS_ITERATIVE_SOLVER_H_INCLUDED



// System includes 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "linear_solver.h"


namespace Kratos
{

  
  /// Base class for all the iterative solvers in Kratos.
  /** This class define the general interface for the iterative solvers in Kratos.
      iterative solver is a template class with this parameter:
      - TSparseSpaceType which specify type
        of the unknowns, coefficients, sparse matrix, vector of
	unknowns, right hand side vector and their respective operators.
      - TDenseMatrixType which specify type of the
        matrices used as temporary matrices or multi solve unknowns and
	right hand sides and their operators.  
      - TPreconditionerType  which specify type of the preconditioner to be used.
      - TStopCriteriaType for specifying type of the object which control the stop criteria for iteration loop.
      - TReordererType which specify type of the Orderer that performs the reordering of matrix to optimize the solution.
  */
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
    class IterativeSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of IterativeSolver
      typedef boost::shared_ptr<IterativeSolver> Pointer;
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
     
      
      /// Default constructor.
      IterativeSolver() : mpPreconditioner(new TPreconditionerType()), mTolerance(0), mResidualNorm(0), mMaxIterationsNumber(0), 
	mIterationsNumber(0) {}

      IterativeSolver(double NewTolerance) : mpPreconditioner(new TPreconditionerType()), 
	mTolerance(NewTolerance), mResidualNorm(0), mMaxIterationsNumber(0), mIterationsNumber(0) {}

      IterativeSolver(double NewTolerance, unsigned int NewMaxIterationsNumber) : mpPreconditioner(new TPreconditionerType()),
	mTolerance(NewTolerance), mResidualNorm(0), 
	mMaxIterationsNumber(NewMaxIterationsNumber), 
	mIterationsNumber(0) {}

      IterativeSolver(double NewTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : 
	mpPreconditioner(pNewPreconditioner),
	mTolerance(NewTolerance), mResidualNorm(0), 
	mMaxIterationsNumber(NewMaxIterationsNumber), 
	mIterationsNumber(0) {}

      /// Destructor.
      virtual ~IterativeSolver(){}
      


      
      virtual typename TPreconditionerType::Pointer GetPreconditioner(void)
      {
	      return mpPreconditioner;
      }

      virtual const typename TPreconditionerType::Pointer GetPreconditioner(void) const
      {
	      return mpPreconditioner;
      }

      virtual void SetPreconditioner(typename TPreconditionerType::Pointer pNewPreconditioner)
      {
	      mpPreconditioner = pNewPreconditioner;
      }
      


      virtual void SetMaxIterationsNumber(unsigned int NewMaxIterationsNumber)
      {
	      mMaxIterationsNumber = NewMaxIterationsNumber;
      }


      virtual unsigned int GetMaxIterationsNumber()
      {
	      return mMaxIterationsNumber;
      }


      virtual void SetIterationsNumber(unsigned int NewIterationNumber)
      {
	      mIterationsNumber = NewIterationNumber;
      }


      virtual unsigned int GetIterationsNumber()
      {
	     return mIterationsNumber;
      }


      virtual void SetTolerance(double NewTolerance)
      {
	     mTolerance = NewTolerance;
      }


      virtual double GetTolerance()
      {
	      return mTolerance;
      }


      virtual void SetResidualNorm(double NewResidualNorm)
      {
		  if(mIterationsNumber == 1) mFirstResidualNorm = NewResidualNorm;

	      mResidualNorm = NewResidualNorm;
      }

      virtual double GetResidualNorm()
      {
	      return mResidualNorm;
      }
      
     

      virtual bool IterationNeeded()
      {
         return ((mIterationsNumber < mMaxIterationsNumber) && (mResidualNorm > mTolerance * mBNorm)); 
      }

      virtual bool IsConverged()
      {
	      return (mResidualNorm <= mTolerance * mBNorm);	  	  
      }
      
     
      
      String Info() const
      {
	      String buffer;
	      buffer << "IterativeSolver" << GetPreconditioner()->Info();
	      return  buffer;
      }
      
      void PrintInfo(std::ostream& OStream) const
      {
	      OStream << "IterativeSolver";
      }
      
      void PrintData(std::ostream& OStream) const
      {
        OStream << std::endl;
       
        if (mResidualNorm > mTolerance * mBNorm)
        {
            OStream << "!!!!!!!!!!! NO CONVERGENCE !!!!!!!!!!!" << std::endl << std::endl;
        }
        else
        {
            OStream << "--------> CONVERGENCE ACHIEVED <-------- " << std::endl << std::endl;
        }

        OStream << "||b||          = " << mBNorm                 << std::endl;
        OStream << "||Ax-b||       = " << mResidualNorm          << std::endl;
        OStream << "||Ax-b||/||b|| = " << (mResidualNorm/mBNorm) << std::endl;

        OStream << std::endl;

        OStream << "Max. iterations: " << mMaxIterationsNumber   << std::endl;
        OStream << "Tolerance      : " << mTolerance             << std::endl;
        OStream << "Iterations     : " << mIterationsNumber      << std::endl;
      }
      
      
     
      
    protected:
     
        

      double mResidualNorm;

      double mFirstResidualNorm;

      unsigned int mIterationsNumber;

      double mBNorm;

      double mTolerance;

      unsigned int mMaxIterationsNumber;
        
     

      void PreconditionedMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
      {
	      GetPreconditioner()->Mult(rA, rX, rY);
      }

      void PreconditionedTransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
      {
	      GetPreconditioner()->TransposeMult(rA, rX, rY);
      }
        
        
     
      
    private:
      
        
      /// A counted pointer to the preconditioner object.
      typename TPreconditionerType::Pointer mpPreconditioner;
             
      /// Assignment operator.
      IterativeSolver& operator=(const IterativeSolver& Other);

      /// Copy constructor.
      IterativeSolver(const IterativeSolver& Other);
  
        
    }; // Class IterativeSolver 

   
        
 
    /// input stream function
  template<class TSparseSpaceType, 
           class TDenseSpaceType, 
           class TPreconditionerType, 
           class TReordererType>
    inline std::istream& operator >> (std::istream& IStream, 
				      IterativeSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {    
    }

    /// output stream function
  template<class TSparseSpaceType, 
           class TDenseSpaceType, 
           class TPreconditionerType, 
           class TReordererType>
   inline std::ostream& operator << (std::ostream& OStream, 
				      const IterativeSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }

    ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_ITERATIVE_SOLVER_H_INCLUDED  defined 


