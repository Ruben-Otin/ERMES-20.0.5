//********************************************************************************   
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: Rubén $
//*   Date:                $Date: 2003/06/04 11:38:33 $
//*   Revision:            $Revision: 1.4 $
//********************************************************************************


#if !defined(KRATOS_CG_SOLVER_JIN_H_INCLUDED )
#define  KRATOS_CG_SOLVER_JIN_H_INCLUDED


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "iterative_solver.h"


namespace Kratos
{

    template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType, class TReordererType>
    class CGSolverJin : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
      public:
      
      /// Counted pointer of CGSolverJin
      typedef boost::shared_ptr<CGSolverJin> Pointer;

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
      
      /// Default constructor.
      CGSolverJin(){}

      CGSolverJin(double NewTolerance) : BaseType(NewTolerance){}

      CGSolverJin(double NewTolerance, 
                  unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber){}

      CGSolverJin(double NewMaxTolerance, 
                  unsigned int NewMaxIterationsNumber, 
                  typename TPreconditionerType::Pointer pNewPreconditioner) :  BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Destructor.
      virtual ~CGSolverJin(){}
      
     
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.*/
      bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	      if(IsNotConsistent(rA, rX, rB)) return false;

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
 	  @param rB. Right hand side vector.*/
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
	      buffer << "Conjugate gradient (from Jin) linear solver with " << GetPreconditioner()->Info();
	      return  buffer;
      }
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
      {
	      OStream << "Conjugate gradient (from Jin) linear solver with ";
	      GetPreconditioner()->PrintInfo(OStream);
      }

      /// Print object's data.
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

            double alpha;

            mBNorm = TSparseSpaceType::TwoNorm(rB);

            //r = rB;
            VectorType r(rB);
          
            //v_aux = rA*r
            VectorType v_aux(size);
            //TSparseSpaceType::Mult(rA,r,v_aux);
            PreconditionedMult(rA,r,v_aux);


            //beta = 1/norm(v_aux)^2;
            double beta = 1.00 / TSparseSpaceType::Dot(v_aux,v_aux);

            //p = beta*v_aux
            VectorType p(v_aux);
            TSparseSpaceType::Scale(beta,p);

            do
            {
                //v_aux = rA*p;
                //TSparseSpaceType::Mult(rA,p,v_aux);
                PreconditionedMult(rA,p,v_aux);

                //alpha = 1/norm(v_aux)^2;
                alpha = 1.00 / TSparseSpaceType::Dot(v_aux,v_aux);

                //rX = rX + alpha*p;
                TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);

                //r = r - alpha*v_aux;
                TSparseSpaceType::ScaleAndAdd(-alpha, v_aux, 1.00, r);

                //v_aux = rA*r
                //TSparseSpaceType::Mult(rA,r,v_aux);
                PreconditionedMult(rA,r,v_aux);

                //beta = 1/norm(v_aux)^2;
                beta = 1.00 / TSparseSpaceType::Dot(v_aux,v_aux);

                //p = p + beta*v_aux;
                TSparseSpaceType::ScaleAndAdd(beta, v_aux, 1.00, p);

                mResidualNorm =TSparseSpaceType::TwoNorm(r);

                mIterationsNumber++;

                //Print iteration info
                info_count++;
                if (info_count==100)
                {
                   std::cout<<"it = "<<mIterationsNumber<<"   res = "<<(mResidualNorm/mBNorm)<<std::endl;
                   info_count=0;
                }

            } while(IterationNeeded());   
                       
            std::cout<<std::endl;

            return IsConverged();
        }




        /// Assignment operator.
        CGSolverJin& operator=(const CGSolverJin& Other);

        /// Copy constructor.
        CGSolverJin(const CGSolverJin& Other);

                 
    }; //********** End Class CGSolverJin ***************** 

 
        
 
   /// input stream function
   template<class TSparseSpaceType, 
            class TDenseSpaceType, 
            class TPreconditionerType, 
            class TReordererType>
   inline std::istream& operator >> (std::istream& IStream, 
				      CGSolverJin<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
    }

    /// output stream function
    template<class TSparseSpaceType, 
             class TDenseSpaceType, 
             class TPreconditionerType, 
             class TReordererType>
    inline std::ostream& operator << (std::ostream& OStream, 
				      const CGSolverJin<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
        rThis.PrintInfo(OStream);
        OStream << std::endl;
        rThis.PrintData(OStream);

        return OStream;
    }

    
}  // namespace Kratos.

#endif // KRATOS_CG_SOLVER_JIN_H_INCLUDED  defined 
