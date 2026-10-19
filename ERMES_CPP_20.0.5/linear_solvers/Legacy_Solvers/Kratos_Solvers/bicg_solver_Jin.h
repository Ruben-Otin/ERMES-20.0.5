//*********************************************************************   
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: ruben $
//*   Date:                $Date: 2003/06/04 11:38:32 $
//*   Revision:            $Revision: 1.4 $
//*********************************************************************

#if !defined(KRATOS_BICG_SOLVER_JIN_H_INCLUDED )
#define  KRATOS_BICG_SOLVER_JIN_H_INCLUDED



// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "iterative_solver.h"


namespace Kratos
{


   template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType, class TReordererType>
   class BiCGSolverJin : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
   {
    public:
     
      //Counted pointer of BiCGSolverJin
      typedef boost::shared_ptr<BiCGSolverJin> Pointer;

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      
      //Default constructor.
      BiCGSolverJin(){}

      BiCGSolverJin(double NewTolerance) : BaseType(NewTolerance){}

      BiCGSolverJin(double NewTolerance, 
                    unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber){}

      BiCGSolverJin(double NewMaxTolerance, 
                    unsigned int NewMaxIterationsNumber, 
                    typename TPreconditionerType::Pointer pNewPreconditioner) : 
                	BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      //Destructor.
      virtual ~BiCGSolverJin(){}
      
     
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
 	  @param rB. Right hand side vector*/
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

     
      //Return information about this object.
      String Info() const
      {
	      String buffer;
	      buffer << "Biconjugate gradient (from Jin) linear solver with " << GetPreconditioner()->Info();
	      return  buffer;
      }
      
      //Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
      {
	      OStream << "Biconjugate gradient (from Jin) linear solver with ";
	      GetPreconditioner()->PrintInfo(OStream);
      }

      //Print object's data.
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

            double alpha, gamma, rnorm;

            mBNorm = TSparseSpaceType::TwoNorm(rB);

            //v_aux 
            VectorType v_aux(size);
           
            //r = rB;
            VectorType r(rB);

            //p = r;
            VectorType p(r);

            //rnorm = (r,r)
            rnorm = TSparseSpaceType::Dot(r,r);

            do
            {
                //v_aux = rA*p;
                //TSparseSpaceType::Mult(rA,p,v_aux);
                PreconditionedMult(rA,p,v_aux);

                //alpha = rnorm / (v_aux,p)
                alpha = rnorm / TSparseSpaceType::Dot(v_aux, p);

                //rX = alpha*p + rX;
                TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);

                 //r = (-alpha)*v_aux + r;
                TSparseSpaceType::ScaleAndAdd(-alpha, v_aux, 1.00, r);

                //gamma = (r,r)/rnorm
                gamma = TSparseSpaceType::Dot(r,r) / rnorm;

                //p = r + gamma*p;
                TSparseSpaceType::ScaleAndAdd(1.00, r, gamma, p);

                //rnorm = gamma*rnorm
                rnorm = gamma*rnorm;
                
                //Cheking residual norm
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



        // Assignment operator.
        BiCGSolverJin& operator=(const BiCGSolverJin& Other);

        // Copy constructor.
        BiCGSolverJin(const BiCGSolverJin& Other);
 
        
    }; //********** Class BiCGSolverJin *************


 
   /// input stream function
   template<class TSparseSpaceType, 
            class TDenseSpaceType, 
            class TPreconditionerType, 
            class TReordererType>
    inline std::istream& operator >> (std::istream& IStream, 
				      BiCGSolverJin<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
    }

    /// output stream function
    template<class TSparseSpaceType, 
             class TDenseSpaceType, 
             class TPreconditionerType, 
             class TReordererType>
    inline std::ostream& operator << (std::ostream& OStream, 
				    const BiCGSolverJin<TSparseSpaceType, TDenseSpaceType, 
				    TPreconditionerType, TReordererType>& rThis)
    {
        rThis.PrintInfo(OStream);
        OStream << std::endl;
        rThis.PrintData(OStream);
        return OStream;
    }
  
  
}  // namespace Kratos.

#endif // KRATOS_BICG_SOLVER_JIN_H_INCLUDED  defined 

