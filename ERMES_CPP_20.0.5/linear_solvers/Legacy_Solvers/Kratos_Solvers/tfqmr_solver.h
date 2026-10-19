//*********************************************************************************   
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: Rubén $
//*   Date:                $Date: 2003/06/04 11:38:33 $
//*   Revision:            $Revision: 1.4 $
//*********************************************************************************

#if !defined(KRATOS_TFQMR_SOLVER_H_INCLUDED )
#define  KRATOS_TFQMR_SOLVER_H_INCLUDED


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"

// Project includes
#include "../includes/define.h"
#include "iterative_solver.h"

namespace Kratos
{

    template<class TSparseSpaceType, 
             class TDenseSpaceType,
             class TPreconditionerType,
             class TReordererType>
    class TFQMRSolver : public IterativeSolver<TSparseSpaceType, 
                                               TDenseSpaceType, 
                                               TPreconditionerType, 
                                               TReordererType>
    {
    public:
      
      /// Counted pointer of TFQMRSolver
      typedef boost::shared_ptr<TFQMRSolver> Pointer;

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;



      
      /// Default constructor.
      TFQMRSolver(){}

      TFQMRSolver(double NewTolerance) : BaseType(NewTolerance){}

      TFQMRSolver(double NewTolerance, 
                  unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber){}

      TFQMRSolver(double NewMaxTolerance, 
                  unsigned int NewMaxIterationsNumber, 
                  typename TPreconditionerType::Pointer pNewPreconditioner) : 
	              BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Destructor.
      virtual ~TFQMRSolver(){}
      
      
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector*/
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


      /// Return information about this object.
      String Info() const
      {
	      String buffer;
	      buffer << "Tranpose-free QMR linear solver with " << GetPreconditioner()->Info();
	      return  buffer;
      }
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
      {
	      OStream << "Tranpose-free QMR linear solver with ";
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

            mBNorm        = TSparseSpaceType::TwoNorm(rB);

            double errtol = mBNorm*mTolerance;

            VectorType r(rB);
            VectorType w(r);

            VectorType y1(r);
            VectorType y2(size);
            for(VectorType::iterator ity2 = y2.begin(); ity2 != y2.end(); ity2++) *ity2 = 0.00;

            VectorType d(size);
            for(VectorType::iterator itd = d.begin(); itd != d.end(); itd++) *itd = 0.00;

            VectorType v(size);
            PreconditionedMult(rA, y1, v);

            VectorType u1(v);
            VectorType u2(size);
            for(VectorType::iterator itu2 = u2.begin(); itu2 != u2.end(); itu2++) *itu2 = 0.00;

            double theta = 0.00; 
            double eta   = 0.00;
            double sigma = 0.00;
            double alpha = 0.00;
            double c     = 0.00;

            double tau   = TSparseSpaceType::TwoNorm(r); 
            double rho   = tau*tau;
            double rhon  = 0.00;
            double beta  = 0.00;

            int m = 0;


            while(mIterationsNumber < mMaxIterationsNumber)
            { 

                mIterationsNumber++;

                sigma = TSparseSpaceType::Dot(r,v);
            
                if (sigma == 0.00) break;
             
                alpha = rho/sigma;

           
                ////////////////even iteration///////////////////

                m = 2*mIterationsNumber-1;  

                //w=w-alpha*u1;
                TSparseSpaceType::ScaleAndAdd(-alpha, u1, 1.00, w);  

                //d=y1+(theta*theta*eta/alpha)*d; 
                TSparseSpaceType::ScaleAndAdd(1.00, y1, (theta*theta*eta/alpha), d);  
                
                theta = TSparseSpaceType::TwoNorm(w)/tau; 
                c     = 1.00/sqrt(1+theta*theta);
                tau   = tau*theta*c;
                eta   = c*c*alpha;

                //x=x+eta*d;
                TSparseSpaceType::ScaleAndAdd(eta, d, 1.00, rX); 
                 
                mResidualNorm = tau*sqrt(m+1);
                 
                if (mResidualNorm <= errtol) break;


                ////////////////odd iteration///////////////////

                m = 2*mIterationsNumber;

                //y2=y1-alpha*v;
                TSparseSpaceType::ScaleAndAdd(1.00, y1, -alpha, v, y2);  

                //u2=A*y2;
                PreconditionedMult(rA,y2,u2);    
                 
                //w=w-alpha*u2;
                TSparseSpaceType::ScaleAndAdd(-alpha, u2, 1.00, w);  
         
                //d=y2+(theta*theta*eta/alpha)*d; 
                TSparseSpaceType::ScaleAndAdd(1.00, y2, (theta*theta*eta/alpha), d);  

                theta = TSparseSpaceType::TwoNorm(w)/tau; 
                c     = 1.00/sqrt(1+theta*theta);
                tau   = tau*theta*c;
                eta   = c*c*alpha;

                //x=x+eta*d;
                TSparseSpaceType::ScaleAndAdd(eta, d, 1.00, rX); 

                mResidualNorm = tau*sqrt(m+1);
                 
                if (mResidualNorm <= errtol) break;
       

                ///////////////////////////////////////////////


                if (rho == 0.00) break;
            
                rhon = TSparseSpaceType::Dot(r,w);
                beta = rhon/rho; 
                rho  = rhon;

                //y1= w + beta*y2;
                TSparseSpaceType::ScaleAndAdd(1.00, w, beta, y2, y1);  

                //u1 = A*y1;
                PreconditionedMult(rA,y1,u1); 
                
                //v=u1+beta*(u2+beta*v);
                TSparseSpaceType::ScaleAndAdd(1.00, u2, beta, v);  
                TSparseSpaceType::ScaleAndAdd(1.00, u1, beta, v);  

                 //Print iteration info
                info_count++;
                if (info_count==100)
                {
                   std::cout<<"it = "<<mIterationsNumber<<"   res = "<<(mResidualNorm/mBNorm)<<std::endl;               
                   info_count=0;
                }

            }

            std::cout<<std::endl;
            return IsConverged();
        }
      

        /// Assignment operator.
        TFQMRSolver& operator=(const TFQMRSolver& Other);

        /// Copy constructor.
        TFQMRSolver(const TFQMRSolver& Other);
        
    }; // ********** Class TFQMRSolver ***************** 
 
        
 
  /// input stream function
   template<class TSparseSpaceType, 
            class TDenseSpaceType, 
            class TPreconditionerType, 
            class TReordererType>
    inline std::istream& operator >> (std::istream& IStream, TFQMRSolver<TSparseSpaceType, 
                                                                         TDenseSpaceType, 
				                                                         TPreconditionerType, 
                                                                         TReordererType>& rThis)
    {
    }

   /// output stream function
    template<class TSparseSpaceType, 
             class TDenseSpaceType, 
             class TPreconditionerType, 
             class TReordererType>
     inline std::ostream& operator << (std::ostream& OStream,const TFQMRSolver<TSparseSpaceType, 
                                                                               TDenseSpaceType, 
				                                                               TPreconditionerType, 
                                                                               TReordererType>& rThis)
    {
        rThis.PrintInfo(OStream);
        OStream << std::endl;
        rThis.PrintData(OStream);
        return OStream;
    } 
    
}  // namespace Kratos.

#endif // KRATOS_TFQMR_SOLVER_H_INCLUDED  defined 