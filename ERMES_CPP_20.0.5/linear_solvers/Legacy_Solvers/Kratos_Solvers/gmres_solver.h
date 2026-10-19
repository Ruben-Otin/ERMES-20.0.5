

#if !defined(KRATOS_GMRES_SOLVER_H_INCLUDED )
#define  KRATOS_GMRES_SOLVER_H_INCLUDED



// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"

// Project includes
#include "../includes/define.h"
#include "iterative_solver.h"


namespace Kratos
{


  
    template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
    class GMRESSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
    public:
    
      
      /// Counted pointer of GMRESSolver
      typedef boost::shared_ptr<GMRESSolver> Pointer;

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
   
      GMRESSolver(){}

      GMRESSolver(double NewTolerance) : BaseType(NewTolerance){}

      GMRESSolver(double NewTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber){}

      GMRESSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : 
	BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Destructor.
      virtual ~GMRESSolver(){}
      

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
	  buffer << "linear solver with " << GetPreconditioner()->Info();
	  return  buffer;
	}
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
	{
	  OStream << "Generalized minimum resdual linear solver with ";
	  GetPreconditioner()->PrintInfo(OStream);
	}

      /// Print object's data.
      void  PrintData(std::ostream& OStream) const 
	{
	  BaseType::PrintData(OStream);
	}
      
            
  
    private:
      

      unsigned int mRestartValue;
        
        
      
      bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	      const int size = TSparseSpaceType::Size(rX);
	      

	      mRestartValue = 25;

	      mIterationsNumber = 0;
          int info_count    = 0;

	      mBNorm = TSparseSpaceType::TwoNorm(rB);
    
	      int i,j,k,l;
    
	      VectorType r(size,double());
	      VectorType x0(rX);
	      VectorType w(size);
    
	      DenseMatrixType h(mRestartValue + 1,mRestartValue);
	      DenseMatrixType v(mRestartValue + 1,size);
	      VectorType g(mRestartValue + 1,double());
	      VectorType y(mRestartValue,double());
	      VectorType s(mRestartValue,double());
	      VectorType c(mRestartValue,double());
    
	      double beta;
    
	      do
	        {
	          PreconditionedMult(rA,rX,r);
	          TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);
	          
	          beta = TSparseSpaceType::TwoNorm(r);

	          if(beta == 0.00)
		    KRATOS_ERROR(std::runtime_error, "GMRESSolver::IterativeSolve" , "Unstable system", "")
    
		      g[0] = beta;
      
	          for(i = 0; i < size ; i++)
		    v(0,i) = r[i] / beta;
                    
	          h = double();
    
	          for(j = 0 ; j < mRestartValue ; j++)
		    {
		      PreconditionedMult(rA,v[j],w);
		      for(i = 0 ; i <= j ; i++)
		        {
		          h(i,j) = TSparseSpaceType::Dot(w,v[i]);
		          TSparseSpaceType::ScaleAndAdd(-h(i,j),v[i], 1.00, w);
		        }
		      h(j+1,j) = TSparseSpaceType::TwoNorm(w);
		      if(h(j+1,j) == 0.00)
		        {
		          j++;
		          break;
		        }
          
		      for(i = 0; i < size ; i++)
		        v(j+1,i) = w[i] / h(j+1,j);
              
		      for(i = 0 ; i < j ; i++)
		        {
		          double raux = c[i] * h(i,j) - s[i] * h(i+1, j); 
		          h(i+1, j) = s[i] * h(i,j) + c[i] * h(i+1, j);
		          h(i,j) = raux;
		        }
        
		      double square_root = sqrt(h(j,j)*h(j,j) + h(j+1,j)*h(j+1,j));
		      c[j] = h(j,j) / square_root;
		      s[j] = -h(j+1,j) / square_root;
        
		      h(j,j) = square_root;
		      h(j+1,j) = double();
      
          
		      g[j+1] = s[j] *g[j];
		      g[j] = c[j] * g[j];
        
		      SetResidualNorm(fabs(g[j+1]));
		    }
      
	          j--;

	          for(l = j ; l >= 0 ; l--)
		    {
		      y[l] = g[l];
		      for(k = l + 1 ; k <= j ; k++)
		        y[l] -= h(l,k) * y[k];
          
		      y[l] = y[l] / h(l,l); 
		    }
    
	          double temp;
      
	          for(l = 0 ; l < size ; l++)
		    {
		      temp = rX[l];
		      for(k = 0 ; k <= j ; k++)
		        temp += v(k,l)*y[k];
          
		      rX[l] = temp;
		    }
     
	          mIterationsNumber++;


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
      GMRESSolver& operator=(const GMRESSolver& Other);

      /// Copy constructor.
      GMRESSolver(const GMRESSolver& Other);

        
     
        
    }; // Class GMRESSolver 

  
 
  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::istream& operator >> (std::istream& IStream, 
				      GMRESSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const GMRESSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }

  
}  // namespace Kratos.

#endif // KRATOS_GMRES_SOLVER_H_INCLUDED defined 


