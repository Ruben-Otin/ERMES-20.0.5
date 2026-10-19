

#if !defined(KRATOS_UPPER_TRIANGULAR_SOLVER_H_INCLUDED )
#define  KRATOS_UPPER_TRIANGULAR_SOLVER_H_INCLUDED



#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"

#include "../includes/define.h"
#include "direct_solver.h"


namespace Kratos
{

  /// A very simple solver for a system with upper triangular matrix.
  /** This class solve a upper triangular system of equation by  
      backsubstitution method. upper triangular solver is a template 
      class with this parameter: 
      - TSparseSpaceType which specify type of the unknowns,
        coefficients, sparse matrix, vector of unknowns, right hand side
        vector and their respective operators.  
      - TDenseMatrixType which
        specify type of the matrices used as temporary matrices or multi
	solve unknowns and right hand sides and their operators.  
      -	TReordererType which specify type of the Orderer that performs
	the reordering of matrix to optimize the solution.
  */
    template<class TSparseSpaceType, class TDenseSpaceType, 
    class TReordererType>
    class UpperTriangularSolver : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
    {
    public:
     
      
      /// Counted pointer of UpperTriangularSolver
      typedef boost::shared_ptr<UpperTriangularSolver> Pointer;

      typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
 
      UpperTriangularSolver(){}

      virtual ~UpperTriangularSolver(){}
      
      
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	      if(IsNotConsistent(rA, rX, rB))
	        return false;

	      const int size = TSparseSpaceType::Size(rX);

	      TSparseSpaceType::Set(rX, 0.00);

	      double a;

	      for(int i = size - 1 ; i >= 0 ; i--)
	        if((a = rA(i,i)) == 0.00)
		    KRATOS_ERROR(std::runtime_error, "UpperTriangularSolver::Solve" , "Zero founded in diagonal, Unstable system", "")
		      else
	          rX[i] = (rB[i] - TSparseSpaceType::RowDot(i, rA, rX)) / a;

	      return true;
      }
      
      /** Multi solve method for solving a set of linear systems with same coefficient matrix.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
      {
	      const int size1 = TDenseSpaceType::Size1(rX);
	      const int size2 = TDenseSpaceType::Size2(rX);

 	      bool is_solved = true;

	      VectorType x(size1);
	      VectorType b(size1);

	      for(int i = 0 ; i < size2 ; i++)
	      {
	          TDenseSpaceType::GetColumn(i,rX, x);
	          TDenseSpaceType::GetColumn(i,rB, b);
	          
	          is_solved &= Solve(rA,x,b);

	          TDenseSpaceType::SetColumn(i,rX, x);
	          TDenseSpaceType::SetColumn(i,rB, b);
	      }

	      return is_solved;
      }

     
      
      virtual String Info() const
      {
	     return "Upper triangular solver";
      }

      /// Print information about this object.
      void  PrintInfo(std::ostream& rOStream) const
      {
	     rOStream << "Upper triangular solver";
      }

      /// Print object's data.
      void  PrintData(std::ostream& rOStream) const 
      {
      }
      
             
    private:
     
      
      /// Assignment operator.
      UpperTriangularSolver& operator=(const UpperTriangularSolver& Other);

      /// Copy constructor.
      UpperTriangularSolver(const UpperTriangularSolver& Other);

        
     
        
    }; // Class UpperTriangularSolver 

  

  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TReordererType>
  inline std::istream& operator >> (std::istream& rIStream, 
				      UpperTriangularSolver<TSparseSpaceType, TDenseSpaceType, 
				      TReordererType>& rThis)
    {
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const UpperTriangularSolver<TSparseSpaceType, TDenseSpaceType, 
				      TReordererType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

  
  
}  // namespace Kratos.

#endif // KRATOS_UPPER_TRIANGULAR_SOLVER_H_INCLUDED  defined 


