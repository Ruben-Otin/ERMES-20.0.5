
#if !defined(KRATOS_DIAGONAL_SOLVER_H_INCLUDED )
#define  KRATOS_DIAGONAL_SOLVER_H_INCLUDED



// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "direct_solver.h"


namespace Kratos
{
  /// A very simple solver for a system with diagonal matrix.
  /** This class solve a diagonal system of equation just by dividing
      the right handside by diagonal member of matrix.  diagonal
      solver is a template class with this parameter: 
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
    class DiagonalSolver : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
    {
    public:
   
      typedef boost::shared_ptr<DiagonalSolver> Pointer;

      typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      /// Default constructor.
      DiagonalSolver(){}

      /// Destructor.
      virtual ~DiagonalSolver(){}
      
      
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

	  double a;

	  for(int i = 0 ; i < size ; i++)
	    if((a = rA(i,i)) == 0.00)
		KRATOS_ERROR(std::runtime_error, "DiagonalSolver::Solve" , "Zero founded in diagonal, Unstable system", "")
	    else
	      rX[i] = rB[i] / a;

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

	  double a;

	  for(int i = 0 ; i < size1 ; i++)
	    if((a = rA(i,i)) == 0.00)
		KRATOS_ERROR(std::runtime_error, "DiagonalSolver::Solve" , "Zero founded in diagonal, Unstable system", "")
	    else
	      for(int j = 0 ; j < size2 ; j++)
		rX(i, j) = rB(i, j) / a;

	  return true;
	}
    
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "Diagonal solver";
	}

      /// Print object's data.
      void  PrintData(std::ostream& rOStream) const 
	{
	}
  
    private:
     
      /// Assignment operator.
      DiagonalSolver& operator=(const DiagonalSolver& Other);

      /// Copy constructor.
      DiagonalSolver(const DiagonalSolver& Other);
 
        
    }; // Class DiagonalSolver 

 
  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TReordererType>
  inline std::istream& operator >> (std::istream& rIStream, 
				      DiagonalSolver<TSparseSpaceType, TDenseSpaceType, 
				      TReordererType>& rThis)
    {
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const DiagonalSolver<TSparseSpaceType, TDenseSpaceType, 
				      TReordererType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

  
}  // namespace Kratos.

#endif // KRATOS_DIAGONAL_SOLVER_H_INCLUDED  defined 


