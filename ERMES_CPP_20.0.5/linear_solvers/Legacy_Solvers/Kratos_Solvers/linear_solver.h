

#if !defined(KRATOS_LINEAR_SOLVER_H_INCLUDED )
#define  KRATOS_LINEAR_SOLVER_H_INCLUDED


#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"
#include "../includes/define.h"



namespace Kratos
{
  
   // Base class for all the linear solvers in Kratos.
   /* This class define the general interface for the linear solvers in Kratos.
      There is three template parameter: 
      - TSparseSpaceType which specify type
        of the unknowns, coefficients, sparse matrix, vector of
	    unknowns, right hand side vector and their respective operators.
      - TDenseMatrixType which specify type of the
        matrices used as temporary matrices or multi solve unknowns and
	    right hand sides and their operators.  
      - TReordererType which specify type of the Orderer that performs the reordering of matrix to optimize the solution.
   */
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    class LinearSolver
    {
        public:
       
          typedef boost::shared_ptr<LinearSolver> Pointer;
  
          typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
          typedef typename TSparseSpaceType::VectorType VectorType;
  
          typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
      

          LinearSolver() : mpReorderer(new TReordererType()){}

          virtual ~LinearSolver(){}
          
          /** Normal solve method.
	      Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	      rVectorx is also th initial guess for iterative methods.
	      @param rA. System matrix
	      @param rX. Solution vector. it's also the initial 
	      guess for iterative linear solvers.
 	      @param rB. Right hand side vector.
          */
          virtual bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
          {
	          return false;
          }
      
      
          /** Multi solve method for solving a set of linear systems with same coefficient matrix.
	      Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	      rVectorx is also th initial guess for iterative methods.
	      @param rA. System matrix
	      @param rX. Solution vector. it's also the initial 
	      guess for iterative linear solvers.
 	      @param rB. Right hand side vector.
          */
          virtual bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
          {
	          return false;
          }
      
      
          virtual typename TReordererType::Pointer GetReorderer(void)
          {
	          return mpReorderer;
          }

          virtual void SetReorderer(typename TReordererType::Pointer pNewReorderer)
          {
	          mpReorderer = pNewReorderer;
          }
      

          virtual bool IsConsistent(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
          {
	          const int size = TSparseSpaceType::Size1(rA);

	          return ((size ==  TSparseSpaceType::Size2(rA)) && 
		              (size ==  TSparseSpaceType::Size(rX))  &&
		              (size ==  TSparseSpaceType::Size(rB))); 
          }
      
          virtual bool IsConsistent(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
          {
	          const int size = TSparseSpaceType::Size1(rA);

	          return ((size ==  TSparseSpaceType::Size2(rA)) && 
		              (size ==  TDenseSpaceType::Size1(rX)) &&
		              (size ==  TDenseSpaceType::Size1(rB)) && 
		              (TDenseSpaceType::Size2(rX) == TDenseSpaceType::Size2(rB)));
          }
      
      
          virtual bool IsNotConsistent(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
          {
	          return (!IsConsistent(rA, rX, rB));
          }

          virtual bool IsNotConsistent(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
          {
	          return (!IsConsistent(rA, rX, rB));
          }

   

          virtual String Info() const
          {
	         return "Linear solver";
          }

          virtual void PrintInfo(std::ostream& OStream) const
          {
	         OStream << "Linear solver";
          }
      
          virtual void PrintData(std::ostream& OStream) const
          {
          }
      
        
        private:
         
          typename TReordererType::Pointer mpReorderer;

          LinearSolver& operator=(const LinearSolver& Other);

          LinearSolver(const LinearSolver& Other);
     
    }; // Class LinearSolver 

  
        
 
    // input stream function
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::istream& operator >> (std::istream& IStream, 
				      LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
    {
    }

    // output stream function
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
				      const LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
    {
          rThis.PrintInfo(rOStream);
          rOStream << std::endl;
          rThis.PrintData(rOStream);

          return rOStream;
    }

}  // namespace Kratos.

#endif // KRATOS_LINEAR_SOLVER_H_INCLUDED  defined 


