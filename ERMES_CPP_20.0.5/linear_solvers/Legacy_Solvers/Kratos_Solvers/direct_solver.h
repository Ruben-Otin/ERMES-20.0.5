//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2003/09/18 14:08:42 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_DIRECT_SOLVER_H_INCLUDED )
#define  KRATOS_DIRECT_SOLVER_H_INCLUDED



// System includes 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "linear_solver.h"


namespace Kratos
{

  
  
  // Base class for all direct solvers in Kratos.
  /* This class define the general interface for direct solvers in Kratos.
     direct solver is a template class with this parameter:
     - TSparseSpaceType which specify type
       of the unknowns, coefficients, sparse matrix, vector of
	   unknowns, right hand side vector and their respective operators.
     - TDenseMatrixType which specify type of the
       matrices used as temporary matrices or multi solve unknowns and
	   right hand sides and their operators.  
     - TReordererType which specify type of the Orderer that performs the reordering of matrix to optimize the solution.
  */
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    class DirectSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
    {
      public:
  
          /// Counted pointer of DirectSolver
          typedef boost::shared_ptr<DirectSolver> Pointer;

          typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
  
          typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
          typedef typename TSparseSpaceType::VectorType VectorType;
  
          typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      
          /// Default constructor.
          DirectSolver(){}

          /// Destructor.
          virtual ~DirectSolver(){}
      
          /// Print information about this object.
          void  PrintInfo(std::ostream& rOStream) const
          {
	          rOStream << "Direct solver";
          }

          /// Print object's data.
          void  PrintData(std::ostream& rOStream) const 
          {
          }
      
    private:
     
          /// Assignment operator.
          DirectSolver& operator=(const DirectSolver& Other);

          /// Copy constructor.
          DirectSolver(const DirectSolver& Other);
        
    }; // Class DirectSolver 

 
        
     /// input stream function
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, 
				                      DirectSolver<TSparseSpaceType, TDenseSpaceType,TReordererType>& rThis)
    {
    }

    /// output stream function
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream,
                                      const DirectSolver<TSparseSpaceType,TDenseSpaceType,TReordererType>& rThis)
    {
         rThis.PrintInfo(rOStream);
         rOStream << std::endl;
         rThis.PrintData(rOStream);

         return rOStream;
    }
   
}  // namespace Kratos.

#endif // KRATOS_DIRECT_SOLVER_H_INCLUDED  defined 


