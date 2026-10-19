//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2003/05/14 12:24:31 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_REORDERER_H_INCLUDED )
#define  KRATOS_REORDERER_H_INCLUDED



// System includes 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"


namespace Kratos
{

  /// Base class for all reorderer objects in Kratos used in linear solvers.
  /** Reorderer is the base class for all reorderer objects in Kratos used in linear solvers.
      two template parameter: 
      - TSparseSpaceType which specify type
        of the unknowns, coefficients, sparse matrix, vector of
	unknowns, right hand side vector and their respective operators.
      - TDenseMatrixType which specify type of the
        matrices used as temporary matrices or multi solve unknowns and
	right hand sides and their operators.  
  */
  template<class TSparseSpaceType, class TDenseSpaceType>
  class Reorderer
  {
    public:
     
      
      /// Counted pointer of Reorderer
      typedef boost::shared_ptr<Reorderer> Pointer;

      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      /// Default constructor.
      Reorderer(){}

      /// Destructor.
      virtual ~Reorderer(){}    


      virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB){}

      virtual void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB){}


      virtual void Reorder(SparseMatrixType& rA, VectorType& rX, VectorType& rB){}
      
      virtual void Reorder(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB){}


      virtual void InverseReorder(SparseMatrixType& rA, VectorType& rX, VectorType& rB){}

      virtual void InverseReorder(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB){}


      virtual void PrintInfo(std::ostream& OStream){}

      /// Print object's data.
      virtual void PrintData(std::ostream& OStream){}
     
    private:    
      
      /// Assignment operator.
      Reorderer& operator=(const Reorderer& Other);

      /// Copy constructor.
      Reorderer(const Reorderer& Other);

        
    }; // Class Reorderer 

 
 
  /// input stream function
//   inline std::istream& operator >> (std::istream& IStream, 
// 				    Reorderer& rThis);
  
  /// output stream function
//   inline std::ostream& operator << (std::ostream& OStream, 
// 				    const Reorderer& rThis);
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_REORDERER_H_INCLUDED  defined 


