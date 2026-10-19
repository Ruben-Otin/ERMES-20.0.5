//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2003/05/22 11:55:45 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_IDENTITY_REORDERER_H_INCLUDED )
#define  KRATOS_IDENTITY_REORDERER_H_INCLUDED



// System includes 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "reorderer.h"


namespace Kratos
{

  /// Identity reorderer object which leave the equations unchanged.
  /** Identity reorderer object won't reorder the equations.  There is
      two template parameter: 
      - TSparseSpaceType which specify type
        of the unknowns, coefficients, sparse matrix, vector of
	unknowns, right hand side vector and their respective operators.
      - TDenseMatrixType which specify type of the
        matrices used as temporary matrices or multi solve unknowns and
	right hand sides and their operators.  
  */
  template<class TSparseSpaceType, class TDenseSpaceType>
  class IdentityReorderer : public Reorderer<TSparseSpaceType, TDenseSpaceType>
  {
  public:
      
    /// Counted pointer of IdentityReorderer
    typedef boost::shared_ptr<IdentityReorderer> Pointer;

    typedef Reorderer<TSparseSpaceType, TDenseSpaceType> BaseType;
  
    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
    typedef typename TSparseSpaceType::VectorType VectorType;
  
    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
    /// Default constructor.
    IdentityReorderer(){}

    /// Destructor.
    virtual ~IdentityReorderer(){}
    
    void Reorder(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }
      
    void Reorder(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& OStream)
    {
      OStream << "IdentityReorderer" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& OStream){}
                      
  private:
       
    /// Assignment operator.
    IdentityReorderer& operator=(const IdentityReorderer& Other);

    /// Copy constructor.
    IdentityReorderer(const IdentityReorderer& Other);

        
    ///@}    
        
  }; // Class IdentityReorderer 

 
  /// input stream function
//   inline std::istream& operator >> (std::istream& IStream, 
// 				    IdentityReorderer& rThis);
  
  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const IdentityReorderer<TSparseSpaceType, TDenseSpaceType>& rThis)
    {
      rThis.PrintInfo(OStream);
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_IDENTITY_REORDERER_H_INCLUDED  defined 


