//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2003/06/04 11:38:33 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_IDENTITY_PRECONDITIONER__H_INCLUDED )
#define  KRATOS_IDENTITY_PRECONDITIONER_H_INCLUDED



// System includes 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "preconditioner.h"


namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
  
  ///@} 
  ///@name  Enum's
  ///@{
      
  ///@}
  ///@name  Functions 
  ///@{
      
  ///@}
  ///@name Kratos Classes
  ///@{
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<class TSparseSpaceType, class TDenseSpaceType>
  class IdentityPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of IdentityPreconditioner
      typedef boost::shared_ptr<IdentityPreconditioner> Pointer;

      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      IdentityPreconditioner(){}

      /// Destructor.
      virtual ~IdentityPreconditioner(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) {}
      
      void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) {}
      
      void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) {TSparseSpaceType::Mult(rA, rX, rY);}

      void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) {TSparseSpaceType::TransposeMult(rA, rX, rY);}

      VectorType& ApplyLeft(VectorType& rX){return rX;}
      
      VectorType& ApplyRight(VectorType& rX){return rX;}
      
      VectorType& ApplyTransposeLeft(VectorType& rX){return rX;}
      
      VectorType& ApplyTransposeRight(VectorType& rX){return rX;}
      
      VectorType& ApplyInverseRight(VectorType& rX){return rX;}
      
      VectorType& Finalize(VectorType& rX){return rX;}
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{
      
      /// Return information about this object.
      String Info() const
	{
	  return  "Identity preconditioner";
	}
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
	{
	  OStream << "Identity preconditioner";
	}

            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
    protected:
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operations
      ///@{ 
        
        
      ///@} 
      ///@name Protected  Access 
      ///@{ 
        
        
      ///@}      
      ///@name Protected Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Protected LifeCycle 
      ///@{ 
      
            
      ///@}
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      IdentityPreconditioner& operator=(const IdentityPreconditioner& Other);

      /// Copy constructor.
      IdentityPreconditioner(const IdentityPreconditioner& Other);

      
      ///@}    
        
    }; // Class IdentityPreconditioner 

    ///@} 
  
    ///@name Type Definitions       
    ///@{ 
  
  
    ///@} 
    ///@name Input and output 
    ///@{ 
        
 
    /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
    inline std::istream& operator >> (std::istream& IStream, 
 				      IdentityPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
    {
    }

    /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
    inline std::ostream& operator << (std::ostream& OStream, 
 				      const IdentityPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
   }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_IDENTITY_PRECONDITIONER_H_INCLUDED  defined 


