//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2004/02/26 12:18:15 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_ILU_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_ILU_PRECONDITIONER_H_INCLUDED



// System includes 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "../includes/kratos_string.h"


namespace Kratos
{

    template<class TSparseSpaceType, class TDenseSpaceType>
    class ILUPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
    {
    public:
     
      /// Counted pointer of ILUPreconditioner
      typedef boost::shared_ptr<ILUPreconditioner> Pointer;

      typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
   
      /// Default constructor.
      ILUPreconditioner()
      {
 
         L  = NULL;
         iL = NULL;
         jL = NULL;
         U  = NULL;
         iU = NULL;
         jU = NULL;
      }

      /// Copy constructor.
      ILUPreconditioner(const ILUPreconditioner& Other){}

      /// Destructor.
      virtual ~ILUPreconditioner()
      {	 
	        if ( L!=NULL) delete[]  L;
	        if (iL!=NULL) delete[] iL;
	        if (jL!=NULL) delete[] jL;
	        if ( U!=NULL) delete[]  U;
	        if (iU!=NULL) delete[] iU;
	        if (jU!=NULL) delete[] jU;
      }
      

    
      
      /// Assignment operator.
      ILUPreconditioner& operator=(const ILUPreconditioner& Other)
	{
	  mILUSize = Other.mILUSize;
	  unsigned int size = Other.iL[mILUSize];
	  L = new double[size];
	  U = new double[size];
	  iL = new int[mILUSize+1];
	  jL = new int[size];
	  iU = new int[mILUSize+1];
	  jU = new int[size];

	  std::copy(Other.L, Other.L+size, L); 
	  std::copy(Other.U, Other.U+size, U); 
	  std::copy(Other.iL, Other.iL+mILUSize+1, iL); 
	  std::copy(Other.jL, Other.jL+size, jL); 
	  std::copy(Other.iU, Other.iU+mILUSize+1, iU); 
	  std::copy(Other.jU, Other.jU+size, jU); 

 	  return *this;
	}

      
    
      

      virtual void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
      {
	VectorType z = rX;
	TSparseSpaceType::Mult(rA,z, rY);
	ApplyLeft(rY);
      }
      
      virtual void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
      {
	VectorType z = rX;
	ApplyTransposeLeft(z);
	TSparseSpaceType::TransposeMult(rA,z, rY);
      }
      
      /** multiply first rX by L^-1 and store result in temp
	  then multiply temp by U^-1 and store result in rX 
	  @param rX  Unknows of preconditioner suystem
      */
      virtual VectorType& ApplyLeft(VectorType& rX)
      {
	const int size = TSparseSpaceType::Size(rX);
	VectorType temp(size);
	double sum;
	int i, indexj;
	for (i=0; i<size; i++) {
		sum=rX[i];
		for (indexj=iL[i]; indexj<iL[i+1]; indexj++) {
			sum=sum-L[indexj]*temp[jL[indexj]];
		}
		temp[i]=sum;
	}
	for (i=size-1; i>=0; i--) {
		sum=temp[i];
		for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++) {
			sum=sum-U[indexj]*rX[jU[indexj]];
		}
		rX[i]=sum/U[iU[i]];
	}
	return rX;
      }
      
      /** Multiply first rX by U^-T and store result in temp
	  then multiply temp by L^-T and store result in rX
	  @param rX  Unknows of preconditioner suystem
      */    
      virtual VectorType& ApplyTransposeLeft(VectorType& rX)
      {
	const int size = TSparseSpaceType::Size(rX);
	VectorType temp(size);
	int i, indexj;
	double tempi, rxi;
	for (i=0; i<size; i++) temp[i]=rX[i];
	for (i=0; i<size; i++) {
		temp[i]=temp[i]/U[iU[i]];
		tempi=temp[i];
		for (indexj=iU[i]+1; indexj<iU[i+1]; indexj++) {
			temp[jU[indexj]]=temp[jU[indexj]]-tempi*U[indexj];
		}
	}
	for (i=0; i<size; i++) rX[i]=temp[i];
	for (i=size-1; i>=0; i--) {
		rxi=rX[i];
		for (indexj=iL[i]; indexj<iL[i+1]; indexj++) {
			rX[jL[indexj]]=rX[jL[indexj]]-rxi*L[indexj];
		}
	}
	return rX;
      }
      

      
      /// Return information about this object.
      virtual String Info() const
      {
	     return "ILUPreconditioner";
      }

      /// Print information about this object.
      virtual void  PrintInfo(std::ostream& OStream) const
      {
	     OStream << "ILUPreconditioner";
      }

      virtual void PrintData(std::ostream& OStream) const
      {
      }

            
      
    protected:
        
      unsigned int mILUSize;
      int *iL, *jL, *iU, *jU;
      double *L, *U;
                
    }; // Class ILUPreconditioner 
  
       
   /// input stream function
   template<class TSparseSpaceType, class TDenseSpaceType>
   inline std::istream& operator >> (std::istream& IStream, 
 				    ILUPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
   }

   /// output stream function
   template<class TSparseSpaceType, class TDenseSpaceType>
   inline std::ostream& operator << (std::ostream& OStream, 
 				    const ILUPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
   }
 
  
  
} 

#endif // KRATOS_ILU_PRECONDITIONER_H_INCLUDED  defined 

