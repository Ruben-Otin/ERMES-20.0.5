//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2004/02/26 12:18:15 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_MILU0_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_MILU0_PRECONDITIONER_H_INCLUDED



// System includes
#include <algorithm> 


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "../includes/kratos_string.h"


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
  
  ///@name  Preconditioners 
  ///@{ 
  
  /// MILU0Preconditioner class. 
  /**   */
  template<class TSparseSpaceType, class TDenseSpaceType>
    class MILU0Preconditioner : public ILUPreconditioner<TSparseSpaceType, TDenseSpaceType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of MILU0Preconditioner
      typedef boost::shared_ptr<MILU0Preconditioner> Pointer;

      typedef ILUPreconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      MILU0Preconditioner(){}

      /// Copy constructor.
      MILU0Preconditioner(const MILU0Preconditioner& Other){}

      /// Destructor.
      virtual ~MILU0Preconditioner()
      {	
      }
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      MILU0Preconditioner& operator=(const MILU0Preconditioner& Other)
	{
	  BaseType::operator=(Other);
 	  return *this;
	}


      void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) 
      {
          BaseType::Initialize(rA, rX, rB);
      }

      
      ///@}
      ///@name Operations
      ///@{
      
      /** MILU0Preconditioner Initialize
	  Initialize preconditioner for linear system rA*rX=rB
	  @param rA  system matrix.
	  @param rX Unknows vector
	  @param rB Right side linear system of equations.
      */
      virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) 
	{
	  // ILU(0) preconditioner
	  // Incomplete LU factorization with same sparcity pattern as original matrix.
	  // Modified version; see pg. 285 Iterative Methods, Saad. This should be a little
	  // bit better than ILU0 for matrices resulting from FEM.	  
        // The diagonal is included in U. L has non-written 1's on its diagonal.
	  // See pg 274 Iterative methods for linear systems, Yousef Saad 
	  // We assume that, within a row, the entries in A are sorted by increasing j
	  const int size = TSparseSpaceType::Size(rX);
	  int i, j, indexj, oldCountL, oldCountU, newCountL, newCountU, fillL, fillU;
	  int k, indexk, indexkj, indexkjlim, jkj, indexjlim;
	  double aik, r;
	  bool diagFound;

	  mILUSize= size;
	  int n = size;
	  // in case a preconditioner is mistakenly initialized twice;
	  if ( L!=NULL) delete[]  L;
	  if (iL!=NULL) delete[] iL;
	  if (jL!=NULL) delete[] jL;
	  if ( U!=NULL) delete[]  U;
	  if (iU!=NULL) delete[] iU;
	  if (jU!=NULL) delete[] jU;

	  for(i = 0 ; i < size ; i++)
	    std::sort(rA[i].begin() , rA[i].end());

	  // Create copy of matrix split in its L and U parts
	  // Traverse matrix to count elements in rows of L, U
	  // If there is no element in the diagonal, make room for a zero
	  iL=new int[n+1];
	  iU=new int[n+1];
	  for (i=0; i<n+1; i++) {iL[i]=0; iU[i]=0;}
	  for (i=0; i<n; i++) {
		diagFound=false;
		for (typename SparseMatrixType::RowIteratorType row_iterator = rA[i].begin() ; row_iterator != rA[i].end() ; ++row_iterator) {
		      j= row_iterator->first; // Changed by Pooyan!!!
			if (i<=j) iU[i]=iU[i]+1;
			if (i>j)  iL[i]=iL[i]+1;
			if (i==j) diagFound=true;
		}
		if (!diagFound) iU[i]=iU[i]+1;
	}
	// iL[i] is now the number of nonzero entries in the 
	// i-th row of L. Transform to CSR indexes
	oldCountL=0; oldCountU=0; 
	for (i=0; i<n+1; i++) {
		newCountL=oldCountL+iL[i];
		iL[i]=oldCountL;
		oldCountL=newCountL;
		newCountU=oldCountU+iU[i];
		iU[i]=oldCountU;
		oldCountU=newCountU;
	}
	// Traverse again the matrix copying its entries
	L =new double[iL[n]];
	jL=new int   [iL[n]];
	U =new double[iU[n]];
	jU=new int   [iU[n]];
	fillL=0;
	fillU=0;
	for (i=0; i<n; i++) {
		diagFound=false;
		for (typename SparseMatrixType::RowIteratorType row_iterator = rA[i].begin() ; row_iterator != rA[i].end() ; ++row_iterator) {
		      j= row_iterator->first; // Changed by Pooyan!!!
			if (i==j) diagFound=true;
			if ( (j>i) && (!diagFound) ) {
				// This row does not have a diagonal entry. Make it and put a zero.
				jU[fillU]=i;
				U[fillU]=0.00;
				fillU++;
				diagFound=true;
			}
			if (i<=j) {
				jU[fillU]=j; 
				U[fillU]= row_iterator->second; 
				fillU++;
			}
			if (i>j) {
				jL[fillL]=j; 
				L[fillL]= row_iterator->second;
				fillL++;
			}
		}
	}

	// Now comes the real factorization:
	// for i=2, ... ,n
	//    for k=1, ... ,i-1    and (i,k) in nonzero pattern
	//       a[i,k]=a[i,k]/a[k,k]
	//       for j=k+1, ... ,n     and (i,j) in nonzero pattern
	//          a[i,j]=a[i,j]-a[i,k]*a[k,j]
	//       end do
	//    end do
	// end do

	for (i=1; i<n; i++) {
		r=0.00; // residual
		for (indexk=iL[i]; indexk<iL[i+1]; indexk++) {
			k=jL[indexk];
			
			L[indexk]=L[indexk]/U[iU[k]];
			aik=L[indexk];

			indexkj   = iU[k  ];  // traverses row k of U
			indexkjlim= iU[k+1];
			jkj       = jU[indexkj];
			while (jkj<=k) jkj=jU[++indexkj]; // add rows only for j>k

			indexj    = indexk+1; // traverses row i of L beyond k
			indexjlim = iL[i+1];
			j         = jL[indexj];
			while ( (indexkj<indexkjlim) && (indexj<indexjlim) ) {
				if (j==jkj) {
					L[indexj]= L[indexj]-aik*U[indexkj];
					indexj++;
					j=jL[indexj];
					indexkj++;
					jkj=jU[indexkj];
				} else {
					if (j<jkj) {
						indexj++;
						j=jL[indexj];
					} else {
						// These are the entries that go to the residual; because 
						// U(k,j) != 0,  and (i,j) is in the zero pattern, otherwise
						// we would have gone under the (j==jkj) condition.
						r=r+aik*U[indexkj];
						indexkj++;
						jkj=jU[indexkj];
					}
				}
			}

			indexj    = iU[i]; // traverses row i of U
			indexjlim = iU[i+1];
			j         = jU[indexj];
			while ( (indexkj<indexkjlim) && (indexj<indexjlim) ) {
				if (j==jkj) {
					U[indexj]= U[indexj]-aik*U[indexkj];
					indexj++;
					j=jU[indexj];
					indexkj++;
					jkj=jU[indexkj];
				} else {
					if (j<jkj) {
						indexj++;
						j=jU[indexj];
					} else {
						// These are the entries that go to the residual; because 
						// U(k,j) != 0,  and (i,j) is in the zero pattern, otherwise
						// we would have gone under the (j==jkj) condition.
						r=r+aik*U[indexkj];
						indexkj++;
						jkj=jU[indexkj];
					}
				}
			}
		}
	}

	// and that's it...
	// for debugging only; write out matrices L, D, U
/* 	std::cout << "L, U:\n"; */
/* 	for (i=0; i<n; i++) { */
/* 		for (indexj=iL[i]; indexj<iL[i+1]; indexj++) { */
/* 			j=jL[indexj]; */
/* 			std::cout << "L[" << i << "," << j << "]=" << L[indexj] << "\n"; */
/* 		} */
/* 	} */
/* 	for (i=0; i<n; i++) { */
/* 		for (indexj=iU[i]; indexj<iU[i+1]; indexj++) { */
/* 			j=jU[indexj]; */
/* 			std::cout << "U[" << i << "," << j << "]=" << U[indexj] << "\n"; */
/* 		} */
/* 	} */


  for (i=0; i<n; i++) if (U[iU[i]]==0.00) {KRATOS_ERROR(std::runtime_error, "MILU0Preconditioner::Initialize", "Zero in U diagonal found!!", "")}

	}



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
      virtual String Info() const
	{
	  return "MILU0Preconditioner";
	}

      /// Print information about this object.
      virtual void  PrintInfo(std::ostream& OStream) const
	{
	  OStream << "MILU0Preconditioner";
	}

      virtual void PrintData(std::ostream& OStream) const
      {
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
      
        
      ///@}    
        
    }; // Class MILU0Preconditioner 
  
  ///@} 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
  inline std::istream& operator >> (std::istream& IStream, 
 				    MILU0Preconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
   }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
  inline std::ostream& operator << (std::ostream& OStream, 
 				    const MILU0Preconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
   }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_MILU0_PRECONDITIONER_H_INCLUDED  defined 

