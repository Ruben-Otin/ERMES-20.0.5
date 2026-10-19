//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2004/02/26 12:18:15 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_MILUN_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_MILUN_PRECONDITIONER_H_INCLUDED



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
  
  /// MILUNPreconditioner class. 
  /**   */
  template<class TSparseSpaceType, class TDenseSpaceType>
    class MILUNPreconditioner : public ILUPreconditioner<TSparseSpaceType, TDenseSpaceType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of MILUNPreconditioner
      typedef boost::shared_ptr<MILUNPreconditioner> Pointer;

      typedef ILUPreconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      MILUNPreconditioner(int NewLevelNeighbors = 0) : levelNeighbors(NewLevelNeighbors) {}

      /// Copy constructor.
      MILUNPreconditioner(const MILUNPreconditioner& Other) : levelNeighbors(Other.levelNeighbors){}

      /// Destructor.
      virtual ~MILUNPreconditioner()
      {	
      }
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      MILUNPreconditioner& operator=(const MILUNPreconditioner& Other)
	{
	  BaseType::operator=(Other);
	  levelNeighbors = Other.levelNeighbors;
 	  return *this;
	}


      void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) 
      {
          BaseType::Initialize(rA, rX, rB);
      }

      
      ///@}
      ///@name Operations
      ///@{
      
      /** MILUNPreconditioner Initialize
	  Initialize preconditioner for linear system rA*rX=rB
	  @param rA  system matrix.
	  @param rX Unknows vector
	  @param rB Right side linear system of equations.
      */
      virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) 
	{
	// MILUN(levelNeighbors) preconditioner
	// Santiago Egido, 2004
	// Incomplete LU factorization with sparcity pattern given by neighbors graph.
	// Modified version; see pg. 285 Iterative Methods, Saad. This should be a little
	// bit better than ILU0 for matrices resulting from FEM.	
	// The matrix is assumed to have symmetric nonzero structure.
	// Zeros can be added on the diagonal.
	// The diagonal of the factorization is included in U. 
	// L has non-written 1's on its diagonal.
	int i, j, indexj, oldCountL, oldCountU, newCountL, newCountU, fillL, fillU;
	int k, indexk, indexkj, indexkjlim, jkj, indexjlim, passes, ii, nnz, nloop;
	int nextia, nextja;
	double aik, r;
	bool diagFound, sorted;
	int *nz, *linz, *iNZ, *jNZ;

	mILUSize=TSparseSpaceType::Size(rX);

	int n = mILUSize;
	// in case a preconditioner is mistakenly initialized twice;
	if ( L!=NULL) delete[]  L;
	if (iL!=NULL) delete[] iL;
	if (jL!=NULL) delete[] jL;
	if ( U!=NULL) delete[]  U;
	if (iU!=NULL) delete[] iU;
	if (jU!=NULL) delete[] jU;

	// make sure the rows of A are sorted by increasing j
	  for(i = 0 ; i < n ; i++)
	    std::sort(rA[i].begin() , rA[i].end());

	// Compute nonzero structure
	// The nonzero structure of the factorization is that of (A+I)^(p+1), so we will do:
	// 1) L=Identity
	// 2) do p+1 times: 
	//        U=A*L
	//        L=U
	
	// create Identity structure in L
	iL=new int[n+1];
	jL=new int[n];
	for (i=0; i<n; i++) {iL[i]=i; jL[i]=i;};
	iL[n]=n;

	nz=new int[n];   // used to mark nonzero elements in row
	for (ii=0; ii<n; ii++) nz[ii]=0;
	linz=new int[n]; // list of indexes of nonzero elements in row
	nnz=0;           // number of nonzeros in row
	for (nloop=0; nloop<=levelNeighbors; nloop++) {

		// U=A*L
		// start counting how many entries will be required in each row of U
		iU=new int[n+1];
		for (ii=n-1; ii>0; ii--) iU[0]=0;
		for (i=0; i<n; i++) {
			// make sure the diagonal appears
			nnz=1;
			linz[0]=i;
			nz[i]=1;
                  for (typename SparseMatrixType::RowIteratorType row_iterator = rA[i].begin() ; row_iterator != rA[i].end() ; ++row_iterator) {
		            j= row_iterator->first; 
				// j has to appear
				if (nz[j]==0) {
					nz[j]=1;
					linz[nnz++]=j;
				}
				for (indexk=iL[j]; indexk<iL[j+1]; indexk++) {
					k=jL[indexk];
					// A(i,j) != 0   and   L(j,k) != 0 , so k has to appear
					if (nz[k]==0) {
						nz[k]=1;
						linz[nnz++]=k;
					}
				}
			}
			// store in iU[i+1] the number of nonzeros in row i of product
			iU[i+1]=nnz;
			// and restore linz, nz, nnz for the next iteration
			for (ii=0; ii<nnz; ii++) nz[linz[ii]]=0;
		}

		// Ok, now we know how large each row of U should be.
		// now convert the iU counts into CSR indexes
		for (ii=0; ii<n; ii++) iU[ii+1]=iU[ii]+iU[ii+1];
		// allocate jL2
		jU=new int[iU[n]];
		// and repeat the previus loop, storing the j's into jL2 now that we 
		// know how many of them there will be. Use iL2 as pointers to storage.
		for (i=0; i<n; i++) {
			// make sure the diagonal appears
			nnz=1;
			linz[0]=i;
			nz[i]=1;
                  for (typename SparseMatrixType::RowIteratorType row_iterator = rA[i].begin() ; row_iterator != rA[i].end() ; ++row_iterator) {
		            j= row_iterator->first; 
				// j has to appear
				if (nz[j]==0) {
					nz[j]=1;
					linz[nnz++]=j;
				}
				for (indexk=iL[j]; indexk<iL[j+1]; indexk++) {
					k=jL[indexk];
					// A(i,j) != 0   and   L(j,k) != 0 , so k has to appear
					if (nz[k]==0) {
						nz[k]=1;
						linz[nnz++]=k;
					}
				}
			}
			// store column indexes into jU
			for (ii=0; ii<nnz; ii++) {
				jU[iU[i]+ii]=linz[ii];
				nz[linz[ii]]=0;
			}
		}

		// L=U
		delete[] iL;
		delete[] jL;
		iL=iU;
		jL=jU;
	} // end for nloop=0; nloop<=levelNeighbors

	// store nonzero pattern in NZ and make sure rows are sorted by increasing j
	iNZ=iL;
	jNZ=jL;
	for (i=0; i<n; i++) {
		std::sort( jNZ+iNZ[i] , jNZ+iNZ[i+1] );
		/*
		sorted=false;
		passes=0;
		while (!sorted) {
			sorted=true;
			for (indexj=iNZ[i]; indexj<iNZ[i+1]-1-passes; indexj++) {
				if (jNZ[indexj]>jNZ[indexj+1]) {
					sorted=false;
					k=jNZ[indexj];
					jNZ[indexj]=jNZ[indexj+1];
					jNZ[indexj+1]=k;
				}
			}
			passes=passes+1; // biggest j has moved to the end, don't sort it again
		}
		*/
	}

	// split NZ into its L and U parts
	// Traverse matrix NZ to count elements in rows of L, U
	iL=new int[n+1];
	iU=new int[n+1];
	for (i=0; i<n+1; i++) {
		iL[i]=0; 
		iU[i]=0;
	}
	for (i=0; i<n; i++) {
		diagFound=false;
		for (indexj=iNZ[i]; indexj<iNZ[i+1]; indexj++) {
			j=jNZ[indexj];
			if (i<=j) iU[i]=iU[i]+1;
			if (i>j)  iL[i]=iL[i]+1;
			if (i==j) diagFound=true;
		}
		if (!diagFound) iU[i]=iU[i]+1;
	}
	// iL[i] is now the number of nonzero entries in the 
	// i-th row of L. Transform to CSR indexes
	oldCountL=0; 
	oldCountU=0; 
	for (i=0; i<n+1; i++) {
		newCountL=oldCountL+iL[i];
		iL[i]=oldCountL;
		oldCountL=newCountL;
		newCountU=oldCountU+iU[i];
		iU[i]=oldCountU;
		oldCountU=newCountU;
	}
	// Traverse again the matrix copying its entries. There is a difference
	// now with respect to ILU0: the matrix NX has no entries; therefore, 
	// typically we will fill in L and U with zeros, until we are "lucky" and
	// hit an entry in the matrix A.

    std::cout << "Memory required by ILUN preconditioner (Mb): " << (int)((iL[n]+iU[n])*12/1048576.00) << std::endl;

	L =new double[iL[n]];
	jL=new int   [iL[n]];
	U =new double[iU[n]];
    jU=new int   [iU[n]];

	fillL=0;
	fillU=0;
	nextia=0; // it will be a bad day when the first entry in the matrix is in the second row...
      typename SparseMatrixType::RowIteratorType nextpa = rA[0].begin(); // pointer to next entry of rA
      typename SparseMatrixType::RowIteratorType nextpalim = rA[0].end(); 
	nextja=nextpa->first; // j corresponding to next entry of a
	for (i=0; i<n; i++) {
		diagFound=false;
		for (indexj=iNZ[i]; indexj<iNZ[i+1]; indexj++) {
			j=jNZ[indexj];
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
				U[fillU]=0.00;
				if ((i==nextia)&&(j==nextja)) { // this corresponds to an entry of a
					U[fillU]=nextpa->second;
					nextpa++;
					if ( nextpa==nextpalim ) { // move to next row
						if (++nextia==n) {
							nextia=-1; // end of matrix reached
						} else {
							nextpa = rA[nextia].begin(); 
							nextpalim = rA[nextia].end(); 
						}
					}
					nextja=nextpa->first;
				}
				fillU++;
			}
			if (i>j) {
				jL[fillL]=j; 
				L[fillL]=0.00;
				if ((i==nextia)&&(j==nextja)) { // this corresponds to an entry of a
					L[fillL]=nextpa->second;
					nextpa++;
					if ( nextpa==nextpalim ) { // move to next row
						if (++nextia==n) {
							nextia=-1; // end of matrix reached
						} else {
							nextpa = rA[nextia].begin(); 
							nextpalim = rA[nextia].end(); 
						}
					}
					nextja=nextpa->first;
				}
				fillL++;
			}
		}
	}

	// at this point, NZ is not needed any more
	delete[] iNZ;
	delete[] jNZ;

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
	// for debugging only; write out matrices L, U
	/*
	std::cout << "L, U:\n";
	for (i=0; i<n; i++) {
		for (indexj=iL[i]; indexj<iL[i+1]; indexj++) {
			j=jL[indexj];
			std::cout << "L[" << i << "," << j << "]=" << L[indexj] << "\n";
		}
	}
	for (i=0; i<n; i++) {
		for (indexj=iU[i]; indexj<iU[i+1]; indexj++) {
			j=jU[indexj];
			std::cout << "U[" << i << "," << j << "]=" << U[indexj] << "\n";
		}
	}
	*/

  for (i=0; i<n; i++) if (U[iU[i]]==0.00) {KRATOS_ERROR(std::runtime_error, "MILUNPreconditioner::Initialize", "Zero in U diagonal found!!", "")}

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
	  return "MILUNPreconditioner";
	}

      /// Print information about this object.
      virtual void  PrintInfo(std::ostream& OStream) const
	{
	  OStream << "MILUNPreconditioner";
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

      int levelNeighbors;
       
        
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
        
    }; // Class MILUNPreconditioner 
  
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
 				    MILUNPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
   }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
  inline std::ostream& operator << (std::ostream& OStream, 
 				    const MILUNPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
   }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_ILU0_PRECONDITIONER_H_INCLUDED  defined 

