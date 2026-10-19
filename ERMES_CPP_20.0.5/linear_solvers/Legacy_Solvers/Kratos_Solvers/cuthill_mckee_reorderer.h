//************************************************************   
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: pooyan $
//*   Date:                $Date: 2004/05/21 10:15:46 $
//*   Revision:            $Revision: 1.4 $
//************************************************************


#if !defined(KRATOS_CUTHILL_MCKEE_REORDERER_H_INCLUDED )
#define  KRATOS_CUTHILL_MCKEE_REORDERER_H_INCLUDED



// System includes 
#include <vector>

// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"


// Project includes
#include "../includes/define.h"
#include "../vectorial_spaces/vector.h"


namespace Kratos
{
 
  /// 
  /** .
      two template parameter: 
      - TSparseSpaceType which specify type
        of the unknowns, coefficients, sparse matrix, vector of
	unknowns, right hand side vector and their respective operators.
      - TDenseMatrixType which specify type of the
        matrices used as temporary matrices or multi solve unknowns and
	right hand sides and their operators.  
  */
template<class TSparseSpaceType, class TDenseSpaceType>
class CuthillMcKeeReorderer : public Reorderer<TSparseSpaceType, TDenseSpaceType>
{
    public:
      
      /// Counted pointer of CuthillMcKeeReorderer
      typedef boost::shared_ptr<CuthillMcKeeReorderer> Pointer;

      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

      typedef unsigned int IndexType;

      typedef Vector<IndexType> IndexVectorType;
  
      /// Default constructor.
      CuthillMcKeeReorderer(){}

       /// Copy constructor.
      CuthillMcKeeReorderer(const CuthillMcKeeReorderer& Other) : mIndexPermutation(Other.mIndexPermutation){}

     /// Destructor.
      virtual ~CuthillMcKeeReorderer(){}
      
      virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	      const unsigned int size = TSparseSpaceType::Size1(rA);
          unsigned int i;

	      mIndexPermutation.resize(size);

	      for(i = 0 ; i < size ; i++) mIndexPermutation[i] = size -i - 1;

 	      CalculateIndexPermutation(rA);
      }




      virtual void InverseReorder(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	      const unsigned int size = TSparseSpaceType::Size1(rA);
          unsigned int i;
	      VectorType temp(size);

	      for(i = 0 ; i < size ; i++) temp[mIndexPermutation[i]] = rX[i];
	      
	      rX = temp;

	      for(i = 0 ; i < size ; i++) temp[mIndexPermutation[i]] = rB[i];
	          
	      rB = temp;
      }
      

      virtual void Reorder(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	      const unsigned int size = TSparseSpaceType::Size1(rA);
          unsigned int i;
	      VectorType temp(size);

	      TSparseSpaceType::ReorderingMatrix(rA, mIndexPermutation);

	      for(i = 0 ; i < size ; i++) temp[i] = rX[mIndexPermutation[i]];
  
	      rX = temp;

	      for(i = 0 ; i < size ; i++) temp[i] = rB[mIndexPermutation[i]];
   
	      rB = temp;
      }
      


      virtual IndexVectorType& CalculateIndexPermutation(SparseMatrixType& rA, IndexType InitialIndex = IndexType())
      {
        const unsigned int size = TSparseSpaceType::Size1(rA);

        mIndexPermutation.resize(size);

        Vector<bool> is_marked(size, false);
        IndexVectorType level_set(1, InitialIndex);
        IndexVectorType conectivity;
        IndexVectorType next_level_set;


        mIndexPermutation[0] = InitialIndex;
        is_marked[InitialIndex] = true;


        unsigned int next = 1;

        while(next < size)
        {
            for(IndexVectorType::iterator i = level_set.begin() ; i != level_set.end() ; ++i)
            {
                TSparseSpaceType::GraphNeighbors(*i, rA, conectivity);

                for(IndexVectorType::iterator j = conectivity.begin() ; j != conectivity.end() ; ++j)
                if(is_marked[*j] == false)
                {
                    mIndexPermutation[next++] = *j;
                    is_marked[*j] = true;
                    next_level_set.push_back(*j);
                }

            }

            level_set = next_level_set;
            next_level_set.clear();
            if(level_set.empty() && (next < size)) // No conected graph
            {
                int k = 0;
                while(is_marked[k])
                k++;
                level_set.resize(1);
                level_set[0] = k;
                mIndexPermutation[next] = k;
                is_marked[k] = true;
            }
        }
        return mIndexPermutation;
      }
      
     

      const IndexVectorType& GetIndexPermutationVector() const
      {
	      return mIndexPermutation;
      }
                
    private:
     

      IndexVectorType mIndexPermutation;
                 
      /// Assignment operator.
      CuthillMcKeeReorderer& operator=(const CuthillMcKeeReorderer& Other);
       
    }; // Class CuthillMcKeeReorderer 
  
}  // namespace Kratos.

#endif // KRATOS_CUTHILL_MCKEE_REORDERER_H_INCLUDED  defined 


