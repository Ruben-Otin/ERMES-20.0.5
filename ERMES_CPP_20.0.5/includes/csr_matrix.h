
#if !defined(KRATOS_CSR_MATRIX_H_INCLUDED)
#define KRATOS_CSR_MATRIX_H_INCLUDED

#include <memory>
#include <algorithm>
#include <vector>

#include "../includes/define.h"
#include "vector.h"

namespace Kratos
{
    template<class TDataType>
    class CSRMatrix
    {
      public:
    
        // Counted pointer of CSRMatrix
        typedef std::shared_ptr<CSRMatrix> Pointer;
    
        /** Index type used in this matrix */
        typedef unsigned int IndexType;
    
        /** Defining a row element as a pair of column's number and value */
        typedef std::pair<IndexType, TDataType> DataSetType;
    
        /** internal data array type which is a vector of vectors of each rows */
        typedef Vector< std::vector<DataSetType> > DataArrayType;
    
        //only for direct solvers
        typedef std::pair<IndexType, double>              DataSetTypeDouble  ;
        typedef Vector< std::vector<DataSetTypeDouble> >  DataArrayTypeDouble;
        typedef DataArrayTypeDouble::value_type::iterator RowIteratorType    ;
    
        // Default constructor.
        CSRMatrix() : mData()
        {
        }
    
        CSRMatrix( IndexType size ) : mData( size )
        {
        }
    
	    // Copy constructor.
        CSRMatrix( const CSRMatrix& Other ) : mData( Other.mData )
        {
        }
    
        // Destructor.
        virtual ~CSRMatrix()
        {
        }
    
        class IndexCheck
        {
	        IndexType mI;
    
          public:
    
	        IndexCheck( int I ) : mI( I )
            {
            }
            
	        bool operator()( const DataSetType& I )
            {
                return I.first == mI;
            }
        };
    
        TDataType& operator()( int I, int J )
        {
	        typename DataArrayType::value_type::iterator member;
    
	        if ( ( member = std::find_if( mData[I].begin(), mData[I].end(), IndexCheck(J) ) ) != mData[I].end() )
            {
	            return member->second;
            }
    
	        return mData[I].insert( mData[I].end(), DataSetType( J, 0.00 ) )->second;
        }
    
        TDataType operator()( int I, int J ) const
        {
            typename DataArrayType::value_type::const_iterator member;
    
            if ( ( member = std::find_if( mData[I].begin(), mData[I].end(), IndexCheck(J) ) ) != mData[I].end() )
            {
                return member->second;
            }
    
            return TDataType();
        }
    
        TDataType getValue( int I, int J )
        {
            typename DataArrayType::value_type::const_iterator member;
    
            if ( ( member = std::find_if( mData[I].begin(), mData[I].end(), IndexCheck(J) ) ) != mData[I].end() )
            {
                return member->second;
            }
    
            return TDataType();
        }
    
        std::vector<DataSetType>& operator[]( int I )
        {
            return mData[I];
        }
    
        void Resize( IndexType Size )
        {
	        if( Size != mData.size() ) 
            {
                mData.resize( Size );
            }
        }
    
	    DataArrayType& Data()
        {
            return mData;
        }
    
	    IndexType RowsNumber() const
        {
            return mData.size();
        }
    
        // Free memory
        void FreeData()
        {
            for ( int i=0; i<mData.size(); i++ )
            {
                mData[i].clear();

                std::vector<DataSetType>().swap( mData[i] );
            }
    
            mData.FreeData();
        }
    
        // Print information about this object.
        virtual void PrintInfo( std::ostream& OStream ) const
        {
	        OStream << "Compressed row sparse matrix( " << RowsNumber() << " * " << RowsNumber() << " )";
        }
    
        // Print object's data.
        virtual void PrintData( std::ostream& OStream ) const
        {
	        OStream.setf ( std::ios::scientific, std::ios::floatfield );

	        OStream.width( 5 );

	        for ( typename CSRMatrix<TDataType>::IndexType i=0; i<RowsNumber(); i++ )
	        {
	            for( int j=0; j<RowsNumber(); j++ ) OStream << operator()(i,j) << " , ";

	            OStream << std::endl;
	        }

	        OStream.width (                      );
	        OStream.unsetf( std::ios::scientific );
	        OStream.unsetf( std::ios::floatfield );
        }
    
        private:
    
          DataArrayType mData;
    
          // Assignment operator.
          CSRMatrix& operator=( const CSRMatrix& Other );
    }; 

    // input stream function
    template<class TDataType>
    inline std::istream& operator >> ( std::istream& IStream, CSRMatrix<TDataType>& rThis );

    // output stream function
    template<class TDataType>
    inline std::ostream& operator << ( std::ostream& OStream, const CSRMatrix<TDataType>& rThis )
    {
          rThis.PrintInfo( OStream );
          OStream << std::endl;
          rThis.PrintData( OStream );

	      return OStream;
    }
}

#endif 


