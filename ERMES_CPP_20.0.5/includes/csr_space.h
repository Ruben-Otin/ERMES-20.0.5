
#if !defined(KRATOS_CSR_SPACE_H_INCLUDED)
#define KRATOS_CSR_SPACE_H_INCLUDED

#include <numeric>

#include <memory>

#include "../includes/define.h"
#include "../includes/csr_matrix.h"
#include "../includes/vector.h"

namespace Kratos
{
    template<class TDataType>
    class CSRSpace
    {
      public:
    
        typedef std::shared_ptr<CSRSpace> Pointer   ;
        typedef TDataType                   DataType  ;
        typedef CSRMatrix<TDataType>        MatrixType;
        typedef Vector<TDataType>           VectorType;
        typedef unsigned int                IndexType ;
    
        CSRSpace()
        {
        }
    
        virtual ~CSRSpace()
        {
        }

        // return size of vector rV
        static IndexType Size( VectorType& rV )
        {
            return rV.size();
        }
    
        // return number of rows of rM
        static IndexType Size1( MatrixType& rM )
        {
            return rM.RowsNumber();
        }
    
        // return number of columns of rM (ERMES CSR matrices are always square)
        static IndexType Size2( MatrixType& rM )
        {
            return rM.RowsNumber();
        }
    
        // rXi = rMij
        static void GetColumn( int j, MatrixType& rM, VectorType& rX )
        {
            for( int i=0; i<rM.RowsNumber(); i++ ) rX[i] = rM(i,j);
        }
    
        // rY = rX
        static void Copy( VectorType& rX, VectorType& rY )
        {
            rY = rX;
        }
    
        // rX * rY
        static TDataType Dot( VectorType& rX, VectorType& rY )
        {
            return std::inner_product( rX.begin(), rX.end(), rY.begin(), TDataType() );
        }
    
        // ||rX|| (sqrt(rX*rX))
        static double TwoNorm( VectorType& rX )
        {
            return sqrt( std::inner_product( rX.begin(), rX.end(), rX.begin(), double() ) );
        }
    
        // rY = rA*rX
        static void Mult( MatrixType& rA, VectorType& rX, VectorType& rY )
        {
	        typename MatrixType::DataArrayType::iterator row_iterator = rA.Data().begin();
            typename MatrixType::DataArrayType::iterator end_iterator = rA.Data().end  ();

	        typename MatrixType::DataArrayType::value_type::iterator data_iterator    ;
	        typename MatrixType::DataArrayType::value_type::iterator end_data_iterator;
	        
	        typename VectorType::iterator c_iterator = rY.begin();
    
	        double result;
    
	        while(row_iterator != end_iterator)
	        {
	            data_iterator = row_iterator->begin();
	            end_data_iterator = row_iterator->end();
    
	            result = 0;
    
	            while(data_iterator != end_data_iterator)
                {
	  	            result += data_iterator->second * rX[data_iterator->first];
	  	            data_iterator++;
                }
    
	            *c_iterator = result;
    
	            row_iterator++;
	            c_iterator++;
	        }
        }
    
        // rY = rAT * rX
        static void TransposeMult( MatrixType& rA, VectorType& rX, VectorType& rY )
        {
	        typename MatrixType::DataArrayType::iterator row_iterator = rA.Data().begin();
	        typename MatrixType::DataArrayType::iterator end_iterator = rA.Data().end  ();

	        typename MatrixType::DataArrayType::value_type::iterator data_iterator    ;
	        typename MatrixType::DataArrayType::value_type::iterator end_data_iterator;

	        std::fill( rY.begin(), rY.end(), TDataType() );
    
	        typename VectorType::iterator c_iterator = rX.begin();
    
	        while( row_iterator != end_iterator )
	        {
	            data_iterator     = row_iterator->begin();
	            end_data_iterator = row_iterator->end  ();
    
	            while( data_iterator != end_data_iterator )
                {
	  	            rY[data_iterator->first] += data_iterator->second * *c_iterator;
	  	            data_iterator++;
                }
    
	            row_iterator++;
	            c_iterator++;
	        }
        }
    
        // rZ = (A * rX) + (B * rY)
        static void ScaleAndAdd( double A, VectorType& rX, double B, VectorType& rY, VectorType& rZ )
        {
	        typename VectorType::iterator x_iterator   = rX.begin();
	        typename VectorType::iterator y_iterator   = rY.begin();
	        typename VectorType::iterator z_iterator   = rZ.begin();
	        typename VectorType::iterator end_iterator = rX.end  ();
    
	        while( x_iterator != end_iterator )
            {
	            *z_iterator++ = ( A * *x_iterator++ ) + ( B * *y_iterator++ );
            }
        }
    
        // rY = (A * rX) + (B * rY)
        static void ScaleAndAdd( double A, VectorType& rX, double B, VectorType& rY )
        {
            typename VectorType::iterator x_iterator = rX.begin();
	        typename VectorType::iterator y_iterator = rY.begin();
    
            while( x_iterator != rX.end() )
	        {
	            *y_iterator = ( A * (*x_iterator) ) + ( B * (*y_iterator) );
    
                y_iterator++;
                x_iterator++;
	        }
        }
    
        // rY = (B * rY)
        static void Scale( double B, VectorType& rY )
        {
	        typename VectorType::iterator y_iterator = rY.begin();
    
            for( y_iterator = rY.begin(); y_iterator != rY.end(); y_iterator++ )
	        {
	            *y_iterator = ( B * (*y_iterator) );
	        }
        }
    
        // rA[i] * rX
        static double RowDot( unsigned int i, MatrixType& rA, VectorType& rX )
        {
	        typename MatrixType::DataArrayType::value_type::iterator data_iterator     = rA.Data()[i].begin();
	        typename MatrixType::DataArrayType::value_type::iterator end_data_iterator = rA.Data()[i].end  ();
    
	        TDataType result = TDataType();
    
	        while( data_iterator != end_data_iterator )
	        {
	            result += data_iterator->second * rX[data_iterator->first];
	            data_iterator++;
	        }
    
	        return result;
        }
    
        // rX = A
        static void Set( VectorType& rX, TDataType A )
        {
            std::fill( rX.begin(), rX.end(), A );
        }
    
        static void GraphNeighbors( unsigned int i, MatrixType& rA, Vector<IndexType>& rConectivity )
        {
	        typename MatrixType::DataArrayType::value_type::iterator data_iterator     = rA.Data()[i].begin();
	        typename MatrixType::DataArrayType::value_type::iterator end_data_iterator = rA.Data()[i].end  ();
    
	        rConectivity.clear();
    
	        while( data_iterator != end_data_iterator )
	        {
	            rConectivity.push_back( data_iterator->first );
	  	        data_iterator++;
	        }
        }
    
        static unsigned int GraphDegree( unsigned int i, MatrixType& rA )
        {
	        return rA.Data()[i].size();
        }
    
	    static void ReorderingMatrix( MatrixType& rA, Vector<unsigned int>& rPermutation )
	    {
	        const int size = Size1( rA );
	        unsigned int i;
    
	        Vector<unsigned int> inverse_permutation( size );
    
	        for( i=0; i<size; i++ ) 
            {
                inverse_permutation[rPermutation[i]] = i;
            }
    
	        MatrixType temp_matrix( size );
    
	        typename MatrixType::DataArrayType::iterator row_iterator = rA.Data().begin();
	        typename MatrixType::DataArrayType::iterator end_iterator = rA.Data().end  ();

	        typename MatrixType::DataArrayType::value_type::iterator data_iterator    ;
	        typename MatrixType::DataArrayType::value_type::iterator end_data_iterator;

	        i = 0;

	        while( row_iterator != end_iterator )
	        {
	            data_iterator     = row_iterator->begin();
	            end_data_iterator = row_iterator->end  ();
    
	            while( data_iterator != end_data_iterator )
                {
	  	            temp_matrix( inverse_permutation[i], inverse_permutation[data_iterator->first] ) = data_iterator->second;
	  	            data_iterator++;
                }
    
	            row_iterator->clear();
                MatrixType::DataArrayType().swap( *row_iterator );
    
	            row_iterator++;
	            i++;
	        }
    
	        rA.Data().swap( temp_matrix.Data() );
	    }
    
	    static void ReorderingVector( Vector<double>& rB, Vector<unsigned int>& rPermutation )
	    {
	  	    const int size = Size( rB );
    
	  	    VectorType temp_vector( rB );
    
	  	    for( int i=0; i<size; i++ ) rB[i] = temp_vector[rPermutation[i]];
	    }
    
	    static void RecoverReorderingVector( Vector<double>& rB, Vector<unsigned int>& rPermutation )
	    {
	  	    const int size = Size( rB );
    
	  	    VectorType temp_vector( rB );
    
	  	    for( int i=0; i<size; i++ ) rB[rPermutation[i]] = temp_vector[i];
	    }
    
        // Print information about this object.
        virtual void PrintInfo( std::ostream& OStream ) const
        {
	        OStream << "Compressed sparse row space";
        }
    
        // Print object's data.
        virtual void PrintData( std::ostream& OStream ) const
        {
        }
    
      private:
    
        // Assignment operator.
        CSRSpace& operator=( const CSRSpace& Other );
    
        // Copy constructor.
        CSRSpace( const CSRSpace& Other );
    }; 
    
    // input stream function
    template<class TDataType>
    inline std::istream& operator >> ( std::istream& IStream, CSRSpace<TDataType>& rThis );
    
    // output stream function
    template<class TDataType>
    inline std::ostream& operator << ( std::ostream& OStream, const CSRSpace<TDataType>& rThis )
    {
        rThis.PrintInfo( OStream );
        OStream << std::endl;
        rThis.PrintData( OStream );
    
	    return OStream;
    }
} 

#endif 


