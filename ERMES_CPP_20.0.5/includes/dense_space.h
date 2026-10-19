
#if !defined(KRATOS_DENSE_SPACE_H_INCLUDED)
#define KRATOS_DENSE_SPACE_H_INCLUDED

#include <memory>

#include "../includes/define.h"
#include "../includes/csr_matrix.h"
#include "../includes/vector.h"

namespace Kratos
{
    template<class TDataType>
    class DenseSpace
    {
      public:
       
        typedef std::shared_ptr<DenseSpace>  Pointer   ;
        typedef TDataType                      DataType  ;
        typedef Matrix<TDataType>              MatrixType;
        typedef Vector<TDataType>              VectorType;
        typedef typename MatrixType::IndexType IndexType ;
    
        DenseSpace()
        {
        }
    
        virtual ~DenseSpace()
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
    
        // return number of columns of rM
        static IndexType Size2( MatrixType& rM )
        {
            return rM.ColumnsNumber();
        } 
    
        // rXi = rMij
        static void GetColumn( int j, MatrixType& rM, VectorType& rX )
        {
            for( int i=0; i<rM.RowsNumber(); i++ ) 
            {
                rX[i] = rM(i,j);
            }
        } 
    
        // rMij = rXi
        static void SetColumn( int j, MatrixType& rM, VectorType& rX )
        {
            for( int i=0; i<rM.RowsNumber(); i++ ) 
            {
                rM(i,j) = rX[i];
            }
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
    
        // ||rX||2
        static double TwoNorm( VectorType& rX )
        {
            return sqrt( std::inner_product( rX.begin(), rX.end(), rX.begin(), double() ) );
        } 
    
        static void Mult( MatrixType& rA, VectorType& rX, VectorType& rY )
        {
	        IndexType size = Size1( rA );

	        for( int i=0; i<size; i++ ) 
            {
                rY[i] = Dot(rA[i], rX);
            }
        }
    
        // rY = rAT * rX
        static void TransposeMult( MatrixType& rA, VectorType& rX, VectorType& rY )
        {
	        IndexType size = Size1( rA );
    
	        rY = TDataType();
    
	        for( int i=0; i<size; i++ )
            {
	            for( int j=0; j<size; j++ )
                {
	                rY[i] += Dot( rA(j,i), rX[j] );
                }
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
	            *z_iterator++ = (A * *x_iterator++) + (B * *y_iterator++);
            }
        } 

        // rY = (A * rX) + (B * rY) 
        static void ScaleAndAdd( double A, VectorType& rX, double B, VectorType& rY ) 
        {
	        typename VectorType::iterator x_iterator   = rX.begin();
	        typename VectorType::iterator y_iterator   = rY.begin();
	        typename VectorType::iterator end_iterator = rX.end  ();
	        
	        double c = B - double(1.00);
    
	        while( x_iterator != end_iterator )
	        {
	            *y_iterator += (A * *x_iterator++) + (c * *y_iterator);
	            y_iterator++;
	        }
        } 
      
        // Print information about this object.
        virtual void PrintInfo( std::ostream& OStream ) const
        {
	        OStream << "Dense space";
        }

        // Print object's data.
        virtual void PrintData( std::ostream& OStream ) const
        {
        }
                  
      private:
           
        // Assignment operator.
        DenseSpace& operator=( const DenseSpace& Other );
    
        // Copy constructor.
        DenseSpace( const DenseSpace& Other );
         
    }; 
    
    // input stream function
    template<class TDataType>
    inline std::istream& operator >> ( std::istream& IStream, DenseSpace<TDataType>& rThis );
    
    // output stream function
    template<class TDataType>
    inline std::ostream& operator << ( std::ostream& OStream, const DenseSpace<TDataType>& rThis )
    {
        rThis.PrintInfo( OStream );
        OStream << std::endl;
        rThis.PrintData( OStream );
    
	    return OStream;
    }
}  

#endif


