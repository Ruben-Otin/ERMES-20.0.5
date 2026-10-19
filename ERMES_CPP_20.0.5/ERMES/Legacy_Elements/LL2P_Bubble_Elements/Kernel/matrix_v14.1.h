//**********************************************************
//*   Author:    $Author: $
//*   Date:      $Date: 2003/09/18 13:58:31 $
//*   Revision:  $Revision: 1.4 $
//**********************************************************

#if !defined(DEFINE_KRATOS_MATRIX )
#define  DEFINE_KRATOS_MATRIX

#include <vector>
#include <numeric>
#include <iostream>
#include <complex>

#include "vector.h"
#include "../includes/exception.h"

namespace Kratos
{
    template<class TDataType>
    class Matrix
    {
      public:

        typedef Matrix<TDataType> MatrixType;
        typedef Vector<TDataType> VectorType;
        
        typedef std::vector<VectorType> ContainerType;
        typedef unsigned int IndexType;
    
        Matrix(): mRowsNumber      (IndexType()),
                  mColumnsNumber   (IndexType()),
                  mRowsIndexBase   (IndexType()),
                  mColumnsIndexBase(IndexType()),
                  mData            ()
        {
        }
    
        Matrix(IndexType NewRowsNumber   ,
               IndexType NewColumnsNumber,
               IndexType NewRowsIndexBase    = IndexType(),
               IndexType NewColumnsIndexBase = IndexType()): mRowsNumber      (NewRowsNumber)      ,
                                                             mColumnsNumber   (NewColumnsNumber)   ,
                                                             mRowsIndexBase   (NewRowsIndexBase)   ,
                                                             mColumnsIndexBase(NewColumnsIndexBase),
                                                             mData            (NewRowsNumber, Vector<TDataType>(NewColumnsNumber, TDataType()))
        {
        }
    
        Matrix(const Vector<TDataType>& Other): mRowsNumber      (1)           ,
                                                mColumnsNumber   (Other.size()),
                                                mRowsIndexBase   (IndexType()) ,
                                                mColumnsIndexBase(IndexType()) ,
                                                mData            (1,Other)
        {
        }
    
        Matrix(const ContainerType& Other)
        {
            mRowsIndexBase    = IndexType();
            mColumnsIndexBase = IndexType();
            mRowsNumber       = Other.size();
    
            if(mRowsNumber == 0)
            {
                mColumnsNumber = 0;
                return;
            }
    
            mColumnsNumber = Other[1].size();
    
            for(typename ContainerType::const_iterator i = Other.begin(); i != Other.end(); ++i)
            if(i->size() != mColumnsNumber)
            KRATOS_ERROR(std::invalid_argument,
            "Matrix::Matrix(const Vector<Vector<TDataType> >& Other)",
            "Column's number are not equal for all the rows", *this)
    
            mData = Other;
        }
    
        const TDataType& operator()(IndexType I, IndexType J) const
        {
            KRATOS_TRY
                return mData[I-mRowsIndexBase][J-mColumnsIndexBase];
            KRATOS_CATCH("Matrix::operator() const", *this)
        }
    
        TDataType& operator()(IndexType I, IndexType J)
        {
            KRATOS_TRY
                return mData[I-mRowsIndexBase][J-mColumnsIndexBase];
            KRATOS_CATCH("Matrix::operator()", *this)
        }
    
        Vector<TDataType>& operator[](IndexType I)
        {
            return mData[I];
        }
    
        const Vector<TDataType>& operator[](IndexType I) const
        {
            return mData[I-mRowsIndexBase];
        }
    
        Matrix& operator*=(double ScaleValue)
        {
            for(typename ContainerType::iterator i = mData.begin() ; i != mData.end() ; ++i)
            for(IndexType j = 0 ; j < mColumnsNumber ; j++)
            (*i)[j] *= ScaleValue;
    
            return *this;
        }
    
        Matrix& operator*=(std::complex<double> ScaleValue)
        {
            for(typename ContainerType::iterator i = mData.begin() ; i != mData.end() ; ++i)
            for(IndexType j = 0 ; j < mColumnsNumber ; j++)
            (*i)[j] *= ScaleValue;
    
            return *this;
        }
    
        Matrix& operator/=(double ScaleValue)
        {
            for(typename ContainerType::iterator i = mData.begin() ; i != mData.end() ; ++i)
            for(IndexType j = 0 ; j < mColumnsNumber ; j++)
            (*i)[j] /= ScaleValue;
            return *this;
        }
    
        Matrix& operator/=(std::complex<double> ScaleValue)
        {
            for(typename ContainerType::iterator i = mData.begin() ; i != mData.end() ; ++i)
            for(IndexType j = 0 ; j < mColumnsNumber ; j++)
            (*i)[j] /= ScaleValue;
            return *this;
        }
    
        MatrixType& operator+=(const MatrixType& rOtherMatrix)
        {
            if( (mRowsNumber != rOtherMatrix.mRowsNumber) || (mColumnsNumber != rOtherMatrix.mColumnsNumber) )
            {
                String buffer;
                buffer << "Non consistance matrices with matrix" <<
                rOtherMatrix.mRowsNumber << "*" << rOtherMatrix.mColumnsNumber;
                KRATOS_ERROR(std::invalid_argument, "Matrix operator+=", buffer, *this)
            }
    
            for(IndexType i = 0 ; i < mRowsNumber ; i++)
            transform(mData[i].begin(),
            mData[i].end(),
            rOtherMatrix.mData[i].begin(),
            mData[i].begin(),
            std::plus<TDataType>());
            return *this;
        }
    
        MatrixType& operator-=(const MatrixType& rOtherMatrix)
        {
            if( (mRowsNumber != rOtherMatrix.mRowsNumber) || (mColumnsNumber != rOtherMatrix.mColumnsNumber) )
            {
                String buffer;
                buffer << "Non consistance matrices with matrix" <<
                rOtherMatrix.mRowsNumber << "*" << rOtherMatrix.mColumnsNumber;
                KRATOS_ERROR(std::invalid_argument, "Matrix operator-=", buffer, *this)
            }
    
            for(IndexType i = 0 ; i < mRowsNumber ; i++)
            transform(mData[i].begin(),
            mData[i].end(),
            rOtherMatrix.mData[i].begin(),
            mData[i].begin(),
            std::minus<TDataType>());
            return *this;
        }
    
        Matrix<TDataType> operator*(const MatrixType& rOtherMatrix)
        {
            if( mColumnsNumber != rOtherMatrix.mRowsNumber )
            KRATOS_ERROR(std::invalid_argument, "Matrix operator*", "Non consistance matrices", *this)
    
            Matrix<TDataType> result(mRowsNumber, rOtherMatrix.mColumnsNumber);
            for(IndexType i = 0 ; i < mRowsNumber ; ++i)
            for(IndexType j = 0 ; j < rOtherMatrix.mColumnsNumber ; j++)
            for(IndexType k = 0 ; k < rOtherMatrix.mRowsNumber ; k++)
            result(i,j) += mData[i][k] * rOtherMatrix.mData[k][j];
            return result;
        }
    
        VectorType operator*(const VectorType& rVector)
        {
            if( mColumnsNumber != rVector.size() )
            KRATOS_ERROR(std::invalid_argument, "Matrix-Vector operator*", "Non consistance matrix-vector", *this)
            VectorType result(mRowsNumber);
            for(IndexType i = 0 ; i < mRowsNumber ; ++i)
            result[i] = std::inner_product(mData[i].begin(),mData[i].end(),rVector.begin(),TDataType());
            return result;
        }
    
        Matrix<TDataType>& operator=(const MatrixType& rOtherMatrix)
        {
            mData=rOtherMatrix.mData;
            mRowsNumber = rOtherMatrix.mRowsNumber;
            mColumnsNumber = rOtherMatrix.mColumnsNumber;
            return *this;
        }
    
        Matrix<TDataType>& operator=(const TDataType& rValue)
        {
            for(int i = 0 ; i < mRowsNumber ; i++)
            std::fill(mData[i].begin(), mData[i].end(), rValue);
            return *this;
        }
    
        void Resize(IndexType NewRowsNumber, IndexType NewColumnsNumber)
        {
            if( (mRowsNumber == NewRowsNumber) && (mColumnsNumber == NewColumnsNumber) )
            return;
    
            mData.resize(NewRowsNumber);
            for(typename ContainerType::iterator i = mData.begin() ; i != mData.end() ; ++i)
            i->resize(NewColumnsNumber);
            mRowsNumber    = NewRowsNumber;
            mColumnsNumber = NewColumnsNumber;
        }
    
        void Resize(IndexType NewRowsNumber, IndexType NewColumnsNumber, const TDataType& rValue)
        {
            Resize(NewRowsNumber, NewColumnsNumber);
            *this = rValue;
        }
    
        void FreeData()
        {
            mData.clear();
            std::vector<VectorType>().swap(mData);
        }
    
        void SqMatrixIni(int size)
        {
            mData.clear();
            mData.resize(size);
    
            typename ContainerType::iterator it;
    
            for(it=mData.begin(); it!=mData.end(); ++it) it->resize(size);
    
            mRowsNumber    = size;
            mColumnsNumber = size;
        }
    
        ~Matrix()
        {
        }
    
        Matrix<TDataType>& Transpose()
        {
            if(mRowsNumber != mColumnsNumber)
            {
                Matrix<TDataType> temp_matrix(*this);
                Resize(mColumnsNumber, mRowsNumber);
                for(IndexType i = 0 ; i < mRowsNumber ; i++)
                for(IndexType j = 0 ; j < mColumnsNumber ; j++)
                mData[i][j] = temp_matrix(j,i);
            }
            else
            {
                for(IndexType i = 0 ; i < mRowsNumber ; i++)
                {
                    for(IndexType j = i + 1 ; j < mColumnsNumber ; j++)
                    {
                        TDataType temp = mData[i][j];
                        mData[i][j] = mData[j][i];
                        mData[j][i] = temp;
                    }
                }
            }
    
            return *this;
        }
    
        IndexType RowsNumber() const
        {
            return mRowsNumber;
        }
    
        IndexType ColumnsNumber() const
        {
            return mColumnsNumber;
        }
    
      private:
    
        IndexType mRowsNumber;
        IndexType mColumnsNumber;
        IndexType mRowsIndexBase;
        IndexType mColumnsIndexBase;
    
        std::vector<VectorType> mData;
    };
    
    // input stream function
    template<class TDataType> inline std::istream&
    operator >> (std::istream& IStream, Matrix<TDataType>& rThisMatrix)
    {
        for (typename Matrix<TDataType>::IndexType i=0; i < rThisMatrix.RowsNumber(); i++)
          for (typename Matrix<TDataType>::IndexType j=0; j < rThisMatrix.ColumnsNumber(); j++)
             IStream >> rThisMatrix(i,j);
        return IStream;
    }
    
    // output stream function
    template<class TDataType> inline std::ostream&
    operator << (std::ostream& OStream, const Matrix<TDataType>& rThisMatrix)
    {
        OStream << "Matrix( " << rThisMatrix.RowsNumber() << " * " << rThisMatrix.ColumnsNumber() << " ) :" << std::endl;
        OStream.setf(std::ios::scientific, std::ios::floatfield);
        OStream.width(5);
    
        for (typename Matrix<TDataType>::IndexType i=0; i < rThisMatrix.RowsNumber(); i++)
        {
            std::copy(rThisMatrix[i].begin(), rThisMatrix[i].end(), std::ostream_iterator<TDataType>(OStream, "\t"));
            OStream << std::endl;
        }
    
        OStream.width();
        OStream.unsetf(std::ios::scientific);
        OStream.unsetf(std::ios::floatfield);
    
        return OStream;
    }
}

#endif /* DEFINE_KRATOS_MATRIX  defined */

