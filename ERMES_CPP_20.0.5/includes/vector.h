
#if !defined(KRATOS_VECTOR)
#define KRATOS_VECTOR

#include <ios>
#include <iterator>
#include <vector>
#include <algorithm>
#include <functional>

namespace Kratos
{
	template<class TDataType>
	class Vector : public std::vector<TDataType>
    {
      public:

	    typedef std::vector<TDataType> BaseType  ;

		typedef Vector<TDataType>      VectorType;

		Vector()
        {
        }

		Vector( int NewSize, const TDataType& NewValue = TDataType() ) : BaseType( NewSize, NewValue )
        {
        }

		Vector( const BaseType& OtherVector ) : BaseType( OtherVector )
        {
        }

		Vector( const VectorType& OtherVector ) : BaseType( OtherVector )
        {
        }

		virtual ~Vector()
        {
        }

        void FreeData()
        {
            this->clear();
            std::vector<TDataType>().swap( *this );
        }

		VectorType& operator=( const TDataType& OtherValue )
		{
			std::fill( this->begin(), this->end(), OtherValue );
			return *this;
		}

		VectorType& operator=( const VectorType& OtherVector )
		{
			BaseType::operator = ( OtherVector );
			return *this;
		}

		VectorType& operator=( const std::vector<TDataType>& OtherVector )
		{
			BaseType::operator = ( OtherVector );
			return *this;
		}

        VectorType& operator+=( const VectorType& OtherVector )
        {
            std::transform( this->begin(), this->end(), OtherVector.begin(), this->begin(), std::plus<TDataType>() );
            return *this;
        }

        VectorType& operator-=( const VectorType& OtherVector )
        {
            std::transform( this->begin(), this->end(), OtherVector.begin(), this->begin(), std::minus<TDataType>() );
            return *this;
        }

        VectorType& operator+=( const TDataType& OtherValue )
        {
            // ( std::bind2nd was removed in C++17 and is not available in Visual C++ 2022 )
            for( typename BaseType::iterator it = this->begin(); it != this->end(); ++it ) *it += OtherValue;
            return *this;
        }

        VectorType& operator-=( const TDataType& OtherValue )
        {
            // ( std::bind2nd was removed in C++17 and is not available in Visual C++ 2022 )
            for( typename BaseType::iterator it = this->begin(); it != this->end(); ++it ) *it -= OtherValue;
            return *this;
        }

        VectorType& operator*=( const TDataType& OtherValue )
        {
            // ( std::bind2nd was removed in C++17 and is not available in Visual C++ 2022 )
            for( typename BaseType::iterator it = this->begin(); it != this->end(); ++it ) *it *= OtherValue;
            return *this;
        }

        VectorType& operator/=( const TDataType& OtherValue )
        {
            // ( std::bind2nd was removed in C++17 and is not available in Visual C++ 2022 )
            for( typename BaseType::iterator it = this->begin(); it != this->end(); ++it ) *it /= OtherValue;
            return *this;
        }
    }; 

    template<class TDataType>
    TDataType operator*( const Vector<TDataType>& FirstVector, const Vector<TDataType>& SecondVector )
    {
        return inner_product( FirstVector.begin(), FirstVector.end(), SecondVector.begin(), TDataType() );
    }

    // output stream function
    template<class TStreamType, class TDataType> inline TStreamType&
    operator << ( TStreamType& OStream, const Vector<TDataType>& rThisVector )
    {
        OStream.setf ( std::ios::scientific, std::ios::floatfield );
	    OStream.width( 5 );

	    std::copy( rThisVector.begin(), rThisVector.end(), std::ostream_iterator<TDataType>( OStream, "\t" ) );

	    OStream.width (                      );
        OStream.unsetf( std::ios::scientific );
	    OStream.unsetf( std::ios::floatfield );

        return OStream;
    }
} 

#endif 

