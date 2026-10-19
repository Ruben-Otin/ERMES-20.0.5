
#if !defined(KRATOS_POINT)
#define KRATOS_POINT

#include <valarray>
#include <cmath>

#include <memory>
#include "../includes/exception.h"

namespace Kratos
{
    const int POINT_ZERO_BASE_INDEX = 0;
    const int POINT_ONE_BASE_INDEX  = 1;

    template<class TDataType>
    class Point : public std::valarray<TDataType>
    {
      public:

        typedef std::shared_ptr<Point<double> > Pointer                      ;
        typedef std::valarray<TDataType>          CoordinatesArrayType         ;
        typedef std::indirect_array<TDataType>    CoordinatesArrayReferenceType;

        // Constructor using coordinates of a point. 
        // It could be used as default constructor and it make all the coordinate equal to zero.
        Point( TDataType Xval=0, TDataType Yval=0, TDataType Zval=0 )
        {
            this->resize(3);
            this->operator[](0) = Xval;
            this->operator[](1) = Yval;
            this->operator[](2) = Zval;
        }

        // Copy constructor. Initialize this point with the coordinates of given point.
        Point( const Point<TDataType>& rOtherPoint ) : std::valarray<TDataType>( rOtherPoint )
        {
        }

        // Constructor using coordinates stored in given array. 
        // Initialize this point with the coordinates in the array. 
        Point( const CoordinatesArrayType& rOtherCoordinates ) : std::valarray<TDataType>( rOtherCoordinates )
        {
        }

        virtual ~Point()
        {
        }

        // Returns X coordinate 
        TDataType X() const
        {
			return this->operator[](0);
        }

        // Returns Y coordinate 
        TDataType Y() const
        {
			return this->operator[](1);
        }

        // Returns Z coordinate 
        TDataType Z() const
        {
	        return this->operator[](2);
        }

        // Returns X coordinate 
        TDataType& X()
        {
			return this->operator[](0);
        }

        // Returns Y coordinate 
        TDataType& Y()
        {
			return this->operator[](1);
        }

        // Returns Z coordinate 
        TDataType& Z()
        {
	        return this->operator[](2);
        }

        /** This is an access method to point's coordinate by indices and with
        specified index base. For example with POINT_ONE_BASE_INDEX this
        function return x, y and z coordinate whith 1, 2 and 3 as input
        respectively and with POINT_ZERO_BASE_INDEX return x, y and z with 0,
        1 and 2 as input.
        NOTE: For efficency it's recommended to use direct access
        methods (X, Y,...) instead of these function. It throw an exception
        if index is not in the range*/
        TDataType Coordinate( int CoordinateIndex, int IndexBase = POINT_ONE_BASE_INDEX ) const
        {
            return this->operator[]( CoordinateIndex - IndexBase );
        }

        /** This is an access method to get a reference to point's coordinate by
        indices and with specified index base. For example with POINT_ONE_BASE_INDEX
        this function return references to x, y and z coordinate whith 1, 2
        and 3 as input respectively and with POINT_ZERO_BASE_INDEX return references
        to x, y and z with 0, 1 and 2 as input.
        NOTE: For efficency it's recommended to use direct access
        methods (X, Y,...) instead of these function. It throw an exception
        if index is not in the range */
        TDataType& rCoordinate( int CoordinateIndex, int IndexBase = POINT_ONE_BASE_INDEX )
        {
            return this->operator[]( CoordinateIndex - IndexBase );
        }

		unsigned int Dimension()
        {
            return this->size();
        }

        /** This operator will compare the length of the vector wich
        conected to these points. In the other word it calculate the sum of 
        the squares of the coordinates for each point and compare them. */
        friend bool operator<( const Point<TDataType>& rFirstPoint, const Point<TDataType>& rSecondPoint )
        {
            double d1 = 0;
            double d2 = 0;

            int i;
            int first_size  = rFirstPoint.size ();
            int second_size = rSecondPoint.size();

            for( i = 0 ; i < first_size ; i++ )
            {
                d1 += rFirstPoint[i] * rFirstPoint[i];
            }

            for( i = 0 ; i < second_size ; i++ )
            {
                d2 += rSecondPoint[i] * rSecondPoint[i];
            }

            return d1 < d2;
        }

        /** This operator will compare the length of the vector wich
        conected to these points. In the other word it calculate the sum of 
        the squares of the coordinates for each point and compare them. */
        friend bool operator>( const Point<TDataType>& rFirstPoint, const Point<TDataType>& rSecondPoint )
        {
            TDataType d1 = 0;
            TDataType d2 = 0;

            int i;
            int first_size  = rFirstPoint.size ();
            int second_size = rSecondPoint.size();

            for( i = 0 ; i <  first_size ; i++ )
            {
                d1 += rFirstPoint[i] * rFirstPoint[i];
            }

            for( i = 0 ; i < second_size ; i++ )
            {
                d2 += rSecondPoint[i] * rSecondPoint[i];
            }

            return d1 > d2;
        }
    }; 
}  

#endif 

