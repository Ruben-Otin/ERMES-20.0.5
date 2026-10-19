
#if !defined(KRATOS_POINTS_CONTAINER)
#define KRATOS_POINTS_CONTAINER

#include "../includes/point_data.h"

namespace Kratos
{
    class PointsContainer
    {
      public:
       
        typedef std::vector<PointData::Pointer> ContainerType      ;

        typedef Vector<Point<double>::Pointer>  PointsContainerType;
        
        PointsContainer()
        {
        }
        
        virtual ~PointsContainer()
        {
        }
        
        PointData::Pointer AddPoint( const Point<double>::Pointer pNewPoint )
        {
            PointData::Pointer p_new_point_data( new PointData( pNewPoint, mData.size() ) );
        
            mData.push_back( p_new_point_data );
        
            mPoints.push_back( pNewPoint ); 
        
            return p_new_point_data;
        }
        
	    Point<double>::Pointer GetPoint( const int PointIndex ) 
        {
            return mPoints[PointIndex];
        }
      
      private:
     
        ContainerType mData;

    	PointsContainerType mPoints;     
    }; 
}  

#endif 

