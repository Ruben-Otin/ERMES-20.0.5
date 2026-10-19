
#if !defined(KRATOS_POINT_DATA)
#define KRATOS_POINT_DATA

#include <memory>
#include <memory>

#include "../includes/point.h"
#include "../includes/variables_container.h"

namespace Kratos
{
    class Element;

    class Node;

    class PointData
    {
      public:
       
        typedef std::shared_ptr<PointData> Pointer;
        
        PointData( const Point<double>::Pointer& pThisPoint, unsigned int NewIndex );
       
        virtual ~PointData();
        
        Point<double>::Pointer& pPoint()
        {
            return mpPoint;
        }
        
        unsigned int Index()
        {
            return mIndex;
        }
           
      private:
          
        Point<double>::Pointer mpPoint;
        
        unsigned int mIndex;
        
        PointData();
        
        PointData( const PointData& rOtherPointData );
        
        const PointData& operator=( const PointData& rOtherPointData );                    
    }; 
}  

#endif 

