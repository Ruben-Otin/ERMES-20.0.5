
#include <functional>

#include "../includes/define.h"
#include "../includes/element.h"

namespace Kratos
{
    PointData::PointData( const Point<double>::Pointer& pThisPoint, unsigned int NewIndex ) : mpPoint( pThisPoint ), mIndex( NewIndex )
    {
    }

    PointData::~PointData()
    {    
    }
}
