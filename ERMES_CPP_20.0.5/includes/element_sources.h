
#if !defined(KRATOS_ELEMENT_SOURCES)
#define KRATOS_ELEMENT_SOURCES

#include "../includes/sources.h"

namespace Kratos
{
    enum SourceTarget
    {
        ST_VOLUME,
        
	    ST_SURFACE_1, ST_SURFACE_2, ST_SURFACE_3, ST_SURFACE_4, ST_SURFACE_5, ST_SURFACE_6,
        
	    ST_LINE_1, ST_LINE_2, ST_LINE_3, ST_LINE_4 , ST_LINE_5 , ST_LINE_6 ,
	    ST_LINE_7, ST_LINE_8, ST_LINE_9, ST_LINE_10, ST_LINE_11, ST_LINE_12,
        
	    ST_POINT_1 , ST_POINT_2 , ST_POINT_3 , ST_POINT_4 , ST_POINT_5 , ST_POINT_6 ,
	    ST_POINT_7 , ST_POINT_8 , ST_POINT_9 , ST_POINT_10, ST_POINT_11, ST_POINT_12,
	    ST_POINT_13, ST_POINT_14, ST_POINT_15, ST_POINT_16, ST_POINT_17, ST_POINT_18,
	    ST_POINT_19, ST_POINT_20, ST_POINT_21, ST_POINT_22, ST_POINT_23, ST_POINT_24
    }; 
    
    class ElementSources
    {
      public:
        
        typedef std::shared_ptr<ElementSources> Pointer         ;

        typedef std::vector<Sources::Pointer>     SourcesArrayType;
           
        ElementSources()
        {
        }
    
        virtual ~ElementSources()
        {
        }
        
        std::vector<Sources::Pointer> operator()( const SourceTarget& rTarget )
        {
            if( HasSourceFor( rTarget ) ) 
            {
                return mData[rTarget];
            }
            else
            {
                return std::vector<Sources::Pointer>();
            }
        }
        
        void AddSources( const SourceTarget& rTarget, Sources::Pointer pNewSources )
        {
            mData[rTarget].push_back( pNewSources );
        }

        bool HasSourceFor( const SourceTarget& rTarget )
        {
            return ( mData.find( rTarget ) != mData.end() );
        }
 
      private:
            
        std::map<SourceTarget, std::vector<Sources::Pointer> > mData;
    }; 
}  

#endif 

