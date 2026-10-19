
#if !defined(KRATOS_MODEL_GEOMETRY_CONTAINER)
#define KRATOS_MODEL_GEOMETRY_CONTAINER

#include "../includes/nodes_container.h"
#include "../includes/points_container.h"

namespace Kratos
{
    class ModelGeometryContainer
    {
      public:
     
        typedef std::shared_ptr<ModelGeometryContainer> Pointer            ;
        typedef Vector<Node::Pointer>                     NodesArrayType     ;
        typedef unsigned int                              IndexType          ;
        typedef std::map<IndexType, Properties::Pointer>  PropertiesArrayType;
        typedef std::vector<Sources::Pointer>             SourcesArrayType   ;
        
        ModelGeometryContainer()
        {
        }
        
        virtual ~ModelGeometryContainer()
        {
        }
             
        int NumberOfNodes()
        {
            return mNodes.NumberOfNodes();
        }
        
        Node::Pointer GetNode( int NodeId )
        {
            return mNodes.GetNode( NodeId );
        }
    
        Point<double>::Pointer GetPoint( const int PointIndex )
        {
            return mPoints.GetPoint( PointIndex );
        }
        
        void AddNode( Node::Pointer pNewNode )
        {
            pNewNode->SetData( mPoints.AddPoint( pNewNode ) );
            mNodes.AddNode( pNewNode );
        }
        
        NodesArrayType& GetNodesArray()
        {
            return mNodes.GetNodesArray();
        }
	    
        Properties::Pointer GetProperties( int PropertiesId )
        {
            return mPropertiesList[PropertiesId];
        }
        
        void AddProperties( IndexType PropertiesId, Properties::Pointer pNewProperties )
        {
            mPropertiesList[PropertiesId] = pNewProperties;
        }
        
        void AddSources( Sources::Pointer pNewSources )
        {
            mSourcesList.push_back( pNewSources );
        }
     
      private:
                
        NodesContainer mNodes;
        
        PointsContainer mPoints;
        
        PropertiesArrayType mPropertiesList;
        
        SourcesArrayType mSourcesList;                        
    }; 
} 

#endif 

