
#if !defined(KRATOS_NODES_CONTAINER)
#define KRATOS_NODES_CONTAINER

#include "../includes/node.h"

namespace Kratos
{
    class NodesContainer
    {
      public:
           
        typedef Vector<Node::Pointer> NodesArrayType;
        
        NodesContainer()
        {
        }
        
        virtual ~NodesContainer()
        {
        }
        
        int NumberOfNodes()
        {
            return mData.size();
        }
               
        void AddNode( const Node::Pointer pNewNode )
        {
            mData.push_back( pNewNode );
        }
        
        Node::Pointer GetNode( int Id )
        {
            return mData[Id];
        }
        
        NodesArrayType& GetNodesArray()
        {
            return mData;
        }
      
      private:
              
        NodesArrayType mData;        
    }; 
} 

#endif 

