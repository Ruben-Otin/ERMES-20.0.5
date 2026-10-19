
#if !defined(KRATOS_SOLUTION_STEP_DATA)
#define KRATOS_SOLUTION_STEP_DATA

#include <memory>

#include "../includes/define.h"
#include "../includes/model_geometry_container.h"
#include "../includes/variables_container.h"

namespace Kratos
{
    class SolutionStepData
    {
      public:
            
	    typedef std::shared_ptr<SolutionStepData> Pointer;
        
	    typedef unsigned int IndexType;
        
        SolutionStepData() : mStepNumber( int() ), mpModelGeometry( new ModelGeometryContainer() )
        {
             String Buffer = "Solution step data created for Solution = "; 
             Buffer << mStepNumber;
             KRATOS_TRACE( "SolutionStepData::SolutionStepData()", Buffer );
        }
        
	    SolutionStepData( SolutionStepData& rOther ) : mStepNumber( rOther.StepNumber() ), mpModelGeometry( rOther.pModelGeometry() )
	    {
             std::vector<VariablesContainer::Pointer>::iterator i;
             
             for( i = rOther.mData.begin(); i != rOther.mData.end(); i++ )
             {
	             mData.push_back( VariablesContainer::Pointer( new VariablesContainer(**i) ) );
             }
        }
        
        SolutionStepData( ModelGeometryContainer::Pointer pNewModelGeometry, int NewStepNumber = int() ) : mStepNumber( NewStepNumber ), mpModelGeometry( pNewModelGeometry )
        {
            String Buffer = "Solution step data created for Solution = "; 
            Buffer << mStepNumber;
            KRATOS_TRACE( "SolutionStepData::SolutionStepData()", Buffer );
        }
        
        virtual ~SolutionStepData()
        {
            String Buffer = "Solution step data deleted for Solution = "; 
            Buffer << mStepNumber;
            KRATOS_TRACE( "SolutionStepData::SolutionStepData()", Buffer );
        }
        
        int NumberOfNodes()
        {
            return mpModelGeometry->NumberOfNodes();
        }
        
        int StepNumber( void ) const
        {
            return mStepNumber;
        }
        
        void StepNumber( int NewStepNumber )
        {
            mStepNumber = NewStepNumber;
        }
        
	    ModelGeometryContainer::Pointer pModelGeometry()
	    {
            return mpModelGeometry;
        }
        
        Node::Pointer GetNode( int NodeId )
        {
            return mpModelGeometry->GetNode( NodeId );
        }
        
        Point<double>::Pointer GetPoint( int PointIndex )
        {
            return mpModelGeometry->GetPoint( PointIndex );
        }
        
        void AddNode( Node::Pointer pNewNode )
        {
            mpModelGeometry->AddNode( pNewNode );

            mData.push_back( VariablesContainer::Pointer( new VariablesContainer() ) );
        }
        
        Vector<Node::Pointer>& GetNodesArray()
        {
            return mpModelGeometry->GetNodesArray();
        }
	    
        Properties::Pointer GetProperties( int PropertiesId )
        {
            return mpModelGeometry->GetProperties( PropertiesId );
        }
        
        void AddSources( Sources::Pointer pNewSources )
        {
           mpModelGeometry->AddSources( pNewSources );
        }
        
        void AddProperties( IndexType PropertiesId, Properties::Pointer pNewProperties )
        {
           mpModelGeometry->AddProperties( PropertiesId, pNewProperties );
        }
        
        VariablesContainer::Pointer pVariables( const IndexType Index )
        {
            return mData[Index];
        }
       
      private:
                
        int mStepNumber;
         
        ModelGeometryContainer::Pointer mpModelGeometry;
         
        std::vector<VariablesContainer::Pointer> mData;             
    }; 
}  

#endif 

