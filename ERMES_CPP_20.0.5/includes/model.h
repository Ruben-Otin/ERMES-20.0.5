
#if !defined(KRATOS_MODEL)
#define KRATOS_MODEL

#include <set>

#include <memory>

#include "../includes/element.h"
#include "../includes/node.h"
#include "../includes/process.h"
#include "../includes/exception.h"
#include "../includes/properties.h"

#include "../includes/points_container.h"
#include "../includes/nodes_container.h"
#include "../includes/time_steps_container.h"
#include "../includes/solution_steps_data_container.h"
#include "../includes/model_properties.h"
#include "../includes/process_info.h"
#include "../includes/kratos_vectorial_spaces.h"
#include "../linear_solvers/kratos_linear_solvers.h"

namespace Kratos
{
    class Model
    {
      public:
       
        typedef std::shared_ptr<Model> Pointer;

        typedef Vector<Node::Pointer> NodesArrayType;

        typedef unsigned int IndexType;
        
        Model(void) : mpSolutionStepsData( new SolutionStepsDataContainer( 3 ) ), 
                      mpTimeStepsData    ( new TimeStepsContainer( mpSolutionStepsData->CurrentSolutionStepData(), 3 ) ), 
                      mSolutionStep( 0   ), 
                      mTime        ( 0.0 ),
                      mDeltaTime   ( 0.0 )
        {
        }
        
        virtual ~Model()
        {
        }
        
        template<class TDataType> 
        TDataType& operator()( const Variable<TDataType>& rV, const Node& rThisNode, double ThisTime )
        {
            return mpTimeStepsData->GetNearestStep( ThisTime )->pVariables( rThisNode.Index() )->GetValue( rV );
        }
    
        template<class TDataType> 
        TDataType& operator()( const Variable<TDataType>& rV, const Node& rThisNode )
        {
            return mpSolutionStepsData->CurrentSolutionStepData()->pVariables( rThisNode.Index() )->GetValue( rV );
        }
       
        operator std::shared_ptr<PropertiesDefaultFunction<double> >()
        {
            return PropertiesDefaultFunction<double>::Pointer( new ModelProperties<double>(*this) );
        }
        
        operator std::shared_ptr<PropertiesDefaultFunction< std::vector<double> > >()
        {
            return PropertiesDefaultFunction< std::vector<double> >::Pointer( new ModelProperties<std::vector<double> >(*this) );
        }
    
        operator std::shared_ptr<PropertiesDefaultFunction<std::complex<double> > >()
        {
            return PropertiesDefaultFunction< std::complex<double> >::Pointer( new ModelProperties< std::complex<double> >(*this) );
        }
    
        operator std::shared_ptr<PropertiesDefaultFunction<Vector<double> > >()
        {
            return PropertiesDefaultFunction<Vector<double> >::Pointer( new ModelProperties<Vector<double> >(*this) );
        }
        
        operator std::shared_ptr<PropertiesDefaultFunction< Matrix<double> > >()
        {
            return PropertiesDefaultFunction< Matrix<double> >::Pointer( new ModelProperties< Matrix<double> >(*this) );
        }
        
        template<class TDataType>
        TDataType& Value( const Variable<TDataType>& rVariable )
        {
            return mGlobalVariables.GetValue( rVariable );
        }
    
        template<class TDataType>
        TDataType& Value( const Variable<TDataType>& rVariable, const Dof::Pointer pThisDof )
        {
            return mpSolutionStepsData->CurrentSolutionStepData()->pVariables( pThisDof->NodeIndex() )->GetValue( rVariable );
        }
    
        template<class TDataType> 
        TDataType& Value( const Variable<TDataType>& rV, const Node& rThisNode )
        {
            return mpSolutionStepsData->CurrentSolutionStepData()->pVariables( rThisNode.Index() )->GetValue( rV );
        }
    
	    double CurrentTime()
        {
            return mTime;
        } 
    
	    void SetCurrentTime( double NewTime )
        {
            mDeltaTime = NewTime - mTime;
            mTime      = NewTime;   
        }
    
        double CurrentDeltaTime()
        {
            return mDeltaTime;
        }
    
	    void SetDeltaTime( double NewDeltaTime )
        {
            mDeltaTime = NewDeltaTime;
        }
        
	    void SetCurrentProcessInfo() 
        {
            mCurrentProcessInfo.SetCurrentTime        ( CurrentTime     () );
	  	    mCurrentProcessInfo.SetDeltaTime          ( CurrentDeltaTime() );
	  	    mCurrentProcessInfo.SetCurrentSolutionStep( GetSolutionStep () );
        }

	    void SetNumPrintableNodes( unsigned int NumPrintableNodes )
        {
            mNumPrintableNodes = NumPrintableNodes;
        }

	    ProcessInfo& GetCurrentProcessInfo()	
        {
            return mCurrentProcessInfo;
        } 

        Point<double>::Pointer GetPoint( const int PointIndex )
	    {
	       return mpSolutionStepsData->CurrentSolutionStepData()->GetPoint( PointIndex );
	    }

        Node::Pointer GetNode( int NodeId );

        Properties::Pointer GetProperties( int PropertyId );
        
        void AddNode( Node::Pointer pNewNode );
    
	    void AddProperties( IndexType PropertiesId, Properties::Pointer pNewProperties );
    
	    void AddSources( Sources::Pointer pNewSources );
    
        void IncreaseSolutionStep();
        
        void AddProcess( Process::Pointer pProcess );
        
        unsigned int GetSolutionStep();

        unsigned int GetNodesArraySize();
      
        NodesArrayType& GetNodesArray();

        NodesArrayType& GetPrintableNodesArray();
        
        Process::ProcessArrayType& GetProcessArray();
    
        SolutionStepData::Pointer CreateNewTimeStep( double NewTime ); 

        SolutionStepData::Pointer CreateNewSolutionStep( void );

        SolutionStepData::Pointer CreateWithCurrentSolutionStepData( void );
       
      private:
       
        SolutionStepsDataContainer::Pointer mpSolutionStepsData;
        
        TimeStepsContainer::Pointer mpTimeStepsData;
    
	    Process::ProcessArrayType mProcesses;
             
        VariablesContainer mGlobalVariables;
    
        ProcessInfo mCurrentProcessInfo;
   
        double mTime;
    
	    double mDeltaTime;
        
        unsigned int mSolutionStep; 

        unsigned int mNumPrintableNodes;
    };
}  

#endif 

