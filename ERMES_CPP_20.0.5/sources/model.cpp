
#include <iostream>
#include "../includes/model.h"

namespace Kratos
{
    void Model::AddNode( Node::Pointer pNewNode )
    {    
        mpSolutionStepsData->CurrentSolutionStepData()->AddNode( pNewNode );   
    }
    
    void Model::AddProperties( IndexType PropertiesId, Properties::Pointer pNewProperties )
    {
	    mpSolutionStepsData->CurrentSolutionStepData()->AddProperties( PropertiesId, pNewProperties );
    }
    
    void Model::AddSources( Sources::Pointer pNewSources )
    {
	    mpSolutionStepsData->CurrentSolutionStepData()->AddSources( pNewSources );
    }
    
    Model::NodesArrayType& Model::GetNodesArray()
    {
        return mpSolutionStepsData->CurrentSolutionStepData()->GetNodesArray();  
    }

    Model::NodesArrayType& Model::GetPrintableNodesArray()
    {
        NodesArrayType& PrintableNodes = mpSolutionStepsData->CurrentSolutionStepData()->GetNodesArray();

        PrintableNodes.resize( mNumPrintableNodes );
        
        return PrintableNodes;  
    }
    
    Properties::Pointer Model::GetProperties( int PropertyId )
    {
        return mpSolutionStepsData->CurrentSolutionStepData()->GetProperties( PropertyId );
    }
    
    Node::Pointer Model::GetNode( int NodeId )
    {
        return mpSolutionStepsData->CurrentSolutionStepData()->GetNode( NodeId );
    }
    
    unsigned int Model::GetSolutionStep()
    {
        return mSolutionStep;
    }

    unsigned int Model::GetNodesArraySize()
    {
        return ( mpSolutionStepsData->CurrentSolutionStepData()->GetNodesArray() ).size();  
    }
    
    void Model::IncreaseSolutionStep()
    {
        mSolutionStep++;   
    }

    void Model::AddProcess( Process::Pointer pProcess )
    {
        mProcesses.push_back( pProcess );
    }

    Process::ProcessArrayType& Model::GetProcessArray()
    {
        return mProcesses;
    }

    SolutionStepData::Pointer Model::CreateNewTimeStep( double NewTime )
    {
        return mpTimeStepsData->CreateNewTimeStepData( NewTime, CreateWithCurrentSolutionStepData() );
    }
    
    SolutionStepData::Pointer Model::CreateNewSolutionStep( void )
    {
	    mSolutionStep++;  
        return mpSolutionStepsData->CreateNewSolutionStepData();
    }
    
    SolutionStepData::Pointer Model::CreateWithCurrentSolutionStepData( void )
    {
        mSolutionStep++;  
        return mpSolutionStepsData->CreateWithCurrentSolutionStepData();
    }
}; 
