/* *********************************************************   
 *          
 *   Last modified by:    
 *   Date:                $Date: 2003/06/04 11:38:33 $
 *   Revision:            $Revision: 1.21 $
 *
 * ***********************************************************/



#include <iostream>

#include "../includes/model.h"


namespace Kratos
{

  void Model::AddNode(Node::Pointer pNewNode)
  {    
      mpSolutionStepsData->CurrentSolutionStepData()->AddNode(pNewNode);   
  }

  void Model::AddProperties(IndexType PropertiesId, Properties::Pointer pNewProperties)
  {
	  mpSolutionStepsData->CurrentSolutionStepData()->AddProperties(PropertiesId, pNewProperties);
  }
  
  void Model::AddSources(Sources::Pointer pNewSources)
  {
	  mpSolutionStepsData->CurrentSolutionStepData()->AddSources(pNewSources);
  }
  
  void Model::AddElement(const String& rElementGroupName, Element::Pointer pNewElement)
  {
      mpSolutionStepsData->CurrentSolutionStepData()->AddElement(rElementGroupName, pNewElement);
  }
  
  
  PElementArray& Model::GetElements(const String& rElementGroupName)
  {
      return mpSolutionStepsData->CurrentSolutionStepData()->GetElements(rElementGroupName);
  }
  

  Model::NodesArrayType& Model::GetNodesArray()
  {
      return mpSolutionStepsData->CurrentSolutionStepData()->GetNodesArray();  
  }

  Properties::Pointer Model::GetProperties(int PropertyId)
  {
      return mpSolutionStepsData->CurrentSolutionStepData()->GetProperties(PropertyId);
  }

  
  Node::Pointer  Model::GetNode(int NodeId)
  {
      return mpSolutionStepsData->CurrentSolutionStepData()->GetNode(NodeId);
  }
  

  Element::Pointer  Model::GetElement(int ElementId)
  {
    KRATOS_TRY
      return mpSolutionStepsData->CurrentSolutionStepData()->GetElement(ElementId);
    KRATOS_CATCH("Model::GetElement(int ElementId)", "")
  }


  Model::ElementsArrayType Model::GetAllElements()
  { 
      return mpSolutionStepsData->CurrentSolutionStepData()->GetAllElements();
  }


  void Model::AddProcess(Process::Pointer pProcess)
  {
      mProcesses.push_back(pProcess);
  }

  
  Process::ProcessArrayType& Model::GetProcessArray()
  {
      return mProcesses;
  }


  void Model::DeleteElementGroup(const String& rElementGroupName)
  {
      mpSolutionStepsData->CurrentSolutionStepData()->DeleteElementGroup(rElementGroupName);
  }
  
  
  unsigned int Model::GetSolutionStep()
  {
      return mSolutionStep;
  }

  void Model::IncreaseSolutionStep()
  {
      mSolutionStep++;   
  }

  SolutionStepData::Pointer Model::CreateNewTimeStep(double NewTime)
  {
      return mpTimeStepsData->CreateNewTimeStepData(NewTime, CreateWithCurrentSolutionStepData());
  }

  SolutionStepData::Pointer Model::CreateNewSolutionStep(void)
  {
	  mSolutionStep++;  
      return mpSolutionStepsData->CreateNewSolutionStepData();
  }

  SolutionStepData::Pointer Model::CreateWithCurrentSolutionStepData(void)
  {
      mSolutionStep++;  
      return mpSolutionStepsData->CreateWithCurrentSolutionStepData();
  }
  

} // Namespace Kratos
