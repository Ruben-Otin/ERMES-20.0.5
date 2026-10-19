

#if !defined(KRATOS_MODEL )
#define  KRATOS_MODEL

#include <set>

#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"

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
#include "../vectorial_spaces/kratos_vectorial_spaces.h"
#include "../linear_solvers/kratos_linear_solvers.h"


namespace Kratos
{

  class Model
  {
    public:
     
      /** A counted pointer of Model */
      typedef boost::shared_ptr<Model> Pointer;
      
      /** Type definition of an array of nodes in model */
      typedef Vector<Node::Pointer> NodesArrayType;

      /** Type definition of an array of elements in model */
      typedef std::vector<Element::Pointer> ElementsArrayType;
     
      typedef unsigned int IndexType;
      



      Model(void) : mpSolutionStepsData(new SolutionStepsDataContainer(3)), 
                    mpTimeStepsData(new TimeStepsContainer(mpSolutionStepsData->CurrentSolutionStepData(), 3)), 
                    mSolutionStep(0), 
                    mTime(0.00),
                    mDeltaTime(0.00)
      {
      }
      

      virtual ~Model()
      {
      }
      


      /** This method interpolate value of a variable for given node for given time*/
      template<class TDataType> 
      TDataType& operator()(const Variable<TDataType>& rV, const Node& rThisNode, double ThisTime)
      {
          return mpTimeStepsData->GetNearestStep(ThisTime)->pVariables(rThisNode.Index())->GetValue(rV);
      }

      template<class TDataType> 
      TDataType& operator()(const Variable<TDataType>& rV, const Node& rThisNode)
      {
          return mpSolutionStepsData->CurrentSolutionStepData()->pVariables(rThisNode.Index())->GetValue(rV);
      }
     


      operator boost::shared_ptr<PropertiesDefaultFunction<double> >()
      {
          return PropertiesDefaultFunction<double>::Pointer(new ModelProperties<double>(*this));
      }
      

      operator boost::shared_ptr<PropertiesDefaultFunction<std::vector<double> > >()
      {
          return PropertiesDefaultFunction<std::vector<double> >::Pointer(new ModelProperties<std::vector<double> >(*this));
      }

      operator boost::shared_ptr<PropertiesDefaultFunction<std::complex<double> > >()
      {
          return PropertiesDefaultFunction<std::complex<double> >::Pointer(new ModelProperties<std::complex<double> >(*this));
      }
      

      operator boost::shared_ptr<PropertiesDefaultFunction<Vector<double> > >()
      {
          return PropertiesDefaultFunction<Vector<double> >::Pointer(new ModelProperties<Vector<double> >(*this));
      }
      

      operator boost::shared_ptr<PropertiesDefaultFunction<Matrix<double> > >()
      {
          return PropertiesDefaultFunction<Matrix<double> >::Pointer(new ModelProperties<Matrix<double> >(*this));
      }
      



      template<class TDataType>
      TDataType& Value(const Variable<TDataType>& rVariable)
      {
        return mGlobalVariables.GetValue(rVariable);
      }

      template<class TDataType>
      TDataType& Value(const Variable<TDataType>& rVariable, const Dof::Pointer pThisDof)
      {
           return mpSolutionStepsData->CurrentSolutionStepData()->pVariables(pThisDof->NodeIndex())->GetValue(rVariable);
      }

      template<class TDataType> 
      TDataType& Value(const Variable<TDataType>& rV, const Node& rThisNode)
      {
          return mpSolutionStepsData->CurrentSolutionStepData()->pVariables(rThisNode.Index())->GetValue(rV);
      }



	  
	  double CurrentTime()
      {
          return mTime;
      } 

	  void SetCurrentTime(double NewTime)
      {
          mDeltaTime = NewTime-mTime;
          mTime      = NewTime;   
      }

	  /** Get the value of the Current DeltaTime.
	  The Delta Time is understood as the diference interval between
	  two Time Steps.*/
	  double CurrentDeltaTime()
      {
          return mDeltaTime;
      }

	  void SetDeltaTime(double NewDeltaTime)
      {
          mDeltaTime = NewDeltaTime;
      }
      


	  void SetCurrentProcessInfo() 
      {
          mCurrentProcessInfo.SetCurrentTime(CurrentTime());
		  mCurrentProcessInfo.SetDeltaTime(CurrentDeltaTime());
		  mCurrentProcessInfo.SetCurrentSolutionStep(GetSolutionStep());
      }
      
      
	  ProcessInfo& GetCurrentProcessInfo()	
      {
          return mCurrentProcessInfo;
      } 


      
      /** Returns the Node::Pointer  corresponding to it's identifier */
      Node::Pointer GetNode(int NodeId);


      Point<double>::Pointer GetPoint(const int PointIndex)
	  {
	     return mpSolutionStepsData->CurrentSolutionStepData()->GetPoint(PointIndex);
	  }



	  /** Returns the Element::Pointer  corresponding to it's identifier */
	  Element::Pointer GetElement(int ElementId);


	  /** Returns a pointer to property with given id */
	  Properties::Pointer GetProperties(int PropertyId);
      


      /** Inserts a node in the model
	  Adds a new pointer to a node to the Model. All entities 
	  on Model are stored as smart pointers, preventing it's 
	  destruction if are deleted*/
      void AddNode(Node::Pointer pNewNode);



	  /** Add a Property to the model */
	  void AddProperties(IndexType PropertiesId, Properties::Pointer pNewProperties);



	  /** Add a Property to the model */
	  void AddSources(Sources::Pointer pNewSources);



      /** Increases the analysis step.*/
      void IncreaseSolutionStep();
      


      /** Inserts an element to the model 
	  to add a new element to the model is necessary to specify 
	  a String key asociated to a group of elements. This key 
	  must be used to get all the elements identified by this key.*/
      void AddElement(const String& rElementGroupName, Element::Pointer pNewElement);
      
    
      /** Adds a initial value to the model */
      void AddInitialValue(String InitialValueName,
                           int INode,
			               Variable<double> Variable,
			               double Value);    
      
      

      /**Deletes the choosen Element group from model
	  Nothing happens if not exists.*/
      void DeleteElementGroup(const String& rElementGroupName);
      


      /** Return the elements coresponging to ElelemntGroup. 
	  Return A refference to a PElementArray that contains 
	  all the elements identified by  Group. 
	  throws an exception if ElementGroup not exists in the model*/
      ElementsArrayType& GetElements(const String& rElementGroupName);

      
      /** Return a copy of all the elements in Model*/
	  ElementsArrayType GetAllElements();
	  

      /** Returns the array of all nodes stored in Model*/
      NodesArrayType& GetNodesArray();
      
      /** Return the current analysis Step */
      unsigned int GetSolutionStep();
    
          
      /** Adds a process to the model.Process is the key to flexibility in kratos.*/
      void AddProcess(Process::Pointer pProcess);
      
      /** Return a refference to an array of pointers to processes.
	  The object of this array is to */
      Process::ProcessArrayType& GetProcessArray();

      SolutionStepData::Pointer CreateNewTimeStep(double NewTime); 
      SolutionStepData::Pointer CreateNewSolutionStep(void);
      SolutionStepData::Pointer CreateWithCurrentSolutionStepData(void);

      
         
    private:
     
      SolutionStepsDataContainer::Pointer mpSolutionStepsData;
      
      TimeStepsContainer::Pointer mpTimeStepsData;

	  Process::ProcessArrayType mProcesses;
           
      VariablesContainer mGlobalVariables;

      ProcessInfo mCurrentProcessInfo;



      double mTime;

	  double mDeltaTime;
      
      unsigned int mSolutionStep; 
      
    }; /* Class Model */  
  
}  /* namespace Kratos.*/

#endif /* KRATOS_MODEL  defined */

