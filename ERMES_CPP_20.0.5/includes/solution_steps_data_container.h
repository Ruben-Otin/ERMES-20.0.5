
#if !defined(KRATOS_SOLUTION_STEPS_DATA_CONTAINER)
#define KRATOS_SOLUTION_STEPS_DATA_CONTAINER

#include "../includes/solution_step_data.h"

namespace Kratos
{
    class SolutionStepsDataContainer
    {
      public:
       
        /** A counted pointer to the SolutionStepsDataContainer */
        typedef std::shared_ptr<SolutionStepsDataContainer> Pointer      ;
        typedef std::deque<SolutionStepData::Pointer>         ContainerType;
        typedef ContainerType::size_type                      SizeType     ;
         
        /** Default Constructor. Initialize the current solution step data with
        an empty SolutionStepData object and set number of buffer steps equal to zero.*/
        SolutionStepsDataContainer() : mpCurrentSolutionStepData( new SolutionStepData() ), mStepsNumber( int() )
        {
            mNumberOfBufferSteps = 0;
        }
        
        SolutionStepsDataContainer( SizeType NewNumberOfBufferSteps ) : mpCurrentSolutionStepData( new SolutionStepData() ), mStepsNumber( int() )
        {
            mNumberOfBufferSteps = NewNumberOfBufferSteps;
        }
      
        virtual ~SolutionStepsDataContainer()
        {
        }
             
        /** Creating new solution step and make it current. It also return
        a counted pointer to new and also current solution step data. */
        SolutionStepData::Pointer CreateNewSolutionStepData()
        {
            mData.push_front( mpCurrentSolutionStepData );
      
            if( mData.size() > mNumberOfBufferSteps ) mData.pop_back();
      
            mpCurrentSolutionStepData = SolutionStepData::Pointer( new SolutionStepData() );
      
            return mpCurrentSolutionStepData;
        }
        
        /** Create a copy of current solution step data and put it into the buffer.
        This would be helpfull where model geometry for new step is not changed
        so the new solution step have the same amount of data of previous one. */
        SolutionStepData::Pointer CreateWithCurrentSolutionStepData()
        {
            mData.push_front( mpCurrentSolutionStepData );
	        mpCurrentSolutionStepData=SolutionStepData::Pointer( new SolutionStepData( *mpCurrentSolutionStepData ) );
	        mStepsNumber++;
	        mpCurrentSolutionStepData->StepNumber( mStepsNumber );
            if( mData.size() > mNumberOfBufferSteps ) mData.pop_back();
            return mpCurrentSolutionStepData;
        }
        
        SolutionStepData::Pointer CurrentSolutionStepData()
        {
            return mpCurrentSolutionStepData;
        }
      
	    int GetStepsNumber()
	    {
            return mStepsNumber;
        }
        
        SolutionStepData::Pointer PreviouseSolutionStepData( int BackStepsNumber )
        {   
	  	    if( BackStepsNumber <= 0                    ) return mpCurrentSolutionStepData;
	  	    if( BackStepsNumber <  mNumberOfBufferSteps ) return mData[BackStepsNumber]   ;
	  	  
	  	    KRATOS_ERROR( std::runtime_error, "SolutionStepDataContainer::SolutionStepData", "No Data in buffer!", "" );
	    }
        
        SizeType SetNumberOfBufferSteps( SizeType NewNumberOfBufferSteps )
        {
            mNumberOfBufferSteps = NewNumberOfBufferSteps;
      
            while( mData.size() > mNumberOfBufferSteps ) mData.pop_back();
      
            return mNumberOfBufferSteps;
        }
        
      private:
              
        ContainerType mData;
          
        SolutionStepData::Pointer mpCurrentSolutionStepData;
          
        unsigned int mNumberOfBufferSteps;
      
	  	int mStepsNumber;
    }; 
} 

#endif 

