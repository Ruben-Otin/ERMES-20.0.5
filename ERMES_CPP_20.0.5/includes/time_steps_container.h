
#if !defined(KRATOS_TIME_STEPS_CONTAINER)
#define KRATOS_TIME_STEPS_CONTAINER

#include <list>
#include <queue>
#include <cmath>

#include <memory>

#include "../includes/solution_step_data.h"

namespace Kratos
{
    class TimeStepsContainer
    {
      public:

        typedef std::shared_ptr<TimeStepsContainer>        Pointer                  ;
        typedef std::pair<double, SolutionStepData::Pointer> TagedSolutionStepDataType;
        typedef std::deque<TagedSolutionStepDataType>        ContainerType            ;
        typedef ContainerType::size_type                     SizeType                 ;
        
        TimeStepsContainer( SolutionStepData::Pointer pThisSolutionStepData ) : mpCurrentTimeStepData( pThisSolutionStepData ),
                                                                                mCurrentTime         ( double() ),
                                                                                mCurrentDeltaTime    ( double() )
        {
            mNumberOfBufferSteps = 0;
        }
        
        TimeStepsContainer( SolutionStepData::Pointer pThisSolutionStepData, 
                            SizeType                  NewNumberOfBufferSteps ) : mpCurrentTimeStepData( pThisSolutionStepData ),
                                                                                 mCurrentTime         ( double() ),
                                                                                 mCurrentDeltaTime    ( double() )
        {
            mNumberOfBufferSteps = NewNumberOfBufferSteps;
        }
        
        virtual ~TimeStepsContainer()
        {
        }
        
        // Creating new time step and make it current. 
        // It also return a counted pointer to new and also current time step data. 
        SolutionStepData::Pointer CreateNewTimeStepData( double NewTime, SolutionStepData::Pointer pThisSolutionStepData )
        {
            if ( ( NewTime - mCurrentTime ) <= 0 )
            {
                 KRATOS_ERROR( std::runtime_error, 
                               "TimeStepsContainer::CreateNewTimeStepData", 
                               "Can not create time step older or in the same time as current", "" );
            }
        
            if ( pThisSolutionStepData == mpCurrentTimeStepData )
            {
                 KRATOS_ERROR( std::runtime_error, 
                               "TimeStepsContainer::CreateNewTimeStepData", 
                               "Can not create time step with the same step data as current", "" );
            }
        
            mData.push_front( TagedSolutionStepDataType( mCurrentTime, mpCurrentTimeStepData ) );
        
            mCurrentDeltaTime = NewTime - mCurrentTime;
            mCurrentTime      = NewTime;
        
            if ( mData.size() > mNumberOfBufferSteps )
            {
                mData.pop_back();
            }
        
            mpCurrentTimeStepData = pThisSolutionStepData;
        
            return mpCurrentTimeStepData;
        }
        
        SolutionStepData::Pointer CurrentTimeStepData()
        {
            return mpCurrentTimeStepData;
        }
        
        double CurrentTime()
        {
            return mCurrentTime;
        }
        
	    void SetCurrentTime( double NewTime )
	    {
            mCurrentTime = NewTime;
        }
        
        double CurrentDeltaTime()
        {
            return mCurrentDeltaTime;
        }
        
        void SetDeltaTime( double NewDeltaTime )
        {
            mCurrentDeltaTime = NewDeltaTime;
        }
        
        SizeType SetNumberOfBufferSteps( SizeType NewNumberOfBufferSteps )
        {
            mNumberOfBufferSteps = NewNumberOfBufferSteps;
        
            while( mData.size() > mNumberOfBufferSteps ) 
            {
                mData.pop_back();
            }
        
            return mNumberOfBufferSteps;
        }
        
        SolutionStepData::Pointer GetNearestStep( double ThisTime )
        {
            SolutionStepData::Pointer result = mpCurrentTimeStepData;
        
	        if( ( mCurrentTime - ThisTime ) <= ( mCurrentDeltaTime / 2.00 ) )
            {
            }
	        else
            {
		        double distance = mCurrentTime;

		        for( ContainerType::iterator i = mData.begin(); i != mData.end(); ++i )
                {
		  	        if( fabs( ThisTime - i->first ) < distance )
                    {
		  		        distance = fabs( ThisTime - i->first );
		  		        result   = i->second;
                    }
                }
            }
        
	        return result;
        }

      private:

        ContainerType mData;

        SolutionStepData::Pointer mpCurrentTimeStepData;

        unsigned int mNumberOfBufferSteps;

        double mCurrentTime;

        double mCurrentDeltaTime;

        TimeStepsContainer();

        TimeStepsContainer( const TimeStepsContainer& OtherTimeStepsContainer );

        TimeStepsContainer& operator=( const TimeStepsContainer& OtherTimeStepsContainer );
    }; 
} 

#endif 

