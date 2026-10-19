
#if !defined(KRATOS_PROCESS_INFO)
#define KRATOS_PROCESS_INFO

#include "../includes/variables_container.h"

namespace Kratos
{
    class ProcessInfo : public VariablesContainer
    {
      public:
       
        ProcessInfo()
        {
        }
        
	    ProcessInfo( const ProcessInfo& Other ) : VariablesContainer( Other                ),
	                                              mDeltaTime        ( Other.mDeltaTime     ),
                                                  mCurrentTime      ( Other.mCurrentTime   ), 
                                                  mStepNumber       ( Other.mStepNumber    ), 
                                                  mTimeStepNumber   ( Other.mTimeStepNumber)
        {
        }
        
	    ProcessInfo(double CurrentTime, 
                    double DeltaTime, 
                    int    SolutionStep, 
                    int    CurrentTimeStep=0) : mDeltaTime     ( DeltaTime       ),
                                                mCurrentTime   ( CurrentTime     ), 
                                                mStepNumber    ( SolutionStep    ),
                                                mTimeStepNumber( CurrentTimeStep )
        {
        }
        
        virtual ~ProcessInfo()
        {
        }
          
	    void SetCurrentTime( double CurrentTime )
	    {
            mCurrentTime = CurrentTime;
        }
	    
	    double GetCurrentTime( void ) const 
	    {
            return mCurrentTime;
        }
	    
	    void SetDeltaTime( double NewDeltaTime )
	    {
            mDeltaTime = NewDeltaTime;
        }
	    
	    double GetDeltaTime( void ) const 
	    {
            return mDeltaTime;
        }
        
	    void SetCurrentTimeStep( int NewCurrentTimeStep )
	    {
            mTimeStepNumber = NewCurrentTimeStep;
        }
	    
	    double GetCurrentTimeStep( void ) const 
	    {
            return mTimeStepNumber;
        }
        
	    void SetCurrentSolutionStep( int NewCurrentSolutionStep )
	    {
            mStepNumber = NewCurrentSolutionStep;
        }
	    
	    double GetCurrentSolutionStep( void ) const 
	    {
            return mStepNumber;
        }
    
      private:
         
        double mDeltaTime;
	    double mCurrentTime;
        
	    int mStepNumber;
	    int mTimeStepNumber;
    }; 
} 

#endif 

