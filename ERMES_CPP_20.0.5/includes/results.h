
#if !defined(KRATOS_RESULTS)
#define KRATOS_RESULTS

#include "../includes/kratos_string.h"

namespace Kratos
{
    class Results
    {
      public:
      
        Results();
        
        virtual ~Results();
        
        void AddResult( const String& rResultName, double ResultTime, double Value )
        {
            mData[rResultName]= Value;
        }
        
        double GetResult( const String& rResultName, double ResultTime )
        {
            return mData[rResultName];
        }
       
      private:
        
        std::map<String, double> mData;
    }; 
}  

#endif 

