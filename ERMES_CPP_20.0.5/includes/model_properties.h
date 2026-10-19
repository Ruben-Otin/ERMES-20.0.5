
#if !defined(KRATOS_MODEL_PROPERTIES)
#define KRATOS_MODEL_PROPERTIES

#include <memory>
#include "../includes/properties_default_function.h"

namespace Kratos
{
    class Model;

    template<class TDataType>
    class ModelProperties : public PropertiesDefaultFunction<TDataType>
    {
      public:
      
        typedef std::shared_ptr<ModelProperties<TDataType> > Pointer;
        
        ModelProperties( Model& rModel ) : mrModel( rModel )
        {
        }
        
        virtual ~ModelProperties()
        {
        }
        
        TDataType operator()( const Variable<TDataType>& rVariable, const Node& rThisNode, double Time )
        {
            return mrModel( rVariable, rThisNode, Time );
        }
        
        TDataType operator()( const Variable<TDataType>& rVariable, const Node& rThisNode )
        {
            return mrModel( rVariable, rThisNode );
        }
          
      private:
      
        Model& mrModel; 
    }; 
}  

#endif 

