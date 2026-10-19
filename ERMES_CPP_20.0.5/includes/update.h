
#if !defined(KRATOS_UPDATE)
#define KRATOS_UPDATE

#include <set>

#include <memory>

#include "../includes/model.h"
#include "../includes/compare_pdof.h"
#include "../includes/exception.h"

namespace Kratos
{
    class Update
    {
      public:
      
        typedef std::shared_ptr<Update> Pointer;
    
        Update( Model::Pointer pNewModel );
    
        virtual ~Update();
    
        virtual void UpdateVariables( Vector<double>& rX, std::set<Dof::Pointer,ComparePDof>& rDofSet );
    
      protected:
        
        Model::Pointer mpModel;
         
      private:
              
        Update( const Update& Other );      
    }; 
}  

#endif 

