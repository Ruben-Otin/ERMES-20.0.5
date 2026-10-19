
#if !defined(DEFINE_KRATOS_DOF)
#define DEFINE_KRATOS_DOF

#include <vector>
#include <map>

#include <memory>

#include "../includes/kratos_string.h"
#include "../includes/exception.h"
#include "../includes/variables.h"

namespace Kratos
{
    /** Degree of freedom class. Stores all data related to a degree of freedom */
    class Dof
    {
      public:
         
        /** Counted pointer of degree of freedom */
        typedef std::shared_ptr<Dof>  Pointer;
        
        template<class TDataType>
        Dof( const Variable<TDataType>& rNewVariable, unsigned int NodeIndex )
        {
            mVariableIndex = rNewVariable.getKey();
            mIsFixed       = false;
            mNodeIndex     = NodeIndex;
        }
      
        virtual ~Dof();
        
        /** Returns variable assigned to this degree of freedom */
        unsigned int GetVariableKey();
        
        template<class TDataType>
        void SetVariable( const Variable<TDataType>& rThisVariable )
        {
            mVariableIndex = rThisVariable.getKey();
        }
    
        /** Returns true if the degree of freedom is fixed */
        bool IsFixed() const;
    
        /** Return the Equation Id related to this degree eof freedom */
        unsigned int EquationId();
        
        /** returns the id corresponding to the node related to this DoF */
        unsigned int NodeIndex() const;
        
        /** Sets the Equation Id to the desired value */
        void SetEquationId( unsigned int Id );
              
        /** Fixes the Dof. */
        void FixDof();
        
        /** Frees the degree of freedom */
        void FreeDof();
        
      private:
        
        /** Default constructor private.
        Prevents the dof to be created without reffernce to it's node */
        Dof();
        
        /** Private copy constructor prevent object from being copied */
        Dof( const  Dof& rNewDof );
        
        /** Assignement opertor private */
        const Dof& operator=(const Dof& SourceDof);
        
        /** True is is fixed */
        bool mIsFixed;
        
        /** Equation identificator of the degree of freedom */
        unsigned int mEquationId;
        
        /** Identification of the node that correspond to this dof */
        unsigned int mNodeIndex;
        
        /** Variable of the degree of freedom */
        unsigned short int mVariableIndex;
    }; 
}  

#endif 

