
#include "../includes/dof.h"

namespace Kratos
{
    Dof::~Dof()
    {
    }
    
    unsigned int Dof::GetVariableKey()
    {
        return mVariableIndex;
    }
    
    bool  Dof::IsFixed() const
    {
        return mIsFixed;
    }
    
    unsigned int Dof::EquationId()
    {
        return mEquationId;
    }
    
    void  Dof::SetEquationId( unsigned int NewEquationId )
    {
        mEquationId = NewEquationId;
    }
    
    void Dof::FixDof()
    {
        mIsFixed = true;
    }
    
    void  Dof::FreeDof()
    {
        mIsFixed = false;
    }
    
    unsigned int Dof::NodeIndex() const
    {
        return mNodeIndex;
    }
} 
