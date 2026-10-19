
#if !defined(KRATOS_GENERIC_ROBIN_1ST_ELECTROSTATIC)
#define KRATOS_GENERIC_ROBIN_1ST_ELECTROSTATIC

#include "../includes/element.h"

namespace Kratos
{
    class GenericRobin_1st_Electrostatic : public Element
    {
      public:
      
        GenericRobin_1st_Electrostatic( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 3;
            mNumDofs  = 3;

            mArea = Calculate_Area();
        }

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetStiffnessMatrix( Matrix<double>& StiffMatrix );

        void GetResidualVector( Vector<double>& ResidualVector );

      private:
      
        int mNumNodes;
        int mNumDofs ;

        double mArea;
        
        double Calculate_Area();
    }; 
} 

#endif 