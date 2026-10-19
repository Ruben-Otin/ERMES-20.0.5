
#if !defined(KRATOS_JSOURCE_2ND)
#define  KRATOS_JSOURCE_2ND

#include "../includes/element.h"

namespace Kratos
{
    class JSource_2nd : public Element
    {
        public:

            JSource_2nd(NodesArrayType& rNodes, 
                        Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}
              
	        void GetResidualVector  (Vector<std::complex<double> >& ResidualVector);
              
            void GetEquationIdVector(Vector<int>& EquationId);

                 
        private:

			double X(int i){ return (mNodes[i-1]->X()); }
            double Y(int i){ return (mNodes[i-1]->Y()); }
            double Z(int i){ return (mNodes[i-1]->Z()); }

            double Calculate_Volume();     
    }; 
} 

#endif 
