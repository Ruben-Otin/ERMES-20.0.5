
#if !defined(KRATOS_JSOURCE_4TH)
#define  KRATOS_JSOURCE_4TH

#include "../includes/element.h"

namespace Kratos
{
    class JSource_4th : public Element
    {
        public:

            JSource_4th(NodesArrayType& rNodes, 
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
