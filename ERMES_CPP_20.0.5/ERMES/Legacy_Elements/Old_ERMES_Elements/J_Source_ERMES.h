

#if !defined(KRATOS_J_SOURCE_ERMES)
#define  KRATOS_J_SOURCE_ERMES 


#include "../includes/element.h"


namespace Kratos
{

    class J_Source_ERMES : public Element
    {
        public:

            J_Source_ERMES(NodesArrayType& rNodes, 
                           Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}

            void GetStiffMatrix12x12(static std::complex<double> StiffMatrix_12x12[12][12]);

            void GetResidualVector(static std::complex<double>* ResidualVector_12); 

            void GetMagnitudeVector(std::complex<double>* MagnitudVector);

            void GetEquationIdVector(static int* EqIdVector_12);

        
        private:

            double X(int i){ return (mNodes[i-1]->X()); }

            double Y(int i){ return (mNodes[i-1]->Y()); }

            double Z(int i){ return (mNodes[i-1]->Z()); }

            double Calculate_Volume();     
    }; 
 
} 


#endif 
