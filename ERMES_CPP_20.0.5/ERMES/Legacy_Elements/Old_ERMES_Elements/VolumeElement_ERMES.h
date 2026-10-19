

#if !defined(KRATOS_VOLUME_ELEMENT_ERMES)
#define  KRATOS_VOLUME_ELEMENT_ERMES 


#include "../includes/element.h"


namespace Kratos
{

    class VolumeElement_ERMES : public Element
    {
        public:

            VolumeElement_ERMES(NodesArrayType& rNodes, 
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

            void Calculate_DN(double DN[3][4]);

            double Calculate_Volume();     
    }; 
 
} 


#endif 
