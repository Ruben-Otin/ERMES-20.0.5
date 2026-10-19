
#if !defined(KRATOS_JSOURCE_2NDQ)
#define  KRATOS_JSOURCE_2NDQ

#include "../includes/element.h"

namespace Kratos
{
    class JSource_2ndQ : public Element
    {
        public:

            JSource_2ndQ(NodesArrayType& rNodes, 
                         Properties::Pointer rProperties): Element(rNodes, rProperties)
            {}
              
	        void GetResidualVector  (Vector<std::complex<double> >& ResidualVector);
              
            void GetEquationIdVector(Vector<int>& EquationId);
                 
        private:

			void Calculate_detJ(std::vector<double>& detJ,
							    std::vector<double>& cX, 
							    std::vector<double>& cY, 
							    std::vector<double>& cZ);

			void NaturalDerivatives3D_2ndQ(Matrix<double>& dNda,
										   Matrix<double>& dNdb,
										   Matrix<double>& dNdu,
										   std::vector<double>& cX, 
										   std::vector<double>& cY, 
	  								       std::vector<double>& cZ);

    }; 
} 

#endif 
