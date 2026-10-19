
#if !defined(KRATOS_VOLUME_ELEMENT_2NDQ)
#define KRATOS_VOLUME_ELEMENT_2NDQ

#include "../includes/element.h"

namespace Kratos
{

    class VolumeElement_2ndQ : public Element
    {
        public:

            VolumeElement_2ndQ(NodesArrayType& rNodes, 
                               Properties::Pointer rProperties) : Element(rNodes, rProperties)
            {};

            void GetStiffMatrix (Matrix<std::complex<double> >& StiffMatrix); 
              
            void GetEquationIdVector(Vector<int>& EquationId);

            void SetPeso(double Peso);

			void CalculateDerivatives(Vector<std::complex<double> >& ncHx,
		                              Vector<std::complex<double> >& ncHy,
									  Vector<std::complex<double> >& ncHz);
        
        private:
            
            double mPeso;

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

			void LagrangeDerivatives3D_2ndQ(Matrix<double>& dNdx, 
				                            Matrix<double>& dNdy, 
										    Matrix<double>& dNdz,
										    std::vector<double>& cX, 
										    std::vector<double>& cY, 
										    std::vector<double>& cZ);

			void DerivativesInNodes(Matrix<double>& ndNdx, 
		                            Matrix<double>& ndNdy, 
									Matrix<double>& ndNdz);

    }; 
} 

#endif