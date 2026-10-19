
#if !defined(KRATOS_VOLUME_ELEMENT_3TH)
#define KRATOS_VOLUME_ELEMENT_3TH 

#include "../includes/element.h"

namespace Kratos
{

    class VolumeElement_3th : public Element
    {
        public:

            VolumeElement_3th(NodesArrayType& rNodes, 
                              Properties::Pointer rProperties):Element(rNodes, rProperties)
            {}

            void GetStiffMatrix (Matrix<std::complex<double> >& StiffMatrix); 
              
            void GetEquationIdVector(Vector<int>& EquationId);

            void SetPeso(double Peso);

			void CalculateDerivatives(Vector<std::complex<double> >& ncHx,
		                              Vector<std::complex<double> >& ncHy,
									  Vector<std::complex<double> >& ncHz);
        
        private:
            
            double mPeso;

			double RtoS(double alpha, double beta, double gamma);

			double X(int i){ return (mNodes[i-1]->X()); }
            double Y(int i){ return (mNodes[i-1]->Y()); }
            double Z(int i){ return (mNodes[i-1]->Z()); }

            double Calculate_Volume(); 

			void Calculate_DN(double DN[3][4]);

			void LagrangeDerivatives3D_3th(Matrix<double>& dNdx, 
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