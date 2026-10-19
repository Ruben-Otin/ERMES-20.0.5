
#if !defined(KRATOS_VOLUME_ELEMENT_2NDQ_DIVD_PLASMA)
#define KRATOS_VOLUME_ELEMENT_2NDQ_DIVD_PLASMA

#include "../includes/element.h"
#include "../ERMES/ColdPlasma.h"

namespace Kratos
{

    class VolumeElement_2ndQ_DivD_Plasma : public Element
    {
        public:

            VolumeElement_2ndQ_DivD_Plasma( NodesArrayType&     rNodes, 
				                            Properties::Pointer rProperties ) : Element(rNodes, rProperties)
			{
				numNodes = 10;
				numDofs  = numNodes * 3;
                freq     = (*mProperties)(FREQUENCY);
				eo       = 8.8541878176e-12;
				pi       = 3.1415926535897932384626433832795;
				mo       = (4.00e-7) * pi;      
                Is_IHL   = false;
			}

			void Set_Plasma_Model  ( ColdPlasma::Pointer pColdPlasma ) { mpColdPlasma = pColdPlasma; }
            void SetPeso           ( double Peso                     ) { mPeso        = Peso;        }

            void Set_Plasma_or_IHL();     
            
            void GetEquationIdVector ( Vector<int>& EquationId );
            
            void GetStiffMatrix      ( Matrix<std::complex<double> >& StiffMatrix ); 

			void CalculateDerivatives( Vector<std::complex<double> >& ncHx,
		                               Vector<std::complex<double> >& ncHy,
									   Vector<std::complex<double> >& ncHz );
        
        private:

            ColdPlasma::Pointer mpColdPlasma;

            double mPeso;
            double freq;

			double eo;
			double pi;
			double mo;

            int numNodes;
			int numDofs;

            bool Is_IHL;
			
			void Generate_UpperDiagonal_CurlCurl( Matrix<std::complex<double> >& StiffMatrix );
            void Generate_UpperDiagonal_CurlDiv ( Matrix<std::complex<double> >& StiffMatrix );
			void Generate_LowerDiagonal         ( Matrix<std::complex<double> >& StiffMatrix );

            void Calculate_detJ            ( Vector<double>& detJ, 
                                             Vector<double>&   cX, Vector<double>&     cY, Vector<double>&    cZ );

            void DerivativesInNodes        ( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz );

			void NaturalDerivatives3D_2ndQ ( Matrix<double>&  dNda, Matrix<double>&  dNdb, Matrix<double>&  dNdu,
                                             Vector<double>&    cX, Vector<double>&    cY, Vector<double>&    cZ );

			void LagrangeDerivatives3D_2ndQ( Matrix<double>&  dNdx, Matrix<double>&  dNdy, Matrix<double>&  dNdz,
                                             Vector<double>&    cX, Vector<double>&    cY, Vector<double>&    cZ );

            void Rotation_Matrix                 ( Matrix<std::complex<double> >& T, double Bx, double By, double Bz );

            void PermittivityTensor_OnNodes      ( Vector<Matrix<std::complex<double> > >& TEn );

            void PermittivityTensor_OnGaussPoints( Vector<Matrix<std::complex<double> > >& TEgp,  
                                                   Vector<Matrix<std::complex<double> > >& TEn , 
                                                   Matrix<double>& N );

            void DivD_a_DivD_OnGaussPoints( Vector<Matrix<std::complex<double> > >& Div_ce_N, 
                                            Vector<Matrix<std::complex<double> > >& aDiv_e_E, 
                                            Vector<Matrix<std::complex<double> > >& TEn ,
                                            Vector<Matrix<std::complex<double> > >& TEgp,  
                                            Matrix<double>& N   , 
                                            Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz );

            void PermittivityTensor_Derivatives( Vector<Matrix<std::complex<double> > >& dTEdx, 
                                                 Vector<Matrix<std::complex<double> > >& dTEdy,
                                                 Vector<Matrix<std::complex<double> > >& dTEdz,
                                                 Vector<Matrix<std::complex<double> > >& TEn  ,
                                                 Matrix<double>& dNdx, 
                                                 Matrix<double>& dNdy, 
                                                 Matrix<double>& dNdz );
    }; 
} 

#endif