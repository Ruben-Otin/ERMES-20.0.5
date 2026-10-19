
#if !defined(KRATOS_FAR_FIELD_2NDQ_PLASMA)
#define  KRATOS_FAR_FIELD_2NDQ_PLASMA

#include "../includes/element.h"

namespace Kratos
{
    class FarField_2ndQ_Plasma: public Element
    {
        public:

            FarField_2ndQ_Plasma( NodesArrayType& rNodes ) : Element( rNodes )
            {
                pi = 3.1415926535897932384626433832795;
                mu = pi*4.00e-7;
                eo = 8.8541878176e-12;

                numNodes = 6;
                numDofs  = numNodes*3; 
            }

            void Set_Plasma_Model( ColdPlasma::Pointer pColdPlasma ) { mpColdPlasma = pColdPlasma ; }
            void SetFrequency    ( double NewFrequency )             { freq         = NewFrequency; }

            void GetStiffMatrix     ( Matrix<std::complex<double> >& StiffMatrix ); 
            void GetEquationIdVector( Vector<int>& EquationId );

        private:

            ColdPlasma::Pointer mpColdPlasma;

            double freq;

            double pi;
            double mu;
            double eo;

            int numNodes;
            int numDofs;

			void Calculate_detJ ( Vector<double>& detJ, Vector<double>& cX, Vector<double>& cY );

            void Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz );

            void Tensor_K_OnNodes      ( Vector<Matrix<std::complex<double> > >& TKn );

            void Tensor_K_OnGaussPoints( Vector<Matrix<std::complex<double> > >& TKgp, Matrix<double>& N );
    }; 
} 

#endif 