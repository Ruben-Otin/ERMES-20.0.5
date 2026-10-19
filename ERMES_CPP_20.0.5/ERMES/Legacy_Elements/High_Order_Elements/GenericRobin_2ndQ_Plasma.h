
#if !defined(KRATOS_GENERIC_ROBIN_2NDQ_PLASMA)
#define  KRATOS_GENERIC_ROBIN_2NDQ_PLASMA

#include "../includes/element.h"

namespace Kratos
{
     class GenericRobin_2ndQ_Plasma : public Element
     {
          public:
          
              GenericRobin_2ndQ_Plasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
              {
                  pi = 3.1415926535897932384626433832795;
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
              double pi  ;
              double eo  ;
              
              int numNodes;
              int numDofs;
              
              double Calculate_Area(); 

              void Calculate_detJ ( Vector<double>& detJ, Vector<double>& cX, Vector<double>& cY );
              
              void Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz );
              
              void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, Matrix<double>& N );      
     }; 
} 

#endif 