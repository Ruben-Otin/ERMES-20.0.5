
#if !defined(KRATOS_GENERIC_ROBIN_3SB_COLDPLASMA)
#define  KRATOS_GENERIC_ROBIN_3SB_COLDPLASMA

#include "../includes/element.h"

namespace Kratos
{
    class GenericRobin_3sb_ColdPlasma : public Element
    {
      public:
         
        GenericRobin_3sb_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            numNodes = 6;
            numDofs  = numNodes * 3;               
        }
        
        void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ) 
        { 
            mpColdPlasma = pColdPlasma; 
        }
        
        void SetFrequency( double NewFrequency )             
        { 
            freq = NewFrequency; 
        }
        
        void GetStiffnessMatrix( Matrix<std::complex<double> >& StiffMatrix ); 
        
        void GetEquationIdVector( Vector<int>& EquationId );
         
      private:
         
        ColdPlasma::Pointer mpColdPlasma;
        
        double freq;
        
        int numNodes;
        int numDofs;
        
        double Calculate_Area(); 
        
        void Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz );
        
        void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, Matrix<double>& N );      
    }; 
} 

#endif 