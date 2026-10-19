
#if !defined(KRATOS_FAR_FIELD_1ST_COLDPLASMA)
#define KRATOS_FAR_FIELD_1ST_COLDPLASMA

#include "../includes/element.h"

namespace Kratos
{
    class FarField_1st_ColdPlasma : public Element
    {
      public:
      
        FarField_1st_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 3;
            mNumDofs  = 9; 

            mPotentialsOn = false;

            mPeso = 1.0; 
            mArea = Calculate_Area();
        }
        
        void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ) 
        { 
            mpColdPlasma = pColdPlasma; 
        }

        void SetFrequency( double NewFrequency )             
        { 
            mFreq = NewFrequency; 
        }

        void SetPeso( double Peso ) 
        { 
            mPeso = Peso;        
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 12;
                mPotentialsOn = true; 
            }
        }

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ); 
      
      private:
      
        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;

        double mFreq;
        double mPeso;
        double mArea;

        ColdPlasma::Pointer mpColdPlasma;
        
        String Get_Wave_Type();
        
        double Calculate_Area(); 
        
        void Calculate_ExtNormal( Vector<double>& extN );

        void jwK_OnGaussPoints( Vector< std::complex<double> >& jwKgp, Matrix<double>& N );

        void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, Matrix<double>& N );

        void Get_invTau_OnGaussPoints( Vector< std::complex<double> >& invTau, Vector< Matrix< std::complex<double> > >& TEgp );

        void Get_Non_Regularized_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Get_Stiffness_Matrix_Ef( Matrix< std::complex<double> >& StiffMatrix );

        void Get_Stiffness_Matrix_AV( Matrix< std::complex<double> >& StiffMatrix );
    }; 
} 

#endif 