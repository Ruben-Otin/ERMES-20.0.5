
#if !defined(KRATOS_FAR_FIELD_3SB_COLDPLASMA)
#define KRATOS_FAR_FIELD_3SB_COLDPLASMA

#include "../includes/element.h"

namespace Kratos
{
    class FarField_3sb_ColdPlasma : public Element
    {
      public:

        FarField_3sb_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 6;
            mNumDofs  = 18; 

            mPotentialsOn = false;

            mArea = Calculate_Area();
        }

        void SetPeso( double Peso )
        { 
            mPeso = Peso;        
        }

        void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ) 
        { 
            mpColdPlasma = pColdPlasma; 
        }
        
        void SetFrequency( double NewFrequency )             
        { 
            mFreq = NewFrequency; 
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 21;
                mPotentialsOn = true; 
            }
        } 

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ); 

      private:

        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;

        double mPeso;
        double mFreq;
        double mArea;

        ColdPlasma::Pointer mpColdPlasma;
        
        String Get_Wave_Type();
        
        double Calculate_Area(); 
        
        void Calculate_ExtNormal( Vector<double>& extN );

        void jwK_OnGaussPoints( Vector< std::complex<double> >& jwKgp, Matrix<double>& N );

        void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, Matrix<double>& N );

        void Get_invTau_InGaussPoints( Vector< std::complex<double> >& invTau, Vector< Matrix< std::complex<double> > >& TEgp );

        void Get_Non_Regularized_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Get_Stiffness_Matrix_Ef( Matrix< std::complex<double> >& StiffMatrix );

        void Get_Stiffness_Matrix_AV( Matrix< std::complex<double> >& StiffMatrix );
    }; 
} 

#endif 