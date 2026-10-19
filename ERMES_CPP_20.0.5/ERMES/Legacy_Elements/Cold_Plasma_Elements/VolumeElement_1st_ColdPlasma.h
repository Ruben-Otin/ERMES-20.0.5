
#if !defined(KRATOS_VOLUME_ELEMENT_1ST_COLDPLASMA)
#define KRATOS_VOLUME_ELEMENT_1ST_COLDPLASMA

#include "../includes/element.h"
#include "../ERMES/ColdPlasma.h"

namespace Kratos
{
    class VolumeElement_1st_ColdPlasma : public Element
    {
      public:

        VolumeElement_1st_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
	    {
            mNumNodes     = 4    ;
            mNumDofs      = 12   ; 
            mPeso         = 1.00 ; 
            mPotentialsOn = false;
            mVolume       = Calculate_Volume();	    

            SetPlasmaIHL();
	    }

	    void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ) 
        { 
            mpColdPlasma = pColdPlasma; 
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs      = 16  ;
                mPotentialsOn = true; 
            }
        }

        void SetPeso( double Peso ) 
        { 
            mPeso = Peso;        
        }

        bool GetIsIHL()
        {
            return mIs_IHL;
        }
        
        void GetEquationIdVector( Vector< int >& EquationId );

        void GetStiffnessMatrix ( Matrix< std::complex<double> >& StiffMatrix ); 

        void GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue    , 
                                          Matrix  <               std::complex<double> >& EleStiffMatrix,  
                                          Vector  <               std::complex<double> >& ResidualVector );

      private:

        int    mNumNodes    ;
        int    mNumDofs     ;

        bool   mPotentialsOn;
        bool   mIs_IHL      ;

        double mPeso        ;
        double mVolume      ;        
          
        ColdPlasma::Pointer mpColdPlasma;

	    double X( int i ){ return mNodes[ i - 1 ]->X(); }
        double Y( int i ){ return mNodes[ i - 1 ]->Y(); }
        double Z( int i ){ return mNodes[ i - 1 ]->Z(); }

        double Calculate_Volume(); 

        void SetPlasmaIHL();  

        void Calculate_DN( double DN[3][4] );

        void Rotation_Matrix( Matrix< std::complex<double> >& T, double Bx, double By, double Bz );

        void DerivativesOnPoints( Matrix< double >& dNdx, 
                                  Matrix< double >& dNdy, 
                                  Matrix< double >& dNdz,
         						  Vector< double >&   cX, 
                                  Vector< double >&   cY, 
                                  Vector< double >&   cZ );

	    void PermittivityTensor_OnNodes      ( Vector< Matrix< std::complex<double> > >& TEn );
        
        void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp,  
                                               Vector< Matrix< std::complex<double> > >& TEn , 
                                               Matrix< double >& N );
        
        void PermittivityTensor_Derivatives  ( Vector< Matrix< std::complex<double> > >& dTEdx, 
                                               Vector< Matrix< std::complex<double> > >& dTEdy,
                                               Vector< Matrix< std::complex<double> > >& dTEdz,
                                               Vector< Matrix< std::complex<double> > >& TEn,
                                               Matrix< double >& dNdx, 
                                               Matrix< double >& dNdy, 
                                               Matrix< double >& dNdz );

        void Get_invTau_InGaussPoints ( Vector<         std::complex<double>   >& invTau, 
                                        Vector< Matrix< std::complex<double> > >& TEgp   );

        void DivD_a_DivD_OnGaussPoints( Vector< Matrix< std::complex<double> > >& Div_ce_N, 
                                        Vector< Matrix< std::complex<double> > >& aDiv_e_E, 
                                        Vector< Matrix< std::complex<double> > >& TEn ,
                                        Vector< Matrix< std::complex<double> > >& TEgp,  
                                        Matrix< double >& N   , 
                                        Matrix< double >& dNdx, 
                                        Matrix< double >& dNdy, 
                                        Matrix< double >& dNdz );

        void Get_NonRegularized_Matrix( Matrix< std::complex<double> >& StiffMatrix );

        void Get_Regularized_Ef_Matrix( Matrix< std::complex<double> >& StiffMatrix );

        void Get_Regularized_AV_Matrix( Matrix< std::complex<double> >& StiffMatrix );

        void Generate_Lower_Diagonal  ( Matrix< std::complex<double> >& StiffMatrix );
    }; 
} 

#endif