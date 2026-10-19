
#if !defined(KRATOS_VOLUME_ELEMENT_2ND_COLDPLASMA)
#define KRATOS_VOLUME_ELEMENT_2ND_COLDPLASMA

#include "../includes/element.h"
#include "../ERMES/ColdPlasma.h"

namespace Kratos
{
    class VolumeElement_2nd_ColdPlasma : public Element
    {
      public:

        VolumeElement_2nd_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
	    {
	    	mNumVerts =  4;
            mNumNodes = 10;
	    	mNumDofs  = 30;

            mPotentialsOn = false;

            mPeso   = 1.0; 
            mVolume = Calculate_Volume();	    

            SetPlasmaIHL();
	    }
        
	    void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ) 
        { 
            mpColdPlasma = pColdPlasma; 
        }

        void SetPeso( double Peso ) 
        { 
            mPeso = Peso;        
        }

        bool GetIsIHL()
        {
            return mIs_IHL;
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 34;
                mPotentialsOn = true; 
            }
        }
        
        void GetEquationIdVector( Vector< int >& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, 
                                          Matrix< std::complex<double> >& EleStiffMatrix, 
                                          Vector< std::complex<double> >& ResidualVector );

      private:

        int mNumVerts;
        int mNumNodes;
        int mNumDofs;

        bool mIs_IHL;
        bool mPotentialsOn;

        double mPeso;
        double mVolume;        

        ColdPlasma::Pointer mpColdPlasma;

        double Calculate_Volume(); 

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void SetPlasmaIHL();  
        
        void Calculate_DN( double DN[ 3 ][ 4 ] );

	    void PermittivityTensor_OnNodes( Vector< Matrix< std::complex<double> > >& TEvt );

        void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, 
                                               Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, 
                                               Vector< Matrix< std::complex<double> > >& TEvt, 
                                               Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void PermittivityTensor_Derivatives( Vector< Matrix< std::complex<double> > >& dTEdx, 
                                             Vector< Matrix< std::complex<double> > >& dTEdy, 
                                             Vector< Matrix< std::complex<double> > >& dTEdz,
                                             Vector< Matrix< std::complex<double> > >& TEvt , 
                                             int numGaussPoints );

        void DivD_a_DivD_OnGaussPoints( Vector< Matrix< std::complex<double> > >& Div_ce_N, 
                                        Vector< Matrix< std::complex<double> > >& aDiv_e_E, 
                                        Vector< Matrix< std::complex<double> > >& TEvt,
                                        Vector< Matrix< std::complex<double> > >& TEgp,  
                                        Matrix< double >& N, 
                                        Matrix< double >& dNdx, Matrix< double >& dNdy, Matrix< double >& dNdz );

        void Derivatives_On_Points( Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz, 
                                    Vector<double>& cX  , Vector<double>& cY  , Vector<double>& cZ );

        void Get_invTau_OnGaussPoints( Vector< std::complex<double> >& invTau, Vector< Matrix< std::complex<double> > >& TEgp );

        void Add_AV_VA_VV_Matrix( Matrix< std::complex<double> >& StiffMatrix, Vector< Matrix< std::complex<double> > >& TEgp );

        void Get_Non_Regularized_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Get_Regularized_Matrix_Ef( Matrix< std::complex<double> >& StiffMatrix ); 

        void Get_Regularized_Matrix_AV( Matrix< std::complex<double> >& StiffMatrix ); 
    }; 
} 

#endif