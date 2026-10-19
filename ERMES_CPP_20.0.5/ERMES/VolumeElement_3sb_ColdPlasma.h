
#if !defined(KRATOS_VOLUME_ELEMENT_3SB_COLDPLASMA)
#define KRATOS_VOLUME_ELEMENT_3SB_COLDPLASMA

#include "../includes/element.h"
#include "../ERMES/ColdPlasma.h"

namespace Kratos
{
    class VolumeElement_3sb_ColdPlasma : public Element
    {
      public:

        VolumeElement_3sb_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
	    {
	    	mNumNodes = 17;
	    	mNumDofs  = 51;

            mPotentialsOn = false;

            mVolume = Calculate_Volume();	    

            SetPlasmaIHL();
	    }

        void SetNormalType( std::map< unsigned int, char >& Type_Of_BCNormal )
        {
	        gpNormalType = &Type_Of_BCNormal;
	    }

        void SetPeso( double Peso )
        { 
            mPeso = Peso;        
        }

        void SetFactor( double Factor )
	    { 
	        mhkFactor = Factor;
	    }

	    void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ) 
        { 
            mpColdPlasma = pColdPlasma; 
        }

        bool GetIsIHL()
        {
            return mIs_IHL;
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if ( PotentialsOn ) 
            {
                mNumDofs = 55;
                mPotentialsOn = true; 
            }
        }
        
        void GetEquationIdVector( Vector< int >& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, Matrix< std::complex<double> >& EleStiffMatrix, Vector< std::complex<double> >& ResidualVector );

      private:

        int mNumNodes;
        int mNumDofs;

        bool mIs_IHL;
        bool mPotentialsOn;

        double mPeso;
        double mVolume;
        double mhkFactor;

        ColdPlasma::Pointer mpColdPlasma;

        std::map< unsigned int, char >* gpNormalType;

        double Calculate_Volume(); 

        double Calculate_Area( Vector<int>& face );  

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void SetPlasmaIHL();  

        void Calculate_DN( double DN[ 3 ][ 4 ] );

        void Add_C_NxN_Matrix( Matrix< std::complex<double> >& C_cur ); 

        void Calculate_ExtNormal( Vector<double>& extN, Vector<int>& face );

        void Derivatives_On_Points( Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

	    void PermittivityTensor_OnNodes( Vector< Matrix< std::complex<double> > >& TEn );

        void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, Vector< Matrix< std::complex<double> > >& TEn, Matrix< double >& N );

        void PermittivityTensor_Derivatives( Vector< Matrix< std::complex<double> > >& dTEdx, 
                                             Vector< Matrix< std::complex<double> > >& dTEdy, 
                                             Vector< Matrix< std::complex<double> > >& dTEdz,
                                             Vector< Matrix< std::complex<double> > >& TEn, 
                                             Matrix< double >& dNdx, 
                                             Matrix< double >& dNdy, 
                                             Matrix< double >& dNdz );

        void Get_invTau_OnGaussPoints( Vector< std::complex<double> >& invTau, Vector< Matrix< std::complex<double> > >& TEgp );

        void Add_CurlDiv_Stab_Matrix_EA( Matrix< std::complex<double> >& StiffMatrix ); 

        void Get_CurlDiv_LL2P_Matrix_Ef( Matrix< std::complex<double> >& StiffMatrix ); 

        void Get_CurlDiv_LL2P_Matrix_AV( Matrix< std::complex<double> >& StiffMatrix ); 
    }; 
} 

#endif