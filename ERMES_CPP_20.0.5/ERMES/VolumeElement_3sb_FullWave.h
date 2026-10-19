
#if !defined(KRATOS_VOLUME_ELEMENT_3SB_FULLWAVE)
#define KRATOS_VOLUME_ELEMENT_3SB_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class VolumeElement_3sb_FullWave : public Element
    {
      public:

        VolumeElement_3sb_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 17;
            mNumDofs  = 51;

            mPotentialsOn = false;
            mComplxFreqOn = false;

            mVolume = Calculate_Volume();
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

        void SetPotentials( bool PotentialsOn )
        {   
            if ( PotentialsOn ) 
            {
                mNumDofs = 55;
                mPotentialsOn = true; 
            }
        }

        void SetFrequency( std::complex<double> NewFrequency )
        {
            mComplxFreqOn = true;
            mComplxFreq   = NewFrequency;
        }

        int GetInnerGiDGaussPoints( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, int numResultsOnGPs );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix );    
                                 
        void GetEquationIdVector( Vector<int>& EquationId ); 

        void GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, Matrix< std::complex<double> >& EleStiffMatrix, Vector< std::complex<double> >& ResidualVector );
                                 
	    void Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes );

        void Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs );

        void Calculate_E_field_OnNodes( cVector2Type& Ef_OnNodes );

        void Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs );

      private:
          
        int mNumNodes;
        int mNumDofs;

        bool mPotentialsOn;
        bool mComplxFreqOn;

        double mPeso;
        double mVolume;
        double mhkFactor;

        std::complex<double> mComplxFreq;

        std::map< unsigned int, char >* gpNormalType;

        double Calculate_Volume();

        double Calculate_Area( Vector<int>& face );   

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void Calculate_ExtNormal( Vector<double>& extN, Vector<int>& face );

        void Calculate_DN( double DN[ 3 ][ 4 ] );

	    void Derivatives_On_Points( Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void Calculate_Rot_P1_Projection( cVectorType& cRotP1x, cVectorType& cRotP1y, cVectorType& cRotP1z );

        void Add_C_NxN_Matrix( Matrix< std::complex<double> >& C_cur ); 

        void Add_CurlDiv_LL2P_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_CurlDiv_Stab_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_AV_VA_VV_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 
    }; 
} 

#endif