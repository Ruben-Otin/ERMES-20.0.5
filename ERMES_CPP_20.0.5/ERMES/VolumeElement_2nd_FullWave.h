
#if !defined(KRATOS_VOLUME_ELEMENT_2NDQ_FULLWAVE)
#define KRATOS_VOLUME_ELEMENT_2NDQ_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class VolumeElement_2nd_FullWave : public Element
    {
      public:

        VolumeElement_2nd_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 10;
            mNumDofs  = 30; 

            mPotentialsOn = false;
            mComplxFreqOn = false;

            mPeso = 1.0; 
            mVolume = Calculate_Volume();
        };

        void SetPeso( double Peso )
        { 
            mPeso = Peso;        
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 34;
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

        std::complex<double> mComplxFreq;

        double Calculate_Volume();

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void Calculate_DN( double DN[ 3 ][ 4 ] );

        void Derivatives_On_Points( Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void Get_Curl_Div_Matrix ( Matrix< std::complex<double> >& StiffMatrix ); 

        void Get_Curl_Curl_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Get_AV_VA_VV_Matrix ( Matrix< std::complex<double> >& StiffMatrix ); 
    }; 
} 

#endif