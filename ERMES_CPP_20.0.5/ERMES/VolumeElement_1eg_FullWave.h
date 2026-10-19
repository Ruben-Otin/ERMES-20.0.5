
#if !defined(KRATOS_VOLUME_ELEMENT_1EG_FULLWAVE)
#define KRATOS_VOLUME_ELEMENT_1EG_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class VolumeElement_1eg_FullWave : public Element
    {
      public:

        VolumeElement_1eg_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumEdges = 6;
            mNumVerts = 4;

            mNumNodes = 10;
            mNumDofs  = 30; 

            mPotentialsOn = false;
            mStabilizerOn = false;
            mComplxFreqOn = false;

            mVolume = Calculate_Volume();
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 34;
                mPotentialsOn = true; 
            }
        }

        void SetStabilizer( bool StabilizerOn )
        {
            mStabilizerOn = StabilizerOn;
        }

        void SetFrequency( std::complex<double> NewFrequency )
        {
            mComplxFreqOn = true;
            mComplxFreq   = NewFrequency;
        }

        int GetInnerGiDGaussPoints( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, int numResultsOnGPs );

        void GetEquationIdVector( Vector<int>& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix );   

        void GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, Matrix< std::complex<double> >& EleStiffMatrix, Vector< std::complex<double> >& ResidualVector );

	    void Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes );

        void Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs );

        void Calculate_E_field_OnNodes( cVector2Type& Ef_OnNodes );

        void Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs );

      private:

        int mNumEdges;
        int mNumVerts;

        int mNumNodes;
        int mNumDofs;
        
        bool mPotentialsOn;
        bool mStabilizerOn;
        bool mComplxFreqOn;

        double mVolume;

        std::complex<double> mComplxFreq;

        double Calculate_Volume();

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void Calculate_DN( double DN[ 3 ][ 4 ] );

        void Calculate_Edge_Lengths( Vector<double>& eLength );

        void Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void Derivatives_On_Points( Vector<double>& dNx_dX, Vector<double>& dNx_dY, Vector<double>& dNx_dZ, Vector<double>& dNy_dX, Vector<double>& dNy_dY, Vector<double>& dNy_dZ, Vector<double>& dNz_dX, Vector<double>& dNz_dY, Vector<double>& dNz_dZ );

        void Get_CurlCurl_NiNj_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_Stabilization_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_Potentials_AV_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 
    }; 
} 

#endif