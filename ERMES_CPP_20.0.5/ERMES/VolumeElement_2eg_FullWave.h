
#if !defined(KRATOS_VOLUME_ELEMENT_2EG_FULLWAVE)
#define KRATOS_VOLUME_ELEMENT_2EG_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class VolumeElement_2eg_FullWave : public Element
    {
      public:

        VolumeElement_2eg_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumFaces = 4;            
            mNumEdges = 6;
            mNumVerts = 4;

            mNumBases = 2*mNumEdges + 2*mNumFaces;

            mNumNodes = 20;
            mNumDofs  = 60; 

            mPotentialsOn = false;
            mStabilizerOn = false;
            mComplxFreqOn = false;

            mVolume = Calculate_Volume();
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 64;
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

        int mNumFaces;
        int mNumEdges;
        int mNumVerts;
        int mNumBases;

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

        void Order_Face_Bases( int n[ 8 ], int e[ 8 ] );

        void Calculate_Edge_Lengths( Vector<double>& eLength );

        void Face_Base( int edgeBase[ 2 ], int nodeBase[ 2 ], int nodeID[ 3 ], int edgeID[ 3 ], int oppnID[ 3 ] );

        void Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void Whitney3D_Ni_2nd( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void Derivatives_Whitney_1st( Vector<double>& dNx_dX, Vector<double>& dNx_dY, Vector<double>& dNx_dZ, Vector<double>& dNy_dX, Vector<double>& dNy_dY, Vector<double>& dNy_dZ, Vector<double>& dNz_dX, Vector<double>& dNz_dY, Vector<double>& dNz_dZ );

        void Derivatives_On_Points( Matrix<double>& dNx_dX, Matrix<double>& dNx_dY, Matrix<double>& dNx_dZ, Matrix<double>& dNy_dX, Matrix<double>& dNy_dY, Matrix<double>& dNy_dZ, Matrix<double>& dNz_dX, Matrix<double>& dNz_dY, Matrix<double>& dNz_dZ, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void Get_Nodes_Natural_Coordinates( Vector<double>& nX, Vector<double>& nY, Vector<double>& nZ ); 

        void Get_CurlCurl_NiNj_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_Stabilization_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_Potentials_AV_Matrix( Matrix< std::complex<double> >& StiffMatrix ); 
    }; 
} 

#endif