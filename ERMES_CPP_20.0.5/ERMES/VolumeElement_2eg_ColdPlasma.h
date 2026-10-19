
#if !defined(KRATOS_VOLUME_ELEMENT_2EG_COLDPLASMA)
#define KRATOS_VOLUME_ELEMENT_2EG_COLDPLASMA

#include "../includes/element.h"
#include "../ERMES/ColdPlasma.h"

namespace Kratos
{
    class VolumeElement_2eg_ColdPlasma : public Element
    {
      public:

        VolumeElement_2eg_ColdPlasma( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
	    {
            mNumFaces = 4;            
            mNumEdges = 6;
            mNumVerts = 4;

            mNumBases = 2*mNumEdges + 2*mNumFaces;

            mNumNodes = 20;
            mNumDofs  = 60; 

            mPotentialsOn = false;
            mStabilizerOn = false;

            mVolume = Calculate_Volume();    

            SetPlasmaIHL();
	    }
        
	    void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ) 
        { 
            mpColdPlasma = pColdPlasma; 
        }

        void SetStabilizer( bool StabilizerOn )
        {
            mStabilizerOn = StabilizerOn;
        }

        bool GetIsIHL()
        {
            return mIs_IHL;
        }

        void SetPotentials( bool PotentialsOn )
        {   
            if( PotentialsOn ) 
            {
                mNumDofs = 64;
                mPotentialsOn = true; 
            }
        }
        
        void GetEquationIdVector( Vector< int >& EquationId );

        void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ); 

        void GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, Matrix< std::complex<double> >& EleStiffMatrix, Vector< std::complex<double> >& ResidualVector );

      private:

        int mNumFaces;
        int mNumEdges;
        int mNumVerts;
        int mNumBases;

        int mNumNodes;
        int mNumDofs;
        
        bool mIs_IHL;
        bool mPotentialsOn;
        bool mStabilizerOn;

        double mVolume;

        ColdPlasma::Pointer mpColdPlasma;

        double Calculate_Volume(); 

	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void SetPlasmaIHL();  
        
        void Calculate_DN( double DN[ 3 ][ 4 ] );

        void Order_Face_Bases( int n[ 8 ], int e[ 8 ] );

        void Calculate_Edge_Lengths( Vector<double>& eLength );

        void Face_Base( int edgeBase[ 2 ], int nodeBase[ 2 ], int nodeID[ 3 ], int edgeID[ 3 ], int oppnID[ 3 ] );

        void Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void Whitney3D_Ni_2nd( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

        void Derivatives_Whitney_1st( Vector<double>& dNx_dX, Vector<double>& dNx_dY, Vector<double>& dNx_dZ, Vector<double>& dNy_dX, Vector<double>& dNy_dY, Vector<double>& dNy_dZ, Vector<double>& dNz_dX, Vector<double>& dNz_dY, Vector<double>& dNz_dZ );

        void Derivatives_On_Points( Matrix<double>& dNx_dX, Matrix<double>& dNx_dY, Matrix<double>& dNx_dZ, Matrix<double>& dNy_dX, Matrix<double>& dNy_dY, Matrix<double>& dNy_dZ, Matrix<double>& dNz_dX, Matrix<double>& dNz_dY, Matrix<double>& dNz_dZ, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

	    void PermittivityTensor_OnNodes( Vector< Matrix< std::complex<double> > >& TEn );
        
        void PermittivityTensor_OnGaussPoints( Vector< Matrix< std::complex<double> > >& TEgp, Matrix< double >& N );

        void Get_Edge_Formulation_Matrix_EA( Matrix< std::complex<double> >& StiffMatrix );

        void Add_Edge_Stabilizatn_Matrix_EA( Matrix< std::complex<double> >& StiffMatrix, Vector< Matrix< std::complex<double> > >& TEgp, Vector<double>& W, Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz );

        void Add_Edge_Potentials_Matrix_AVV( Matrix< std::complex<double> >& StiffMatrix, Vector< Matrix< std::complex<double> > >& TEgp, Vector<double>& W, Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz );
    }; 
} 

#endif