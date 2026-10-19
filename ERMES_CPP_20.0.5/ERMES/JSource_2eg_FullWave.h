
#if !defined(KRATOS_JSOURCE_2EG_FULLWAVE)
#define KRATOS_JSOURCE_2EG_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class JSource_2eg_FullWave : public Element
    {
      public:

        JSource_2eg_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumFaces = 4;            
            mNumEdges = 6;
            mNumVerts = 4;

            mNumBases = 2*mNumEdges + 2*mNumFaces;

            mNumNodes = 20;
            mNumDofs  = 60; 

            mPotentialsOn = false;

            mVolume = Calculate_Volume();
        }
              
        void SetPotentials( bool PotentialsOn )
        {   
            mPotentialsOn = PotentialsOn; 
        }

        void Set_Cartesian_J_Nodal( Vector< Vector< std::complex<double> > >& Nodal_J )
        {
            mNodal_J = Nodal_J;
        }

        void GetEquationIdVector( Vector<int>& EquationId );

	    void GetResidualVector( Vector< std::complex<double> >& ResidualVector );

        void Get_J_Imposed_OnNodes( std::map< unsigned int, Vector<double> >& rJimp, std::map< unsigned int, Vector<double> >& iJimp );
          
      private:

        int mNumFaces;
        int mNumEdges;
        int mNumVerts;
        int mNumBases;

        int mNumNodes;
        int mNumDofs ;

        bool mPotentialsOn;

        double mVolume;

        Vector< Vector< std::complex<double> > > mNodal_J;

        double Calculate_Volume(); 
        
	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }

        void Calculate_DN( double DN[ 3 ][ 4 ] );

        void Calculate_Edge_Lengths( Vector<double>& eLength );

        void Order_Face_Bases( int n[ 8 ], int e[ 8 ] );

        void Face_Base( int edgeBase[ 2 ], int nodeBase[ 2 ], int nodeID[ 3 ], int edgeID[ 3 ], int oppnID[ 3 ] );

        void Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
        void Whitney3D_Ni_2nd( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
        
        void Get_CartesianXYZ_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
        void Get_Axisymmetric_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
        void Get_PlasmaModeJs_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
    }; 
} 

#endif 
