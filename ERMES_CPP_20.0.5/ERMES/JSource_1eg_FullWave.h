
#if !defined(KRATOS_JSOURCE_1EG_FULLWAVE)
#define KRATOS_JSOURCE_1EG_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class JSource_1eg_FullWave : public Element
    {
      public:

        JSource_1eg_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumEdges = 6;
            mNumVerts = 4;

            mNumNodes = 10;
            mNumDofs  = 30; 

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

        int mNumEdges;
        int mNumVerts;

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

        void Whitney3D_Ni_1st( Matrix<double>& Nx, Matrix<double>& Ny, Matrix<double>& Nz, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
        
        void Get_CartesianXYZ_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
        void Get_Axisymmetric_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
        void Get_PlasmaModeJs_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
    }; 
} 

#endif 
