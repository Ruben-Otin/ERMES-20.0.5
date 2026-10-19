
#if !defined(KRATOS_JSOURCE_1ST_FULLWAVE)
#define KRATOS_JSOURCE_1ST_FULLWAVE

#include "../includes/element.h"

namespace Kratos
{
    class JSource_1st_FullWave : public Element
    {
      public:

        JSource_1st_FullWave( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 4 ;
            mNumDofs  = 12; 

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

        int mNumNodes;
        int mNumDofs ;

        bool mPotentialsOn;

        double mVolume;

        Vector< Vector< std::complex<double> > > mNodal_J;

        double Calculate_Volume(); 
        
	    double X( int i ){ return mNodes[ i-1 ]->X(); }
        double Y( int i ){ return mNodes[ i-1 ]->Y(); }
        double Z( int i ){ return mNodes[ i-1 ]->Z(); }
        
        void Get_CartesianXYZ_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
        void Get_Axisymmetric_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
        void Get_PlasmaModeJs_Sources( Vector< std::complex<double> >& cJx, Vector< std::complex<double> >& cJy, Vector< std::complex<double> >& cJz );
    }; 
} 

#endif 
