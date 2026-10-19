
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
            mNumDofs  = mNumNodes * 3;
            mVolume   = Calculate_Volume();
        }

        void SetNormalType( std::map< unsigned int, char >& mType_Of_BCNormal )
        {
	        gpNormalType = &mType_Of_BCNormal;
	    }

        void SetFactor( double Factor )
	    { 
	        mhkFactor = Factor;
	    }

        int Get_InnerGiD_GaussPoints( std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, int numResultsOnGPs );

        void GetStiffMatrix( Matrix< std::complex<double> >& StiffMatrix );    
                                 
        void GetEquationIdVector( Vector<int>& EquationId ); 
                                 
	    void Calculate_Rotational_OnNodes( cVector2Type& rot_OnNodes );

        void Calculate_Rotational_OnGaussPoints( cVector2Type& rot_OnGPs, int numResultsOnGPs );

        void Calculate_E_field_OnGaussPoints( cVector2Type& Ef_OnGPs, int numResultsOnGPs );

      private:
          
        int mNumNodes;
        int mNumDofs;

		double mVolume;
        double mhkFactor;

        std::map< unsigned int, char >* gpNormalType;
        
        double X( int i ){ return mNodes[i-1]->X(); }
        double Y( int i ){ return mNodes[i-1]->Y(); }
        double Z( int i ){ return mNodes[i-1]->Z(); }

        double Calculate_Volume(); 

        double Calculate_Area( Vector<int>& face ); 

        double Calculate_Element_Diameter(); 
       
        void Calculate_ExtNormal( Vector<double>& extN, Vector<int>& face );
        
        void GetStiffMatrix_LL2P_CurlDivE( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_hk2CurlDiv_Stabilization( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_C_NxN_Matrix( Matrix< std::complex<double> >& C_cur ); 

        void Add_D_NxN_Matrix( Matrix< std::complex<double> >& D_div ); 
        
        void Calculate_Ni_Derivatives( Matrix<double>& DN );

        void Calculate_Rot_P1_Projection( cVectorType& cRotP1x, cVectorType& cRotP1y, cVectorType& cRotP1z );

	    void LagrangeDerivatives3D_3sb( Matrix<double>&    dNdx, 
                                        Matrix<double>&    dNdy, 
                                        Matrix<double>&    dNdz, 
                                        std::vector<double>& cX, 
                                        std::vector<double>& cY, 
                                        std::vector<double>& cZ );
    }; 
} 

#endif