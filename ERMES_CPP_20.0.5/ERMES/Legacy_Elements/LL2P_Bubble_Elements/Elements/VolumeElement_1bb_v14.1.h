
#if !defined(KRATOS_VOLUME_ELEMENT_1BB)
#define KRATOS_VOLUME_ELEMENT_1BB

#include "../includes/element.h"
#include <list>

namespace Kratos
{
    class VolumeElement_1bb : public Element
    {
      public:

        VolumeElement_1bb( NodesArrayType& rNodes, Properties::Pointer rProperties ) : Element( rNodes, rProperties )
        {
            mNumNodes = 5;
            mNumDofs  = mNumNodes*3; 
            mFreq     = (*mProperties)(FREQUENCY);
		    mEo       = 8.8541878176e-12;
		    mPi       = 3.1415926535897932384626433832795;
			mMo       = (4.00e-7) * mPi;      
            mVolume   = Calculate_Volume();
            mZero     = std::complex<double>( 0.00, 0.00 );
        }

        void GetStiffMatrix         ( Matrix< std::complex<double> >& StiffMatrix );   
                                    
        void GetEquationIdVector    ( Vector<int>& EquationId );
                                    
        void SetPeso                ( double Peso );
                                  
        void Set_Mesh_Connectivities( std::map< unsigned int, Properties::Pointer >& PropertiesMap    , 
                                      Vector  < Vector<unsigned int>              >& NodesConnectivity, 
                                      Vector  < Vector<unsigned int>              >& ElementsInfo     ,
                                      std::map< unsigned int, char                >& mType_Of_BCNormal );

       void Set_Global_Derivatives  ( Vector<        double  >& Element_Volumes, 
                                      Vector< Vector<double> >& Element_dNk_dXs, 
                                      Vector< Vector<double> >& Element_dNk_dYs, 
                                      Vector< Vector<double> >& Element_dNk_dZs );

	    void CalculateDerivatives   ( Vector< std::complex<double> >& ncHx, 
                                      Vector< std::complex<double> >& ncHy, 
                                      Vector< std::complex<double> >& ncHz );

      private:
          
        std::map< unsigned int, Properties::Pointer >* gpPropertiesMap;

        std::map< unsigned int, char >* gpNormalType;
        
        Vector< Vector<unsigned int> >* mpNodesConnectivity; 
        Vector< Vector<unsigned int> >* mpElementsInfo     ;     
                
        Vector<        double  >* gpVolumes;
        Vector< Vector<double> >* gp_dNk_dX;
        Vector< Vector<double> >* gp_dNk_dY;
        Vector< Vector<double> >* gp_dNk_dZ; 

        std::complex<double> mZero;
          
        int mNumNodes;
        int mNumDofs ;
        
        double mPeso  ;
        double mVolume;

        double mFreq;
		double mEo;
		double mPi;
		double mMo;

        double Calculate_Volume(); 

        double Calculate_Element_Diameter(); 

        double RtoS( double alpha, double beta, double gamma );

	    double X( int i ){ return ( mNodes[i-1]->X() ); }
        double Y( int i ){ return ( mNodes[i-1]->Y() ); }
        double Z( int i ){ return ( mNodes[i-1]->Z() ); }

        void GetStiffMatrix_Plain_Regular( Matrix< std::complex<double> >& StiffMatrix ); 
        
        void GetStiffMatrix_PLL2P_SimpTau( Matrix< std::complex<double> >& StiffMatrix );   
        void GetStiffMatrix_PLL2P_FullTau( Matrix< std::complex<double> >& StiffMatrix );

        void GetStiffMatrix_LL2P_VolumAvg( Matrix< std::complex<double> >& StiffMatrix ); 
        void GetStiffMatrix_LL2P_CurlDivE( Matrix< std::complex<double> >& StiffMatrix ); 
        void GetStiffMatrix_LL2P_CurlDivQ( Matrix< std::complex<double> >& StiffMatrix ); 
        void GetStiffMatrix_LL2P_EdgeElem( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_hk2CurlDiv_Stabilization( Matrix< std::complex<double> >& StiffMatrix ); 
        void Add_CurlDivGBk_Stabilization( Matrix< std::complex<double> >& StiffMatrix ); 

        void EdgeElementBase_3D( Matrix<double>&   Lx, Matrix<double>&   Ly, Matrix<double>&   Lz, Matrix<double>& N, 
                                 Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz                   ); 

        void Calculate_EdgeElem_Matrices( Matrix< std::complex<double> >& C_cur, Matrix< std::complex<double> >& iA_cur ); 

        void EdgeElementBase_Derivatives_3D( Matrix<double>& dLxdx, Matrix<double>& dLxdy, Matrix<double>& dLxdz,
		                                     Matrix<double>& dLydx, Matrix<double>& dLydy, Matrix<double>& dLydz,
		                                     Matrix<double>& dLzdx, Matrix<double>& dLzdy, Matrix<double>& dLzdz,
		                                     Matrix<double>& N, 
		                                     Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz   ); 

        void Invert_Matrix    ( Matrix<double>& a );

        int  NumberOfNeigbours( int i, int j );

        void Get_k_Nodes_List ( int i, int j, Vector<unsigned int>& k_Nodes_List );

        void Calculate_Ni_Derivatives ( Matrix<double>& DN );
                                                      
        void DerivativesInNodes       ( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz );

	    void LagrangeDerivatives3D_1bb( Matrix     <double>& dNdx, Matrix     <double>& dNdy, Matrix     <double>& dNdz, 
                                        std::vector<double>&   cX, std::vector<double>&   cY, std::vector<double>&   cZ );

        void Calculate_NiNj_matrix    ( Matrix< std::complex<double> >& w2_Ni_ep_Nj );
        void Add_C_NxN_Matrix         ( Matrix< std::complex<double> >& C_cur       );
        void Add_D_NxN_Matrix         ( Matrix< std::complex<double> >& D_div       );

        void   Calculate_ExtNormal    ( Vector<double>& extN, Vector<int>& face );
        double Calculate_Area         ( Vector<int>& face ); 

        void Get_Elements_List        ( Vector  < unsigned int >& Elements_List );

        void Get_Elements_Properties  ( Vector  < unsigned int >& Elements_List, 
                                        std::map< unsigned int, std::complex<double> >& eEp,
		                                std::map< unsigned int, std::complex<double> >& eMu );

        void Get_alpha_tbeta_kk       ( std::map< unsigned int, std::complex<double> >& alpha_kk,
									    std::map< unsigned int, std::complex<double> >& tbeta_kk );

        void Calculate_L2P_CurlDiv_Terms    ( Vector< std::map< unsigned int, double> >& Ci_dNkdX, 
		                                      Vector< std::map< unsigned int, double> >& Ci_dNkdY,    
		                                      Vector< std::map< unsigned int, double> >& Ci_dNkdZ );
                                        
        void Calculate_L2P_Curl_ep_Div_Terms( Vector< std::map< unsigned int, std::complex<double> > >& Ci_ep_dNkdX, 
		                                      Vector< std::map< unsigned int, std::complex<double> > >& Ci_ep_dNkdY,    
		                                      Vector< std::map< unsigned int, std::complex<double> > >& Ci_ep_dNkdZ );

    }; 
} 

#endif