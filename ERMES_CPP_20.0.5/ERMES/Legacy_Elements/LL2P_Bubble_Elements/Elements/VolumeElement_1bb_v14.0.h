
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
            mNumNodes   = 17;
            mNumDofs    = mNumNodes*3; 
            
            ctcNumNodes = mNumNodes;
            ctcNumDofs  = mNumDofs; 

            mFreq = (*mProperties)(FREQUENCY);

		    mPi = 3.1415926535897932384626433832795;
            mEo = 8.8541878176e-12;
			mMo = (4.00e-7) * mPi;   

            mVolume = Calculate_Volume();
            mZero   = std::complex<double>( 0.00, 0.00 );
        }

        void GetStiffMatrix         ( Matrix< std::complex<double> >& StiffMatrix );   
                                    
        void GetEquationIdVector    ( Vector<int>& EquationId );
                                    
        void SetPeso                ( double Peso );
                                  
        void Set_Mesh_Connectivities( std::map< unsigned int, Properties::Pointer >& PropertiesMap    , 
                                      Vector  < Vector<unsigned int>              >& NodesConnectivity, 
                                      Vector  < Vector<unsigned int>              >& ElementsInfo     ,
                                      std::map< unsigned int, char                >& mType_Of_BCNormal );

        void Set_Global_Derivatives ( Vector<        double  >& Element_Volumes, 
                                      Vector< Vector<double> >& Element_dNk_dXs, 
                                      Vector< Vector<double> >& Element_dNk_dYs, 
                                      Vector< Vector<double> >& Element_dNk_dZs );

        void Set_Contact_Info       ( Vector< Vector<Node::Pointer> >& CtC_pNodes, 
                                      Vector< std::complex<double>  >& CtC_Properties );
        
	    void CalculateDerivatives   ( Vector<std::complex<double> >& ncHx, Vector<std::complex<double> >& ncHy, Vector<std::complex<double> >& ncHz );

      private:
          
        std::map< unsigned int, Properties::Pointer >* gpPropertiesMap;

        std::map< unsigned int, char >* gpNormalType;
        
        Vector< Vector<unsigned int> >* mpNodesConnectivity; 
        Vector< Vector<unsigned int> >* mpElementsInfo     ;     
                
        Vector<        double  >* gpVolumes;
        Vector< Vector<double> >* gp_dNk_dX;
        Vector< Vector<double> >* gp_dNk_dY;
        Vector< Vector<double> >* gp_dNk_dZ; 

        Vector< Vector< Node::Pointer > > mCtC_pNodes;
        Vector<  std::complex<double>   > mCtC_Properties;
        Vector< bool > mCTC_FaceUse; 
        Vector< Vector<int> > mCTC_DuplicateFace; 

        std::complex<double> mZero;

        int mThereIsACtCFace;
		int mFaceId;
          
        int mNumNodes;
        int mNumDofs ;

        int ctcNumNodes;
        int ctcNumDofs ;
        
        double mPeso  ;
        double mVolume;

        double mFreq;
		double mEo;
		double mPi;
		double mMo;

        double Calculate_Volume(); 

        double RtoS( double alpha, double beta, double gamma );

	    double X( int i ){ return ( mNodes[i-1]->X() ); }
        double Y( int i ){ return ( mNodes[i-1]->Y() ); }
        double Z( int i ){ return ( mNodes[i-1]->Z() ); }

        void TangencialCoordinates( std::vector<double>& n, std::vector<double>& t, std::vector<double>& b );

        void GetStiffMatrix_PLL2P_SimpTau( Matrix< std::complex<double> >& StiffMatrix );   
        void GetStiffMatrix_InHomogeneous( Matrix< std::complex<double> >& StiffMatrix ); 
        void GetStiffMatrix_Plain_Regular( Matrix< std::complex<double> >& StiffMatrix ); 

        void GetStiffMatrix_LL2P_CurlDivE    ( Matrix< std::complex<double> >& StiffMatrix ); 
        //void GetStiffMatrix_LL2P_CurlDivE_ctc( Matrix< std::complex<double> >& StiffMatrix ); 
        void GetStiffMatrix_LL2P_CurlDivE_CtC( Matrix< std::complex<double> >& StiffMatrix ); 

        void GetStiffMatrix_LL2P_VolumAvg( Matrix< std::complex<double> >& StiffMatrix ); 
        void GetStiffMatrix_LL2P_CurlDivQ( Matrix< std::complex<double> >& StiffMatrix ); 
        void GetStiffMatrix_LL2P_EdgeElem( Matrix< std::complex<double> >& StiffMatrix ); 

        void GetStiffMatrix_LL2P_CtC     ( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_hk2CurlDiv_Stabilization( Matrix< std::complex<double> >& StiffMatrix ); 

        void Add_CurGBk_Stabilization   ( Matrix< std::complex<double> >& Sh_cur );  
        void Add_DivGBk_Stabilization   ( Matrix< std::complex<double> >& Sh_div );
        void Add_DivCurGBk_Stabilization( Matrix< std::complex<double> >& Sh_div , Matrix< std::complex<double> >& Sh_cur );
        void Add_CtCDiv_Stabilization   ( Matrix< std::complex<double> >& CtC_Sh_div ); 
        void Add_nxnCur_Stabilization   ( Matrix< std::complex<double> >& nxn_Sh_cur ); 

        double Calculate_Element_Diameter(); 

        void EdgeElementBase_3D( Matrix<double>&   Lx, Matrix<double>&   Ly, Matrix<double>&   Lz, Matrix<double>& N, 
                                 Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz                   ); 

        void EdgeElementBase_2D( Matrix<double>& Lx, Matrix<double>& Ly, Matrix<double>& Lz, Matrix<double>&  N, Vector<int>& face ); 

        void Calculate_EdgeElem_Matrices( Matrix< std::complex<double> >& C_cur, Matrix< std::complex<double> >& iA_cur ); 

        void EdgeElementBase_Derivatives_3D( Matrix<double>& dLxdx, Matrix<double>& dLxdy, Matrix<double>& dLxdz,
		                                     Matrix<double>& dLydx, Matrix<double>& dLydy, Matrix<double>& dLydz,
		                                     Matrix<double>& dLzdx, Matrix<double>& dLzdy, Matrix<double>& dLzdz,
		                                     Matrix<double>& N, 
		                                     Matrix<double>& dNdx, Matrix<double>& dNdy, Matrix<double>& dNdz   ); 

        void Invert_Matrix ( Matrix<double>& a );

        int    svdcmp      (double **a, int nRows, int nCols, double *w, double **v);
        double pythag      (double   a, double b);
        void   printMatrix (double **a, int rows, int cols); 
        void   printVector (double  *v, int size);
        
        int  NumberOfNeigbours        ( int locNodeId_i, int locNodeId_j );

        void Get_k_Nodes_List         ( int locNodeId_i, int locNodeId_j, Vector<unsigned int>& k_Nodes_List );

        void Get_Node_Connectivity    ( int locNodeId  , Vector<unsigned int>& Node_Connectivity );

        void Calculate_Ni_Derivatives ( Matrix<double>& DN );
                                 
        void DerivativesInNodes1bb    ( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz );
        void DerivativesInNodes2bb    ( Matrix<double>& ndNdx, Matrix<double>& ndNdy, Matrix<double>& ndNdz );

	    void LagrangeDerivatives3D_1bb( Matrix<double>&  dNdx, Matrix<double>&  dNdy, Matrix<double>&  dNdz, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ );
        void LagrangeDerivatives3D_2bb( Matrix<double>&  dNdx, Matrix<double>&  dNdy, Matrix<double>&  dNdz, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ );
        void LagrangeDerivatives3D_1st( Matrix<double>&  dNdx, Matrix<double>&  dNdy, Matrix<double>&  dNdz, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ );
        void LagrangeDerivatives3D_2nd( Matrix<double>&  dNdx, Matrix<double>&  dNdy, Matrix<double>&  dNdz, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ );
        void LagrangeDerivatives3D_3th( Matrix<double>&  dNdx, Matrix<double>&  dNdy, Matrix<double>&  dNdz, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ );
        void LagrangeDerivatives3D_4th( Matrix<double>&  dNdx, Matrix<double>&  dNdy, Matrix<double>&  dNdz, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ );
     
        void Calculate_NiNj_matrix    ( Matrix< std::complex<double> >& w2_Ni_ep_Nj );
        void Add_C_NxN_Matrix         ( Matrix< std::complex<double> >& C_cur       );
        void Add_C_NxN_Matrix_ctc     ( Matrix< std::complex<double> >& C_cur       );
        void Add_C_NxN_Matrix_edge    ( Matrix< std::complex<double> >& C_cur       );
        void Add_D_NxN_Matrix         ( Matrix< std::complex<double> >& D_div       );
        void Add_CtC_D_NxN_Matrix     ( Matrix< std::complex<double> >& D_div       ); 

        void   Calculate_ExtNormal( Vector<double>& extN, Vector<int>& face );
        double Calculate_Area     ( Vector<int>& face ); 

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