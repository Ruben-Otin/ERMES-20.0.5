
#if !defined(KRATOS_COLD_PLASMA)
#define KRATOS_COLD_PLASMA

#include "../includes/kratos_string.h"
#include "../includes/node.h"
#include "../includes/vector.h"
#include <memory>

namespace Kratos
{
	class ColdPlasma 
	{
	  public:
      
	    typedef std::shared_ptr<ColdPlasma> Pointer;

        typedef Vector<std::complex<double> > ComplexVector;

        ColdPlasma( String base_file_name )
        {
            mBase_FileName = base_file_name;

            mNumSpecies = 0;
        };

	  	~ColdPlasma();

	    int Get_Material_Id()
	    {
		    return mMaterialId;
	    }

        double Get_Kij_Tol() 
        { 
            return mKij_Tol; 
        }

        double Get_Epar_Tol() 
        { 
            return mEpar_Tol; 
        }

        bool Is_Epar_Tol_Off() 
        { 
            if( mEpar_Tol == 0.0 )
            {
                return true;
            }
            else
            {
                return false; 
            }
        }

        bool Is_Full_Matrix() 
        { 
            if( mMatrixStorageFormat == 0 )
            {
                return true;
            }
            else
            {
                return false; 
            }
        }

        bool Is_HermFull_Matrix() 
        { 
            if( mMatrixStorageFormat == 1 )
            {
                return true;
            }
            else
            {
                return false; 
            }
        }

        bool Is_HermSymm_Matrix() 
        { 
            if( mMatrixStorageFormat == 2 )
            {
                return true;
            }
            else
            {
                return false; 
            }
        }

	  	void Load_Data( String data_file_name );

        void Resize_Profile_3D_Vectors( unsigned int NewSize );

        void Linear_Profile_Interpolate( Vector<int> HONodesId );

        void Get_Bext_SDPRL_OnNode( Vector<double>& Bext, ComplexVector& SDPRL, Node::Pointer pNode, double wFreq );

	  	void Get_PermittivityTensor_OnNode( Matrix<std::complex<double> >& TEn, Node::Pointer pNode, double wFreq );

	  	void Get_ConductivityTensor_OnNode( Matrix<std::complex<double> >& TCn, Node::Pointer pNode, double wFreq );

        void Get_AllPlasmaParameters_OnNode( double& eDensity, Vector<double>& Bext, ComplexVector& SDPRL, Node::Pointer pNode, double wFreq );
      
	  	double Get_eDensity_OnNode( Node::Pointer pNode );

	  	Vector<double> Get_Bext_OnNode( Node::Pointer pNode );

        ComplexVector Get_SDPRL_OnNode( Node::Pointer pNode, double wFreq );

	  private:

	    // Cold plasma material Id
	    int mMaterialId;
        
	    // Plasma geometry extrusion type:
        // 0 = Flat-1D
        // 1 = Curved
        // 2 = Full-3D
        // 3 = Tensor
        int mGeometryExtType;

	    // Tokamak cartesian coordinate system
	    Vector<double> mTK_O; 
	    Vector<double> mTK_n;
        Vector<double> mTK_t;
        Vector<double> mTK_b;

        // First point of the 1D measurement line
        Vector<double> mL1_O; 
        
	    // Poloidal curvatures
	    double mPoloidalCurv_ne;
        double mPoloidalCurv_Be;
        
	    // Hydrogen ion isotopes relative density
	    double mH1_div_ne; 
        double mH2_div_ne; 
        double mH3_div_ne;
        
        // Helium ion isotopes relative density
        double mHe3_q1_div_ne; 
        double mHe3_q2_div_ne; 
        double mHe4_q1_div_ne; 
        double mHe4_q2_div_ne;
        
        // Berylium ion isotopes relative density
        double mBe9_q1_div_ne; 
        double mBe9_q2_div_ne; 
        double mBe9_q3_div_ne; 
        double mBe9_q4_div_ne;
        
        // Heavy ion isotopes relative density
        double mNe20_q10_div_ne; 
        double mNi60_q20_div_ne;
        double mW184_q25_div_ne;
        
        // Matrix storage format: 
        // 0 = Full_matrix
        // 1 = Herm-Full
        // 2 = Herm-Symm 
        int mMatrixStorageFormat;

        // Damping mode:
        // 0 = No_damping 
        // 1 = Background_ctc
        // 2 = Background_den
        // 3 = ComplexFreq_all
        // 4 = ComplexFreq_ele 
        int mDamping_mode;

        // Damping value
        double mDamping_value;

        // Minimum value for matrix elements
        double mKij_Tol;

        // E parallel tolerance
        double mEpar_Tol;

	    // 1D line start point 
	    double mR_start; 
        double mZ_start;

        // Problem file name
        String mBase_FileName;

	    // 1D measurement points of the electron density file (strictly ascending)
	    Vector<double> mLo;

	    // 1D measurement points of the B field file (strictly ascending, may differ from mLo)
	    Vector<double> mLo_B;

	    // Measured electron density 
	    Vector<double> mNe;
        
	    // 1D measured B field
	    Vector<double> mB_mod, mB_ang;

	    // 3D measured B field
	    Vector<double> mBx, mBy, mBz;

        // Imported TE tensor
        Vector< Matrix<std::complex<double> > > mImp_K_Tensor; 

        // Species table (charge [C], mass [kg], density relative to ne), built once in Load_Data().
        // Only species with non-zero density are stored, so absent species cost nothing in the
        // per-node permittivity evaluation and cannot trigger spurious resonance errors.
        static const int mMaxSpecies = 15;
        int    mNumSpecies;
        double mSp_Q    [ mMaxSpecies ];
        double mSp_M    [ mMaxSpecies ];
        double mSp_Ratio[ mMaxSpecies ];

	    // Private functions:
	    void Set_Tokamak_ntb();

        void Set_Species_Table();

        String Read_Path_Line( std::istream& rStream );

        void Check_Matrix_Storage_vs_Damping();

        void Check_1D_Grid( Vector<double>& Grid, String FileName );

        bool Is_Data_Line( const std::string& FileLine );

        void Set_3D_Vector_Ne_Bxyz();

        void Read_eDensity_File_1D( String FileName );
	    void Read_sBfields_File_1D( String FileName );
            
        void Read_eDensity_File_3D( String FileName );
        void Read_sBfields_File_3D( String FileName );

        void Read_K_Tensor_File_3D( String FileName );

	    void Send_Error_Msg( String ErrorMsg, String CheckMsg );

        // Same as Send_Error_Msg for messages built with std::string concatenation
	    void Send_Error_Msg_Str( const std::string& ErrorMsg, const std::string& CheckMsg )
        {
            Send_Error_Msg( String( ErrorMsg.c_str() ), String( CheckMsg.c_str() ) );
        }

        void Add_PoloidalCurv_to_TKBfield( Node::Pointer pNode, Vector<double>& B_ntb );

        void Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz );

        void Get_eDensity_Bext_OnNode( double& eDensity, Vector<double>& Bext, Node::Pointer pNode );

        double Get_eDensity_1D( double incR, Vector<int>& index );
        
        double Node_XYZ_to_incR( Node::Pointer pNode, double PoloidalCurv );

        Vector<int> Get_Indexes_In_Lo( double incR );

        Vector<int> Get_Indexes_In_Grid( Vector<double>& Grid, double x );
        
        Vector<double> Node_XYZ_to_TK_XYZ( Node::Pointer pNode );

        Vector<double> TK_Bfield_to_XYZ( Node::Pointer pNode, double modB, double angB );

        Vector<double> Get_Bext_1D( double incR, Vector<int>& index, Node::Pointer pNode );

        ComplexVector Calculate_PlasmaPermittivity_SDPRL( double eDensity, Vector<double>& Bext, double wFreq );
	};
}

#endif