
#if !defined(KRATOS_COLD_PLASMA)
#define KRATOS_COLD_PLASMA

#include "../includes/kratos_string.h"
#include "../includes/node.h"
#include "../vectorial_spaces/vector.h"
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"

namespace Kratos
{
	class ColdPlasma 
	{
		public:

			typedef boost::shared_ptr<ColdPlasma> Pointer;

			ColdPlasma ();
			~ColdPlasma();

			int Get_Material_Id();

			void Load_Data( String data_file_name );

			double         Get_ElectronDensity_InNode         ( Node::Pointer pNode );
			Vector<double> Get_Bexternal_InNode               ( Node::Pointer pNode );
			Vector<double> Get_PlasmaPermittivity_SDPRL_InNode( Node::Pointer pNode, double wFreq );
			Vector<double> Get_AllPlasmaParameters_InNode     ( Node::Pointer pNode, double wFreq );

            bool Is_dEdx_zero() { return mCSym_dEdx; }
            bool Is_dEdy_zero() { return mCSym_dEdy; }
            bool Is_dEdz_zero() { return mCSym_dEdz; }

            bool Is_Matrix_Symmetric () { return mMatrix_Is_Symmetric; }

            double Get_Kij_Tol () { return mKij_Tol; }

            char Far_Field_Wave_Type() { return mFar_field_wave; }

            String Get_Eparallel_Tolerance() { return mEpar_tolerance; }

		private:

			// Cold plasma material Id
			int mMaterialId;

			// Geometry type
			bool mFlatGeometry;
			bool mCurvedGeometry;

			// Origin of the Tokamak cartesian coordinate system
			Vector<double> mTK_OO; 

			// Unit vectors of the Tokamak cartesian coordinate system
			Vector<double> mTK_n, mTK_t, mTK_b;

			// Vertical radius of curvature
			double mPoloidalCurv;

			// Hydrogen ion isotopes relative density
			double mH1_div_ne, mH2_div_ne, mH3_div_ne;

            // Helium ion isotopes relative density
            double mHe4_q1_div_ne, mHe4_q2_div_ne;

            // Berylium ion isotopes relative density
            double mBe9_q1_div_ne, mBe9_q2_div_ne;

            // Wave type (R, L, P) applied to the far field boundary condition 
            char mFar_field_wave;

            // E parallel tolerance
            String mEpar_tolerance;

            // Cylindrical symmetry (dE/dx, dE/dy, dE/dz = 0 in the plasma volume)
            bool mCSym_dEdx, mCSym_dEdy, mCSym_dEdz;

            // Format of the sistem matrix
            bool mMatrix_Is_Symmetric;

            // Tolerance value for the matrix elements
            double mKij_Tol;

			// Starting point of the Ne, aB measurement line
			double mR_start, mZ_start;

			// Measurement points
			Vector<double> mLo;
			
			// Measured electron density 
			Vector<double> mNe;

			// Measured applied B field
			Vector<double> mB_mod, mB_ang;

			// Private functions
			void Set_Tokamak_ntb();
			
			void Add_PoloidalCurv_to_TKBfield( Node::Pointer pNode, Vector<double>& B_ntb );

			void Send_Error_Msg( String ErrorMsg, String CheckMsg );

			void Read_eDensity_File( String FileName );
			void Read_aBfield_File ( String FileName );

            double Node_XYZ_to_incR( Node::Pointer pNode );
            
			double         Get_ElectronDensity_FromNe( double incR, Vector<int>& index );
			Vector<double> Get_Bexternal_FromBmodBang( double incR, Vector<int>& index, Node::Pointer pNode );
			Vector<double> Calculate_PlasmaPermittivity_SDPRL( double eDensity, Vector<double>& Bext, double wFreq );

			Vector<int> Get_Indexes_In_Lo ( double incR );

			Vector<double> Node_XYZ_to_TK_XYZ( Node::Pointer pNode );
			Vector<double> TK_Bfield_to_XYZ  ( Node::Pointer pNode, double modB, double angB);
	};
}

#endif