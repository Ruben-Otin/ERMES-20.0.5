
#include "ColdPlasma.h"
#include "Utils/PhysicalConstants.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cctype>
#include <algorithm>

namespace Kratos
{
	//***************************************************************************************************************************
	//* - Load plasma data from files.
	//***************************************************************************************************************************
	//***************************************************************************************************************************
	//* - Next non-empty line of the plasma data file, without leading/trailing blanks ( file paths ).
	//***************************************************************************************************************************
	String ColdPlasma::Read_Path_Line( std::istream& rStream )
	{
		std::string line;

		while( std::getline( rStream, line ) )
		{
			std::size_t b = line.find_first_not_of( " \t\r" );

			if( b == std::string::npos ) continue;

			std::size_t e = line.find_last_not_of( " \t\r" );

			line = line.substr( b, e - b + 1 );

			return String( line.c_str() );
		}

		return String( "" );
	}

	//***************************************************************************************************************************
	//* - Load plasma data from files ( written by 12-Cold_Plasma.bas ).
	//***************************************************************************************************************************
	void ColdPlasma::Load_Data( String data_file_name )
	{
		// Open cold plasma properties data file in reading mode
		std::fstream dataFile( data_file_name, std::fstream::in );

		if( !dataFile.is_open() )
		{
			Send_Error_Msg_Str( "Error opening cold plasma data file " + data_file_name, "Please, check the plasma input files." );
		}

		// Plasma material Id
		dataFile >> mMaterialId;

		if( dataFile.eof() )
		{
			Send_Error_Msg( "No cold plasma applied", "Please, apply plasma material to volumes." );
		}

		// Plasma geometry type:
		// 0 = Flat1D
		// 1 = Curved
		// 2 = Full3D 
		// 3 = Tensor
		dataFile >> mGeometryExtType;

		// Poloidal curvatures
		dataFile >> mPoloidalCurv_ne;
		dataFile >> mPoloidalCurv_Be;

		// First point of the 1D measurement line
		mL1_O.resize( 3 );

		dataFile >> mL1_O[ 0 ];
		dataFile >> mL1_O[ 1 ];
		dataFile >> mL1_O[ 2 ];

		// Coordinates of the center of the tokamak
		mTK_O.resize( 3 );	

		dataFile >> mTK_O[ 0 ];
		dataFile >> mTK_O[ 1 ];
		dataFile >> mTK_O[ 2 ];

		// Central axis of the tokamak
		mTK_n.resize( 3 );

		dataFile >> mTK_n[ 0 ];
		dataFile >> mTK_n[ 1 ];
		dataFile >> mTK_n[ 2 ];

		// Electron density data file
		String eDensityFile;

		// File paths are read as whole lines: they may contain spaces ( e.g. Windows folders )
		eDensityFile = Read_Path_Line( dataFile );

		// External static B fields data file
		String sBfieldsFile;

		sBfieldsFile = Read_Path_Line( dataFile );

		// Densities relative to the electron density
		dataFile >> mH1_div_ne;
        dataFile >> mH2_div_ne;
        dataFile >> mH3_div_ne;

        dataFile >> mHe3_q1_div_ne;
        dataFile >> mHe3_q2_div_ne;

        dataFile >> mHe4_q1_div_ne;
        dataFile >> mHe4_q2_div_ne;

        dataFile >> mBe9_q1_div_ne;
        dataFile >> mBe9_q2_div_ne;
        dataFile >> mBe9_q3_div_ne;
        dataFile >> mBe9_q4_div_ne;

        dataFile >> mNe20_q10_div_ne; 
        dataFile >> mNi60_q20_div_ne; 
        dataFile >> mW184_q25_div_ne;

        // Matrix storage format: 
		// 0 = Full_matrix
		// 1 = Herm-Full
		// 2 = Herm-Symm 
        dataFile >> mMatrixStorageFormat;

        // Elements matrix tolerance
        dataFile >> mKij_Tol;

        // E parallel tolerance
        dataFile >> mEpar_Tol;

        // Damping mode:
		// 0 = No_damping
		// 1 = Background_ctc
		// 2 = Background_den
		// 3 = ComplexFreq_all
		// 4 = ComplexFreq_ele 
        dataFile >> mDamping_mode;

		// Damping value
		dataFile >> mDamping_value;

		// Closing cold plasma data file
        dataFile.close();

        // Build the species table once ( it is used in every permittivity evaluation )
        Set_Species_Table();

        // Hermitian matrix storage is only valid for a lossless ( Hermitian ) plasma tensor
        Check_Matrix_Storage_vs_Damping();

		// Complete cold plasma set-up
		if( mGeometryExtType < 2 )
		{
			Set_Tokamak_ntb();	
			Read_eDensity_File_1D( eDensityFile );
		    Read_sBfields_File_1D( sBfieldsFile );

			// lower_bound() based interpolation needs strictly ascending grids
			Check_1D_Grid( mLo  , eDensityFile );
			Check_1D_Grid( mLo_B, sBfieldsFile );
		}
		else
		{
			Set_3D_Vector_Ne_Bxyz();
			Read_eDensity_File_3D( eDensityFile );
		    Read_sBfields_File_3D( sBfieldsFile );		
		}

		if( mGeometryExtType == 3 )
		{
		    Read_K_Tensor_File_3D( eDensityFile );
		}
	}

	//***************************************************************************************************************************
	//* - Resize 3D vectors plasma profiles.
	//***************************************************************************************************************************
	void ColdPlasma::Resize_Profile_3D_Vectors( unsigned int NewSize )
	{
	    mNe.resize( NewSize, 0.0 );
		mBx.resize( NewSize, 0.0 );
		mBy.resize( NewSize, 0.0 );
		mBz.resize( NewSize, 0.0 );

		if( mGeometryExtType == 3 ) 
		{
			int OldSize = mImp_K_Tensor.size();
			
			mImp_K_Tensor.resize( NewSize );

			for( unsigned int ki=OldSize; ki<NewSize; ki++ ) 
		    { 
			    mImp_K_Tensor[ ki ].Resize( 3, 3, std::complex<double>( 0.0, 0.0 ) );
		    }
		}
	}

	//***************************************************************************************************************************
	//* - Tokamak coordinates system with origen mTK_O and main axis mTK_n.
	//***************************************************************************************************************************
	void ColdPlasma::Set_Tokamak_ntb()
	{
		// Tokamak central axis unitary vector n
		double n_norm = std::sqrt( mTK_n[ 0 ]*mTK_n[ 0 ] + mTK_n[ 1 ]*mTK_n[ 1 ] + mTK_n[ 2 ]*mTK_n[ 2 ] );

		if( n_norm > 0.0 ) 
		{
			mTK_n *= ( 1.0 / n_norm );
		}
		else 
		{
			Send_Error_Msg( "Wrong tokamak vertical axis", "Please, check plasma input data." );
		}

	    // Position vector of 1D line start point from mTK_O
		Vector<double> R1D( 3 );

		R1D[ 0 ] = mL1_O[ 0 ] - mTK_O[ 0 ];
		R1D[ 1 ] = mL1_O[ 1 ] - mTK_O[ 1 ];
		R1D[ 2 ] = mL1_O[ 2 ] - mTK_O[ 2 ];

		// Unitary vector b = ( n x R1D ) / || n x R1D ||
		mTK_b.resize( 3 );

		mTK_b[ 0 ] = mTK_n[ 1 ] * R1D[ 2 ] - mTK_n[ 2 ] * R1D[ 1 ];
		mTK_b[ 1 ] = mTK_n[ 2 ] * R1D[ 0 ] - mTK_n[ 0 ] * R1D[ 2 ];
		mTK_b[ 2 ] = mTK_n[ 0 ] * R1D[ 1 ] - mTK_n[ 1 ] * R1D[ 0 ];

		double b_norm = std::sqrt( mTK_b[ 0 ]*mTK_b[ 0 ] + mTK_b[ 1 ]*mTK_b[ 1 ] + mTK_b[ 2 ]*mTK_b[ 2 ] );

		if( b_norm > 0.0 ) 
		{
			mTK_b *= ( 1.0 / b_norm );
		}
		else 
		{
			Send_Error_Msg( "Wrong L1_O_XYZ point coordinates", "Please, check plasma input data." );
		}

		// Unitary vector t = ( b x n ) / || b x n ||
		mTK_t.resize( 3 );

		mTK_t[ 0 ] = mTK_b[ 1 ] * mTK_n[ 2 ] - mTK_b[ 2 ] * mTK_n[ 1 ];
		mTK_t[ 1 ] = mTK_b[ 2 ] * mTK_n[ 0 ] - mTK_b[ 0 ] * mTK_n[ 2 ];
		mTK_t[ 2 ] = mTK_b[ 0 ] * mTK_n[ 1 ] - mTK_b[ 1 ] * mTK_n[ 0 ];

		mTK_t *= ( 1.0 / std::sqrt( mTK_t[ 0 ]*mTK_t[ 0 ] + mTK_t[ 1 ]*mTK_t[ 1 ] + mTK_t[ 2 ]*mTK_t[ 2 ] ) );

		// 1D line start point in cylindrical mTK system
		Vector<double> RZ_O( 3 );

		RZ_O[ 0 ] = R1D[ 0 ] * mTK_t[ 0 ] + R1D[ 1 ] * mTK_t[ 1 ] + R1D[ 2 ] * mTK_t[ 2 ];
		RZ_O[ 1 ] = R1D[ 0 ] * mTK_b[ 0 ] + R1D[ 1 ] * mTK_b[ 1 ] + R1D[ 2 ] * mTK_b[ 2 ];
		RZ_O[ 2 ] = R1D[ 0 ] * mTK_n[ 0 ] + R1D[ 1 ] * mTK_n[ 1 ] + R1D[ 2 ] * mTK_n[ 2 ];

		mR_start = std::sqrt( RZ_O[ 0 ]*RZ_O[ 0 ] + RZ_O[ 1 ]*RZ_O[ 1 ]);

		mZ_start = RZ_O[ 2 ];
	}

	//***************************************************************************************************************************
	//* - Transform node cartesian coordinates to tokamak coordinates.
	//***************************************************************************************************************************
	Vector<double> ColdPlasma::Node_XYZ_to_TK_XYZ( Node::Pointer pNode )
	{
		// Node cartesian coordinates 
		Vector<double> OG_XYZ( 3 );

		OG_XYZ[ 0 ] = pNode->X();
		OG_XYZ[ 1 ] = pNode->Y();
		OG_XYZ[ 2 ] = pNode->Z();

		// Position vector of node from mTK_O
		Vector<double> Ptk( 3 );

		Ptk[ 0 ] = OG_XYZ[ 0 ] - mTK_O[ 0 ];
		Ptk[ 1 ] = OG_XYZ[ 1 ] - mTK_O[ 1 ];
		Ptk[ 2 ] = OG_XYZ[ 2 ] - mTK_O[ 2 ];

		// TK coordinates of the node
		Vector<double> TK_XYZ( 3 );

		TK_XYZ[ 0 ] = Ptk[ 0 ] * mTK_t[ 0 ] + Ptk[ 1 ] * mTK_t[ 1 ] + Ptk[ 2 ] * mTK_t[ 2 ];
		TK_XYZ[ 1 ] = Ptk[ 0 ] * mTK_b[ 0 ] + Ptk[ 1 ] * mTK_b[ 1 ] + Ptk[ 2 ] * mTK_b[ 2 ];
		TK_XYZ[ 2 ] = Ptk[ 0 ] * mTK_n[ 0 ] + Ptk[ 1 ] * mTK_n[ 1 ] + Ptk[ 2 ] * mTK_n[ 2 ];

		return TK_XYZ;
	}

	//***************************************************************************************************************************
	//* - Node cartesian coordinates to incR ( distance from Lo[ 0 ] ).
	//***************************************************************************************************************************
	double ColdPlasma::Node_XYZ_to_incR( Node::Pointer pNode, double PoloidalCurv )
	{
		// Node coordinates in tokamak coordinate system
		Vector<double> TK_XYZ = Node_XYZ_to_TK_XYZ( pNode );

		// Node distance to the tokamak axis
		double R;
		
		if( mGeometryExtType == 0 )
		{
			R = std::abs( TK_XYZ[ 0 ] );
		}
		else
		{
			R = std::sqrt( TK_XYZ[ 0 ]*TK_XYZ[ 0 ] + TK_XYZ[ 1 ]*TK_XYZ[ 1 ] );
		}

		// Add poloidal curvature
		if( PoloidalCurv > 0.0 )
		{
			double Z = TK_XYZ[ 2 ];

			double Ro = mR_start - PoloidalCurv;
			double Zo = mZ_start;

			double dist = std::sqrt( ( R - Ro )*( R - Ro ) + ( Z - Zo )*( Z - Zo ) );

			R = Ro + dist;
		}

        return ( mR_start - R );
	}

	//***************************************************************************************************************************
	//* - B field in XYZ cartesian coordinate system from [ modB, angB ] 1D format on pNode.
	//***************************************************************************************************************************
	Vector<double> ColdPlasma::TK_Bfield_to_XYZ( Node::Pointer pNode, double modB, double angB )
	{
		// B in ntb coordinates and xyz coordinates
		Vector<double> B_ntb( 3 );
		Vector<double> B_xyz( 3 );

		// Degrees to radians
		double pi       = 3.1415926535898;
		double angB_rad = angB * pi / 180;

		// B field in curved geometry mode
		if( mGeometryExtType == 1 )
		{
			Vector<double> N_tbn = Node_XYZ_to_TK_XYZ( pNode );

			double modN_tb = std::sqrt( N_tbn[ 0 ]*N_tbn[ 0 ] + N_tbn[ 1 ]*N_tbn[ 1 ] );

            double modB_tb = modB * std::cos( angB_rad );

			double inv_mod;

			if( modN_tb != 0.0 ) 
			{
				inv_mod = 1.0 / modN_tb;
			}
			else
			{
				inv_mod = 0.0;
			}

			B_ntb[ 0 ] = + modB_tb * N_tbn[ 1 ] * inv_mod;
			B_ntb[ 1 ] = - modB_tb * N_tbn[ 0 ] * inv_mod;
			B_ntb[ 2 ] = + modB    * std::sin( angB_rad );
		}
		// B field in flat geometry mode
		else
		{
			B_ntb[ 0 ] = + 0.0;
			B_ntb[ 1 ] = - modB * std::cos( angB_rad );
			B_ntb[ 2 ] = + modB * std::sin( angB_rad );
		}

		// Add poloidal curvature to Be field
		if( mPoloidalCurv_Be > 0.0 ) 
		{
			Add_PoloidalCurv_to_TKBfield( pNode, B_ntb );
		}

		// B in ntb CS to xyz CS
		B_xyz[ 0 ] = B_ntb[ 0 ] * mTK_t[ 0 ] + B_ntb[ 1 ] * mTK_b[ 0 ] + B_ntb[ 2 ] * mTK_n[ 0 ];
		B_xyz[ 1 ] = B_ntb[ 0 ] * mTK_t[ 1 ] + B_ntb[ 1 ] * mTK_b[ 1 ] + B_ntb[ 2 ] * mTK_n[ 1 ];
		B_xyz[ 2 ] = B_ntb[ 0 ] * mTK_t[ 2 ] + B_ntb[ 1 ] * mTK_b[ 2 ] + B_ntb[ 2 ] * mTK_n[ 2 ];

		return B_xyz;
	}

	//***************************************************************************************************************************
	//* - Add poloidal curvature to B ext.
	//***************************************************************************************************************************
	void ColdPlasma::Add_PoloidalCurv_to_TKBfield( Node::Pointer pNode, Vector<double>& B_ntb )
	{
		// Node coordinates in TK system
		Vector<double> rNode = Node_XYZ_to_TK_XYZ( pNode );

		// Curvature center coordinates in TK sytem
		Vector<double> rCurC( 3 );

		// Curved geometry mode
		if( mGeometryExtType == 1 )
		{
		    double modN_tb = std::sqrt( rNode[ 0 ]*rNode[ 0 ] + rNode[ 1 ]*rNode[ 1 ] );

		    double inv_modN_tb;

			if( modN_tb != 0.0 ) 
			{
				inv_modN_tb = 1.0 / modN_tb;
			}
			else  
			{
				inv_modN_tb = 0.0;
			}

			rCurC[ 0 ] = ( mR_start - mPoloidalCurv_Be ) * rNode[ 0 ] * inv_modN_tb;
			rCurC[ 1 ] = ( mR_start - mPoloidalCurv_Be ) * rNode[ 1 ] * inv_modN_tb;
			rCurC[ 2 ] = ( mZ_start );
		}
		// Flat geometry mode
		else
		{
			rNode[ 1 ] = 0.0;

			rCurC[ 0 ] = mR_start - mPoloidalCurv_Be;
			rCurC[ 1 ] = 0.0;
			rCurC[ 2 ] = mZ_start;
		}

		// Node position vector from curvature center
		Vector<double> rN_CC( 3 );

		rN_CC[ 0 ] = rNode[ 0 ] - rCurC[ 0 ];
		rN_CC[ 1 ] = rNode[ 1 ] - rCurC[ 1 ];
		rN_CC[ 2 ] = rNode[ 2 ] - rCurC[ 2 ];

		// Vector perpendicular to rN_CC and mTK_n axis ( P_rN_n = rN_CC x mTK_n ) 
		Vector<double> P_rN_n( 3 );

		P_rN_n[ 0 ] = + rN_CC[ 1 ];
		P_rN_n[ 1 ] = - rN_CC[ 0 ];
		P_rN_n[ 2 ] = 0.0;

		// Vector perpendicular to P_rN_n and rN_CC ( P_rN_CC = P_rN_n x rN_CC )
		Vector<double> P_rN_CC( 3 );

		P_rN_CC[ 0 ] = P_rN_n[ 1 ] * rN_CC[ 2 ] - P_rN_n[ 2 ] * rN_CC[ 1 ];
		P_rN_CC[ 1 ] = P_rN_n[ 2 ] * rN_CC[ 0 ] - P_rN_n[ 0 ] * rN_CC[ 2 ];
		P_rN_CC[ 2 ] = P_rN_n[ 0 ] * rN_CC[ 1 ] - P_rN_n[ 1 ] * rN_CC[ 0 ];

		// Change direction after crossing mTK_n axis
		if( rN_CC[ 0 ] < 0 )
		{
		    P_rN_CC *= -1;
		}

		// Change direction after crossing mTK_b axis
		if( rCurC[ 0 ] < 0 )
		{
		    P_rN_CC *= -1;
		}

		// Redefine P_rN_CC if rN_CC is parallel to mTK_n
		if( ( rN_CC[ 0 ] == 0.0 ) && ( rN_CC[ 1 ] == 0.0 ) )
		{
		    P_rN_CC = rCurC;
			
			if( rN_CC[ 2 ] > 0.0 ) 
			{
				P_rN_CC *= -1;
			}
		}
		
		// P_rN_CC vector normalization
		double modPrNCC = std::sqrt( P_rN_CC[ 0 ]*P_rN_CC[ 0 ] + P_rN_CC[ 1 ]*P_rN_CC[ 1 ] + P_rN_CC[ 2 ]*P_rN_CC[ 2 ] );

	    if( modPrNCC > 0.0 ) 
		{
			P_rN_CC *= ( 1.0 / modPrNCC );
		}
		
		// Projection of Bz on the normalized vector P_rN_CC
		B_ntb[ 0 ] += B_ntb[ 2 ] * P_rN_CC[ 0 ];
		B_ntb[ 1 ] += B_ntb[ 2 ] * P_rN_CC[ 1 ];
		B_ntb[ 2 ]  = B_ntb[ 2 ] * P_rN_CC[ 2 ];
	}

	//***************************************************************************************************************************
	//* - Lower and upper indexes in mLo containing incR:
	//*   index_in_Lo[ 0 ] = lower_index
	//*   index_in_Lo[ 1 ] = upper_index
	//***************************************************************************************************************************
	Vector<int> ColdPlasma::Get_Indexes_In_Lo( double incR )
	{
		return Get_Indexes_In_Grid( mLo, incR );
	}

	//***************************************************************************************************************************
	//* - Lower and upper indexes of an ascending grid containing x:
	//*   index[ 0 ] = lower_index
	//*   index[ 1 ] = upper_index
	//*   Values outside the grid are clamped to the first / last point ( index[ 0 ] == index[ 1 ] ).
	//***************************************************************************************************************************
	Vector<int> ColdPlasma::Get_Indexes_In_Grid( Vector<double>& Grid, double x )
	{
		Vector<double>::iterator it_G = std::lower_bound( Grid.begin(), Grid.end(), x );

		Vector<int> index( 2 );

		if( it_G == Grid.begin() )
		{
			index[ 0 ] = 0;
			index[ 1 ] = 0;
		}
		else if( it_G == Grid.end() )
		{
			index[ 0 ] = Grid.size() - 1;
			index[ 1 ] = Grid.size() - 1;
		}
		else
		{
			index[ 0 ] = it_G - Grid.begin() - 1;
			index[ 1 ] = it_G - Grid.begin();
		}

		return index;
	}

	//***************************************************************************************************************************
	//* - Bext in XYZ cartesian components from vectors mB_mod and mB_ang by linear interpolation.
	//***************************************************************************************************************************
	Vector<double> ColdPlasma::Get_Bext_1D( double incR, Vector<int>& index, Node::Pointer pNode )
	{
		// index refers to the B field grid mLo_B ( see Get_Indexes_In_Grid ), which may differ from the density grid
		double modB, angB;

		if( index[ 0 ] == index[ 1 ] )
		{
			modB = mB_mod[ index[ 0 ] ];
			angB = mB_ang[ index[ 0 ] ];
		}
		else
		{
			double alpha = ( incR - mLo_B[ index[ 0 ] ] ) / ( mLo_B[ index[ 1 ] ] - mLo_B[ index[ 0 ] ] );

			modB = alpha * mB_mod[ index[ 1 ] ] + ( 1.0 - alpha ) * mB_mod[ index[ 0 ] ];
			angB = alpha * mB_ang[ index[ 1 ] ] + ( 1.0 - alpha ) * mB_ang[ index[ 0 ] ];
		}

		return TK_Bfield_to_XYZ( pNode, modB, angB );
	}

	//****************************************************************************************************************************
	//* - Electron density from vector mNe by linear interpolation.
	//****************************************************************************************************************************
	double ColdPlasma::Get_eDensity_1D( double incR, Vector<int>& index )
	{
		if( index[ 0 ] == index[ 1 ] )
		{
			return mNe[ index[ 0 ] ];
		}

		double alpha = ( incR - mLo[ index[ 0 ] ] ) / ( mLo[ index[ 1 ] ] - mLo[ index[ 0 ] ] );

		return alpha * mNe[ index[ 1 ] ] + ( 1.0 - alpha ) * mNe[ index[ 0 ] ];
	}

	//***************************************************************************************************************************
	//* - eDensity on pNode.
	//***************************************************************************************************************************
	double ColdPlasma::Get_eDensity_OnNode( Node::Pointer pNode )
	{
		double eDensity = 0.0;

	    if( mGeometryExtType < 2 )
	    {
		    double      incRe = Node_XYZ_to_incR ( pNode, mPoloidalCurv_ne );
		    Vector<int> indxe = Get_Indexes_In_Lo( incRe );
		    
		    eDensity = Get_eDensity_1D( incRe, indxe );
		}
		else
		{
			int indxp = pNode->Id() - 1;

			if( ( indxp >= 0 ) && ( indxp < (int)mNe.size() ) )
			{
                eDensity = mNe[ indxp ];
			}			
		}

		return eDensity;
	}

	//***************************************************************************************************************************
	//* - Bext on pNode.
	//***************************************************************************************************************************
	Vector<double> ColdPlasma::Get_Bext_OnNode( Node::Pointer pNode )
	{
		Vector<double> Bext( 3, 0.0 );

	    if( mGeometryExtType < 2 )
	    {
		    double      incRb = Node_XYZ_to_incR ( pNode, mPoloidalCurv_Be );
		    Vector<int> indxb = Get_Indexes_In_Grid( mLo_B, incRb );
		    
		    Bext = Get_Bext_1D( incRb, indxb, pNode );
		}
		else
		{
			int indxp = pNode->Id() - 1;

			if( ( indxp >= 0 ) && ( indxp < (int)mBx.size() ) )
			{
                Bext[ 0 ] = mBx[ indxp ];
			    Bext[ 1 ] = mBy[ indxp ];
			    Bext[ 2 ] = mBz[ indxp ];
			}			
		}

		return Bext;
	}

    //***************************************************************************************************************************
    // - eDensity and Bext on pNode.
    //***************************************************************************************************************************
    void ColdPlasma::Get_eDensity_Bext_OnNode( double& eDensity, Vector<double>& Bext, Node::Pointer pNode )
    {
	    eDensity = 0.0;
			
		Bext.resize( 3, 0.0 );	    
		
		if( mGeometryExtType < 2 )
	    {
		    double      incRe = Node_XYZ_to_incR ( pNode, mPoloidalCurv_ne );
		    Vector<int> indxe = Get_Indexes_In_Lo( incRe );

            eDensity = Get_eDensity_1D( incRe, indxe );
		    
		    double      incRb = Node_XYZ_to_incR ( pNode, mPoloidalCurv_Be );
		    Vector<int> indxb = Get_Indexes_In_Grid( mLo_B, incRb );

		    Bext = Get_Bext_1D( incRb, indxb, pNode );
		}
		else
		{
			int indxp = pNode->Id() - 1;

			if( ( indxp >= 0 ) && ( indxp < (int)mNe.size() ) )
			{
                eDensity = mNe[ indxp ];
			}

			if( ( indxp >= 0 ) && ( indxp < (int)mBx.size() ) )
			{
                Bext[ 0 ] = mBx[ indxp ];
			    Bext[ 1 ] = mBy[ indxp ];
			    Bext[ 2 ] = mBz[ indxp ];
			}			
		}
    }

    //***************************************************************************************************************************
	// - eDensity, Bext and [ S, D, P, R, L ] components of the permittivity tensor on pNode.
    //***************************************************************************************************************************
	void ColdPlasma::Get_AllPlasmaParameters_OnNode( double& eDensity, Vector<double>& Bext, ComplexVector& SDPRL, Node::Pointer pNode, double wFreq )
	{
		Get_eDensity_Bext_OnNode( eDensity, Bext, pNode ); 

        SDPRL = Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
	}

    //***************************************************************************************************************************
    // - Bext and [ S, D, P, R, L ] components of the permittivity tensor on pNode.
    //***************************************************************************************************************************
    void ColdPlasma::Get_Bext_SDPRL_OnNode( Vector<double>& Bext, ComplexVector& SDPRL, Node::Pointer pNode, double wFreq )
    {
        double eDensity; 

		Get_eDensity_Bext_OnNode( eDensity, Bext, pNode ); 

        SDPRL = Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
    }

	//***************************************************************************************************************************
	//* - [ S, D, P, R, L ] components of the permittivity tensor on pNode.
	//***************************************************************************************************************************
    ColdPlasma::ComplexVector ColdPlasma::Get_SDPRL_OnNode( Node::Pointer pNode, double wFreq )
	{
		double eDensity;

		Vector<double> Bext;

		Get_eDensity_Bext_OnNode( eDensity, Bext, pNode ); 

        return Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
	}

	//***************************************************************************************************************************
	//* - [ S, D, P, R, L ] components of the permittivity tensor.
	//***************************************************************************************************************************
    ColdPlasma::ComplexVector ColdPlasma::Calculate_PlasmaPermittivity_SDPRL( double eDensity, Vector<double>& Bext, double wFreq )
	{
		// Stix cold plasma dielectric tensor ( time convention exp(-i*w*t) ):
		//
		//   R = 1 - sum_s wps^2 / ( w ( w + i*nu_s + Wcs ) )
		//   L = 1 - sum_s wps^2 / ( w ( w + i*nu_s - Wcs ) )
		//   P = 1 - sum_s wps^2 / ( w ( w + i*nu_s ) )
		//
		// with signed cyclotron frequencies Wcs = q_s |B| / m_s ( negative for electrons ).
		// Collisions enter as w -> w + i*nu ( damping modes 3 and 4 ); modes 1 and 2 add i*delta to R, L and P.
		// The species table ( charge, mass, density ratio ) is built once in Set_Species_Table().

		const double modBext = std::sqrt( Bext[ 0 ]*Bext[ 0 ] + Bext[ 1 ]*Bext[ 1 ] + Bext[ 2 ]*Bext[ 2 ] );
		const double Wf      = wFreq;
		const double eo      = PhysConst::Eps0;

		const std::complex<double> cUnit( 0.0, 1.0 );

		std::complex<double> R( 1.0, 0.0 );
		std::complex<double> L( 1.0, 0.0 );
		std::complex<double> P( 1.0, 0.0 );

		for( int s=0; s<mNumSpecies; s++ )
		{
			const double Qs = mSp_Q[ s ];
			const double Ms = mSp_M[ s ];
			const double Ns = eDensity * mSp_Ratio[ s ];

			if( Ns == 0.0 ) continue;

			// Plasma and ( signed ) cyclotron frequencies
			const double W2_ps   = Qs * Qs * Ns / ( eo * Ms );
			const double Omega_s = Qs * modBext / Ms;

			// Complex frequency w + i*nu ( collisional damping )
			std::complex<double> W_col( Wf, 0.0 );

			if( ( mDamping_mode == 3 ) || ( ( mDamping_mode == 4 ) && ( s == 0 ) ) )
			{
				W_col = std::complex<double>( Wf, mDamping_value );
			}

			const std::complex<double> den_R = Wf * ( W_col + Omega_s );
			const std::complex<double> den_L = Wf * ( W_col - Omega_s );
			const std::complex<double> den_P = Wf * ( W_col );

			if( ( den_R == 0.0 ) || ( den_L == 0.0 ) || ( den_P == 0.0 ) )
			{
				Send_Error_Msg( "Exact cyclotron resonance or zero frequency in cold plasma tensor ( R, L or P infinite )", "Please, check input data or add some damping." );
			}

			R -= W2_ps / den_R;
			L -= W2_ps / den_L;
			P -= W2_ps / den_P;
		}

		// Adding damping in 'Background_const' damping mode
		if( mDamping_mode == 1 )
		{
		    R += cUnit * mDamping_value;
		    L += cUnit * mDamping_value;
		    P += cUnit * mDamping_value;
		}

		// Adding damping in 'Background_denst' damping mode
		if( ( mDamping_mode == 2 ) && ( eDensity > 0.0 ) )
		{
	        R += cUnit * mDamping_value;
	        L += cUnit * mDamping_value;
	        P += cUnit * mDamping_value;
		}

		// Components S, D and P of the plasma permittivity tensor (Stix notation)
		// Components R, L and P from the diagonal plasma permittivity tensor
        ComplexVector SDPRL( 5 );

		SDPRL[ 0 ] = 0.5 * ( R + L );
	    SDPRL[ 1 ] = 0.5 * ( R - L );
	    SDPRL[ 2 ] = P;
		SDPRL[ 3 ] = R;
		SDPRL[ 4 ] = L;

        return SDPRL;
	}

	//***************************************************************************************************************************
	//* - Species table: charge, mass and density ratio of every species with non-zero density.
	//*   Index 0 is always the electron ( used by damping mode 4 ). Ion masses are neutral atomic masses minus the mass of
	//*   the removed electrons ( previous versions used neutral masses, ~0.05% too heavy for H+ ).
	//***************************************************************************************************************************
	void ColdPlasma::Set_Species_Table()
	{
		using namespace PhysConst;

		struct SpeciesData { double ChargeState; double AtomicMass_u; double Ratio; };

		const SpeciesData Table[ mMaxSpecies ] =
		{
			{  0.0, 0.0     , 1.0                }, // electron ( mass/charge set below )
			{  1.0, M_H1    , mH1_div_ne         }, // H1   +1
			{  1.0, M_H2    , mH2_div_ne         }, // H2   +1
			{  1.0, M_H3    , mH3_div_ne         }, // H3   +1
			{  1.0, M_He3   , mHe3_q1_div_ne     }, // He3  +1
			{  2.0, M_He3   , mHe3_q2_div_ne     }, // He3  +2
			{  1.0, M_He4   , mHe4_q1_div_ne     }, // He4  +1
			{  2.0, M_He4   , mHe4_q2_div_ne     }, // He4  +2
			{  1.0, M_Be9   , mBe9_q1_div_ne     }, // Be9  +1
			{  2.0, M_Be9   , mBe9_q2_div_ne     }, // Be9  +2
			{  3.0, M_Be9   , mBe9_q3_div_ne     }, // Be9  +3
			{  4.0, M_Be9   , mBe9_q4_div_ne     }, // Be9  +4
			{ 10.0, M_Ne20  , mNe20_q10_div_ne   }, // Ne20 +10
			{ 20.0, M_Ni60  , mNi60_q20_div_ne   }, // Ni60 +20
			{ 25.0, M_W184  , mW184_q25_div_ne   }  // W184 +25
		};

		// Electron is always the first species
		mNumSpecies = 0;
		mSp_Q    [ mNumSpecies ] = -Qe;
		mSp_M    [ mNumSpecies ] =  Me;
		mSp_Ratio[ mNumSpecies ] =  1.0;
		mNumSpecies++;

		double ChargeNeutrality = 0.0;

		for( int s=1; s<mMaxSpecies; s++ )
		{
			if( Table[ s ].Ratio < 0.0 )
			{
				Send_Error_Msg( "Negative ion density ratio in cold plasma data", "Please, check plasma species densities." );
			}

			ChargeNeutrality += Table[ s ].ChargeState * Table[ s ].Ratio;

			if( Table[ s ].Ratio == 0.0 ) continue;

			mSp_Q    [ mNumSpecies ] = Table[ s ].ChargeState * Qe;
			mSp_M    [ mNumSpecies ] = IonMass( Table[ s ].AtomicMass_u, Table[ s ].ChargeState );
			mSp_Ratio[ mNumSpecies ] = Table[ s ].Ratio;
			mNumSpecies++;
		}

		// Quasi-neutrality check ( sum Z_i n_i / n_e should be 1 ). Only a warning: some studies use non-neutral mixtures.
		if( std::abs( ChargeNeutrality - 1.0 ) > 1.0e-3 )
		{
			std::cout << "Warning: cold plasma species are not quasi-neutral, sum( Z_i n_i ) / n_e = " << ChargeNeutrality << std::endl;
		}
	}

	//***************************************************************************************************************************
	//* - Hermitian storage ( Herm-Full / Herm-Symm ) keeps only the upper triangle of the plasma block and rebuilds the lower
	//*   triangle as its complex conjugate. That is exact only for a lossless ( Hermitian ) tensor. With damping the result
	//*   would be wrong, so full storage is selected instead.
	//***************************************************************************************************************************
	void ColdPlasma::Check_Matrix_Storage_vs_Damping()
	{
		const bool HermitianStorage = ( mMatrixStorageFormat == 1 ) || ( mMatrixStorageFormat == 2 );
		const bool DampingOn        = ( mDamping_mode != 0 ) && ( mDamping_value != 0.0 );

		if( HermitianStorage && DampingOn )
		{
			std::cout << std::endl;
			std::cout << "Warning: Hermitian matrix storage requires a lossless plasma tensor, but damping is active." << std::endl;
			std::cout << "         Switching to Full_matrix storage." << std::endl;
			std::cout << std::endl;

			mMatrixStorageFormat = 0;
		}
	}

	//***************************************************************************************************************************
	//* - Checks that a 1D profile grid is not empty and strictly ascending ( required by Get_Indexes_In_Grid ).
	//***************************************************************************************************************************
	void ColdPlasma::Check_1D_Grid( Vector<double>& Grid, String FileName )
	{
		if( Grid.size() == 0 )
		{
			Send_Error_Msg_Str( "No data read from 1D profile file " + FileName, "Please, check file format." );
		}

		for( unsigned int i=1; i<Grid.size(); i++ )
		{
			if( !( Grid[ i ] > Grid[ i - 1 ] ) )
			{
				Send_Error_Msg_Str( "Positions not strictly ascending in 1D profile file " + FileName, "Please, sort the file by its first column." );
			}
		}
	}

	//***************************************************************************************************************************
	//* - True if a profile file line holds numerical data ( first non-blank character starts a number ).
	//*   Blank lines and text/comment lines are skipped. Previous versions only accepted lines starting with a digit, which
	//*   silently dropped lines with leading blanks, signs ( "-0.01" ) or a leading decimal point ( ".5" ).
	//***************************************************************************************************************************
	bool ColdPlasma::Is_Data_Line( const std::string& FileLine )
	{
		std::size_t pos = FileLine.find_first_not_of( " \t\r" );

		if( pos == std::string::npos ) return false;

		char c = FileLine[ pos ];

		if( std::isdigit( (unsigned char)c ) ) return true;

		if( ( c == '-' ) || ( c == '+' ) || ( c == '.' ) )
		{
			return ( pos + 1 < FileLine.size() ) && ( std::isdigit( (unsigned char)FileLine[ pos + 1 ] ) || ( FileLine[ pos + 1 ] == '.' ) );
		}

		return false;
	}

    //***************************************************************************************************************************
    //* - Rotation matrix.
    //***************************************************************************************************************************
    void ColdPlasma::Rotation_Matrix( Matrix<std::complex<double> >& T, double Bx, double By, double Bz )
    {
         T.Resize( 3, 3 );

         double Bnorm = std::sqrt( Bx*Bx + By*By + Bz*Bz );

         // With B = 0 the cold plasma tensor is isotropic ( D = 0, S = P ), so any frame is valid.
         // Return the identity instead of dividing by zero ( previous versions produced NaN ).
         if( Bnorm == 0.0 )
         {
             for( int row=0; row<3; row++ )
             {
                 for( int col=0; col<3; col++ )
                 {
                     T[ row ][ col ] = std::complex<double>( ( row == col ) ? 1.0 : 0.0, 0.0 );
                 }
             }

             return;
         }

         // Unit vectors
         Vector<double> n( 3 ), t( 3 ), b( 3 );

         // b = B / |B|
         b[ 0 ] = Bx / Bnorm;
         b[ 1 ] = By / Bnorm;
         b[ 2 ] = Bz / Bnorm;

         // n _|_ b
         double nnorm = std::sqrt( b[ 0 ]*b[ 0 ] + b[ 1 ]*b[ 1 ] );

         if( nnorm > 0.0 )
         {
             n[ 0 ] = b[ 1 ] / nnorm;
             n[ 1 ] =-b[ 0 ] / nnorm;
             n[ 2 ] = 0.0;
         }
         else
         {
             n[ 0 ] = b[ 2 ];
             n[ 1 ] = 0.0;
             n[ 2 ] = 0.0;
         }

         // t = b x n  ( ( n, t, b ) is right-handed, as required by the Stix tensor with K_xy = -i*D )
         t[ 0 ] = b[ 1 ]*n[ 2 ] - b[ 2 ]*n[ 1 ];
         t[ 1 ] = b[ 2 ]*n[ 0 ] - b[ 0 ]*n[ 2 ];
         t[ 2 ] = b[ 0 ]*n[ 1 ] - b[ 1 ]*n[ 0 ];

         // Rotation matrix
         // T = [ nx tx bx
         //       ny ty by
         //       nz tz bz ]
         for( int row=0; row<3; row++ )
         {
             T[ row ][ 0 ] = std::complex<double>( n[ row ] );
             T[ row ][ 1 ] = std::complex<double>( t[ row ] );
             T[ row ][ 2 ] = std::complex<double>( b[ row ] );
         }
    }

	//***************************************************************************************************************************
	// - Permittivity tensor on pNode.
	//***************************************************************************************************************************
	void ColdPlasma::Get_PermittivityTensor_OnNode( Matrix<std::complex<double> >& TEn, Node::Pointer pNode, double wFreq )
	{
		// Permittivity tensor in a node
		TEn.Resize( 3, 3 );

		// Relative permittivity tensor in b = B/|B| coordinates
		Matrix<std::complex<double> > K( 3, 3 );

		// Bext on node
		Vector<double> Bext = Get_Bext_OnNode( pNode );

		// Constants
        double eo = PhysConst::Eps0;

		// Get K relative in b = B/|B| coordinates
		if( mGeometryExtType == 3 )
		{
		    int indxp = pNode->Id() - 1;

		    if( ( indxp < 0 ) || ( indxp >= (int)mImp_K_Tensor.size() ) )
		    {
		        Send_Error_Msg( "Node without imported K tensor", "Please, check Tensor_K_CPGE.bin and the 3D density file." );
		    }

		    K = mImp_K_Tensor[ indxp ];
		}
		else
		{
		    double eDensity = Get_eDensity_OnNode( pNode );

		    ColdPlasma::ComplexVector SDPRL;

            SDPRL = Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
			
			std::complex<double> S = SDPRL[ 0 ];
			std::complex<double> D = SDPRL[ 1 ];
			std::complex<double> P = SDPRL[ 2 ];

			std::complex<double> cZero( 0.0, 0.0 );
			std::complex<double> cUnit( 0.0, 1.0 );

		    K[ 0 ][ 0 ] = S        ; K[ 0 ][ 1 ] =-cUnit * D; K[ 0 ][ 2 ] = cZero;
			K[ 1 ][ 0 ] = cUnit * D; K[ 1 ][ 1 ] = S        ; K[ 1 ][ 2 ] = cZero;
			K[ 2 ][ 0 ] = cZero    ; K[ 2 ][ 1 ] = cZero    ; K[ 2 ][ 2 ] = P    ;
	    }

		// K absolute
		K *= eo;

		// Rotation matrix
		Matrix< std::complex<double> > T;

		Rotation_Matrix( T, Bext[0], Bext[1], Bext[2] );

		// TEn = T * K * Tt
		// Be carefull!!, T.Transpose() changes T permanently
		TEn = T   * K;
		TEn = TEn * T.Transpose();
	}

    //***************************************************************************************************************************
	// - Conductivity tensor on pNode.
	//***************************************************************************************************************************
	void ColdPlasma::Get_ConductivityTensor_OnNode( Matrix<std::complex<double> >& TCn, Node::Pointer pNode, double wFreq )
	{
		// Coductivity tensor in a node
		TCn.Resize( 3, 3 );

        // Permittivity tensor in a node
        Matrix< std::complex<double> > TEn;

        Get_PermittivityTensor_OnNode( TEn, pNode, wFreq );

        // Constants
        double eo = PhysConst::Eps0;

        std::complex<double> iw( 0.0, wFreq );

        TCn[ 0 ][ 0 ] = -iw * ( TEn[ 0 ][ 0 ] - eo );
		TCn[ 1 ][ 0 ] = -iw * ( TEn[ 1 ][ 0 ]      );
		TCn[ 2 ][ 0 ] = -iw * ( TEn[ 2 ][ 0 ]      );

        TCn[ 0 ][ 1 ] = -iw * ( TEn[ 0 ][ 1 ]      );
		TCn[ 1 ][ 1 ] = -iw * ( TEn[ 1 ][ 1 ] - eo );
		TCn[ 2 ][ 1 ] = -iw * ( TEn[ 2 ][ 1 ]      );

        TCn[ 0 ][ 2 ] = -iw * ( TEn[ 0 ][ 2 ]      );
		TCn[ 1 ][ 2 ] = -iw * ( TEn[ 1 ][ 2 ]      );
		TCn[ 2 ][ 2 ] = -iw * ( TEn[ 2 ][ 2 ] - eo );
	}

	//***************************************************************************************************************************
	//* - Read electron density file in 1D format.
	//***************************************************************************************************************************
	void ColdPlasma::Read_eDensity_File_1D( String FileName )
	{
		// Format: one "position density" pair per line. Text / blank lines are ignored.
		std::fstream eDensityFile( FileName, std::fstream::in );

		if( !eDensityFile.is_open() )
		{
			Send_Error_Msg_Str( "Error opening 1D electron density file " + FileName, "Please, check file path." );
		}

        String FileLine;

        while( std::getline( eDensityFile, FileLine ) ) 
		{ 
		    if( !Is_Data_Line( FileLine ) ) continue;

		    std::stringstream LineStrm( FileLine );

			double Lo, Ne;

		    if( !( LineStrm >> Lo >> Ne ) )
			{
				Send_Error_Msg_Str( "Wrong line in 1D electron density file " + FileName + ": " + FileLine, "Expected: position density" );
			}

            mLo.push_back( Lo );
            mNe.push_back( Ne );	
        } 

		eDensityFile.close();
	}

	//***************************************************************************************************************************
	//* - Set 3D Ne and Bxyz vectors size.
	//***************************************************************************************************************************
	void ColdPlasma::Set_3D_Vector_Ne_Bxyz()
	{
        int MaxNodeId = 0;

        std::fstream Nodes_File( mBase_FileName + "-1.dat", std::fstream::in );  

		if( !Nodes_File.is_open() )
		{
			Send_Error_Msg_Str( "Error opening nodes file " + mBase_FileName + "-1.dat", "Please, check the problem files." );
		}

        // Line from nodes file
		std::string line;

        while( std::getline( Nodes_File, line ) )
        {
            // Ignore comments
            if( line.substr( 0, 2 ) == "//" ) continue;

            // Read formatted line ( lines that do not match are ignored )
            int NodeId = 0; double X, Y, Z;

            if( std::sscanf( line.c_str(), "No[%d] = p(%lf,%lf,%lf);", &NodeId, &X, &Y, &Z ) != 4 ) continue;

			// Update MaxNodeId
            if( NodeId > MaxNodeId ) MaxNodeId = NodeId;
        }

        Nodes_File.close();

		Resize_Profile_3D_Vectors( MaxNodeId );
	}

	//***************************************************************************************************************************
	//* - Read electron density file in 3D format.
	//***************************************************************************************************************************
	void ColdPlasma::Read_eDensity_File_3D( String FileName )
	{
		// Format: one "nodeId density" pair per line ( nodeId starts at 1 ). Text / blank lines are ignored.
		std::fstream eDensityFile( FileName, std::fstream::in );

		if( !eDensityFile.is_open() )
		{
			Send_Error_Msg_Str( "Error opening 3D electron density file " + FileName, "Please, check file path." );
		}

		String FileLine;

		while( std::getline( eDensityFile, FileLine ) ) 
		{ 
		    if( !Is_Data_Line( FileLine ) ) continue;

		    std::stringstream LineStrm( FileLine );

			long nId; double nNe;

		    if( !( LineStrm >> nId >> nNe ) )
			{
				Send_Error_Msg_Str( "Wrong line in 3D electron density file " + FileName + ": " + FileLine, "Expected: nodeId density" );
			}

			if( ( nId < 1 ) || ( nId > (long)mNe.size() ) )
			{
				Send_Error_Msg_Str( "Node index out of range in density file " + FileName, "Please, check density file." );
			}

            mNe[ nId - 1 ] = nNe;	
        } 

		eDensityFile.close();
	}

	//***************************************************************************************************************************
	//* - Read external static B fields file in 1D format.
	//***************************************************************************************************************************
	void ColdPlasma::Read_sBfields_File_1D( String FileName )
	{
		// Format: one "position |B| angle(deg)" triplet per line. Text / blank lines are ignored.
		// The positions are stored in mLo_B: the B grid does not need to match the density grid
		// ( previous versions discarded them and indexed B with the density grid ).
		std::fstream sBfieldsFile( FileName, std::fstream::in );

		if( !sBfieldsFile.is_open() )
		{
			Send_Error_Msg_Str( "Error opening 1D external static B fields file " + FileName, "Please, check file path." );
		}

        String FileLine;

        while( std::getline( sBfieldsFile, FileLine ) )  
		{ 
		    if( !Is_Data_Line( FileLine ) ) continue;

		    std::stringstream LineStrm( FileLine );

			double Lo, modB, angB;

		    if( !( LineStrm >> Lo >> modB >> angB ) )
			{
				Send_Error_Msg_Str( "Wrong line in 1D B field file " + FileName + ": " + FileLine, "Expected: position modB angB" );
			}

            mLo_B .push_back( Lo   );
            mB_mod.push_back( modB );
            mB_ang.push_back( angB );	
        } 

		sBfieldsFile.close();
	}

	//***************************************************************************************************************************
	//* - Read external static B fields file in 3D format.
	//***************************************************************************************************************************
	void ColdPlasma::Read_sBfields_File_3D( String FileName )
	{
		// Format: one "nodeId Bx By Bz" line per node ( nodeId starts at 1 ). Text / blank lines are ignored.
		std::fstream sBfieldsFile( FileName, std::fstream::in );

		if( !sBfieldsFile.is_open() )
		{
			Send_Error_Msg_Str( "Error opening 3D external static B fields file " + FileName, "Please, check file path." );
		}

		String FileLine;

        while( std::getline( sBfieldsFile, FileLine ) )  
		{ 
		    if( !Is_Data_Line( FileLine ) ) continue;

		    std::stringstream LineStrm( FileLine );

			long nId; double nBx, nBy, nBz;

		    if( !( LineStrm >> nId >> nBx >> nBy >> nBz ) )
			{
				Send_Error_Msg_Str( "Wrong line in 3D B field file " + FileName + ": " + FileLine, "Expected: nodeId Bx By Bz" );
			}

			if( ( nId < 1 ) || ( nId > (long)mBx.size() ) )
			{
				Send_Error_Msg_Str( "Node index out of range in Be field file " + FileName, "Please, check Be field file." );
			}

			mBx[ nId - 1 ] = nBx;
            mBy[ nId - 1 ] = nBy;
			mBz[ nId - 1 ] = nBz;	
        } 

		sBfieldsFile.close();
	}

	//***************************************************************************************************************************
	//* - Read imported K tensor from file.
	//***************************************************************************************************************************
	void ColdPlasma::Read_K_Tensor_File_3D( String FileName )
	{
		// Binary file with one 3x3 complex<double> tensor ( row major ) per data line of the 3D density file, in the
		// same order. The file name is fixed ( "Tensor_K_CPGE.bin", in the working directory ).
		std::ifstream K_Tensor_file( "Tensor_K_CPGE.bin", std::ios::binary );

		if( !K_Tensor_file.is_open() )
		{
			Send_Error_Msg( "Error opening imported K tensor file Tensor_K_CPGE.bin", "Please, check file location." );
		}

		// 3D density file gives the node index of each tensor 
		std::fstream Ne_file( FileName, std::fstream::in );

		if( !Ne_file.is_open() )
		{
			Send_Error_Msg_Str( "Error opening 3D electron density file " + FileName, "Please, check file path." );
		}

		// Binary matrix size
        const int VectorSize         = 3; 
        const int MatrixSize         = VectorSize * VectorSize; 
        const int Chunck_Size_Matrix = MatrixSize * sizeof( std::complex<double> );

		// Resize K tensors vector
		mImp_K_Tensor.resize( mNe.size() );

		for( unsigned int ki=0; ki<mNe.size(); ki++ ) 
		{ 
			mImp_K_Tensor[ ki ].Resize( VectorSize, VectorSize, std::complex<double>( 0.0, 0.0 ) );
		}

		// Read imported K tensors
		std::string line;	

		Vector< std::complex<double> > Kb( MatrixSize, std::complex<double>( 0.0, 0.0 ) );

        while( std::getline( Ne_file, line ) ) 
        {
			// Only data lines have a tensor in the binary file ( previous versions also consumed a tensor for
			// header lines and used an undefined node index for them )
			if( !Is_Data_Line( line ) ) continue;

			std::stringstream LineStrm( line );

			long nId; double nNe;

			if( !( LineStrm >> nId >> nNe ) )
			{
				Send_Error_Msg_Str( "Wrong line in 3D electron density file " + FileName + ": " + line, "Expected: nodeId density" );
			}

			if( ( nId < 1 ) || ( nId > (long)mImp_K_Tensor.size() ) )
			{
				Send_Error_Msg_Str( "Node index out of range in density file " + FileName, "Please, check density file." );
			}

            K_Tensor_file.read( reinterpret_cast<char*>( Kb.data() ), Chunck_Size_Matrix );

            if( K_Tensor_file.gcount() != Chunck_Size_Matrix )
            {
				Send_Error_Msg( "Imported K tensor file Tensor_K_CPGE.bin is shorter than the density file", "Please, check that both files correspond to the same mesh." );
            }

			Matrix< std::complex<double> >& K_Tensor = mImp_K_Tensor[ nId - 1 ];

            for( int i=0; i<VectorSize; i++ )
            {
                for( int j=0; j<VectorSize; j++ )
                {
                    K_Tensor[ i ][ j ] = Kb[ i * VectorSize + j ];
                }
            }
        }

        Ne_file      .close();
        K_Tensor_file.close();
	}

	//***************************************************************************************************************************
	//* - Linear interpolation of plasma profiles to high-order nodes.
	//***************************************************************************************************************************
	void ColdPlasma::Linear_Profile_Interpolate( Vector<int> HONodesId )
	{
		// Interpolate plasma profile to 2nd order middle nodes
		if( HONodesId.size() == 10 )
		{
		    // 2nd order edge node Ids
		    int iniNId, midNId, endNId; 
		    
		    // Edge[ 0 , 1 ]
		    iniNId = HONodesId[ 0 ] - 1; 
		    midNId = HONodesId[ 4 ] - 1; 
		    endNId = HONodesId[ 1 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}

		    // Edge[ 1 , 2 ]
		    iniNId = HONodesId[ 1 ] - 1; 
		    midNId = HONodesId[ 5 ] - 1; 
		    endNId = HONodesId[ 2 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}

		    // Edge[ 2 , 0 ]
		    iniNId = HONodesId[ 2 ] - 1; 
		    midNId = HONodesId[ 6 ] - 1; 
		    endNId = HONodesId[ 0 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}
		    
		    // Edge[ 0 , 3 ]
		    iniNId = HONodesId[ 0 ] - 1; 
		    midNId = HONodesId[ 7 ] - 1; 
		    endNId = HONodesId[ 3 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}
		    
		    // Edge[ 1 , 3 ]
		    iniNId = HONodesId[ 1 ] - 1; 
		    midNId = HONodesId[ 8 ] - 1; 
		    endNId = HONodesId[ 3 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}
		    
		    // Edge[ 2 , 3 ]
		    iniNId = HONodesId[ 2 ] - 1; 
		    midNId = HONodesId[ 9 ] - 1; 
		    endNId = HONodesId[ 3 ] - 1;
		    
		    mNe[ midNId ] = 0.5 * ( mNe[ iniNId ] + mNe[ endNId ] );
		    mBx[ midNId ] = 0.5 * ( mBx[ iniNId ] + mBx[ endNId ] );
		    mBy[ midNId ] = 0.5 * ( mBy[ iniNId ] + mBy[ endNId ] );
		    mBz[ midNId ] = 0.5 * ( mBz[ iniNId ] + mBz[ endNId ] );

			if( mGeometryExtType == 3 ) 
			{
				mImp_K_Tensor[ midNId ]  = mImp_K_Tensor[ iniNId ];
				mImp_K_Tensor[ midNId ] += mImp_K_Tensor[ endNId ];
				mImp_K_Tensor[ midNId ] *= 0.5;
			}
		}
		// Interpolate plasma profile to 3sb surface middle nodes
		else if ( HONodesId.size() == 17 )
		{
		    // 3sb surface vertices nodes Ids
		    int vNId0, vNId1, vNId2; 

			// 3sb surface inner nodes Ids
			int sNId0, sNId1, sNId2;  
		    
		    // Face[ 0 , 1 , 2 ]
		    vNId0 = HONodesId[ 0 ] - 1; 
		    vNId1 = HONodesId[ 1 ] - 1; 
		    vNId2 = HONodesId[ 2 ] - 1;

			sNId0 = HONodesId[ 4 ] - 1; 
		    sNId1 = HONodesId[ 5 ] - 1; 
		    sNId2 = HONodesId[ 6 ] - 1;

		    mNe[ sNId0 ] = 0.50 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId0 ] = 0.50 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId0 ] = 0.50 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId0 ] = 0.50 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId1 ] = 0.25 * mNe[ vNId0 ] + 0.50 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId1 ] = 0.25 * mBx[ vNId0 ] + 0.50 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId1 ] = 0.25 * mBy[ vNId0 ] + 0.50 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId1 ] = 0.25 * mBz[ vNId0 ] + 0.50 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId2 ] = 0.25 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.50 * mNe[ vNId2 ];
		    mBx[ sNId2 ] = 0.25 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.50 * mBx[ vNId2 ];
		    mBy[ sNId2 ] = 0.25 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.50 * mBy[ vNId2 ];
		    mBz[ sNId2 ] = 0.25 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.50 * mBz[ vNId2 ];

			if( mGeometryExtType == 3 ) 
			{
				Matrix< std::complex<double> > K0, K1, K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.50; mImp_K_Tensor[ sNId0 ]  = K0; 
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId0 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId0 ] += K2;
				
				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId1 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.50; mImp_K_Tensor[ sNId1 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId1 ] += K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId2 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId2 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.50; mImp_K_Tensor[ sNId2 ] += K2;
			}

			// Face[ 0 , 1 , 3 ]
		    vNId0 = HONodesId[ 0 ] - 1; 
		    vNId1 = HONodesId[ 1 ] - 1; 
		    vNId2 = HONodesId[ 3 ] - 1;

			sNId0 = HONodesId[ 7 ] - 1; 
		    sNId1 = HONodesId[ 8 ] - 1; 
		    sNId2 = HONodesId[ 9 ] - 1;

		    mNe[ sNId0 ] = 0.50 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId0 ] = 0.50 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId0 ] = 0.50 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId0 ] = 0.50 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId1 ] = 0.25 * mNe[ vNId0 ] + 0.50 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId1 ] = 0.25 * mBx[ vNId0 ] + 0.50 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId1 ] = 0.25 * mBy[ vNId0 ] + 0.50 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId1 ] = 0.25 * mBz[ vNId0 ] + 0.50 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId2 ] = 0.25 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.50 * mNe[ vNId2 ];
		    mBx[ sNId2 ] = 0.25 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.50 * mBx[ vNId2 ];
		    mBy[ sNId2 ] = 0.25 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.50 * mBy[ vNId2 ];
		    mBz[ sNId2 ] = 0.25 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.50 * mBz[ vNId2 ];

			if( mGeometryExtType == 3 ) 
			{
				Matrix< std::complex<double> > K0, K1, K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.50; mImp_K_Tensor[ sNId0 ]  = K0; 
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId0 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId0 ] += K2;
				
				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId1 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.50; mImp_K_Tensor[ sNId1 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId1 ] += K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId2 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId2 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.50; mImp_K_Tensor[ sNId2 ] += K2;
			}

			// Face[ 1 , 2 , 3 ]
		    vNId0 = HONodesId[ 1 ] - 1; 
		    vNId1 = HONodesId[ 2 ] - 1; 
		    vNId2 = HONodesId[ 3 ] - 1;

			sNId0 = HONodesId[ 10 ] - 1; 
		    sNId1 = HONodesId[ 11 ] - 1; 
		    sNId2 = HONodesId[ 12 ] - 1;

		    mNe[ sNId0 ] = 0.50 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId0 ] = 0.50 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId0 ] = 0.50 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId0 ] = 0.50 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId1 ] = 0.25 * mNe[ vNId0 ] + 0.50 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId1 ] = 0.25 * mBx[ vNId0 ] + 0.50 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId1 ] = 0.25 * mBy[ vNId0 ] + 0.50 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId1 ] = 0.25 * mBz[ vNId0 ] + 0.50 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId2 ] = 0.25 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.50 * mNe[ vNId2 ];
		    mBx[ sNId2 ] = 0.25 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.50 * mBx[ vNId2 ];
		    mBy[ sNId2 ] = 0.25 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.50 * mBy[ vNId2 ];
		    mBz[ sNId2 ] = 0.25 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.50 * mBz[ vNId2 ];

			if( mGeometryExtType == 3 ) 
			{
				Matrix< std::complex<double> > K0, K1, K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.50; mImp_K_Tensor[ sNId0 ]  = K0; 
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId0 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId0 ] += K2;
				
				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId1 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.50; mImp_K_Tensor[ sNId1 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId1 ] += K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId2 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId2 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.50; mImp_K_Tensor[ sNId2 ] += K2;
			}

			// Face[ 0 , 2 , 3 ]
		    vNId0 = HONodesId[ 0 ] - 1; 
		    vNId1 = HONodesId[ 2 ] - 1; 
		    vNId2 = HONodesId[ 3 ] - 1;

			sNId0 = HONodesId[ 13 ] - 1; 
		    sNId1 = HONodesId[ 14 ] - 1; 
		    sNId2 = HONodesId[ 15 ] - 1;

		    mNe[ sNId0 ] = 0.50 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId0 ] = 0.50 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId0 ] = 0.50 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId0 ] = 0.50 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId1 ] = 0.25 * mNe[ vNId0 ] + 0.50 * mNe[ vNId1 ] + 0.25 * mNe[ vNId2 ];
		    mBx[ sNId1 ] = 0.25 * mBx[ vNId0 ] + 0.50 * mBx[ vNId1 ] + 0.25 * mBx[ vNId2 ];
		    mBy[ sNId1 ] = 0.25 * mBy[ vNId0 ] + 0.50 * mBy[ vNId1 ] + 0.25 * mBy[ vNId2 ];
		    mBz[ sNId1 ] = 0.25 * mBz[ vNId0 ] + 0.50 * mBz[ vNId1 ] + 0.25 * mBz[ vNId2 ];

		    mNe[ sNId2 ] = 0.25 * mNe[ vNId0 ] + 0.25 * mNe[ vNId1 ] + 0.50 * mNe[ vNId2 ];
		    mBx[ sNId2 ] = 0.25 * mBx[ vNId0 ] + 0.25 * mBx[ vNId1 ] + 0.50 * mBx[ vNId2 ];
		    mBy[ sNId2 ] = 0.25 * mBy[ vNId0 ] + 0.25 * mBy[ vNId1 ] + 0.50 * mBy[ vNId2 ];
		    mBz[ sNId2 ] = 0.25 * mBz[ vNId0 ] + 0.25 * mBz[ vNId1 ] + 0.50 * mBz[ vNId2 ];

			if( mGeometryExtType == 3 ) 
			{
				Matrix< std::complex<double> > K0, K1, K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.50; mImp_K_Tensor[ sNId0 ]  = K0; 
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId0 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId0 ] += K2;
				
				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId1 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.50; mImp_K_Tensor[ sNId1 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.25; mImp_K_Tensor[ sNId1 ] += K2;

				K0 = mImp_K_Tensor[ vNId0 ]; K0 *= 0.25; mImp_K_Tensor[ sNId2 ]  = K0;
				K1 = mImp_K_Tensor[ vNId1 ]; K1 *= 0.25; mImp_K_Tensor[ sNId2 ] += K1;
				K2 = mImp_K_Tensor[ vNId2 ]; K2 *= 0.50; mImp_K_Tensor[ sNId2 ] += K2;
			}
		}
	}

	//***************************************************************************************************************************
	//* - Print error messages and exit ERMES.
	//***************************************************************************************************************************
	void ColdPlasma::Send_Error_Msg( String ErrorMsg, String CheckMsg )
	{
		std::cout << std::endl;

		std::cout << "!!!!!!!!!! " << ErrorMsg << " !!!!!!!!!!" << std::endl;

		std::cout << std::endl;

		std::cout << CheckMsg << std::endl;

		std::cout << std::endl;

		std::cout << "ERMES analysis finished." << std::endl;

		exit( EXIT_SUCCESS );
	}

	//***************************************************************************************************************************
	//* - Clear cold plasma properties.
	//***************************************************************************************************************************
	ColdPlasma::~ColdPlasma()
	{
	    mLo.clear(); Vector<double>().swap( mLo );
		mNe.clear(); Vector<double>().swap( mNe );
	    
	    mBx.clear(); Vector<double>().swap( mBx );
	    mBy.clear(); Vector<double>().swap( mBy );
	    mBz.clear(); Vector<double>().swap( mBz );

		mL1_O.clear(); Vector<double>().swap( mL1_O );
	    mTK_O.clear(); Vector<double>().swap( mTK_O );	

		mTK_n.clear(); Vector<double>().swap( mTK_n );
	    mTK_t.clear(); Vector<double>().swap( mTK_t );
	    mTK_b.clear(); Vector<double>().swap( mTK_b );

	    mLo_B .clear(); Vector<double>().swap( mLo_B  );
	    mB_mod.clear(); Vector<double>().swap( mB_mod );
	    mB_ang.clear(); Vector<double>().swap( mB_ang );

		mImp_K_Tensor.clear(); Vector< Matrix<std::complex<double> > >().swap( mImp_K_Tensor ); 
	}
} 
