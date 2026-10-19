
#include "ColdPlasma.h"
#include <iostream>
#include <fstream>

namespace Kratos
{
	//***********************************************
	//* - Initiating Cold plasma properties
	//***********************************************
	ColdPlasma::ColdPlasma()
	{
		// Resizing Tokamak coordinate system
		mTK_n.resize(3);
		mTK_t.resize(3);
		mTK_b.resize(3);

		// Resize origen of Tokamak coordinate system
		mTK_OO.resize(3);

		// Setting Flat geometry as default
		mFlatGeometry   = true ;
		mCurvedGeometry = false;
		mPoloidalCurv   = 0.0;
	}

	//***********************************************
	//* - Reads data from file
	//***********************************************
	void ColdPlasma::Load_Data( String data_file_name )
	{
		// Open to read cold plasma properties data file
		std::fstream dataFile(data_file_name, std::fstream::in);

		// Material Id
		dataFile >> mMaterialId;

		if (dataFile.eof())
		{
			Send_Error_Msg("No cold plasma applied", "Please, apply plasma material to volumes.");
		}

		// Sets the radial density profile to flat or curved 
		bool GeometryOption; 
		dataFile >> GeometryOption;

		if (GeometryOption == 0) 
		{ 
			mFlatGeometry   = true ; 
			mCurvedGeometry = false; 
		}
		else                     
		{ 
			mFlatGeometry   = false; 
			mCurvedGeometry = true ; 
		}

		// Vertical radius of curvature
		dataFile >> mPoloidalCurv;

		// Coordinates of the center of the Tokamak.
		dataFile >> mTK_OO[0]; 
		dataFile >> mTK_OO[1]; 
		dataFile >> mTK_OO[2];

		// Unitary vectors of the main central axis of the Tokamak.
		dataFile >> mTK_n[0]; 
		dataFile >> mTK_n[1]; 
		dataFile >> mTK_n[2];

		// Normal to the planes of constant electron density
		dataFile >> mTK_t[0];
		dataFile >> mTK_t[1];
		dataFile >> mTK_t[2];

		// Setting Tokamak Cartesian coordinate system
		Set_Tokamak_ntb();

		// Electron density data files 
		String eDensityFile;
		dataFile >> eDensityFile;
		Read_eDensity_File( eDensityFile );

		// Applied B field data files
		String aBfieldFile;
		dataFile >> aBfieldFile;
		Read_aBfield_File( aBfieldFile );		

		// Species density relative to the electron density 
		dataFile >> mH1_div_ne ; 
        dataFile >> mH2_div_ne ; 
        dataFile >> mH3_div_ne ;

        dataFile >> mHe4_q1_div_ne;
        dataFile >> mHe4_q2_div_ne;

        dataFile >> mBe9_q1_div_ne;
        dataFile >> mBe9_q2_div_ne;

        // Wave type applied to the far field boundary condition
        dataFile >> mFar_field_wave;

        // E parallel tolerance
        dataFile >> mEpar_tolerance;

        // Cylindrical symmetry in the plasma volume
        dataFile >> mCSym_dEdx;
        dataFile >> mCSym_dEdy;
        dataFile >> mCSym_dEdz;

        // System matrix format
        String Matrix_Format;
        
        dataFile >> Matrix_Format;

        if ( Matrix_Format == "Symmetric" ) mMatrix_Is_Symmetric = true ;
        else                                mMatrix_Is_Symmetric = false;

        // Elements matrix tolerance
        dataFile >> mKij_Tol;

		// Closing cold plasma data file
        dataFile.close();
	}

	//********************************************************************************************
	//* - Returns material Id of the palsma
	//********************************************************************************************
	int ColdPlasma::Get_Material_Id()
	{
		return mMaterialId;
	}

	//********************************************************************************************
	//* - Sets a Cartesian coordinate system with the origen in mTK_OO and the main axis mTK_n
	//********************************************************************************************
	void ColdPlasma::Set_Tokamak_ntb()
	{   
		// Setting unit vector n
		double n_norm = std::sqrt( mTK_n[0]*mTK_n[0] + mTK_n[1]*mTK_n[1] + mTK_n[2]*mTK_n[2] );

		if ( n_norm > 0.0 ) mTK_n *= (1.0 / n_norm);
		else Send_Error_Msg("Wrong Tokamak vertical axis", "Please, check plasma input data.");

		// Setting unit vector t
		double t_norm = std::sqrt( mTK_t[0]*mTK_t[0] + mTK_t[1]*mTK_t[1] + mTK_t[2]*mTK_t[2] );

		if ( t_norm > 0.0 ) mTK_t *= (1.0 / t_norm);
		else Send_Error_Msg("Wrong normal plane vector", "Please, check normal plane vector input data.");

		// Checking that t and n are ortogonal
		double dot_n_t = std::abs( mTK_n[0]*mTK_t[0] + mTK_n[1]*mTK_t[1] + mTK_n[2]*mTK_t[2] );

		if ( dot_n_t > 1e-6 ) 
			Send_Error_Msg("Wrong normal plane vector", "Normal plane vector must be perpendicular to Tokamak axis.");

		// Setting a unit vector b ortogonal to t and n (b = n x t)
		mTK_b[0] = mTK_n[1] * mTK_t[2] - mTK_n[2] * mTK_t[1];
		mTK_b[1] = mTK_n[2] * mTK_t[0] - mTK_n[0] * mTK_t[2];
		mTK_b[2] = mTK_n[0] * mTK_t[1] - mTK_n[1] * mTK_t[0];

		mTK_b *= ( 1.0 / std::sqrt( mTK_b[0]*mTK_b[0] + mTK_b[1]*mTK_b[1] + mTK_b[2]*mTK_b[2] ) );
	}

	//***************************************************************************************
	//* - Transform node coordinates to tokamak coordinate system 
	//***************************************************************************************
	Vector<double> ColdPlasma::Node_XYZ_to_TK_XYZ( Node::Pointer pNode )
	{
		// Read coordinates from the original coordinate system (OG_CS)
		Vector<double> OG_XYZ(3);
		OG_XYZ[0] = pNode->X();
		OG_XYZ[1] = pNode->Y();
		OG_XYZ[2] = pNode->Z();

		// Vector from mTK_OO to the point OG_XYZ (in OG coordinates)
		Vector<double> Ptk(3);
		Ptk[0] = OG_XYZ[0] - mTK_OO[0];
		Ptk[1] = OG_XYZ[1] - mTK_OO[1];
		Ptk[2] = OG_XYZ[2] - mTK_OO[2];

		// Projection of Ptk on the axis of the tokamak coordinate system 
		Vector<double> TK_XYZ(3);
		TK_XYZ[0] = Ptk[0] * mTK_t[0] + Ptk[1] * mTK_t[1] + Ptk[2] * mTK_t[2];
		TK_XYZ[1] = Ptk[0] * mTK_b[0] + Ptk[1] * mTK_b[1] + Ptk[2] * mTK_b[2];
		TK_XYZ[2] = Ptk[0] * mTK_n[0] + Ptk[1] * mTK_n[1] + Ptk[2] * mTK_n[2];

		return TK_XYZ;
	}

	//***************************************************************************************
	//* - Transform node coordinates to incR (node distance from Lo[0])
	//***************************************************************************************
	double ColdPlasma::Node_XYZ_to_incR( Node::Pointer pNode )
	{
        // Obtain node coordinates in the Tokamak cartesian coordinate system  
		Vector<double> TK_XYZ = Node_XYZ_to_TK_XYZ( pNode );

		// Distance to the Tokamak axis
		double R;

		// Considering toroidal curvature or not
		if (mFlatGeometry == true)
		{
			R = std::abs(TK_XYZ[0]);
		}
		else
		{
			R = std::sqrt(TK_XYZ[0] * TK_XYZ[0] + TK_XYZ[1] * TK_XYZ[1]);
		}

		// Adding poloidal curvature
		if (mPoloidalCurv > 0.0)
		{
			double Z    = TK_XYZ[2];
			double Ro   = mR_start - mPoloidalCurv;
			double Zo   = mZ_start;
			double dist = std::sqrt( (R - Ro) * (R - Ro) + (Z - Zo) * (Z - Zo) );

			R = Ro + dist;
		}

		// return incR (node distance from Lo[0])
        return ( mR_start - R );
	}

	//****************************************************************
	//* - Get the lower and upper index in Lo containing incR
	//*   index_in_Lo[0] = lower_index
	//*   index_in_Lo[1] = upper_index
	//****************************************************************
	Vector<int> ColdPlasma::Get_Indexes_In_Lo ( double incR )
	{
	    // Looking for the closest values of incR in mLo[i] 
		std::vector<double>::iterator it_mLo;

		// Obtaining the upper index (it_mLo - mLo.begin())
		it_mLo = std::lower_bound (mLo.begin(), mLo.end(), incR );

		// Indexes of mLo containing incR between them
		Vector<int> index_in_Lo(2);

		if (it_mLo == mLo.begin())
		{
			index_in_Lo[0] = 0;
			index_in_Lo[1] = 0;
		}
		else if (it_mLo == mLo.end())
		{
			index_in_Lo[0] = mLo.size() - 1;
			index_in_Lo[1] = mLo.size() - 1;
		}
		else
		{
			index_in_Lo[0] = it_mLo - mLo.begin() - 1;
			index_in_Lo[1] = it_mLo - mLo.begin();
		}	
		
		return index_in_Lo;
	}

	//****************************************************************************************
	// - Returns all the plasma parameters in a node:
	//      PlasmaParameters[0] = Electron density
	//      PlasmaParameters[1] = Bx
    //      PlasmaParameters[2] = By
	//      PlasmaParameters[3] = Bz
	//      PlasmaParameters[4] = S
	//      PlasmaParameters[5] = D
	//      PlasmaParameters[6] = P
	//      PlasmaParameters[7] = R
	//      PlasmaParameters[8] = L
	//*****************************************************************************************
	Vector<double> ColdPlasma::Get_AllPlasmaParameters_InNode( Node::Pointer pNode, double wFreq )
	{
		// Vector containing all the plasma parameters: eDensity, Bext, SDPRL
		Vector<double> PlasmaParameters(9);

		// Calculating incR = node distance from mLo[0]
		double incR = Node_XYZ_to_incR ( pNode );

		// Getting upper and lower indexes from mLo
		Vector<int> index = Get_Indexes_In_Lo( incR  );

		// Getting electron density from vector mNe
		double eDensity = Get_ElectronDensity_FromNe( incR, index );
		PlasmaParameters[0] = eDensity;

		// Gets Bext from vectors mB_mod, mB_ang
		Vector<double> Bext = Get_Bexternal_FromBmodBang( incR, index, pNode );
		PlasmaParameters[1] = Bext[0];
		PlasmaParameters[2] = Bext[1];
		PlasmaParameters[3] = Bext[2];

		// Calculates plasma permittivity tensor components
		Vector<double> SDPRL = Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
		PlasmaParameters[4] = SDPRL[0];
		PlasmaParameters[5] = SDPRL[1];
		PlasmaParameters[6] = SDPRL[2];
		PlasmaParameters[7] = SDPRL[3];
		PlasmaParameters[8] = SDPRL[4];

        return PlasmaParameters;	  
	}

	//*****************************************************************************************************
	//* - Obatins the S, P, D, R, L components of the permittivity tensor in a node
	//*****************************************************************************************************
	Vector<double> ColdPlasma::Get_PlasmaPermittivity_SDPRL_InNode( Node::Pointer pNode, double wFreq )
	{
		// Getting electron density distribution 
		double eDensity = Get_ElectronDensity_InNode( pNode );

		// Getting magnetic field modulus
		Vector<double> Bext = Get_Bexternal_InNode( pNode );

		// Getting S, P, D, R, L components of the permittivity tensor
		return Calculate_PlasmaPermittivity_SDPRL( eDensity, Bext, wFreq );
	}

	//*******************************************************************************************************************
	//* - Calculates the S, P, D, R, L components of the permittivity tensor
	//*******************************************************************************************************************
	Vector<double> ColdPlasma::Calculate_PlasmaPermittivity_SDPRL( double eDensity, Vector<double>& Bext, double wFreq )
	{
		// Module of Bext
		double modBext  = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );
        
		// Useful constants
		double eo = 8.8542e-12;   // Vacuum electric permittivity
		double Qe = 1.6022e-19;   // Electron charge
		double Me = 9.1095e-31;   // Electron mass
        double Wf = wFreq;        // Problem frequency (angular)

        double MH1  = 1.6727e-27; // Hidrogen mass
        double MH2  = 3.3437e-27; // Deuterium mass
        double MH3  = 5.0076e-27; // Tritium mass

        double MHe4 =  6.6465e-27; // Helium 4 mass
        double MBe9 = 14.9651e-27; // Beryllium 9 mass

        // Number of species
        int numSpecies = 8;

		// Defining species properties (charge, mass and density)
		Vector<double> Qs(numSpecies);  
		Vector<double> Ms(numSpecies);  
		Vector<double> Ns(numSpecies);

		// Electron
		Qs[0] = -Qe;  
		Ms[0] =  Me;  
		Ns[0] =  eDensity;   

		// Hydrogen (H1)
		Qs[1] = +Qe;  
		Ms[1] =  MH1;  
		Ns[1] =  eDensity * mH1_div_ne;

        // Deuterium (H2)
		Qs[2] = +Qe;  
		Ms[2] =  MH2;  
		Ns[2] =  eDensity * mH2_div_ne; 

        // Tritium (H3)
		Qs[3] = +Qe;  
		Ms[3] =  MH3;  
		Ns[3] =  eDensity * mH3_div_ne; 

        // Helium 4 (He4), Q=+1
        Qs[4] = +Qe;  
        Ms[4] =  MHe4;  
        Ns[4] =  eDensity * mHe4_q1_div_ne; 

        // Helium 4 (He4), Q=+2
        Qs[5] = +2.0*Qe;  
        Ms[5] =  MHe4;  
        Ns[5] =  eDensity * mHe4_q2_div_ne; 

        // Beryllium 9 (Be9), Q=+1
        Qs[6] = +Qe;  
        Ms[6] =  MBe9;  
        Ns[6] =  eDensity * mBe9_q1_div_ne; 

        // Beryllium 9 (Be9), Q=+2
        Qs[7] = +2.0*Qe;  
        Ms[7] =  MBe9;  
        Ns[7] =  eDensity * mBe9_q2_div_ne; 

        // Computing electric permittivity tensor parameters
		Vector<double> W2_ps(numSpecies); 
		for (int i=0; i<numSpecies; i++) 
		{
			W2_ps[i] = Qs[i] * Qs[i] * Ns[i] / ( eo * Ms[i] );
		}

		Vector<double> Omega_s(numSpecies);
		for (int i=0; i<numSpecies; i++) 
		{
			Omega_s[i] = Qs[i] * modBext / Ms[i];
		}

		double R = 1.0; 
		for (int i=0; i<numSpecies; i++) 
		{
			if ( std::abs( Wf + Omega_s[i] ) > 0.0 ) 
			{
				R -= W2_ps[i] / ( Wf * ( Wf + Omega_s[i] ) );
			}
			else 
			{
				Send_Error_Msg("R plasma parameter infinite", "Please, please check input data.");
			}
		}

		double L = 1.0; 
		for (int i=0; i<numSpecies; i++) 
		{
			if ( std::abs( Wf - Omega_s[i] ) > 0.0 )
			{
				L -= W2_ps[i] / ( Wf * ( Wf - Omega_s[i] ) );
			}
			else 
			{
				Send_Error_Msg("L plasma parameter infinite", "Please, please check input data.");
			}
		}

		double P = 1.0; 
		for (int i=0; i<numSpecies; i++) 
		{
			P -= W2_ps[i] / ( Wf * Wf );
		}

		// Components S, D and P of the plasma permittivity tensor (Stix notation) 
		// Components R, L and P from the diagonal plasma permittivity tensor
		Vector<double> SDPRL(5);
		SDPRL[0] = 0.5 * ( R + L );
	    SDPRL[1] = 0.5 * ( R - L );
	    SDPRL[2] = P;
		SDPRL[3] = R;
		SDPRL[4] = L;

	    return SDPRL;
	}

	//****************************************************************
	//* - Calculates electron density at the given node
	//****************************************************************
	double ColdPlasma::Get_ElectronDensity_InNode( Node::Pointer pNode )
	{
		// Calculating incR = node distance from mLo[0]
		double incR = Node_XYZ_to_incR ( pNode );

		// Getting upper and lower indexes from mLo
		Vector<int> index = Get_Indexes_In_Lo( incR  );

		// Getting electron density from vector mNe 
		return Get_ElectronDensity_FromNe( incR, index );
	}

	//****************************************************************
	//* - Extracts electron density from vector mNe
	//****************************************************************
	double ColdPlasma::Get_ElectronDensity_FromNe( double incR, Vector<int>& index )
	{
		// Calculating electron density from mNe (linear interpolation)
		double eDensityInNode;

		if (index[0] == index[1])
		{
			eDensityInNode = mNe[index[0]];
		}
		else
		{
			double   alpha = ( incR - mLo[index[0]] ) / std::abs( mLo[index[1]] - mLo[index[0]] );
			eDensityInNode = alpha * mNe[index[1]] + (1.0 - alpha) * mNe[index[0]];
		}
			
		return eDensityInNode;
	}

	//****************************************************************
	//* - Calculates applied external B field at the given node
	//****************************************************************
	Vector<double> ColdPlasma::Get_Bexternal_InNode( Node::Pointer pNode )
	{
		// Calculating incR = node distance from mLo[0]
		double incR = Node_XYZ_to_incR( pNode );

		// Getting upper and lower indexes from mLo
		Vector<int> index = Get_Indexes_In_Lo( incR );

		// Getting Bext from vectors mB_mod and mB_ang
		return Get_Bexternal_FromBmodBang( incR, index, pNode );
	}

	//****************************************************************
	//* - Extracts Bext from vectors mB_mod, mB_ang
	//****************************************************************
	Vector<double> ColdPlasma::Get_Bexternal_FromBmodBang( double incR, Vector<int>& index, Node::Pointer pNode )
	{
		// Calculating module and phase of B from mB_mod and mB_ang (linear interpolation)
		double modB, angB;

		if (index[0] == index[1])
		{
			modB = mB_mod[index[0]];
			angB = mB_ang[index[0]];
		}
		else
		{
			double alpha = ( incR - mLo[index[0]] ) / std::abs (mLo[index[1]] - mLo[index[0]] );
			modB = alpha * mB_mod[index[1]] + (1.0 - alpha) * mB_mod[index[0]];
			angB = alpha * mB_ang[index[1]] + (1.0 - alpha) * mB_ang[index[0]];
		}

		return TK_Bfield_to_XYZ( pNode, modB, angB );
	}

	//****************************************************************
	//* - Calculates B field in XYZ coordinate system
	//****************************************************************
	Vector<double> ColdPlasma::TK_Bfield_to_XYZ( Node::Pointer pNode, double modB, double angB )
	{
		// B in two different coordinate systems (TK and XYZ) 
		Vector<double> B_ntb(3), B_xyz(3);

		// Transforming angle in degrees to radians
		double pi = 3.1415926535897932384626433832795;
		double angB_rad = angB * pi / 180;

		// B field in curved geometry mode
		if (mCurvedGeometry == true)
		{
			Vector<double> N_tbn = Node_XYZ_to_TK_XYZ( pNode );

			double modN_tb = std::sqrt( N_tbn[0]*N_tbn[0] + N_tbn[1]*N_tbn[1] );
            double modB_tb = modB * std::cos(angB_rad);

			double inv_mod; 
			if (modN_tb != 0.0) inv_mod = 1.0 / modN_tb;
			else                inv_mod = 0.0;

			B_ntb[0] = + modB_tb * N_tbn[1] * inv_mod;
			B_ntb[1] = - modB_tb * N_tbn[0] * inv_mod;
			B_ntb[2] = + modB * std::sin(angB_rad);
		}
		// B field in flat geometry mode 
		else
		{
			B_ntb[0] = + 0.0;
			B_ntb[1] = - modB * std::cos(angB_rad);
			B_ntb[2] = + modB * std::sin(angB_rad);
		}

		// Adding poloidal curvature
		if (mPoloidalCurv > 0.0) Add_PoloidalCurv_to_TKBfield( pNode, B_ntb );
		
		// Transforming B in ntb CS to xyz CS
		B_xyz[0] = B_ntb[0] * mTK_t[0] + B_ntb[1] * mTK_b[0] + B_ntb[2] * mTK_n[0];
		B_xyz[1] = B_ntb[0] * mTK_t[1] + B_ntb[1] * mTK_b[1] + B_ntb[2] * mTK_n[1];
		B_xyz[2] = B_ntb[0] * mTK_t[2] + B_ntb[1] * mTK_b[2] + B_ntb[2] * mTK_n[2];

		return B_xyz;
	}

	//*********************************************************************
	//* - Adds poloidal curvature to the applied B field
	//*********************************************************************
	void ColdPlasma::Add_PoloidalCurv_to_TKBfield( Node::Pointer pNode, Vector<double>& B_ntb )
	{
		// Node coordinates in TK system
		Vector<double> rNode = Node_XYZ_to_TK_XYZ( pNode );

		// Coordinates of the curvature center in TK sytem
		Vector<double> rCurC(3);

		// Curved geometry mode
		if (mCurvedGeometry == true)
		{
			double modN_tb = std::sqrt(rNode[0]*rNode[0] + rNode[1]*rNode[1]);
			
			double inv_modN_tb;
			if (modN_tb != 0.0) inv_modN_tb = 1.0 / modN_tb;
			else                inv_modN_tb = 0.0;

			rCurC[0] = (mR_start - mPoloidalCurv) * rNode[0] * inv_modN_tb;
			rCurC[1] = (mR_start - mPoloidalCurv) * rNode[1] * inv_modN_tb;
			rCurC[2] = (mZ_start);
		}
		// Flat geometry mode
		else
		{
			rNode[1] = 0.0;

			rCurC[0] = mR_start - mPoloidalCurv;
			rCurC[1] = 0.0;
			rCurC[2] = mZ_start;
		}

		// Position vector from the curvature center to the node
		Vector<double> rN_CC(3);
		rN_CC[0] = rNode[0] - rCurC[0];
		rN_CC[1] = rNode[1] - rCurC[1];
		rN_CC[2] = rNode[2] - rCurC[2];

		// Perpendicular vector to rN_CC vector
		Vector<double> P_rN_CC(3);

        double modrNCC_tb = std::sqrt( rN_CC[0]*rN_CC[0] + rN_CC[1]*rN_CC[1] );

		if ( rN_CC[0] > 0.0 ) 
		{
			P_rN_CC[0] = -rN_CC[2] * rN_CC[0] / modrNCC_tb;
		    P_rN_CC[1] = -rN_CC[2] * rN_CC[1] / modrNCC_tb;
			P_rN_CC[2] = +modrNCC_tb;  
		}
		else if ( rN_CC[0] < 0.0 ) 
		{      
			P_rN_CC[0] = +rN_CC[2] * rN_CC[0] / modrNCC_tb;
		    P_rN_CC[1] = +rN_CC[2] * rN_CC[1] / modrNCC_tb;
			P_rN_CC[2] = -modrNCC_tb;  
		}
		else
		{
			double modrCurC_tb = std::sqrt( rCurC[0]*rCurC[0] + rCurC[1]*rCurC[1] );

			double inv_modrCurC_tb;
			if (modrCurC_tb != 0.0) inv_modrCurC_tb = 1.0 / modrCurC_tb;
			else                    inv_modrCurC_tb = 0.0;

			P_rN_CC[0] = -rCurC[0] * inv_modrCurC_tb;
		    P_rN_CC[1] = -rCurC[1] * inv_modrCurC_tb;
			P_rN_CC[2] =  0.0;  
		}

		// Normalization of the P_rN_CC vector
		double modPNCC = std::sqrt( P_rN_CC[0]*P_rN_CC[0] + P_rN_CC[1]*P_rN_CC[1] + P_rN_CC[2]*P_rN_CC[2] );
		
		double inv_modPNCC;
	    if ( modPNCC > 0.0 ) inv_modPNCC = 1.0 / modPNCC;
		else                 inv_modPNCC = 0.0;

		// Projection of Bz on the normalized vector P_rN_CC
		B_ntb[0]+= B_ntb[2] * P_rN_CC[0] * inv_modPNCC;
		B_ntb[1]+= B_ntb[2] * P_rN_CC[1] * inv_modPNCC;
		B_ntb[2] = B_ntb[2] * P_rN_CC[2] * inv_modPNCC;
	}

	//*********************************************************************
	//* - Reads electron density file 
	//*********************************************************************
	void ColdPlasma::Read_eDensity_File( String FileName )
	{
		std::fstream eDensityFile(FileName, std::fstream::in);

		if ( !eDensityFile.is_open() ) 
		{ 
			Send_Error_Msg("Error opening electron density file", "Please, check file path."); 
		}
	
		String FileInfo;

		// Pulse number and time
		eDensityFile >> FileInfo; eDensityFile >> FileInfo; 
		// R start
		eDensityFile >> FileInfo; eDensityFile >> mR_start;
		// Z start
		eDensityFile >> FileInfo; eDensityFile >> mZ_start;
		// Legend info
		eDensityFile >> FileInfo; eDensityFile >> FileInfo;

		// Reading distance and electron density
		double Lo, Ne;
		while (true)
		{
			eDensityFile >> Lo; 
			eDensityFile >> Ne; 
			
			if ( eDensityFile.eof() ) break;

            mLo.push_back( std::abs(Lo) );
            mNe.push_back( std::abs(Ne) );	
		}

		eDensityFile.close();
	}

	//*********************************************************************
	//* - Reads applied B field file
	//*********************************************************************
	void ColdPlasma::Read_aBfield_File( String FileName )
	{
		std::fstream aBfieldFile(FileName, std::fstream::in);

		if ( !aBfieldFile.is_open() )
		{
			Send_Error_Msg("Error opening applied B field file", "Please, check file path.");
		}
		
		String FileInfo;

		// Pulse number and time
		aBfieldFile >> FileInfo; aBfieldFile >> FileInfo;
		// R start
		aBfieldFile >> FileInfo; aBfieldFile >> mR_start;
		// Z start
		aBfieldFile >> FileInfo; aBfieldFile >> mZ_start;
		// Legend info
		aBfieldFile >> FileInfo; aBfieldFile >> FileInfo; aBfieldFile >> FileInfo;

		// Reading distance and electron density
		double Lo, modB, angB;
		while (true)
		{
			aBfieldFile >> Lo;
			aBfieldFile >> modB;
			aBfieldFile >> angB;

			if ( aBfieldFile.eof() ) break;

			mB_mod.push_back( std::abs(modB) );
			mB_ang.push_back( std::abs(angB) );
		}

		aBfieldFile.close();
	}

	//***********************************************
	//* - Print error messages and exit ERMES
	//***********************************************
	void ColdPlasma::Send_Error_Msg( String ErrorMsg, String CheckMsg )
	{
		std::cout << std::endl;
		std::cout << "!!!!!!!!!! " << ErrorMsg << " !!!!!!!!!!" << std::endl;
		std::cout << std::endl;
		std::cout << CheckMsg << std::endl;
		std::cout << std::endl;
		std::cout << "ERMES analysis finished." << std::endl;
		exit(EXIT_SUCCESS);
	}

	//***********************************************
	//* - Clearing Cold plasma properties
	//***********************************************
	ColdPlasma::~ColdPlasma()
	{
		mTK_OO.clear(); mTK_OO.swap(Vector<double>(mTK_OO));
		
		mTK_n.clear(); mTK_n.swap(Vector<double>(mTK_n));
		mTK_t.clear(); mTK_t.swap(Vector<double>(mTK_t));
		mTK_b.clear(); mTK_b.swap(Vector<double>(mTK_b));

		mLo.clear(); mLo.swap(Vector<double>(mLo));
		mNe.clear(); mNe.swap(Vector<double>(mNe));

		mB_mod.clear(); mB_mod.swap(Vector<double>(mB_mod));
		mB_ang.clear(); mB_ang.swap(Vector<double>(mB_ang));
	}

} /* end namespace Kratos */