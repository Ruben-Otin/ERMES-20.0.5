
#include "JSource_2nd_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{ 
	//*****************************************************************************************************************
    //*******************************************  Cartesian mode  ****************************************************
    //*****************************************************************************************************************  
    // 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = Module of the cartesian current X component.
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Phase  of the cartesian current X component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Module of the cartesian current Y component.
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = Phase  of the cartesian current Y component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Module of the cartesian current Z component.
    // SINUSOIDAL_SURFACE_CURRENT[ 5 ] = Phase  of the cartesian current Z component.
    //  
    //*****************************************************************************************************************
    //*******************************************  Axysimmetric mode  *************************************************
    //*****************************************************************************************************************  
    //
    // COMPLEX_IBC[ 0 ] = Axisymmetric current module.
    // COMPLEX_IBC[ 1 ] = Axisymmetric current phase.
    //
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of curvature.    
    //
    // COMPLEX_NEUMANN_FLOW[ 3 ] = X component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Y component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 5 ] = Z component of the curvature axis vector.
    // 
    //*****************************************************************************************************************
    //*******************************************  Helicoidal mode  ***************************************************
    //*****************************************************************************************************************    
    // 
    // COMPLEX_IBC_2o[ 0 ] = Toroid major radius.
	// 
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of curvature.    
    //
    // COMPLEX_NEUMANN_FLOW[ 3 ] = X component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Y component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 5 ] = Z component of the curvature axis vector.
	// 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = Module of the poloidal current component.
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Phase  of the poloidal current component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Module of the toroidal current component.
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = Phase  of the toroidal current component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 5 ] = 0 
	//                                  
    //*****************************************************************************************************************
    //*******************************************  Plasma mode  *******************************************************
    //*****************************************************************************************************************  
	//  
    // COMPLEX_IBC_2o[ 0 ] = Toroid major radius.
    // 
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of curvature.
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of curvature.    
    //
    // COMPLEX_NEUMANN_FLOW[ 3 ] = X component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Y component of the curvature axis vector.
    // COMPLEX_NEUMANN_FLOW[ 5 ] = Z component of the curvature axis vector.
	// 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = Module of the poloidal current component.
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Phase  of the poloidal current component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Module of the toroidal current component.
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = Phase  of the toroidal current component.
	// 
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Plasma mode index. 
	// SINUSOIDAL_SURFACE_CURRENT[ 5 ] = 1 
	// 
    //*****************************************************************************************************************
	
    //*****************************************************************************************************************
    //* - Global index of the nodes
    //*****************************************************************************************************************	
    void JSource_2nd_FullWave::GetEquationIdVector( Vector<int>& EquationId )
    {
        EquationId.resize( mNumDofs );

		if( mPotentialsOn )
		{
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                EquationId[ i               ] = mNodes[ i ]->pDofcAx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcAy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcAz()->EquationId();
            }				
		}
		else
		{
            for( int i=0; i<mNumNodes; i++ ) 
            { 
                EquationId[ i               ] = mNodes[ i ]->pDofcEx()->EquationId();
                EquationId[ i +   mNumNodes ] = mNodes[ i ]->pDofcEy()->EquationId();
		    	EquationId[ i + 2*mNumNodes ] = mNodes[ i ]->pDofcEz()->EquationId();
            }		
		}
    }

    //*****************************************************************************************************************
    //* - Element volume
    //*****************************************************************************************************************
    double JSource_2nd_FullWave::Calculate_Volume()
    {
        double Det;

        Det = + X(2)*Y(3)*Z(4) + X(4)*Y(2)*Z(3) + X(3)*Y(4)*Z(2) - X(4)*Y(3)*Z(2) - X(2)*Y(4)*Z(3) - X(3)*Y(2)*Z(4) 
              - X(1)*Y(3)*Z(4) - X(4)*Y(1)*Z(3) - X(3)*Y(4)*Z(1) + X(4)*Y(3)*Z(1) + X(1)*Y(4)*Z(3) + X(3)*Y(1)*Z(4) 
              + X(1)*Y(2)*Z(4) + X(4)*Y(1)*Z(2) + X(2)*Y(4)*Z(1) - X(4)*Y(2)*Z(1) - X(1)*Y(4)*Z(2) - X(2)*Y(1)*Z(4) 
              - X(1)*Y(2)*Z(3) - X(3)*Y(1)*Z(2) - X(2)*Y(3)*Z(1) + X(3)*Y(2)*Z(1) + X(1)*Y(3)*Z(2) + X(2)*Y(1)*Z(3);
        
        return fabs( Det / 6.0 );
    }

	//*****************************************************************************************************************
    //* - Cartesian sources
    //*****************************************************************************************************************
    void JSource_2nd_FullWave::Get_CartesianXYZ_Sources( Vector< std::complex<double> >& cJx, 
		                                                 Vector< std::complex<double> >& cJy, 
		                                                 Vector< std::complex<double> >& cJz ) 
    {
		// Cartesian source
        Vector<double> vJ = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT );

		if( ( vJ[ 0 ] != 0.0 ) || ( vJ[ 2 ] != 0.0 ) || ( vJ[ 4 ] != 0.0 ) ) 
		{
			std::complex<double> Element_cJx( vJ[ 0 ] * cos( vJ[ 1 ] ), vJ[ 0 ] * sin( vJ[ 1 ] ) );
			std::complex<double> Element_cJy( vJ[ 2 ] * cos( vJ[ 3 ] ), vJ[ 2 ] * sin( vJ[ 3 ] ) );
			std::complex<double> Element_cJz( vJ[ 4 ] * cos( vJ[ 5 ] ), vJ[ 4 ] * sin( vJ[ 5 ] ) );

			for( int i=0; i<mNumVerts; i++ )
			{
				cJx[ i ] += Element_cJx;
				cJy[ i ] += Element_cJy;
				cJz[ i ] += Element_cJz;
			}
		}

		if( mNodal_J.size() > 0 )
		{
			for( int i=0; i<mNumVerts; i++ )
			{
				cJx[ i ] += mNodal_J[ i ][ 0 ];
				cJy[ i ] += mNodal_J[ i ][ 1 ];
				cJz[ i ] += mNodal_J[ i ][ 2 ];
			}		
		}
	}

	//*****************************************************************************************************************
    //* - Axisymmetric sources
    //*****************************************************************************************************************
    void JSource_2nd_FullWave::Get_Axisymmetric_Sources( Vector< std::complex<double> >& cJx, 
		                                                 Vector< std::complex<double> >& cJy, 
		                                                 Vector< std::complex<double> >& cJz ) 
    {
        // Axisymmetric source
		Vector<double> vJa = (*mProperties)( COMPLEX_IBC );

		if( vJa[ 0 ] == 0.0 )
		{
		    return;
		}

		// Complex axysimmetric current sources
		std::complex<double> cJa( vJa[ 0 ] * cos( vJa[ 1 ] ), vJa[ 0 ] * sin( vJa[ 1 ] ) );

        // Origin of coordinates
		double ox = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ];
		double oy = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ];
		double oz = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 2 ];

		// Rotation axis
		double ax = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ];
		double ay = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];
		double az = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ];

		// Transform sources to cartesian and sent to nodes
		for( int i=0; i<mNumVerts; i++ )
		{
			// Position vector 
			double px = mNodes[ i ]->X() - ox;
			double py = mNodes[ i ]->Y() - oy;
			double pz = mNodes[ i ]->Z() - oz;
               
			// r = a x p 
			double rx = ay*pz - az*py;
			double ry = az*px - ax*pz;
			double rz = ax*py - ay*px;

			double norm = sqrt( rx*rx + ry*ry + rz*rz );

			if( norm > 0.0 )
			{
				rx = rx / norm;
				ry = ry / norm;
				rz = rz / norm;
			}
			else
			{
				rx = 0.0;
				ry = 0.0;
				rz = 0.0;
			}

			cJx[ i ] += cJa * rx;
			cJy[ i ] += cJa * ry;
			cJz[ i ] += cJa * rz;
		}
	}

	//*****************************************************************************************************************
    //* - Plasma modes J sources
    //*****************************************************************************************************************
    void JSource_2nd_FullWave::Get_PlasmaModeJs_Sources( Vector< std::complex<double> >& cJx, 
		                                                 Vector< std::complex<double> >& cJy, 
		                                                 Vector< std::complex<double> >& cJz ) 
    {
		// Torus major radius
		double Ro = (*mProperties)( COMPLEX_IBC_2o )[ 0 ];
		
		if( Ro <= 0.0 ) 
		{
		    return;
		}

		// Plasma mode J current sources 
        Vector<double> vJp = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT );

		if( ( vJp[ 0 ] == 0.0 ) && ( vJp[ 2 ] == 0.0 ) ) 
		{
		    return;
		}

		// Complex poloidal and toroidal current 
        std::complex<double> cJpol( vJp[ 0 ] * cos( vJp[ 1 ] ), vJp[ 0 ] * sin( vJp[ 1 ] ) );
		std::complex<double> cJtor( vJp[ 2 ] * cos( vJp[ 3 ] ), vJp[ 2 ] * sin( vJp[ 3 ] ) );

        // Origin of coordinates
		Vector<double> OO_A( 3 );

		OO_A[ 0 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ];
		OO_A[ 1 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ];
		OO_A[ 2 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 2 ];

	    // Rotation axis
		Vector<double> nAxi( 3 );

		nAxi[ 0 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ];
		nAxi[ 1 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];
		nAxi[ 2 ] = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ];

	    nAxi *= ( 1.0 / sqrt( nAxi[ 0 ]*nAxi[ 0 ] + nAxi[ 1 ]*nAxi[ 1 ] + nAxi[ 2 ]*nAxi[ 2 ] ) );

		// Phi reference axis
		Vector<double> nRef( 3 );

		if( sqrt( nAxi[ 0 ]*nAxi[ 0 ] + nAxi[ 1 ]*nAxi[ 1 ] ) > 0.0 )
        {
            nRef[ 0 ] = nAxi[ 1 ];
		    nRef[ 1 ] =-nAxi[ 0 ];
		    nRef[ 2 ] = 0.0;
        }
        else
        {
            nRef[ 0 ] = nAxi[ 2 ];
            nRef[ 1 ] = 0.0;
            nRef[ 2 ] = 0.0;
        }

		nRef *= ( 1.0 / sqrt( nRef[ 0 ]*nRef[ 0 ] + nRef[ 1 ]*nRef[ 1 ] + nRef[ 2 ]*nRef[ 2 ] ) );

		// Phi reference axis ( ortogonal )
		Vector<double> nOef( 3 );

		nOef[ 0 ] = nAxi[ 1 ] * nRef[ 2 ] - nAxi[ 2 ] * nRef[ 1 ];
		nOef[ 1 ] = nAxi[ 2 ] * nRef[ 0 ] - nAxi[ 0 ] * nRef[ 2 ];
		nOef[ 2 ] = nAxi[ 0 ] * nRef[ 1 ] - nAxi[ 1 ] * nRef[ 0 ];

		nOef *= ( 1.0 / sqrt( nOef[ 0 ]*nOef[ 0 ] + nOef[ 1 ]*nOef[ 1 ] + nOef[ 2 ]*nOef[ 2 ] ) );

		// Transform to cartesian coordinates and sent to nodes
		for( int i=0; i<mNumVerts; i++ )
		{
			// Position vector 
			Vector<double> OO_P( 3 );

			OO_P[ 0 ] = mNodes[ i ]->X() - OO_A[ 0 ];
			OO_P[ 1 ] = mNodes[ i ]->Y() - OO_A[ 1 ];
			OO_P[ 2 ] = mNodes[ i ]->Z() - OO_A[ 2 ];
               
			// Toroidal vector ( nTor = nAxi x OO_P )
			Vector<double> nTor( 3 );

			nTor[ 0 ] = nAxi[ 1 ] * OO_P[ 2 ] - nAxi[ 2 ] * OO_P[ 1 ];
			nTor[ 1 ] = nAxi[ 2 ] * OO_P[ 0 ] - nAxi[ 0 ] * OO_P[ 2 ];
			nTor[ 2 ] = nAxi[ 0 ] * OO_P[ 1 ] - nAxi[ 1 ] * OO_P[ 0 ];

			nTor *= ( 1.0 / sqrt( nTor[ 0 ]*nTor[ 0 ] + nTor[ 1 ]*nTor[ 1 ] + nTor[ 2 ]*nTor[ 2 ] ) );

			// Radial vector ( nRad = nTor x nAxi )
			Vector<double> nRad( 3 );

			nRad[ 0 ] = nTor[ 1 ] * nAxi[ 2 ] - nTor[ 2 ] * nAxi[ 1 ];
			nRad[ 1 ] = nTor[ 2 ] * nAxi[ 0 ] - nTor[ 0 ] * nAxi[ 2 ];
			nRad[ 2 ] = nTor[ 0 ] * nAxi[ 1 ] - nTor[ 1 ] * nAxi[ 0 ];
			 
			nRad *= ( 1.0 / sqrt( nRad[ 0 ]*nRad[ 0 ] + nRad[ 1 ]*nRad[ 1 ] + nRad[ 2 ]*nRad[ 2 ] ) );

			// Inner radial vector ( InnR = OO_P - Ro * nRad )
			Vector<double> InnR( 3 );

			InnR[ 0 ] = OO_P[ 0 ] - Ro * nRad[ 0 ];
			InnR[ 1 ] = OO_P[ 1 ] - Ro * nRad[ 1 ];
			InnR[ 2 ] = OO_P[ 2 ] - Ro * nRad[ 2 ];

			// Poloidal vector ( nPol = nInR x nTor )
			Vector<double> nPol( 3 );

			nPol[ 0 ] = InnR[ 1 ] * nTor[ 2 ] - InnR[ 2 ] * nTor[ 1 ];
			nPol[ 1 ] = InnR[ 2 ] * nTor[ 0 ] - InnR[ 0 ] * nTor[ 2 ];
			nPol[ 2 ] = InnR[ 0 ] * nTor[ 1 ] - InnR[ 1 ] * nTor[ 0 ];
			 
			nPol *= ( 1.0 / sqrt( nPol[ 0 ]*nPol[ 0 ] + nPol[ 1 ]*nPol[ 1 ] + nPol[ 2 ]*nPol[ 2 ] ) );

			// Poloidal ( Theta ) and toroidal ( Phi ) currents
			std::complex<double> J_Theta = cJpol;
			std::complex<double> J_Phi   = cJtor;

		    // Add plasma mode M 
			if( (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 5 ] == 1 )
			{
			    // Plasma mode index
		        int M = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 4 ];

				// Minor radius
			    double rt = sqrt( InnR[ 0 ]*InnR[ 0 ] + InnR[ 1 ]*InnR[ 1 ] + InnR[ 2 ]*InnR[ 2 ] );

			    // Theta angle
				double Y_Theta = InnR[ 0 ] * nAxi[ 0 ] + InnR[ 1 ] * nAxi[ 1 ] + InnR[ 2 ] * nAxi[ 2 ];
				double X_Theta = InnR[ 0 ] * nRad[ 0 ] + InnR[ 1 ] * nRad[ 1 ] + InnR[ 2 ] * nRad[ 2 ]; 
				
				double Theta = atan2( Y_Theta, X_Theta ); 

				// Phi angle
				double XY_Prj = OO_P[ 0 ] * nRad[ 0 ] + OO_P[ 1 ] * nRad[ 1 ] + OO_P[ 2 ] * nRad[ 2 ];

				double Y_Phi = XY_Prj * ( nRad[ 0 ] * nOef[ 0 ] + nRad[ 1 ] * nOef[ 1 ] + nRad[ 2 ] * nOef[ 2 ] );
				double X_Phi = XY_Prj * ( nRad[ 0 ] * nRef[ 0 ] + nRad[ 1 ] * nRef[ 1 ] + nRad[ 2 ] * nRef[ 2 ] ); 
				
				double Phi = atan2( Y_Phi, X_Phi ); 
			    
				// Plasma mode
			    J_Theta *= sin( M * Theta + Phi );
			    J_Phi   *= - M * ( ( Ro / rt ) + cos( Theta ) ) * sin( M * Theta + Phi ) - sin( Theta ) * cos( M * Theta + Phi );
			}

			cJx[ i ] += J_Theta * nPol[ 0 ] + J_Phi * nTor[ 0 ];
			cJy[ i ] += J_Theta * nPol[ 1 ] + J_Phi * nTor[ 1 ];
			cJz[ i ] += J_Theta * nPol[ 2 ] + J_Phi * nTor[ 2 ];
		}
	}

    //*****************************************************************************************************************
    //* - Imposed J currents on element nodes
    //*****************************************************************************************************************
    void JSource_2nd_FullWave::Get_J_Imposed_OnNodes( std::map< unsigned int, Vector<double> >& rJimp, 
		                                              std::map< unsigned int, Vector<double> >& iJimp )
	{
		// J sources on nodes
		Vector< std::complex<double> > cJx( mNumVerts, std::complex<double>( 0.0, 0.0 ) );
		Vector< std::complex<double> > cJy( mNumVerts, std::complex<double>( 0.0, 0.0 ) );
		Vector< std::complex<double> > cJz( mNumVerts, std::complex<double>( 0.0, 0.0 ) );

		// Get cartesian and axisymmetric sources 
		if( (*mProperties)( COMPLEX_IBC_2o )[ 0 ] == 0 ) 
		{
		    Get_CartesianXYZ_Sources( cJx, cJy, cJz );
			Get_Axisymmetric_Sources( cJx, cJy, cJz );
		}
		// Get plasma mode currents sources
		else
		{
		    Get_PlasmaModeJs_Sources( cJx, cJy, cJz );
		}

		// Sent sources to map 
		for( int i=0; i<mNumVerts; i++ )
		{
			Vector<double> real_J_vector( 3 );
			
		    real_J_vector[ 0 ] = std::real( cJx[ i ] );
            real_J_vector[ 1 ] = std::real( cJy[ i ] );
            real_J_vector[ 2 ] = std::real( cJz[ i ] );

			Vector<double> imag_J_vector( 3 );

		    imag_J_vector[ 0 ] = std::imag( cJx[ i ] );
            imag_J_vector[ 1 ] = std::imag( cJy[ i ] );
            imag_J_vector[ 2 ] = std::imag( cJz[ i ] );

			int NodeId = mNodes[ i ]->Id();

			rJimp[ NodeId ] = real_J_vector;
			iJimp[ NodeId ] = imag_J_vector;
		}
	}

    //*****************************************************************************************************************
    //* - Residual vector
    //*****************************************************************************************************************
    void JSource_2nd_FullWave::GetResidualVector( Vector< std::complex<double> >& ResidualVector ) 
    {
		// Initializing residual vector
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

		// J sources on nodes
		Vector< std::complex<double> > cJx( mNumVerts, std::complex<double>( 0.0, 0.0 ) );
		Vector< std::complex<double> > cJy( mNumVerts, std::complex<double>( 0.0, 0.0 ) );
		Vector< std::complex<double> > cJz( mNumVerts, std::complex<double>( 0.0, 0.0 ) );

		// Get cartesian and axisymmetric sources 
		if( (*mProperties)( COMPLEX_IBC_2o )[ 0 ] == 0 )
		{
		    Get_CartesianXYZ_Sources( cJx, cJy, cJz );
			Get_Axisymmetric_Sources( cJx, cJy, cJz );
		}
		// Get plasma mode currents sources
		else
		{
		    Get_PlasmaModeJs_Sources( cJx, cJy, cJz );
		}

		// If potentials are off then the source term is equal to jwJ
		if( mPotentialsOn == false )
		{
		    std::complex<double> jw( 0.0, (*mProperties)( FREQUENCY ) );
		    
			cJx *= jw;
			cJy *= jw; 
			cJz *= jw; 
		}

		// Gaussian points and weights
		Vector<double> cX, cY, cZ, W;

		int numGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W );

		double detJ = 6.0 * mVolume;

		for( int gp=0; gp<numGaussPoints; gp++ ) 
		{
			W[ gp ] *= detJ; 
		}

		// 1st order Lagrangian base functions
		Matrix<double> N1st;

		Lagrange3D_Ni_1st( N1st, cX, cY, cZ );

        // 2nd order Lagrangian base functions
		Matrix<double> N2nd;

		Lagrange3D_Ni_2nd( N2nd, cX, cY, cZ );

		// Add sources to residual vector
		for( int i=0; i<mNumNodes; i++ )
		{
			for( int j=0; j<mNumVerts; j++ )
			{		
				double IntgVol_NiNj = 0.0;

                for( int gp=0; gp<numGaussPoints; gp++ )
				{
					IntgVol_NiNj += W[ gp ] * N2nd[ i ][ gp ] * N1st[ j ][ gp ];
				}
				
				ResidualVector[ i               ] += cJx[ j ] * IntgVol_NiNj;
				ResidualVector[ i +   mNumNodes ] += cJy[ j ] * IntgVol_NiNj;
				ResidualVector[ i + 2*mNumNodes ] += cJz[ j ] * IntgVol_NiNj;
			}
		}
    }
} 