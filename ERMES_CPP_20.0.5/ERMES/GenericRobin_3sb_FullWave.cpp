
#include "GenericRobin_3sb_FullWave.h"
#include "../includes/variables.h"

#include "Utils/GaussIntegrationTables.h"
#include "Utils/LagrangeBases.h"

namespace Kratos
{  	
    //**********************************************************************************************************************
    //*******************************************  Quasi-Static Mode  ******************************************************
    //**********************************************************************************************************************    
    // 
    // COMPLEX_IBC[ 0 ] = Real part of the generic Robin coefficient
    // COMPLEX_IBC[ 1 ] = Imag part of the generic Robin coefficient
    // 
    // COMPLEX_IBC_2o[ 0 ] = Flux module
    // COMPLEX_IBC_2o[ 1 ] = Flux phase
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 0
    //  
    //**********************************************************************************************************************
    //*******************************************  Plane Wave Mode  ********************************************************
    //**********************************************************************************************************************    
    //
    // COMPLEX_IBC[ 0 ] = Real part of the electric conductivity in S/m
    // COMPLEX_IBC[ 1 ] = Imag part of the electric conductivity in S/m
    // 
    // COMPLEX_IBC_2o[ 0 ] = Real part of the relative electric permittivity
    // COMPLEX_IBC_2o[ 1 ] = Imag part of the relative electric permittivity
    //  
    // COMPLEX_NEUMANN_FLOW[ 0 ] = Real part of the relative magnetic permeability
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Imag part of the relative magnetic permeability
    // 
    // COMPLEX_NEUMANN_FLOW[ 3 ] = Module of the E field parallel to the polarization vector
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Phase  of the E field parallel to the polarization vector
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 1 in single   plane wave  mode  
    //                           = 2 in multiple plane waves mode 
    //                            
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = X component of the polarization vector
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Y component of the polarization vector
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Z component of the polarization vector
	//
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = X component of the wave vector
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Y component of the wave vector
    // SINUSOIDAL_SURFACE_CURRENT[ 5 ] = Z component of the wave vector
    // 
    //**********************************************************************************************************************
    //*******************************************  Gaussian Beam Mode  *****************************************************
    //**********************************************************************************************************************    
    // 
    // COMPLEX_IBC[ 0 ] = Module of the E field perpendicular to the polarization and wave vectors at the waist center
    // COMPLEX_IBC[ 1 ] = Phase  of the E field perpendicular to the polarization and wave vectors at the waist center
    // 
    // COMPLEX_IBC_2o[ 0 ] = Gaussian beam waist radius on the direction parallel to the polarization vector
    // COMPLEX_IBC_2o[ 1 ] = Gaussian beam waist radius on the direction perpendicular to the polarization and wave vectors
    //   
    // COMPLEX_NEUMANN_FLOW[ 0 ] = X coordinate of the center of the Gaussian beam waist
    // COMPLEX_NEUMANN_FLOW[ 1 ] = Y coordinate of the center of the Gaussian beam waist
    // COMPLEX_NEUMANN_FLOW[ 2 ] = Z coordinate of the center of the Gaussian beam waist 
    // 
    // COMPLEX_NEUMANN_FLOW[ 3 ] = Module of the E field parallel to the polarization vector at the waist center
    // COMPLEX_NEUMANN_FLOW[ 4 ] = Phase  of the E field parallel to the polarization vector at the waist center
    //
    // COMPLEX_NEUMANN_FLOW[ 5 ] = 3
	// 
    // SINUSOIDAL_SURFACE_CURRENT[ 0 ] = X component of the polarization vector 
    // SINUSOIDAL_SURFACE_CURRENT[ 1 ] = Y component of the polarization vector 
	// SINUSOIDAL_SURFACE_CURRENT[ 2 ] = Z component of the polarization vector
	//
    // SINUSOIDAL_SURFACE_CURRENT[ 3 ] = X component of the wave vector 
	// SINUSOIDAL_SURFACE_CURRENT[ 4 ] = Y component of the wave vector 
    // SINUSOIDAL_SURFACE_CURRENT[ 5 ] = Z component of the wave vector 
    // 
    //**********************************************************************************************************************
    
    //**********************************************************************************************************************
    //* - Global IDs of the DOFs
    //**********************************************************************************************************************	
    void GenericRobin_3sb_FullWave::GetEquationIdVector( Vector<int>& EquationId )
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

            for( int i=0; i<3; i++ ) 
            { 
                EquationId[ i + 3*mNumNodes ] = mNodes[ i ]->pDofcVs()->EquationId();
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

    //**********************************************************************************************************************
    // - Area of the element
    //**********************************************************************************************************************
    double GenericRobin_3sb_FullWave::Calculate_Area()
    {
		double na[ 3 ], v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ 2 ]->X() - mNodes[ 0 ]->X();
        v2[ 1 ] = mNodes[ 2 ]->Y() - mNodes[ 0 ]->Y();
        v2[ 2 ] = mNodes[ 2 ]->Z() - mNodes[ 0 ]->Z();
		    
        v1[ 0 ] = mNodes[ 1 ]->X() - mNodes[ 0 ]->X();
        v1[ 1 ] = mNodes[ 1 ]->Y() - mNodes[ 0 ]->Y();
        v1[ 2 ] = mNodes[ 1 ]->Z() - mNodes[ 0 ]->Z();

        // Area = 0.5 * ||v2 x v1||
        na[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        na[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        na[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

		return ( 0.5 * sqrt( na[ 0 ]*na[ 0 ] + na[ 1 ]*na[ 1 ] + na[ 2 ]*na[ 2 ] ) );
    }

    //**********************************************************************************************************************
    //* - Exterior normal 
    //**********************************************************************************************************************	
    void GenericRobin_3sb_FullWave::Calculate_Exterior_Normal( Vector<double>& Next )
    {
		double v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = mNodes[ 2 ]->X() - mNodes[ 0 ]->X();
        v2[ 1 ] = mNodes[ 2 ]->Y() - mNodes[ 0 ]->Y();
        v2[ 2 ] = mNodes[ 2 ]->Z() - mNodes[ 0 ]->Z();
		    
        v1[ 0 ] = mNodes[ 1 ]->X() - mNodes[ 0 ]->X();
        v1[ 1 ] = mNodes[ 1 ]->Y() - mNodes[ 0 ]->Y();
        v1[ 2 ] = mNodes[ 1 ]->Z() - mNodes[ 0 ]->Z();

        Next.resize( 3 );

        Next[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        Next[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        Next[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

        Next *= ( 1.0 / sqrt( Next[ 0 ]*Next[ 0 ] + Next[ 1 ]*Next[ 1 ] + Next[ 2 ]*Next[ 2 ] ) );
    }

    //**********************************************************************************************************************
    //* - jk / mu = j * w * sqrt ( ep / mu )
    //**********************************************************************************************************************	
    std::complex<double> GenericRobin_3sb_FullWave::Calculate_jK_div_Mu()
    {
        // Gaussian beam 
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;
            
            return std::complex<double>( 0.0, mFreq * sqrt( eo / mo ) );
        }
        // Plane wave 
        else
        {
            if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] == 0.0 ) 
            {
                return std::complex<double>( 0.0, 0.0 );
            }            
            
            double wf = mFreq;
            double w2 = mFreq * mFreq;

            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;            
            
            double sg_real = (*mProperties)( COMPLEX_IBC )[ 0 ];
            double sg_imag = (*mProperties)( COMPLEX_IBC )[ 1 ];

            double ep_real = (*mProperties)( COMPLEX_IBC_2o )[ 0 ] * eo; 
            double ep_imag = (*mProperties)( COMPLEX_IBC_2o )[ 1 ] * eo; 

            double mu_real = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] * mo; 
            double mu_imag = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ] * mo; 

            std::complex<double> cUnit( 0.0                        , 1.0                         );
            std::complex<double> cMu  ( mu_real                    , mu_imag                     );
            std::complex<double> cEpw2( ep_real * w2 - sg_imag * wf, ep_imag * w2 + sg_real * wf );

            return ( cUnit * sqrt( cEpw2 / cMu ) );
        }
    }

    //**********************************************************************************************************************
    //* - jk = j * w * sqrt ( ep * mu )
    //**********************************************************************************************************************	
    std::complex<double> GenericRobin_3sb_FullWave::Calculate_jK()
    {
        // Gaussian beam 
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;
            
            return std::complex<double>( 0.0, mFreq * sqrt( eo * mo ) );
        }
        // Plane wave 
        else
        {
            double wf = mFreq;
            double w2 = mFreq * mFreq;

            double eo = 8.8541878176e-12;
            double mo = 3.14159265358979 * 4.0e-7;            
            
            double sg_real = (*mProperties)( COMPLEX_IBC )[ 0 ];
            double sg_imag = (*mProperties)( COMPLEX_IBC )[ 1 ];

            double ep_real = (*mProperties)( COMPLEX_IBC_2o )[ 0 ] * eo; 
            double ep_imag = (*mProperties)( COMPLEX_IBC_2o )[ 1 ] * eo; 

            double mu_real = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ] * mo; 
            double mu_imag = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ] * mo; 

            std::complex<double> cUnit( 0.0                        , 1.0                         );
            std::complex<double> cMu  ( mu_real                    , mu_imag                     );
            std::complex<double> cEpw2( ep_real * w2 - sg_imag * wf, ep_imag * w2 + sg_real * wf );
            
            return ( cUnit * sqrt( cEpw2 * cMu ) );
        }
    }

    //**********************************************************************************************************************
    //* - NiNj surface integrals matrix
    //**********************************************************************************************************************	
    void GenericRobin_3sb_FullWave::Calculate_NiNj_Matrix( Matrix<double>& NiNj_Matrix )
    {
        // Setting NiNj matrix
		NiNj_Matrix.Resize( mNumNodes, mNumNodes, 0.0 );

		// Gauss points
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		// Basis functions
        Matrix<double> N;

		Lagrange2D_Ni_3sb( N, cX, cY );

		// Jacobian
        double Jacob = 2.0 * mArea;

		// Surface integrals matrix
        for( int i=0; i<mNumNodes; i++ )
        {
        	for( int j=0; j<mNumNodes; j++ )
        	{
                for( int gp=0; gp<numGaussPoints; gp++ ) 
                {
                    NiNj_Matrix[ i ][ j ] += Jacob * W[ gp ] * N[ i ][ gp ] * N[ j ][ gp ];
                }
        	}
        }
	}

    //**********************************************************************************************************************
    //* - Residual vector and quasi-static flux values ( needed for total J sources calculations )
    //**********************************************************************************************************************
    void GenericRobin_3sb_FullWave::GetResidualVector_Flux( std::map< unsigned int, std::complex<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
    {
        ResidualVector.resize( mNumDofs, std::complex<double>( 0.0, 0.0 ) );

        if( ( mPotentialsOn == true ) && ( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 0.0 ) ) 
		{
            double Intg_Ni = mArea / 3.0;
            
            Vector<double> Coef = (*mProperties)( COMPLEX_IBC    );
            Vector<double> Flux = (*mProperties)( COMPLEX_IBC_2o );

            std::complex<double> cFlux( Flux[ 0 ] * cos( Flux[ 1 ] ), Flux[ 0 ] * sin( Flux[ 1 ] ) );

            for( int i=0; i<3; i++ ) 
            { 
                if( ( Coef[ 0 ] == 0.0 ) && ( Coef[ 1 ] == 0.0 ) ) 
                {
                    FluxValue[ mNodes[ i ]->Id() ] = cFlux;
                }
                
                ResidualVector[ i + 3*mNumNodes ] += cFlux * Intg_Ni;
            }
		}
    }

    //**********************************************************************************************************************
    //* - Residual vector for a single plane wave
    //**********************************************************************************************************************
    //***************************************************************************************************************************
    //* - Adds the plane wave source term of the Robin boundary condition ( curl and divergence parts ):
    //*   residual += jK/Mu * Int Ni ( C - D ) E dS,   E = cF exp( jK k.r ),   C = n x n x . - n x k x . ,   D = n ( n.  - k. )
    //*   The plane wave is evaluated at the Gauss points of the face ( previous versions evaluated it at the nodes and
    //*   interpolated it with the element basis, which under-resolves the phase exp( jK k.r ) on coarse meshes ). Same
    //*   treatment as the edge elements and the Gaussian beam source.
    //***************************************************************************************************************************
    void GenericRobin_3sb_FullWave::Add_PlaneWave_Source( Vector< std::complex<double> >& cF, Vector<double>& Kdir, Vector< std::complex<double> >& ResidualVector )
    {
        // Complex values jK and jK/Mu
        std::complex<double> jK        = Calculate_jK();
        std::complex<double> jK_div_Mu = Calculate_jK_div_Mu();

		// Element exterior normal
		Vector<double> Next; 
		
		Calculate_Exterior_Normal( Next );

		// k*n
		double KdirNext = Kdir[ 0 ]*Next[ 0 ] + Kdir[ 1 ]*Next[ 1 ] + Kdir[ 2 ]*Next[ 2 ];

		// Combined coefficient matrix C - D
        double CmD[ 3 ][ 3 ];

        for( int r=0; r<3; r++ )
        {
            for( int c=0; c<3; c++ )
            {
                double Cm = Next[ r ]*Next[ c ] - Kdir[ r ]*Next[ c ] - ( ( r == c ) ? ( 1.0 - KdirNext ) : 0.0 );
                double Dm = Next[ r ]*Next[ c ] - Next[ r ]*Kdir[ c ];

                CmD[ r ][ c ] = Cm - Dm;
            }
        }

        // ( C - D ) * cF ( constant on the face )
        std::complex<double> CmDF[ 3 ];

        for( int r=0; r<3; r++ ) CmDF[ r ] = CmD[ r ][ 0 ] * cF[ 0 ] + CmD[ r ][ 1 ] * cF[ 1 ] + CmD[ r ][ 2 ] * cF[ 2 ];

		// Gauss points of the face, mapped linearly from its three vertices
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.0 * mArea;

		Matrix<double> N;

		Lagrange2D_Ni_3sb( N, cX, cY );

		Vector< std::complex<double> > Wphase( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
			double L2 = cX[ gp ];
			double L3 = cY[ gp ];
			double L1 = 1.0 - L2 - L3;

			double X = L1 * mNodes[ 0 ]->X() + L2 * mNodes[ 1 ]->X() + L3 * mNodes[ 2 ]->X();
			double Y = L1 * mNodes[ 0 ]->Y() + L2 * mNodes[ 1 ]->Y() + L3 * mNodes[ 2 ]->Y();
			double Z = L1 * mNodes[ 0 ]->Z() + L2 * mNodes[ 1 ]->Z() + L3 * mNodes[ 2 ]->Z();

			Wphase[ gp ] = ( W[ gp ] * detJ ) * std::exp( jK * ( Kdir[ 0 ]*X + Kdir[ 1 ]*Y + Kdir[ 2 ]*Z ) );
		}

	    // Residual vector
		for( int i=0; i<mNumNodes; i++ )
        {
            std::complex<double> intg( 0.0, 0.0 );

			for( int gp=0; gp<numGaussPoints; gp++ ) intg += N[ i ][ gp ] * Wphase[ gp ];

            ResidualVector[ i               ] += jK_div_Mu * intg * CmDF[ 0 ];
            ResidualVector[ i +   mNumNodes ] += jK_div_Mu * intg * CmDF[ 1 ];
			ResidualVector[ i + 2*mNumNodes ] += jK_div_Mu * intg * CmDF[ 2 ];
        }  
    }

    void GenericRobin_3sb_FullWave::GetResidualVector_Single_Wave( Vector< std::complex<double> >& ResidualVector )
    {
	    // Module and phase of the primary field ( E/A )
        double modul_F = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ]; 
        double phase_F = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];

		if( modul_F == 0.0 ) return;

		// E/A polarization vector ( normalized )
        Vector<double> Fpol( 3 );

        Fpol[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 0 ];
        Fpol[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 1 ];
        Fpol[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 2 ];

        Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );

        // Complex E/A vector
        Vector< std::complex<double> > cF( 3 );

        cF[ 0 ] = Fpol[ 0 ]; 
        cF[ 1 ] = Fpol[ 1 ]; 
        cF[ 2 ] = Fpol[ 2 ];

        cF *= std::complex<double>( modul_F * cos( phase_F ), modul_F * sin( phase_F ) );

		if( mPotentialsOn ) 
		{
			cF *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		}

        // Wave propagation vector ( normalized )
        Vector<double> Kdir( 3 );

        Kdir[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 3 ];
        Kdir[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 4 ];
        Kdir[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 5 ];

        Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );

        // Source term evaluated at the Gauss points
        Add_PlaneWave_Source( cF, Kdir, ResidualVector );
    }

    //**********************************************************************************************************************
    //* - Residual vector for multiple plane waves
    //**********************************************************************************************************************
    void GenericRobin_3sb_FullWave::GetResidualVector_Multpl_Waves( Vector< Vector<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
	{
        // If plane waves set is empty then return
		if( FluxValue.size() == 0 ) return; 

		// Iterator over plane waves
		Vector< Vector<double> >::iterator pw_it;

        for( pw_it = FluxValue.begin(); pw_it != FluxValue.end(); ++pw_it )
		{
            // Module and phase of the primary field ( E/A )
            double modul_F = (*pw_it)[ 0 ]; 
            double phase_F = (*pw_it)[ 1 ];

			if( modul_F == 0.0 ) continue;

		    // E/A polarization vector ( normalized )
            Vector<double> Fpol( 3 );

            Fpol[ 0 ] = (*pw_it)[ 2 ];
            Fpol[ 1 ] = (*pw_it)[ 3 ];
            Fpol[ 2 ] = (*pw_it)[ 4 ];

            Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );

            // Complex E/A vector
            Vector< std::complex<double> > cF( 3 );

            cF[ 0 ] = Fpol[ 0 ]; 
            cF[ 1 ] = Fpol[ 1 ]; 
            cF[ 2 ] = Fpol[ 2 ];

            cF *= std::complex<double>( modul_F * cos( phase_F ), modul_F * sin( phase_F ) );

		    if( mPotentialsOn ) 
		    {
		    	cF *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		    }

            // Wave propagation vector ( normalized )
            Vector<double> Kdir( 3 );

            Kdir[ 0 ] = (*pw_it)[ 5 ];
            Kdir[ 1 ] = (*pw_it)[ 6 ];
            Kdir[ 2 ] = (*pw_it)[ 7 ];

            Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );

            // Source term evaluated at the Gauss points
            Add_PlaneWave_Source( cF, Kdir, ResidualVector );
		}
	}

    //**********************************************************************************************************************
    //* - Residual vector for Gaussian beams
    //**********************************************************************************************************************
    void GenericRobin_3sb_FullWave::GetResidualVector_Gauss_Beams( Vector< std::complex<double> >& ResidualVector )
    {
       // Gaussian beam waist radius 
        double GBR_pol = (*mProperties)( COMPLEX_IBC_2o )[ 0 ]; 
        double GBR_per = (*mProperties)( COMPLEX_IBC_2o )[ 1 ];

		if( ( GBR_pol <= 0.0 ) || ( GBR_per <= 0.0 ) ) return;

		// Field polarization vector ( normalized )
        Vector<double> Fpol( 3 );

        Fpol[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 0 ];
        Fpol[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 1 ];
        Fpol[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 2 ];

        Fpol *= ( 1.0 / sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] ) );

        // Wave propagation vector ( normalized )
        Vector<double> Kdir( 3 );

        Kdir[ 0 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 3 ];
        Kdir[ 1 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 4 ];
        Kdir[ 2 ] = (*mProperties)( SINUSOIDAL_SURFACE_CURRENT )[ 5 ];

        Kdir *= ( 1.0 / sqrt( Kdir[ 0 ]*Kdir[ 0 ] + Kdir[ 1 ]*Kdir[ 1 ] + Kdir[ 2 ]*Kdir[ 2 ] ) );

        // Remove any component of the polarization vector along the propagation direction (Gram-Schmidt),
        // so the waist-plane axes ( Fpol, Fper ) are orthogonal to Kdir even if the input vectors are not
        // exactly perpendicular. The negated comparison also catches NaN from zero-length input vectors.
        double FpolKdir = Fpol[ 0 ]*Kdir[ 0 ] + Fpol[ 1 ]*Kdir[ 1 ] + Fpol[ 2 ]*Kdir[ 2 ];

        Fpol[ 0 ] -= FpolKdir * Kdir[ 0 ];
        Fpol[ 1 ] -= FpolKdir * Kdir[ 1 ];
        Fpol[ 2 ] -= FpolKdir * Kdir[ 2 ];

        double normFpol = sqrt( Fpol[ 0 ]*Fpol[ 0 ] + Fpol[ 1 ]*Fpol[ 1 ] + Fpol[ 2 ]*Fpol[ 2 ] );

        if( !( normFpol > 1.0e-12 ) ) return;

        Fpol *= ( 1.0 / normFpol );

        // Fper = Kdir X Fpol ( normalized )
        Vector<double> Fper( 3 );

        Fper[ 0 ] = Kdir[ 1 ]*Fpol[ 2 ] - Kdir[ 2 ]*Fpol[ 1 ];
        Fper[ 1 ] = Kdir[ 2 ]*Fpol[ 0 ] - Kdir[ 0 ]*Fpol[ 2 ];
        Fper[ 2 ] = Kdir[ 0 ]*Fpol[ 1 ] - Kdir[ 1 ]*Fpol[ 0 ];

        Fper *= ( 1.0 / sqrt( Fper[ 0 ]*Fper[ 0 ] + Fper[ 1 ]*Fper[ 1 ] + Fper[ 2 ]*Fper[ 2 ] ) );

        // Field parallel to the polarization vector
        Vector< std::complex<double> > Fpar( 3 );

        Fpar[ 0 ] = Fpol[ 0 ]; 
        Fpar[ 1 ] = Fpol[ 1 ]; 
        Fpar[ 2 ] = Fpol[ 2 ];

        double mod_Fpar = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 3 ]; 
        double pha_Fpar = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 4 ];
        
        Fpar *= std::complex<double>( mod_Fpar * cos( pha_Fpar ), mod_Fpar * sin( pha_Fpar ) );

        // Field perpendicular to the polarization and wave vectors
        Vector< std::complex<double> > Fppw( 3 );

        Fppw[ 0 ] = Fper[ 0 ]; 
        Fppw[ 1 ] = Fper[ 1 ]; 
        Fppw[ 2 ] = Fper[ 2 ];

        double mod_Fppw = (*mProperties)( COMPLEX_IBC )[ 0 ]; 
        double pha_Fppw = (*mProperties)( COMPLEX_IBC )[ 1 ];
        
        Fppw *= std::complex<double>( mod_Fppw * cos( pha_Fppw ), mod_Fppw * sin( pha_Fppw ) );

        // Total field on waist central point
        Vector< std::complex<double> > Ftot( 3 );

        Ftot[ 0 ] = Fpar[ 0 ] + Fppw[ 0 ]; 
        Ftot[ 1 ] = Fpar[ 1 ] + Fppw[ 1 ]; 
        Ftot[ 2 ] = Fpar[ 2 ] + Fppw[ 2 ];

		if( mPotentialsOn ) 
		{
			Ftot *= ( 1.0 / std::complex<double>( 0.0, mFreq ) );
		}

		// Complex wave vector jK 
        std::complex<double> jK = Calculate_jK();

        // Rayleight range
        double ZR_pol = 0.5 * std::imag( jK ) * GBR_pol * GBR_pol;
        double ZR_per = 0.5 * std::imag( jK ) * GBR_per * GBR_per;

        // Gaussian beam waist center coordinates
        double Xo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 0 ];			
        double Yo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 1 ];
        double Zo = (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 2 ];

		// Source field evaluated at the Gauss points of the face ( previous versions evaluated it at the nodes and
		// interpolated it with the element basis, which under-resolves the beam phase on coarse port meshes ).
		// The face is mapped linearly from its three vertices, as in Calculate_NiNj_Matrix().
		Vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.0 * mArea;

		Matrix<double> N;

		Lagrange2D_Ni_3sb( N, cX, cY );

		Vector< Vector< std::complex<double> > > cFg( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
			W[ gp ] *= detJ;

			// Gauss point coordinates ( area coordinates L1, L2, L3 on the vertex nodes )
			double L2 = cX[ gp ];
			double L3 = cY[ gp ];
			double L1 = 1.0 - L2 - L3;

            Vector<double> Rn( 3 );

            Rn[ 0 ] = L1 * mNodes[ 0 ]->X() + L2 * mNodes[ 1 ]->X() + L3 * mNodes[ 2 ]->X() - Xo; 
            Rn[ 1 ] = L1 * mNodes[ 0 ]->Y() + L2 * mNodes[ 1 ]->Y() + L3 * mNodes[ 2 ]->Y() - Yo; 
            Rn[ 2 ] = L1 * mNodes[ 0 ]->Z() + L2 * mNodes[ 1 ]->Z() + L3 * mNodes[ 2 ]->Z() - Zo;

            // Distance from waist plane
            double Zd = Kdir[ 0 ]*Rn[ 0 ] + Kdir[ 1 ]*Rn[ 1 ] + Kdir[ 2 ]*Rn[ 2 ];

            // Accumulated phase 
            std::complex<double> AccPhase = jK * Zd;

            // Gouy phase shift ( 1/2 per transverse axis, elliptical beam )
            std::complex<double> GouyPhase( 0.0, 0.5 * ( std::atan( Zd / ZR_pol ) + std::atan( Zd / ZR_per ) ) );

            // Q parameters ( exp(-iwt) convention: q = z - i*zR )
            std::complex<double> Qpol( Zd, -ZR_pol );
            std::complex<double> Qper( Zd, -ZR_per );

            // Projection on waist plane axes
            double Xpol = Fpol[ 0 ]*Rn[ 0 ] + Fpol[ 1 ]*Rn[ 1 ] + Fpol[ 2 ]*Rn[ 2 ];
            double Xper = Fper[ 0 ]*Rn[ 0 ] + Fper[ 1 ]*Rn[ 1 ] + Fper[ 2 ]*Rn[ 2 ];

            // Transverse profile and wavefront curvature
            std::complex<double> QpPhase = jK * 0.5 * ( ( Xpol * Xpol / Qpol ) + ( Xper * Xper / Qper ) );

            // Field normalization: sqrt( w0_pol*w0_per / ( w_pol*w_per ) ) = sqrt( ZR_pol*ZR_per / ( |Qpol|*|Qper| ) )
            double Beta = std::sqrt( ZR_pol * ZR_per / ( std::abs( Qpol ) * std::abs( Qper ) ) );

            cFg[ gp ]  = Ftot;
            cFg[ gp ] *= Beta * std::exp( AccPhase - GouyPhase + QpPhase );
		}

		// Element exterior normal
		Vector<double> Next; 
		
		Calculate_Exterior_Normal( Next );

		// k*n
		double KdirNext = Kdir[ 0 ]*Next[ 0 ] + Kdir[ 1 ]*Next[ 1 ] + Kdir[ 2 ]*Next[ 2 ];

		// Complex value jK/Mu
        std::complex<double> jK_div_Mu = Calculate_jK_div_Mu();

		// Curl and div coefficient matrices: ( n x n x E - n x k x E ) and n ( n . E - k . E ).
		// Their combination C - D is applied in a single pass ( residual += jK/Mu * Int Ni ( C - D ) E ).
        double CmD[ 3 ][ 3 ];

        for( int r=0; r<3; r++ )
        {
            for( int c=0; c<3; c++ )
            {
                double Cm = Next[ r ]*Next[ c ] - Kdir[ r ]*Next[ c ] - ( ( r == c ) ? ( 1.0 - KdirNext ) : 0.0 );
                double Dm = Next[ r ]*Next[ c ] - Next[ r ]*Kdir[ c ];

                CmD[ r ][ c ] = Cm - Dm;
            }
        }

        // ( C - D ) * E on Gauss points, weighted
        Vector< Vector< std::complex<double> > > CEg( numGaussPoints );

		for( int gp=0; gp<numGaussPoints; gp++ )
		{
			CEg[ gp ].resize( 3 );

			for( int r=0; r<3; r++ )
			{
				CEg[ gp ][ r ] = W[ gp ] * ( CmD[ r ][ 0 ] * cFg[ gp ][ 0 ] + CmD[ r ][ 1 ] * cFg[ gp ][ 1 ] + CmD[ r ][ 2 ] * cFg[ gp ][ 2 ] );
			}
		}

	    // Residual vector
		for( int i=0; i<mNumNodes; i++ )
        {
            std::complex<double> intg_x( 0.0, 0.0 );
			std::complex<double> intg_y( 0.0, 0.0 );
			std::complex<double> intg_z( 0.0, 0.0 );

			for( int gp=0; gp<numGaussPoints; gp++ )
			{
				intg_x += N[ i ][ gp ] * CEg[ gp ][ 0 ];
				intg_y += N[ i ][ gp ] * CEg[ gp ][ 1 ];
				intg_z += N[ i ][ gp ] * CEg[ gp ][ 2 ];
			}

            ResidualVector[ i               ] += jK_div_Mu * intg_x;
            ResidualVector[ i +   mNumNodes ] += jK_div_Mu * intg_y;
			ResidualVector[ i + 2*mNumNodes ] += jK_div_Mu * intg_z;
        }  
    }

    //**********************************************************************************************************************
    //* - Residual vector for plane waves and Gaussian beams
    //**********************************************************************************************************************
    void GenericRobin_3sb_FullWave::GetResidualVector_Flux( Vector< Vector<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector )
    {
        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 1.0 )
        {
            GetResidualVector_Single_Wave( ResidualVector );
        }
		else if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 2.0 )
		{
            GetResidualVector_Multpl_Waves( FluxValue, ResidualVector );
		}
        else if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 3.0 )
        {
            GetResidualVector_Gauss_Beams( ResidualVector );
        }
    }

    //**********************************************************************************************************************
    //* - Stiffness matrix
    //**********************************************************************************************************************	
    void GenericRobin_3sb_FullWave::GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ) 
    {
		// Setting element stiffness matrix
        StiffMatrix.Resize( mNumDofs, mNumDofs, std::complex<double>( 0.0, 0.0 ) );
        
        // Robin coefficient
        std::complex<double> cRobinCoeff( 0.0, 0.0 );

        if( (*mProperties)( COMPLEX_NEUMANN_FLOW )[ 5 ] == 0.0 )
        {
		    cRobinCoeff -= std::complex<double>( (*mProperties)( COMPLEX_IBC )[ 0 ], (*mProperties)( COMPLEX_IBC )[ 1 ] );
        }
        else
        {
            cRobinCoeff -= Calculate_jK_div_Mu();
        }

        // NiNj surface integrals matrix
        Matrix<double> NiNj_Matrix;

		Calculate_NiNj_Matrix( NiNj_Matrix );

        // Stiffness matrix
        for( int i=0; i<mNumNodes; i++ )
        {
            for( int j=0; j<mNumNodes; j++ )
            {
                StiffMatrix[ i               ][ j               ] = cRobinCoeff * NiNj_Matrix[ i ][ j ];
                StiffMatrix[ i +   mNumNodes ][ j +   mNumNodes ] = cRobinCoeff * NiNj_Matrix[ i ][ j ];
                StiffMatrix[ i + 2*mNumNodes ][ j + 2*mNumNodes ] = cRobinCoeff * NiNj_Matrix[ i ][ j ];
            }
        }	
    }
} 