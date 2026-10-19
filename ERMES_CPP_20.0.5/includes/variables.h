
#if !defined(KRATOS_VARIABLES)
#define KRATOS_VARIABLES

#include <complex>

#include "../includes/variable.h"
#include "../includes/vector.h"
#include "../includes/matrix.h"

namespace Kratos
{
	// Parameters
	const Variable<double> Sg  ( "Sg", 100 );
	const Variable<double> TIME( "T" , 200 );

	// Current sources
	const Variable< Vector<double> > SINUSOIDAL_VOLUME_CURRENT ( "J", 1000 );
	const Variable< Vector<double> > SINUSOIDAL_SURFACE_CURRENT( "S", 1002 );
	const Variable< Vector<double> > SINUSOIDAL_POINT_CURRENT  ( "I", 1004 );

	// Boundary condition coefficients
	const Variable<double> CONVECTION_COEFFICIENT( "Coef", 2000 );
	const Variable<double> NEUMANN_FLOW          ( "NF"  , 2002 );
														              		      
	const Variable< Vector<double> > COMPLEX_NEUMANN_FLOW( "cNF", 2004 );
	const Variable        <double>   REAL_NEUMANN_FLOW   ( "rNF", 2006 );
	const Variable        <double>   IMAG_NEUMANN_FLOW   ( "iNF", 2008 );
																				   
	const Variable< Vector<double> > COMPLEX_IBC_2o( "cIBC2o", 2010 );
	const Variable< Vector<double> > COMPLEX_IBC   ( "cIBC"  , 2012 );

	const Variable<double> REAL_IMPEDANCE_BOUNDARY_CONDITION( "rIBC", 2014 );
	const Variable<double> IMAG_IMPEDANCE_BOUNDARY_CONDITION( "iIBC", 2016 );

	// Frequency
	const Variable<double> FREQUENCY( "w"    , 3000 );
	const Variable<double> PHASE    ( "phase", 3002 );

	// Material properties
	const Variable<double> REAL_ELECTRIC_CONDUCTIVITY( "sigma_real"  , 4000 );
	const Variable<double> IMAG_ELECTRIC_CONDUCTIVITY( "sigma_imag"  , 4002 );

	const Variable<double> REAL_ELECTRIC_PERMITTIVITY( "epsilon_real", 4004 );
	const Variable<double> IMAG_ELECTRIC_PERMITTIVITY( "epsilon_imag", 4006 );
																		
	const Variable<double> REAL_MAGNETIC_PERMEABILITY( "mu_real"     , 4010 );
	const Variable<double> IMAG_MAGNETIC_PERMEABILITY( "mu_imag"     , 4012 );
																		
	// Complex electric field
    const Variable< std::complex<double> > cEx( "cEx", 5000 );											    
    const Variable< std::complex<double> > cEy( "cEy", 5001 );								    
    const Variable< std::complex<double> > cEz( "cEz", 5002 );

	// Electric field
	const Variable< Vector<double> > E     ( "E(t)"  , 5010 );
	const Variable< Vector<double> > REAL_E( "rE"    , 5012 );
	const Variable< Vector<double> > IMAG_E( "iE"    , 5014 );
	const Variable        <double>   MOD_E ( "mod(E)", 5016 );

	// Magnetic field
	const Variable< Vector<double> > H     ( "H(t)"  , 5020 );
	const Variable< Vector<double> > REAL_H( "rH"    , 5022 );
	const Variable< Vector<double> > IMAG_H( "iH"    , 5024 );
	const Variable        <double>   MOD_H ( "mod(H)", 5026 );

    // Magnetic flux density
	const Variable< Vector<double> > B     ( "B(t)"  , 5030 );
	const Variable< Vector<double> > REAL_B( "rB"    , 5032 );
	const Variable< Vector<double> > IMAG_B( "iB"    , 5034 );
	const Variable        <double>   MOD_B ( "mod(B)", 5036 );

    // Current density	
	const Variable< Vector<double> > J     ( "J(t)"  , 5040 );
	const Variable< Vector<double> > REAL_J( "rJ"    , 5042 );
	const Variable< Vector<double> > IMAG_J( "iJ"    , 5044 );
	const Variable        <double>   MOD_J ( "mod(J)", 5046 );

	// Complex magnetic vector potential
    const Variable< std::complex<double> > cAx( "cAx", 5050 );											    
    const Variable< std::complex<double> > cAy( "cAy", 5051 );								    
    const Variable< std::complex<double> > cAz( "cAz", 5052 );

	// Magnetic vector potential
	const Variable< Vector<double> > A     ( "A(t)"  , 5060 );
	const Variable< Vector<double> > REAL_A( "rA"    , 5062 );
	const Variable< Vector<double> > IMAG_A( "iA"    , 5064 );
	const Variable        <double>   MOD_A ( "mod(A)", 5066 );

	// Complex scalar potential
	const Variable< std::complex<double> > cVs ( "cVs" , 5070 );

	// Scalar potential
	const Variable<double> V     ( "V(t)"  , 5080 );
	const Variable<double> REAL_V( "rV"    , 5082 );
	const Variable<double> IMAG_V( "iV"    , 5084 );
	const Variable<double> MOD_V ( "mod(V)", 5086 );

	// Derived magnitudes
	const Variable        <double>   JOULE_HEATING  ( "Joule_heating"  , 6000 );
    const Variable< Vector<double> > POYNTING_VECTOR( "Poynting_vector", 6002 );
	const Variable< Vector<double> > LORENTZ_FORCE  ( "Lorentz_force"  , 6004 );

	// Boundary and contact normals
	const Variable< Vector<double> > BOUNDARY_NORMALS( "Boundary_normals", 6010 );
	const Variable< Vector<double> > CONTACT_NORMALS ( "Contact_normals" , 6012 );

	// Cold plasma properties
	const Variable<double> ELECTRON_DENSITY     ( "Electron_density"     , 7000 );
	const Variable<double> PLASMA_PERMITTIVITY_S( "Plasma_permittivity_S", 7001 );
	const Variable<double> PLASMA_PERMITTIVITY_D( "Plasma_permittivity_D", 7002 );
	const Variable<double> PLASMA_PERMITTIVITY_P( "Plasma_permittivity_P", 7003 );
	const Variable<double> PLASMA_PERMITTIVITY_R( "Plasma_permittivity_R", 7004 );
	const Variable<double> PLASMA_PERMITTIVITY_L( "Plasma_permittivity_L", 7005 );
    
	// Cold plasma fields
    const Variable< Vector<double> > B_EXT              ( "Bext"               , 7010 );
    const Variable< Vector<double> > E_PARALLEL_T       ( "Eparallel(t)"       , 7011 );
    const Variable        <double>   MOD_E_PARALLEL     ( "mod(Eparallel)"     , 7012 );
    const Variable< Vector<double> > E_PERPENDICULAR_T  ( "Eperpendicular(t)"  , 7013 );
    const Variable        <double>   MOD_E_PERPENDICULAR( "mod(Eperpendicular)", 7014 );
    
	// Electrostatic potential
    const Variable<double> VOLTAGE( "Voltage", 8000 );

	// Electrostatic magnitudes
	const Variable< Vector<double> > ELECTROSTATIC_ELECTRIC_FIELD ( "Electric_field" , 8010 );
    const Variable< Vector<double> > ELECTROSTATIC_CURRENT_DENSITY( "Current_density", 8012 );
    const Variable        <double>   ELECTROSTATIC_JOULE_HEATING  ( "Joule_heating"  , 8014 );
}  

#endif 

