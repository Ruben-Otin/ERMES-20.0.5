
#if !defined(ERMES_PHYSICAL_CONSTANTS_H)
#define ERMES_PHYSICAL_CONSTANTS_H

//**********************************************************************************************************************
//* - Physical constants (CODATA 2018, SI units).
//*
//*   Single place for constants that were previously hard-coded with different precisions in several files
//*   (e.g. eps0 = 8.8542e-12 and 8.8541878176e-12, me = 9.1095e-31). New code should use these values.
//**********************************************************************************************************************
namespace Kratos
{
    namespace PhysConst
    {
        const double Pi       = 3.14159265358979323846;

        const double Eps0     = 8.8541878128e-12;   // Vacuum electric permittivity [F/m]
        const double Mu0      = 1.25663706212e-6;   // Vacuum magnetic permeability [H/m] ( = 4*pi*1e-7 within 1e-9 )
        const double C0       = 299792458.0;        // Speed of light in vacuum [m/s]

        const double Qe       = 1.602176634e-19;    // Elementary charge [C]
        const double Me       = 9.1093837015e-31;   // Electron mass [kg]
        const double AMU      = 1.66053906660e-27;  // Atomic mass unit [kg]

        // Neutral atomic masses [u]. Ion masses are obtained subtracting the mass of the removed electrons.
        const double M_H1     = 1.00782503223;
        const double M_H2     = 2.01410177812;
        const double M_H3     = 3.01604928199;
        const double M_He3    = 3.01602932265;
        const double M_He4    = 4.00260325413;
        const double M_Be9    = 9.012183065;
        const double M_Ne20   = 19.9924401762;
        const double M_Ni60   = 59.93078588;
        const double M_W184   = 183.95093092;

        // Mass [kg] of an ion of neutral atomic mass AtomicMass_u [u] and charge state ChargeState
        inline double IonMass( double AtomicMass_u, double ChargeState )
        {
            return AtomicMass_u * AMU - ChargeState * Me;
        }
    }
}

#endif
