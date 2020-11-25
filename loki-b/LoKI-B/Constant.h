//
// Created by daan on 2-5-19.
//

#ifndef LOKI_CPP_CONSTANT_H
#define LOKI_CPP_CONSTANT_H

#include <cmath>

namespace loki
{
namespace Constant
{
// List of fundamental constants and units as obtained from the NIST database
// (https://physics.nist.gov/cuu/Constants/index.html) on Jule 23rd 2018.

constexpr const double boltzmann = 1.38064852e-23;           // Boltzmann constant in J/K
constexpr const double electronCharge = 1.6021766208e-19;    // Electron charge in C
constexpr const double electronMass = 9.10938356e-31;        // Electron mass in Kg
constexpr const double unifiedAtomicMass = 1.660539040e-27;  // Unified Atomic Mass unit (UAM) in kg
constexpr const double bohrRadius = 5.2917721067e-11;        // Bohr radius in m
constexpr const double vacuumPermittivity = 8.854187817e-12; // Vacuum permittivity in F/m
constexpr const double plank = 6.626070040e-34;              // Plank constant in J s
constexpr const double speedOfLight = 299792458;             // Speed of light in vacuum in m/s
constexpr const double atmosphereInPa = 101325;              // Standard atmosphere in Pa
constexpr const double atmosphereInTorr = 760;               // Standard atmosphere in Torr (not obtained from NIST database)

/// \todo Always define the literal pi value ourselves?
#ifndef M_PI
constexpr const double pi = 3.1415926535897932384626433833;
#else
constexpr const double pi = M_PI;
#endif

constexpr const double kBeV = boltzmann / electronCharge;
constexpr const double plankReduced = plank / (2. * pi);
constexpr const double plankReducedInEv = plank / (2. * pi * electronCharge);

} // namespace Constant

namespace SI {

/** The Townsend is a unit for the reduced field E/N and is defined
 *  to be equal to 1e-17 V*cm^2, or in SI units 1e-21 V*m^2.
 */
constexpr const double Townsend = 1e-21;
/** The constant gamma=sqrt(2e/m_e) appears in many places in the code.
 *  The name gamma for this constant appears in \cite Hagelaar2005 just
 *  below equation (6). Here you find the numerical value of gamma in SI
 *  units, sqrt(C/kg).
 */
constexpr const double gamma = std::sqrt(2*Constant::electronCharge/Constant::electronMass);

}

} // namespace loki

#endif // LOKI_CPP_CONSTANT_H
