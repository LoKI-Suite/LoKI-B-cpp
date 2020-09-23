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

const double boltzmann = 1.38064852e-23;           // Boltzmann constant in J/K
const double electronCharge = 1.6021766208e-19;    // Electron charge in C
const double electronMass = 9.10938356e-31;        // Electron mass in Kg
const double unifiedAtomicMass = 1.660539040e-27;  // Unified Atomic Mass unit (UAM) in kg
const double bohrRadius = 5.2917721067e-11;        // Bohr radius in m
const double vacuumPermittivity = 8.854187817e-12; // Vacuum permittivity in F/m
const double plank = 6.626070040e-34;              // Plank constant in J s
const double speedOfLight = 299792458;             // Speed of light in vacuum in m/s
const double atmosphereInPa = 101325;              // Standard atmosphere in Pa
const double atmosphereInTorr = 760;               // Standard atmosphere in Torr (not obtained from NIST database)

const double pi = M_PI;
const double e = M_E;

const double kBeV = boltzmann / electronCharge;
const double plankReduced = plank / (2. * pi);
const double plankReducedInEv = plank / (2. * pi * electronCharge);
} // namespace Constant
} // namespace loki

#endif // LOKI_CPP_CONSTANT_H
