//
// Created by daan on 24-06-2019.
//

#ifndef LOKI_CPP_POWER_H
#define LOKI_CPP_POWER_H

#include "LoKI-B/Enumeration.h"

namespace loki
{

struct PowerTerm
{
    PowerTerm() : forward{0.0}, backward{0.0} { }
    double net() const { return forward + backward; }
    PowerTerm& operator+=(const PowerTerm& term)
    {
        forward += term.forward;
        backward += term.backward;
        return *this;
    }
    double forward;
    double backward;
};

/** \todo Document that this was changed during the development of the C++ code:
 *  while ionization, attachment only have forward 'inE' parts, we represent these
 *  terms with a full PowerTerm object. That hardly costs anything, and has the
 *  advantage of uniformity and nicer semantics: if somebody is interested in the
 *  net ionization energy losses, he can now say ionization.net(). Of course this
 *  will still be equal to ionization.forward in the absenace of reverse processes.
 */
struct GasPower
{
    GasPower() { }
    PowerTerm excitation;
    PowerTerm vibrational;
    PowerTerm rotational;
    PowerTerm ionization;
    PowerTerm attachment;
    GasPower& operator+=(const GasPower& term)
    {
        excitation += term.excitation;
        vibrational += term.vibrational;
        rotational += term.rotational;
        ionization += term.ionization;
        attachment += term.attachment;
        return *this;
    }
};

struct Power : GasPower
{
    double field{0.};
    double elasticNet{0.}, elasticGain{0.}, elasticLoss{0.};
    double carNet{0.}, carGain{0.}, carLoss{0.};

    double inelastic{0.};
    double superelastic{0.};
    double eDensGrowth{0.};
    double electronElectron{0.};
    double balance{0.};
    double relativeBalance{0.};
    double reference{0.};

    Power& operator+=(const GasPower &gasPower)
    {
	// invoke += on our base class
        GasPower::operator+=(gasPower);
	return *this;
    }
};

} // namespace loki

#endif // LOKI_CPP_POWER_H
