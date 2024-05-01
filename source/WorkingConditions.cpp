//
// Created by daan on 2-5-19.
//

#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Log.h"
#include <limits>

namespace loki
{

// TODO: Comment on WorkingConditions().

WorkingConditions::WorkingConditions(const WorkingConditionsSetup &setup)
    : m_gasPressure(setup.gasPressure),
      m_gasTemperature(setup.gasTemperature),
      m_gasDensity(m_gasPressure / (Constant::boltzmann * m_gasTemperature)),
      m_electronDensity(setup.electronDensity),
      //chamberLength(setup.chamberLength),
      //chamberRadius(setup.chamberRadius),
      m_excitationFrequency(setup.excitationFrequency)
{
    /* set the reducedField and electronTemperature to dummy
     * values. These are set by the JobManager when prepareFirstJobs
     * or nextJob() is called.
     */
    m_reducedField = std::numeric_limits<double>::quiet_NaN();
    // see WorkingConditions.h: electronTemperature does not seem to be used.
    m_electronTemperature = std::numeric_limits<double>::quiet_NaN();

    /// \todo More parameters may need to be made available
    m_argumentMap.emplace("gasTemperature", &m_gasTemperature);
}

WorkingConditions::WorkingConditions(const json_type &cnf)
    : m_gasPressure(cnf.at("gasPressure").at("value").get<double>()),
      m_gasTemperature(cnf.at("gasTemperature").at("value").get<double>()),
      m_gasDensity(m_gasPressure / (Constant::boltzmann * m_gasTemperature)),
      m_electronDensity(cnf.at("electronDensity").at("value").get<double>()),
      //chamberLength(cnf.at("chamberLength").get<double>()),
      //chamberRadius(cnf.at("chamberRadius").get<double>()),
      m_excitationFrequency(cnf.at("excitationFrequency").at("value").get<double>())
{
    /* set the reducedField and electronTemperature to dummy
     * values. These are set by the JobManager when prepareFirstJobs
     * or nextJob() is called.
     */
    m_reducedField = std::numeric_limits<double>::quiet_NaN();
    // see WorkingConditions.h: electronTemperature does not seem to be used.
    m_electronTemperature = std::numeric_limits<double>::quiet_NaN();

    /// \todo More parameters may need to be made available
    m_argumentMap.emplace("gasTemperature", &m_gasTemperature);
}

const double* WorkingConditions::addParameter(const std::string& name, const double* vptr)
{
    if (findParameter(name))
    {
        throw std::runtime_error("Redeclaration of parameter '" + name + "'.");
    }
    if (!vptr)
    {
        throw std::runtime_error("Parameter '" + name + "' must be non-null.");
    }
    m_argumentMap.emplace(name, vptr);
    Log<Message>::Notify("Working conditions: added parameter + '" + name + "'.");
    return vptr;
}

const double* WorkingConditions::findParameter(const std::string& name) const
{
    const ArgumentMap::const_iterator i = m_argumentMap.find(name);
    return i==m_argumentMap.end() ? nullptr : i->second;
}

const double& WorkingConditions::getParameter(const std::string& name) const
{
    const double* vptr = findParameter(name);
    if (!vptr)
    {
        throw std::runtime_error("Parameter '" + name + "' is not available.");
    }
    return *vptr;
}

void WorkingConditions::updateReducedField(double value)
{
    m_reducedField = value;
    updatedReducedField.emit();
}

void WorkingConditions::updateElectronTemperature(double value)
{
    m_electronTemperature = value;
    updatedElectronTemperature.emit();
}

} // namespace loki
