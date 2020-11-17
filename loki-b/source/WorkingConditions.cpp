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
    : m_gasPressure(cnf.at("gasPressure").get<double>()),
      m_gasTemperature(cnf.at("gasTemperature").get<double>()),
      m_gasDensity(m_gasPressure / (Constant::boltzmann * m_gasTemperature)),
      m_electronDensity(cnf.at("electronDensity").get<double>()),
      //chamberLength(cnf.at("chamberLength").get<double>()),
      //chamberRadius(cnf.at("chamberRadius").get<double>()),
      m_excitationFrequency(cnf.at("excitationFrequency").get<double>())
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
