/** \file
 *
 *  Implementation of a class that manages the working conditions of a simulation run.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2025 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
 *  M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
 *  L. L. Alves
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  \author Daan Boer and Jan van Dijk
 *  \date   2 May 2019
 */


#include "LoKI-B/WorkingConditions.h"
#include "LoKI-B/Constant.h"
#include "LoKI-B/Log.h"
#include <limits>

namespace loki
{

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
