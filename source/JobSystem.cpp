/** \file
 *
 *  LoKI-B support for parameterized series of simulations.
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
 *  \author Daan Boer and Jan van Dijk (C++ version)
 */

#include "LoKI-B/JobSystem.h"
#include "LoKI-B/Log.h"

#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace loki
{

/** A Parameter controls one of the parameters of a parameterized model.
 *  It manages a Range object that describes the value(s) of the parameter
 *  for which the model must be run, In addition it manages a callback function,
 *  which is called by the JobMaanager when this parameter value changes:
 *  it must prepare the model to do a run with the new set of values.
 */
class JobManager::Parameter
{
  public:
    /** The callback function must accept a double and returns a void.
     */
    Parameter(const std::string &name, const callback_type callback, const Range *range)
        : m_callback(callback), m_name(name), m_range(range), m_active_ndx(0)
    {
    }
    const std::string &name() const
    {
        return m_name;
    }
    /** m_callback is the function that gets called by the JobManager when a
     *  new value in the range is activated.
     */
    const callback_type m_callback;
    Range::size_type size() const
    {
        return m_range->size();
    }
    double active_value() const
    {
        return m_range->value(m_active_ndx);
    }
    void reset()
    {
        m_active_ndx = 0;
    }
    bool advance()
    {
        return (m_range->size() - 1) > m_active_ndx++;
    }

  private:
    std::string m_name;
    const std::unique_ptr<const Range> m_range;
    Range::size_type m_active_ndx;
};

JobManager::JobManager() : m_jobIndex(0)
{
}

JobManager::JobManager(JobManager&&) = default;

JobManager::~JobManager()
{
}

void JobManager::addParameter(const std::string &name, const callback_type callback, const Range *range)
{
    Log<Message>::Notify("Adding simulation parameter '", name, "'");
    m_parameters.emplace_back(new Parameter{name, callback, range});
}

void JobManager::prepareFirstJob()
{
    for (const auto &p : m_parameters)
    {
        (p->m_callback)(p->active_value());
    }
}

bool JobManager::prepareNextJob()
{
    Parameter &p = *m_parameters[m_jobIndex];

    if (p.advance())
    {
        (p.m_callback)(p.active_value());

        if (m_jobIndex != m_parameters.size() - 1)
        {
            ++m_jobIndex;
        }

        return true;
    }
    else
    {
        if (m_jobIndex == 0)
            return false;

        p.reset();
        --m_jobIndex;

        return prepareNextJob();
    }
}

JobManager::size_type JobManager::njobs() const
{
    if (dimension() == 0)
    {
        return 0;
    }
    size_type n = 1;
    for (const auto &p : m_parameters)
    {
        n *= p->size();
    }
    return n;
}

std::string JobManager::getCurrentJobFolder() const
{
    std::stringstream ss;

    for (const auto &p : m_parameters)
    {
        if (p != m_parameters.front())
        {
            ss << '_';
        }
        ss << p->name() << "_" << p->active_value();
    }
    return ss.str();
}

} // namespace loki
