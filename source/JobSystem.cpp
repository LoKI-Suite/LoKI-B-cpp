/** \file
 *
 *  LoKI-B support for parameterized series of simulations.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2020 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
#include "LoKI-B/Parse.h"

#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <regex>

namespace loki
{

Range::~Range()
{
}

double Range::value(size_type ndx) const
{
    assert(ndx < size());
    return get_value(ndx);
}

/** A degenerate Range implementation that describes a single value.
 */
class RangeSingleValue : public Range
{
  public:
    RangeSingleValue(double value) : Range(1), m_value(value)
    {
        Log<Message>::Notify("Creating single-value range, value = ", value);
    }

  protected:
    virtual double get_value(size_type ndx) const override
    {
        return m_value;
    }

  private:
    const double m_value;
};

/** A range that represents a 'linspace', a series of equidistant
 *  values in a given domain.
 */
class RangeLinSpace : public Range
{
  public:
    RangeLinSpace(double start, double stop, size_type size) : Range(size), start(start), stop(stop)
    {
        Log<Message>::Notify("Creating linspace range"
                             ", start = ",
                             start, ", stop = ", stop, ", size = ", size);
        if (size < 2)
        {
            throw std::runtime_error("Range(linspace): at least two points required.");
        }
    }

  protected:
    virtual double get_value(size_type ndx) const override
    {
        return start + ndx * (stop - start) / (size() - 1);
    }

  private:
    const double start;
    const double stop;
};

/** A range that represents a 'logspace', the values are obtained as
 *  10^v, where v is any of a series of equidistant values in a given
 *  domain.
 */
class RangeLogSpace : public Range
{
  public:
    RangeLogSpace(double log_start, double log_stop, size_type size)
        : Range(size), log_start(log_start), log_stop(log_stop)
    {
        Log<Message>::Notify("Creating logspace range"
                             ", log_start = ",
                             log_start, ", log_stop = ", log_stop, ", size = ", size);
        if (size < 2)
        {
            throw std::runtime_error("Range(logspace): at least two points required.");
        }
    }

  protected:
    virtual double get_value(size_type ndx) const override
    {
        const double log_value = log_start + ndx * (log_stop - log_start) / (size() - 1);
        return std::pow(10., log_value);
    }

  private:
    const double log_start;
    const double log_stop;
};

/** A range that represents an array of predefined values.
 */
class RangeArray : public Range
{
  public:
    RangeArray(const json_type& cnf)
        : Range(cnf.size())
    {
        Log<Message>::Notify("Creating array range size = ", cnf.size());
        for (auto& v : cnf)
        {
            array.push_back(v);
        }
    }
  protected:
    virtual double get_value(size_type ndx) const override
    {
        return array[ndx];
    }

  private:
    std::vector<double> array;
};


static Range *createLinLogRange(const std::string &rangeString)
{
    try
    {
        static const std::regex r(
            R"(\s*((?:logspace)|(?:linspace))\(\s*(-?\d+\.?\d*)\s*,\s*(-?\d+\.?\d*)\s*,\s*(\d+\.?\d*))");
        std::smatch m;

        if (!std::regex_search(rangeString, m, r))
        {
            throw std::runtime_error("Parse error.");
        }

        std::stringstream ss;

        const std::string function = m.str(1);

        double start, stop;
        Range::size_type nvalues;

        ss << m[2];
        ss >> start;
        ss.clear();
        ss << m[3];
        ss >> stop;
        ss.clear();
        ss << m[4];
        ss >> nvalues;

        if (function == "linspace")
            return new RangeLinSpace(start, stop, nvalues);
        else if (function == "logspace")
            return new RangeLogSpace(start, stop, nvalues);
        else
            throw std::runtime_error("Unknown function '" + function + "'.");
    }
    catch (std::exception &exc)
    {
        throw std::runtime_error("Error creating range from string '" + rangeString + ":\n" + std::string(exc.what()));
    }
}

Range *Range::create(const std::string &str)
{
    double value;
    if (Parse::getValue(str, value))
    {
        // if it was a number (getValue succeeded), create a single-value-range
        return new RangeSingleValue{value};
    }
    else
    {
        // otherwise we assume it is a span (linear or logarithmic)
        return createLinLogRange(str);
    }
}

Range *Range::create(const json_type &cnf)
{
    if (cnf.type() == json_type::value_t::string)
    {
        return createLinLogRange(cnf);
    }
    else if (cnf.type() == json_type::value_t::array)
    {
        return new RangeArray(cnf);
    }
    else
    {
        return new RangeSingleValue(cnf.get<double>());
    }
}

/** A Parameter controls one of the parameters of a parametrized model.
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
    /** callback is the function that gets called by the JobManager when a
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
