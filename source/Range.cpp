/** \file
 *
 *  LoKI-B implementation of a value range.
 *
 *  LoKI-B solves a time and space independent form of the two-term
 *  electron Boltzmann equation (EBE), for non-magnetised non-equilibrium
 *  low-temperature plasmas excited by DC/HF electric fields from
 *  different gases or gas mixtures.
 *  Copyright (C) 2018-2024 A. Tejero-del-Caz, V. Guerra, D. Goncalves,
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
 *  \date   May 2019 (initial version)
 */

#include "LoKI-B/Range.h"
#include "LoKI-B/Log.h"

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
    RangeLinSpace(double start, double stop, size_type size) : Range(size), m_start(start), m_stop(stop)
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
        return m_start + ndx * (m_stop - m_start) / (size() - 1);
    }

  private:
    const double m_start;
    const double m_stop;
};

/** A range that represents a 'logspace', the values are obtained as
 *  10^v, where v is any of a series of equidistant values in a given
 *  domain.
 */
class RangeLogSpace : public Range
{
  public:
    RangeLogSpace(double log_start, double log_stop, size_type size)
        : Range(size), m_log_start(log_start), m_log_stop(log_stop)
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
        const double log_value = m_log_start + ndx * (m_log_stop - m_log_start) / (size() - 1);
        return std::pow(10., log_value);
    }

  private:
    const double m_log_start;
    const double m_log_stop;
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

Range *Range::create(const json_type &cnf)
{
    /** \todo handle (or at least check) the unit. Add the target unit as
     *  argument of this function to allow that.
     */
    if (cnf.contains("value"))
    {
        return new RangeSingleValue(cnf.at("value"));
    }
    else if (cnf.contains("range"))
    {
        const std::string function = cnf.at("range");
        if (function == "linspace")
        {
            return new RangeLinSpace(cnf.at("min"), cnf.at("max"), cnf.at("nvalues"));
        }
        else if (function == "logspace")
        {
            return new RangeLogSpace(cnf.at("min"), cnf.at("max"), cnf.at("nvalues"));
        }
        else
        {
            throw std::runtime_error("Unknown range type '" + function +
                "'. Expected 'linspace' or 'logspace'.");
        }
    }
    else if (cnf.type() == json_type::value_t::array)
    {
        /** \todo Define semantics for an array-range. Where to state the unit?
         *  This mode does not exist in MATLAB-LoKI-B, IIRC.
         */
        return new RangeArray(cnf);
    }
    else
    {
        throw std::runtime_error("Cannot create a range from section '" + cnf.dump(2) + "'.");
    }
}

} // namespace loki
