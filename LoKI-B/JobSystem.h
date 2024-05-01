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

#ifndef LOKI_CPP_JOBSYSTEM_H
#define LOKI_CPP_JOBSYSTEM_H

#include "LoKI-B/json.h"
#include "LoKI-B/Exports.h"
#include <cstdlib>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace loki
{

/** A Range class gives access to a sequence of values and has two public
 *  members: size() returns the number of values, value(size_type)
 *  returns the value for a given index. The number of values is passed
 *  to the constructor and kept in a private member, value() delegates
 *  the retrieval of the value to the protected virtual member get_value.
 *  That member must be implemented by derived classes to do the
 *  actual calculation and retrun the value.
 */
class lokib_export Range
{
  public:
    using size_type = std::size_t;

    /** \todo this supports the legacy parser. Remove when the time comes.
     *
     *  Create a new Range object from string \a str. The caller must
     *  assume ownership of the pointer that is returned by this function.
     *
     *  The argument can describe a single value, or a linear of logarithmic
     *  value range. Some sample input and the values it will produce are:
     *  \verbatim
          "42.0"             # 42.0
          "linspan(0,20,3)"  # 0, 10, 20
          "logspan(2,4,3)"   # 1e2, 1e3, 1e4 \endverbatim
     *
     *  \sa RangeSingleValue
     *  \sa RangeLinSpace
     *  \sa RangeLogSpace
     */
    static Range *create(const std::string &str);
    /** Create a new Range object from JSON object \a cnf. The caller must
     *  assume ownership of the pointer that is returned by this function.
     *
     *  \todo Describe the object once we agree on the structure.
     *
     *  \sa RangeSingleValue
     *  \sa RangeLinSpace
     *  \sa RangeLogSpace
     */
    static Range *create(const json_type &cnf);
    virtual ~Range();

    /// returns the number of values of this Range.
    size_type size() const
    {
        return m_size;
    }
    /// Returns the \a ndx'th value. Argument \a ndx must be smaller than size().
    double value(size_type ndx) const;

  protected:
    /// This constructor records \a size, the number of values in the Range
    Range(size_type size) : m_size(size)
    {
    }
    /// Must be overridden to return the value for \a ndx < size().
    virtual double get_value(size_type ndx) const = 0;

  private:
    /// The number of values of this Range
    const size_type m_size;
};

/** Class JobManager manages a collection of parameter definitions and
 *  makes it easy to run a simulation for each combination of parameter
 *  values.
 *
 *  A parameter can be declared by calling member addParameter, passing the
 *  name of the parameter, a pointer to a callback function and a pointer
 *  to a range object. The callback function must accept a double-valued
 *  argument and return void. It defines the action that must be undertaken
 *  by the model when the value of this parameter changes, and takes the
 *  new parameter value as argument. The JobManager assumes ownership of
 *  the Range pointer and will delete it when its destructor is invoked.
 *
 *  After the simulation has been set up, the user should call member
 *  prepareFirstJob() to activate the initial combination of parameter
 *  values. Subsequent parameter combinations can be activated by calling
 *  prepareNextJob(). That member returns true if a new job could be
 *  activated, false if the set of parameter combinations is exhausted.
 *  Member getCurrentJobFolder() returns a case identifier
 *  that contains the names and active values of the parameters, separated
 *  by underscore characters.
 */
class lokib_export JobManager
{
  public:
    using size_type = std::size_t;
    using callback_type = std::function<void(double)>;

    JobManager();
    JobManager(const JobManager&) = delete;
    JobManager(JobManager&&);
    JobManager& operator=(const JobManager&) = delete;
    JobManager& operator=(const JobManager&&) = delete;
    ~JobManager();

    void addParameter(const std::string &name, const callback_type callback, const Range *range);
    void prepareFirstJob();
    bool prepareNextJob();
    std::string getCurrentJobFolder() const;
    size_type dimension() const
    {
        return m_parameters.size();
    }
    size_type njobs() const;

  private:
    class Parameter;
    std::vector<std::unique_ptr<Parameter>> m_parameters;
    size_type m_jobIndex;
};

} // namespace loki

#endif // LOKI_CPP_JOBSYSTEM_H
