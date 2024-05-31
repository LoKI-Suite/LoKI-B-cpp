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

#ifndef LOKI_CPP_RANGE_H
#define LOKI_CPP_RANGE_H

#include "LoKI-B/json.h"
#include "LoKI-B/Exports.h"

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

    /** Create a new Range object from JSON object \a cnf. The caller must
     *  assume ownership of the pointer that is returned by this function.
     *
     *  \todo Describe the possible content of \a cnf once the dust has settled.
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

} // namespace loki

#endif // LOKI_CPP_RANGE_H
