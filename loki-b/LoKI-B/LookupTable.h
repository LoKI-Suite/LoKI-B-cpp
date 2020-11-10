/** \file
 *
 *  A lookup table class for LoKI-B
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
 *  \author Jan van Dijk and Daan Boer
 */

#ifndef LOKI_CPP_LOOKUPTABLE_H
#define LOKI_CPP_LOOKUPTABLE_H

#include "LoKI-B/LinearAlgebra.h"
#include "LoKI-B/json.h"
#include <iostream>
#include <cassert>

namespace loki {

/** A lookup table class that supports extrapolation and clipping.
 *  A lookup table is constructed with two vectors that represent samples of
 *  the abscissa and ordinate of a function. At least table points must
 *  be provided. After a lookup table has been constructed, approximate values
 *  of the function can be obtained for arbitrary points using one of the
 *  members interpolate, interpolate_or_clip or interpolate_or_set. See the
 *  documentation of these members for details.
 *
 *  Lookup table objects can also be created with the static member function
 *  create. Overloads exist for a JSON specification of the data and for the
 *  traditional LXCat format.
 *
 *  \todo Support interpolation for a series of ascending values at once,
 *        for performance reasons (getIndexL can be made faster then). Such
 *        functionality was present in the CrossSection.cpp implementation,
 *        and should be restored.
 *
 *  \author Jan van Dijk and Daan Boer
 *  \date November 2020
 */
class LookupTable
{
public:
    /// the type of the vector that is used for storing table values
    using Vector = loki::Vector;
    /// the type of an index of a Vector object
    using Index = Vector::Index;

    /** A lookup table can be constructed by passing vectors \a x and \a y
     *  of equal lengths to the constructor. These vectors are copied into
     *  private data members for later reference. The elements of \a x must
     *  be in ascending order. The vector \a x represents the values of the
     *  argument for the samples, the vector \a y the corresponding function
     *  values. At least two table points must be specified.
     *
     *  If the vectors have different lengths, the length is smaller than two
     *  of the values of \a x are not in ascending order, a std:::runtime_error
     *  is thrown.
     */
    LookupTable(const Vector& x, const Vector& y);
    /** Create and return a lookup table that is created from the JSON object \a data.
     *  \a data must contain a field "data" that is an array of double pairs, each
     *  representing a table point.
     */
    static LookupTable create(const json_type& data);
    /* Accepts a reference to an input file stream of an LXCat file. This stream should be
     * at a position just after reading a collision description from the LXCat file, since
     * this function searches for the line containing solely dashes, indicating that a
     * cross section follows. The raw cross section is then stored as a vector of pairs of
     * doubles, which the user passes by reference.
     */
    static LookupTable create(std::istream& is);

    /// returns a reference to the tabulated values of the argument
    const Vector& x() const { return m_x; }
    /// returns a reference to the tabulated values of the function
    const Vector& y() const { return m_y; }
    /// returns the number of tabulated points
    Index size() const { return m_x.size(); }
    /// returns the lowest tabulated argument value
    double xMin() const { return m_x[0]; }
    /// returns the highest tabulated argument value
    double xMax() const { return m_x[size()-1]; }
    /** Get the index of the last sample of the abscissa
     *  that is smaller than or equal to \a x. When \a x
     *  exceeds the highest sampled value of the abscissa,
     *  the value size()-2 is returned.
     *  This guarantees that, if the return value is denoted as i,
     *  the samples i and i+1 are always available for interpolating
     *  function values for points inbetween x()[i] and x()[i+1] or,
     *  when i=0 or i==size()-2, for extrapolating function values
     *  for points that are outside the range of the table.
     */
    Index getIndexL(double x) const
    {
        assert(size()>1);
        Index base=1;
        const Index ndxMax = size()-1;
        while(m_x[base]<=x && base!=ndxMax)
        {
            ++base;
        }
        return base-1;
    }
    /** Lookup the indices that contain \a x, use linear interpolation
     *  to interpolate the function value in \a x and return the result.
     *  Note: when \a x is outside of the table range, linear extrapolation
     *  is done, using the lowest (x<xMin()) or highest (x>xMax()) pair of
     *  values.
     */
    double interpolate(double x) const
    {
        const Index il = getIndexL(x);
        const double x_rel = (x-m_x[il])/(m_x[il+1]-m_x[il]);
        return m_y[il]*(1-x_rel)+m_y[il+1]*x_rel;
    }
    /** This function behaves the same way as interpolate, except that the
     *  value is clipped to the lowest tabulated value when x<xMin() and
     *  \a clip_low is true, and to the highest tabulated value when x>xMin()
     *  and \a clip_high is true.
     */
    double interpolate_or_clip(double x, bool clip_low, bool clip_high) const
    {
        if (x<xMin() && clip_low)
        {
            return m_y[0];
        }
        else if (x>xMax() && clip_high)
        {
            return m_y[m_y.size()-1];
        }
        else
        {
            return interpolate(x);
        }
    }
    /** When \a x < xMin(), this function returns \a y_below; when
     *  \a x > xMin() this function returns y_above. For other values
     *  of \a x, which are in the range of the table, the result of
     *  calling interpolate(x) is returned.
     */
    double interpolate_or_set(double x, double y_below, double y_above) const
    {
        return x<xMin() ? y_below : x>xMax() ? y_above : interpolate(x);
    }
private:
    /// the tabulated values of the argument
    const Vector m_x;
    /// the tabulated values of the function
    const Vector m_y;
};

} // namespace loki

#endif // LOKI_CPP_LOOKUPTABLE_H
