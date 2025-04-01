/** \file
 *
 *  Declaration of a function for visualizing the sparsity pattern of a matrix
 *  as an XPM picture.
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
 *  \author Jan van Dijk and Daan Boer
 *  \date November 2020
 */

#ifndef LOKI_CPP_MATRIX2PICTURE_H
#define LOKI_CPP_MATRIX2PICTURE_H

#include "LoKI-B/LinearAlgebra.h"
#include <fstream>
#include <iostream>

namespace loki
{

/** \todo See if there something else that can be created
 *        this easily and more widely supported, instead
 *        of xpm.
 */
inline void writeXPM(const Matrix &m, std::ostream &os)
{
    os << "/* XPM */" << std::endl;
    os << "char* const name[] = {" << std::endl;
    // height, width, number of colors, number of chars per pixel
    os << '"' << m.rows() << ' ' << m.cols() << " 3 1\"" << std::endl;
    // declare color code characters w (white), b (blue) anr r (red)
    os << "\"w c #ffffff\"" << std::endl;
    os << "\"b c #0000ff\"" << std::endl;
    os << "\"r c #ff0000\"" << std::endl;
    for (Matrix::Index r = 0; r != m.rows(); ++r)
    {
        os << '"';
        for (Matrix::Index c = 0; c != m.cols(); ++c)
        {
            // <0: red
            // =0: white
            // >0: blue
            const double v = m(r, c);
            os << (v > 0 ? 'b' : v < 0 ? 'r' : 'w');
        }
        os << '"';
        if (r != m.rows() - 1)
            os << ',';
        os << '\n';
    }
    os << "};\n";
}

inline void writeXPM(const Matrix &m, std::string f)
{
    std::ofstream os(f);
    if (!os)
    {
        throw std::runtime_error("Could not open file '" + f + "' for writing.");
    }
    writeXPM(m, os);
}

} // namespace loki

#endif // LOKI_CPP_MATRIX2PICTURE_H
