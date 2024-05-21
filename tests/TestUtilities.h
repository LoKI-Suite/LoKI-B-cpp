/** \file
 *
 *  A small collection of unit test functions and macros for LoKI-B.
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
 *  This code is based on an implementation of the same idea in PLASIMO
 *  (see https://plasimo.nl) by the author.
 *
 *  \author Jan van Dijk
 *  \date   10 November 2013
 */

#ifndef LOKI_CPP_TEST_TEST_UTILITIES_H
#define LOKI_CPP_TEST_TEST_UTILITIES_H

unsigned ntests=0;
unsigned nerrors=0;

#define test_expr(E) \
do { \
    ++ntests; \
    std::clog << "Testing expression '" << #E << "': "; \
    if (E) \
        std::clog << "OK"; \
    else \
    { \
        std::clog << "***FAILED***"; \
        ++nerrors; \
    } \
    std::clog << std::endl; \
} while (0)

#define test_report \
    std::clog << "Executed " << ntests << " tests, " \
        << nerrors << " failed." << std::endl;

#endif // LOKI_CPP_TEST_TEST_UTILITIES_H
