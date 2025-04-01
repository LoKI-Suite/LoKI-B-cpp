/** \file
 *
 *  Test the lokib::Event<> template.
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
 *  \author Jan van Dijk
 *  \date   11 September 2022
 */

#include "LoKI-B/Event.h"
#include "tests/TestUtilities.h"
#include <iostream>
#include <sstream>

void hello(std::ostream& os)
{
    os << "Hello";
}

struct Writer
{
    void write_separator(std::ostream& os) const
    {
        os << ", ";
    }
    void write_world(std::ostream& os)
    {
        os << "world!";
    }
};

void test1()
{
    /* 1. Create an event type that accepts callbacks that are compatible with
     *    the signature 'void f(std::ostream&)'.
     */
    using OutputEvent = loki::Event<std::ostream&>;
    OutputEvent evt;

    /* 2. Add callbacks to the event. we add free function 'hello', and
     *    member functions 'write_separator' and 'write_world' of the
     *    'Writer' class, both bound to Writer instance w.
     */
    evt.addListener(hello);

    Writer w;
    evt.addListener(&Writer::write_separator,&w);
    evt.addListener(&Writer::write_world,&w);
    /* (Note that Writer::write_world is non-constant. That is OK, because w
     * is a non-constant object, so w.write_world(ss) is valid code. This would
     * not compile if we would have written 'const Writer w;' above.)
     */

    /* 3. Emit the event. This calls all registered callnacks, in the order in
     *    which they were registered. The stringstream ss is passed to each
     *    (a *reference*, as you can see from the definition of OutputEvent).
     *    The registered callbacks will write "Hello", ", " and "world!" to the
     *    stream, so after emit, the string stream should contain the inevitable
     *    message.
     */
    std::stringstream ss;
    evt.emit(ss);
    test_expr(ss.str()=="Hello, world!");
}

int main()
{
    test1();

    test_report;
    return nerrors;
}
