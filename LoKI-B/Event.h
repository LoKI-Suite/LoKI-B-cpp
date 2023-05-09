/** \file
 *
 *  An Event class template for usage with LoKI-B.
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

#ifndef LOKI_CPP_EVENT_H
#define LOKI_CPP_EVENT_H

#include <functional>
#include <vector>
#include <iostream>

namespace loki
{

/** An event class is characterized by a list of types that define the signature
 *  of callback functions, or listeners, that accepts a compatible list of
 *  arguments and returns void. An event object manages a vector of such
 *  callback functions that can be invoked by calling member emit, passing a
 *  compatible argument list. Callback functions can be added with the help of
 *  one of the overloads of member addListener. When emit is called, the callback
 *  functions are guaranteed to be called in the order in which they were
 *  registered with the event.
 *  An event object is default-constructable but non-copyable.
 *
 *  \author Daan Boer
 *  \date   13. May 2019
 */
template <typename... T>
class Event
{
  public:
    /// The default constructor of the Event class is available
    Event() = default;
    /// Event objects are non-copyable, the default copy constructor is deleted.
    Event(const Event &other) = delete;
    /// The default destructor of the Event class is used.
    ~Event() = default;
    /** A callback functions is a std::function that returns a void and can
     *  be called with an argument list that is defined by the list of template
     *  arguments.
     */
    using Callback = std::function<void(T...)>;

    /** Call all registered callback functions with arguments \a args.
     *  The callbacks are called in the order of registration with the
     *  event object.
     */
    void emit(T... args) const
    {
        for (const auto &callback : m_callbacks)
        {
            callback(args...);
        }
    }
    /// Add function \a f to the list of callback functions.
    void addListener(Callback f)
    {
        m_callbacks.emplace_back(f);
    }
    /** This overload adds a callback function that is created from a member
     *  function pointer \a f and the object \a c on which \a f is to be called.
     */
    template <class C>
    void addListener(void (C::*f)(T... Args), C *c)
    {
        m_callbacks.emplace_back([c, f](T... t) -> void { (c->*f)(t...); });
    }
    /// an overload of addListener that accepts a constant object and member function
    template <class C>
    void addListener(void (C::*f)(T... Args) const, const C *c)
    {
        m_callbacks.emplace_back([c, f](T... t) -> void { (c->*f)(t...); });
    }

  private:
    /** The container of registered callbacks.
     */
    std::vector<Callback> m_callbacks;
};

} // namespace loki

#endif // LOKI_CPP_EVENT_H
