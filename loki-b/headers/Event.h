//
// Created by daan on 13-5-19.
//

#ifndef LOKI_CPP_EVENT_H
#define LOKI_CPP_EVENT_H

#include <vector>
#include <functional>
#include <memory>

namespace loki {
    template<typename ...T>
    class Event {
        std::vector<std::function<void(T&...)>> callbacks;

    public:
        Event() = default;

        Event(const Event &other) = delete;

        ~Event() = default;

        void emit(T& ...Args) {
            for (const auto &callback : callbacks)
                callback(Args...);
        }

        void addListener(std::function<void(T&...)> f) {
            callbacks.emplace_back(f);
        }

        template<class C>
        void addListener(void (C::*f)(T&... Args), C *c) {
            callbacks.emplace_back([c, f](T&... t) -> void { (c->*f)(t...); });
        }
    };

    using event = Event<>;
}

#endif //LOKI_CPP_EVENT_H
