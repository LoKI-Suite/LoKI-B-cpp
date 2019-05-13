//
// Created by daan on 12-5-19.
//

#include <iostream>
#include <functional>
#include <vector>
#include <cmath>

template <typename ... T>
class Event {
    std::vector<std::function<void(T...)>> callbacks;

public:
    Event() = default;
    Event(const Event &other) = delete;
    ~Event() = default;

    void emit(T ... Args) {
        for (const auto &callback : callbacks)
            callback(Args...);
    }

    void addListener(std::function<void(T...)> f) {
        callbacks.emplace_back(f);
    }

    template<class C>
    void addListener(void (C::*f)(T... Args), C &c) {
        callbacks.emplace_back([&c, f](T... t)->void{(c.*f)(t...);});
    }
};

void callback_sum(double a, double b) {
    std::cout << "sum = " << a + b << std::endl;
}

class anyClass {
public:
    double a = 0., b = 0.;

    anyClass(double a, double b) : a(a), b(b) {}
    anyClass(const anyClass &other) {a = other.a; b = other.b; std::cout << "I have been coppied" << std::endl;}
    ~anyClass() = default;

    void displayChebyshev(double p) {
        std::cout << pow(pow(a, p) + pow(b, p), 1./p) << std::endl;
    }

    void listen() {
        std::cout << "I'm listening: " << a << std::endl;
    }
};


int main (int argc, char ** argv)
{
    // Example for a static function.
    Event<double,double> event_sum;
    event_sum.addListener(callback_sum);
    event_sum.emit(5.0, 2.5);

    // Example for a member function, using a lambda expression.
    anyClass test_class(3., 4.);
    Event<double> event_double;
    event_double.addListener([&](double p)->void{test_class.displayChebyshev(p);});
    event_double.emit(2.);
    test_class.a = 5.;
    event_double.emit(2.);

    // This will probably be the most used version in loki.
    Event<> loki_event;
    loki_event.addListener(&anyClass::listen, test_class);
//    loki_event.addListener([&]()->void{test_class.listen();});
    loki_event.emit();

    // Since the object is passed by reference, the behaviour of the listening function
    // can change, and thus be controlled, by 
    test_class.a = 3.;
    loki_event.emit();

    return 0;
}