//
// Created by daan on 26-06-2019.
//

#include <cstdint>
#include <iostream>

struct Doubles {
    double one{8.}, two{5.}, three{2.};
};

int main (int argc, char ** argv) {
    Doubles d;

    auto *ptr = (double*)&d;

    for (uint32_t i = 0; i < 3; ++i) {
        std::cout << ptr[i] << std::endl;
    }

    return 0;
}