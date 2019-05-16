//
// Created by daan on 16-5-19.
//

#include <random>
#include <chrono>
#include "GasMixture.h"
#include "Eigen/Dense"

//#include "Working_v1.h"

int main(int argc, char **argv) {
    GasMixture<Boltzmann> mixture;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<uint16_t> distLevel(2, 3);
    std::uniform_int_distribution<uint16_t> dist(0, 200);

    auto begin = std::chrono::high_resolution_clock::now();


    for (int i = 0; i < 10000; ++i) {
        StateEntry rEntry{rotational, "N2", "X", dist(mt), dist(mt), 0};
        mixture.addState(rEntry);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "mus" << std::endl;

    std::cout << mixture.gasses[0].states.size() << std::endl;

//    StateEntry entry{ionic, "N2", "X", 3, 5, 1};
//    StateEntry entry2{rotational, "N2", "X", 3, 8, 0};
//    StateEntry entry3{vibrational, "N2", "A", 9, 0, 0};

//    mixture.addState(entry);
//    mixture.addState(entry2);
//    mixture.addState(entry3);

    return 0;
}