//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_LOG_H
#define LOKI_CPP_LOG_H

#include <cstdint>
#include <string>
#include <iostream>

#define RED     "\e[1;31m"
#define YELLOW  "\e[1;33m"
#define BOLD    "\e[1;37m"
#define NC      "\e[0m"

namespace loki {

    template<typename ErrorType>
    class Log {
    public:
        template<typename ...T>
        inline
        static void Notify(const T &...message) {
            std::cerr << BOLD << "[Notice] ";
            ErrorType::print(message...);
            std::cerr << NC;
        }

        template<typename ...T>
        inline
        static void Warning(const T &...message) {
            std::cerr << YELLOW << "[Warning] ";
            ErrorType::print(message...);
            std::cerr << NC;
        }

        template<typename ...T>
        inline
        static void Error(const T &...message) {
            std::cerr << RED << "[Error] ";
            ErrorType::print(message...);
            std::cerr << NC;
            throw std::runtime_error("Loki Exception");
        }
    };

    struct Message {
        template<typename ...T>
        inline
        static void print(const T &...t) {
            (std::cerr << ... << t) << std::endl;
        }
    };

    struct GasPropertyError {
        template<typename T>
        static void print(const T &t) {
            std::cerr << "Could not parse gas property: " << t << '.' << std::endl;
        }
    };

    struct DoubleCollision {
        template<typename T>
        static void print(const T &t) {
            std::cerr << "Double collision detected: " << t << '.' << std::endl;
        }
    };

    struct LXCatError {
        template<typename T>
        static void print(const T &t) {
            std::cerr << "Could not properly parse " << t << '.' << std::endl;
        }
    };

    struct FileError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Could not find/open the " << t << " file." << std::endl;
        }
    };

    struct ParseFieldError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Could not properly parse the " << t << " field." << std::endl;
        }
    };

    struct ParseSectionError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Could not properly parse the " << t
                      << " section. Please check for missing values and improper "
                         "indentation." << std::endl;
        }
    };

    struct MissingSectionError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "The input file does not contain the " << t << " section." << std::endl;
        }
    };

    struct NumArgumentsError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Wrong number of input elements in the " << t << " property function." << std::endl;
        }
    };

    struct WrongPropertyError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Using " << t << " to set the wrong state property." << std::endl;
        }
    };

    struct PropertyFunctionError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Trying to call nonexistent property function: " << t << "." << std::endl;
        }
    };

    struct PropertyFunctionParseError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Could not parse function name and argument list: " << t << "." << std::endl;
        }
    };

    struct PropertyValueParseError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Could not parse property state value: " << t << "." << std::endl;
        }
    };

    struct PropertyArgumentsError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Cannot extract property function arguments from input file: " << t << "." << std::endl;
        }
    };

    struct PropertyStateError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Trying to assign property to nonexistent state: " << t << "." << std::endl;
        }
    };

    struct ChildrenPopulationError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Populations of children do not add up to 1 for: " << t << "." << std::endl;
        }
    };

    struct ZeroFractionPopulationError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr << "Populations of children do not add up to 0 for gas with fraction 0: " << t << "." << std::endl;
        }
    };

    struct NegativeElastic {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Negative values encountered when extracting elastic cross section from the effective cross section of "
                    << t << ". Negative values have been clipped to 0 and results might be unreliable." << std::endl;
        }
    };

    struct MultipleProductsInReverse {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Detected superelastic collision with multiple products; "
                    << t << "." << std::endl;
        }
    };

    struct MultipleReactantInEedfCol {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Detected EedfCollision with more than one reactant; "
                    << t << "." << std::endl;
        }
    };

    struct SuperElasticForNonReverse {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Super elastic cross section requested for non reverse collision; "
                    << t << "." << std::endl;
        }
    };



    struct NoEnergy {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Energy not found for state: "
                    << t << "." << std::endl;
        }
    };

    struct NoElectricQuadMoment {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Electrical quadrupole moment not found for state: "
                    << t << "." << std::endl;
        }
    };

    struct NoRotationalConstant {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Rotational constant not found for gas: "
                    << t << "." << std::endl;
        }
    };

    struct NoStatWeight {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Statistical weight not found for state: "
                    << t << "." << std::endl;
        }
    };

    struct RotCollisionInCARGas {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Rotational collision found for CAR candidate: "
                    << t << "." << std::endl;
        }
    };

    struct CARForNonExistent {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Trying to enable CAR for non-existent gas: "
                    << t << "." << std::endl;
        }
    };

    struct PowerRatioError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "EEDF in e-e collision routine has converged, but abs(Pee/Pref) = "
                    << t << " > 1e-9." << std::endl;
        }
    };

    struct PowerBalanceError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Relative power balance error larger than "
                    << t << ". Results might be incorrect." << std::endl;
        }
    };

    struct GlobalIterError {
        template<typename T>
        inline
        static void print(const T &t) {
            std::cerr
                    << "Eedf in global cycle, invoked while mixing solutions, did not converge after "
                    << t << " iterations." << std::endl;
        }
    };
}

#endif //LOKI_CPP_LOG_H
