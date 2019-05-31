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

    template <typename ErrorType>
    class Log {
    public:
        template <typename T> inline
        static void Notify(const T &message) {
            std::cerr << BOLD << "[Notice] ";
            ErrorType::print(message);
            std::cerr << NC;
        }

        template<typename T> inline
        static void Warning(const T &message) {
            std::cerr << YELLOW << "[Warning] ";
            ErrorType::print(message);
            std::cerr << NC;
        }

        template<typename T> inline
        static void Error(const T &message) {
            std::cerr << RED << "[Error] ";
            ErrorType::print(message);
            std::cerr << NC;
            throw std::runtime_error("Loki Exception");
        }
    };

    struct Message {
        template<typename T> inline
        static void print(const T &t) {
            std::cerr << t << std::endl;
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
        template<typename T> inline
        static void print(const T &t) {
            std::cerr << "Could not find/open the " << t << " file." << std::endl;
        }
    };

    struct ParseFieldError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Could not properly parse the " << t << " field." << std::endl;
        }
    };

    struct ParseSectionError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Could not properly parse the " << t
                      << " section. Please check for missing values and improper "
                         "indentation." << std::endl;
        }
    };

    struct MissingSectionError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "The input file does not contain the " << t << " section." << std::endl;
        }
    };

    struct NumArgumentsError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Wrong number of input elements in the " << t << " property function." << std::endl;
        }
    };

    struct WrongPropertyError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Using " << t << " to set the wrong state property." << std::endl;
        }
    };

    struct PropertyFunctionError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Trying to call nonexistent property function: " << t << "." << std::endl;
        }
    };

    struct PropertyFunctionParseError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Could not parse function name and argument list: " << t << "." << std::endl;
        }
    };

    struct PropertyValueParseError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Could not parse property state value: " << t << "." << std::endl;
        }
    };

    struct PropertyArgumentsError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Cannot extract property function arguments from input file: " << t << "." << std::endl;
        }
    };

    struct PropertyStateError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Trying to assign property to nonexistent state: " << t << "." << std::endl;
        }
    };

    struct ChildrenPopulationError {
        template <typename T> inline
        static void print(const T &t) {
            std::cerr << "Populations of children do not add up to 1 for: " << t << "." << std::endl;
        }
    };
}

#endif //LOKI_CPP_LOG_H
