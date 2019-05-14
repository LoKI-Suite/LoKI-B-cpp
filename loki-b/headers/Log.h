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
#define NC      "\e[0m"

namespace loki {

    template <typename ErrorType>
    class Log {
    public:
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

    struct GeneralError {
        template<typename T> inline
        static void print(const T &t) {
            std::cerr << t << std::endl;
        }
    };

    struct FileError {
        template<typename T> inline
        static void print(const T &t) {
            std::cerr << t << std::endl;
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

    // Log<FileError>::Warning();
}

#endif //LOKI_CPP_LOG_H
