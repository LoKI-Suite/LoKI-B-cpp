//
// Created by daan on 14-5-19.
//

#ifndef LOKI_CPP_LOG_H
#define LOKI_CPP_LOG_H

#include <cstdint>
#include <string>
#include <iostream>

namespace loki {
    template <typename ErrorType>
    class Log {
    public:
        template<typename T> inline
        static void Warning(const T &message) {
            std::cerr << "[Warning] ";
            ErrorType::print(message);
        }

        template<typename T> inline
        static void Error(const T &message) {
            std::cerr << "[Error] ";
            ErrorType::print(message);
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
