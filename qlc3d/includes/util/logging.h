#ifndef PROJECT_QLC3D_LOGGING_H
#define PROJECT_QLC3D_LOGGING_H
#include <fmt/core.h>
#include <fmt/format.h>
#include <iostream>
#include <string>

class Log {
public:
    template <typename... T>
    static void inline info(fmt::format_string<T...> formatString, T&&... args) {
        std::string message { fmt::vformat(formatString, fmt::make_format_args(args...)) };
        std::cout << "[INFO] " << message << std::endl;
    }

    template <typename... T>
    static void inline error(fmt::format_string<T...> formatString, T&&... args) {
        std::string message { fmt::vformat(formatString, fmt::make_format_args(args...)) };
        std::cerr << "[ERROR] " << message << std::endl;
    }

    template <typename... T>
    static void inline warn(fmt::format_string<T...> formatString, T&&... args) {
        std::string message { fmt::vformat(formatString, fmt::make_format_args(args...)) };
        std::cout << "[WARNING] " << message << std::endl;
    }
};

#endif
