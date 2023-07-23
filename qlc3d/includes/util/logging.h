#ifndef PROJECT_QLC3D_LOGGING_H
#define PROJECT_QLC3D_LOGGING_H
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/color.h>
#include <iostream>
#include <string>
#include <exception>
#include <filesystem>

template <>
class fmt::formatter<std::filesystem::path> {
public:
    constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
    template <typename Context>
    constexpr auto format (std::filesystem::path const& p, Context& ctx) const {
        return format_to(ctx.out(), "{}", p.string());
    }
};

class Log {

    static std::string indent;

    static void fallbackError(const std::string &message) {
        std::cerr << "[ERROR] Could not format log message when formatString=\"" << message << "\"" << std::endl;
        throw std::runtime_error("");
    }

public:
    template <typename... T>
    static void inline info(fmt::format_string<T...> formatString, T&&... args) {
        try {
            std::string message{"[INFO] " + Log::indent + fmt::vformat(formatString, fmt::make_format_args(args...)) + "\n"};
            fmt::print(message);
            fflush(stdout);
        } catch (...) {
            Log::fallbackError(fmt::to_string(formatString));
        }
    }

    template <typename... T>
    static void inline error(fmt::format_string<T...> formatString, T&&... args) {
        try {
            std::string message { "[ERROR] " + fmt::vformat(formatString, fmt::make_format_args(args...)) + "\n"};
            fmt::print(fmt::emphasis::bold | fmt::fg(fmt::color::orange_red), message);
            fflush(stdout);
        } catch (...) {
            Log::fallbackError(fmt::to_string(formatString));
        }
    }

    template <typename... T>
    static void inline warn(fmt::format_string<T...> formatString, T&&... args) {
        try {
            std::string message{"[WARNING] " + fmt::vformat(formatString, fmt::make_format_args(args...)) + "\n"};
            fmt::print(fmt::emphasis::bold, message);
            fflush(stdout);
        } catch (...) {
            Log::fallbackError(fmt::to_string(formatString));
        }
    }

    static void clearIndent() {
        Log::indent.resize(0);
    }

    static void incrementIndent(unsigned int amount = 1) {
        for (int i = 0; i < amount; i++) {
            Log::indent += ' ';
        }
    }

    static void decrementIndent(unsigned int amount = 1) {
        unsigned long newLength = Log::getIndent();
        if (amount < newLength) {
            newLength -= amount;
        } else {
            newLength = 0;
        }
        Log::indent.resize(newLength);
    }

    static unsigned int getIndent() {
        return Log::indent.length();
    }
};

#endif
