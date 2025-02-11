#include <util/stringutil.h>
#include <algorithm>
#include <regex>

std::string StringUtil::toString(const std::vector<std::string> &vector) {
    std::string result = "[";

    std::size_t size = vector.size();
    for (std::size_t i = 0; i < size; ++i) {
        result += vector[i];
        if (i < size - 1) {
            result += ", ";
        }
    }

    result += "]";
    return result;
}

std::string StringUtil::toLowerCase(const std::string in) {
    std::string out = in;
    std::transform(out.begin(), out.end(), out.begin(), ::tolower);
    return out;
}

std::vector<std::string> StringUtil::split(const std::string &str, const std::string &regexDelimiter) {
    std::regex re(regexDelimiter);
    std::sregex_token_iterator it(str.begin(), str.end(), re, -1);
    std::sregex_token_iterator end;
    std::vector<std::string> result(it, end);
    return result;
}