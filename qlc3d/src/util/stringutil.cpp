#include <util/stringutil.h>

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