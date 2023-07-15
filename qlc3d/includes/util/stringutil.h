#ifndef PROJECT_QLC3D_STRINGUTIL_H
#define PROJECT_QLC3D_STRINGUTIL_H

#include <vector>
#include <string>

class StringUtil {
public:
  /** Converts a vector of strings to a string representation of the vector, e.g. ["a", "b", "c"] */
  [[nodiscard]] static std::string toString(const std::vector<std::string> &vector);
};
#endif //PROJECT_QLC3D_STRINGUTIL_H
