#ifndef PROJECT_QLC3D_STRINGUTIL_H
#define PROJECT_QLC3D_STRINGUTIL_H

#include <vector>
#include <string>

class StringUtil {
public:
  /** Converts a vector of strings to a string representation of the vector, e.g. ["a", "b", "c"] */
  [[nodiscard]] static std::string toString(const std::vector<std::string> &vector);

  /** Converts a string to lower case */
  [[nodiscard]] static std::string toLowerCase(const std::string in);
};
#endif //PROJECT_QLC3D_STRINGUTIL_H
