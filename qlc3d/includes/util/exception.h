#ifndef PROJECT_QLC3D_EXCEPTION_H
#define PROJECT_QLC3D_EXCEPTION_H
#include <string>

#define RUNTIME_ERROR(message) throw std::runtime_error((message) \
+ std::string("\nThrown in file ") + __FILE__\
+ ", in function " + __func__\
+ ", on line " + std::to_string(__LINE__));

/**
 * Throw in case encountering some condition that requires functionality that is not yet implemented,
 * but is planned to be in the future.
 */
class NotYetImplementedException : public std::runtime_error {
public:
    explicit NotYetImplementedException(const std::string &message) : std::runtime_error(message) {}
};

#endif //PROJECT_QLC3D_EXCEPTION_H
