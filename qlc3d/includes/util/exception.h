#ifndef PROJECT_QLC3D_EXCEPTION_H
#define PROJECT_QLC3D_EXCEPTION_H
#include <string>

#define RUNTIME_ERROR(message) throw std::runtime_error((message) \
+ std::string("\nThrown in file ") + __FILE__\
+ ", in function " + __func__\
+ ", on line " + std::to_string(__LINE__));

#endif //PROJECT_QLC3D_EXCEPTION_H
