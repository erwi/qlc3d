//
// Created by eero on 07/04/2021.
//
#ifndef PROJECT_QLC3D_TEST_UTIL_H
#define PROJECT_QLC3D_TEST_UTIL_H

#include <string>

namespace TestUtil {

    /**
     * TemporaryFile is a file created in the OS temp directory. It deletes itself when going out of scope.
     */
    class TemporaryFile {
        std::string name_;
        TemporaryFile();
    public:

        ~TemporaryFile();
        const std::string &name() { return name_; }
        static TemporaryFile empty();
        static TemporaryFile withContents(const std::string &fileContents);
    };
}
#endif //PROJECT_QLC3D_TEST_UTIL_H
