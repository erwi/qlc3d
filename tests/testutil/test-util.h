//
// Created by eero on 07/04/2021.
//
#ifndef PROJECT_QLC3D_TEST_UTIL_H
#define PROJECT_QLC3D_TEST_UTIL_H

#include <string>
#include <vector>

namespace TestUtil {
    /**
     * TemporaryFile is a file created in the OS temp directory. It deletes itself when going out of scope.
     */
    class TemporaryFile {
        std::string name_;
        TemporaryFile();
    public:

        ~TemporaryFile();
        [[nodiscard]] const std::string &name() const { return name_; }
        static TemporaryFile empty();
        static TemporaryFile withContents(const std::string &fileContents);

        /**
         * Read the file contents to a vector of strings, each string corresponding to one line
         * of the file.
         * @return the file contents.
         */
        std::vector<std::string> readContentsAsText() const;
    };
}
#endif //PROJECT_QLC3D_TEST_UTIL_H
