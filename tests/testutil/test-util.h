#ifndef PROJECT_QLC3D_TEST_UTIL_H
#define PROJECT_QLC3D_TEST_UTIL_H

#include <filesystem>
#include <string>
#include <vector>

namespace TestUtil {

  const std::string RESOURCE_THIN_GID_MESH = "resources/thin.msh";
  const std::string RESOURCE_SMALL_CUBE_GMSH_MESH = "resources/gmsh-small-cube.msh";

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

  /** TemporaryDirectory is a directory created in the OS temp directory. It deletes itself when going out of scope. */
  class TemporaryDirectory {
    std::filesystem::path path_;
  public:
    ~TemporaryDirectory();
    TemporaryDirectory();
    [[nodiscard]] const std::filesystem::path &path() const { return path_; }
    [[nodiscard]] std::vector<std::filesystem::path> listFiles() const;

  };
}
#endif //PROJECT_QLC3D_TEST_UTIL_H
