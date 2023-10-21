#ifndef PROJECT_QLC3D_TEST_UTIL_H
#define PROJECT_QLC3D_TEST_UTIL_H

#include <filesystem>
#include <string>
#include <vector>

namespace TestUtil {

  const std::string RESOURCE_THIN_GID_MESH = "resources/thin.msh";
  const std::string RESOURCE_SMALL_CUBE_GMSH_MESH = "resources/gmsh-small-cube.msh";
  /** 1 x 0.1 x 1 mesh with left/right Neumann boundaries and front/back periodic boundaries. Electrode1/FixLC1 on top, Electrode2/FixLC2 on bottom */
  const std::string RESOURCE_PSEUDO_2D_NEUMANN_GMSH_MESH = "resources/pseudo-2d-neumann.msh";

  /** 1 x 1 x 1 mesh with left/right/front/back neumann surfaces. Electrode1/FixLC1 on top, Electrode2/FixLC2 on bottom */
  const std::string RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH = "resources/unit-cube-neumann.msh";

  /**
   * TemporaryFile is a file created in the OS temp directory. It deletes itself when going out of scope.
   */
  class TemporaryFile {
    std::filesystem::path name_;
    TemporaryFile();
  public:

    ~TemporaryFile();
    [[nodiscard]] const std::filesystem::path &name() const { return name_; }
    static TemporaryFile empty();
    static TemporaryFile withContents(const std::string &fileContents);

    /**
     * Read the file contents to a vector of strings, each string corresponding to one line
     * of the file.
     * @return the file contents.
     */
    [[nodiscard]] std::vector<std::string> readContentsAsText() const;
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
