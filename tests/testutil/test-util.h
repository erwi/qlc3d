#ifndef PROJECT_QLC3D_TEST_UTIL_H
#define PROJECT_QLC3D_TEST_UTIL_H

#include <filesystem>
#include <string>
#include <vector>

namespace TestUtil {

  const std::string RESOURCE_THIN_GID_MESH = "resources/thin.msh";
  const std::string RESOURCE_SMALL_CUBE_GMSH_MESH = "resources/gmsh-small-cube.msh";
  const std::string RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH = "resources/gmsh-small-cube-quadratic.msh";
  /** 1 x 0.1 x 1 mesh with left/right Neumann boundaries and front/back periodic boundaries. Electrode1/FixLC1 on top, Electrode2/FixLC2 on bottom */
  const std::string RESOURCE_PSEUDO_2D_NEUMANN_GMSH_MESH = "resources/pseudo-2d-neumann.msh";

  /** 1 x 1 x 1 mesh with left/right/front/back neumann surfaces. Electrode1/FixLC1 on top, Electrode2/FixLC2 on bottom */
  const std::string RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH = "resources/unit-cube-neumann.msh";
  /** 1 x 1 x 1 mesh with left/right/front/back neumann surfaces. Electrode1/FixLC1 on top, Electrode2/FixLC2 on bottom. Quadratic elements. */
  const std::string RESOURCE_UNIT_CUBE_NEUMANN_QUADRATIC_GMSH_MESH = "resources/unit-cube-neumann-quadratic.msh";


  /** 1 x 1 x 2 mesh with dielectric bottom half and LC top half. Neumann side boundaries. Electrodes on top and bottom. */
  const std::string RESOURCE_UNIT_CUBE_DIELECTRIC_NEUMAN_GMSH_MESH = "resources/unit-cube-dielectric-neumann.msh";
  const std::string RESOURCE_UNIT_CUBE_DIELECTRIC_NEUMAN_QUADRATIC_GMSH_MESH = "resources/unit-cube-dielectric-neumann-quadratic.msh";


  /** 1 x 1 x 1 mesh with periodic boundaries on all sides. */
  const std::string RESOURCE_UNIT_CUBE_FULLY_PERIODIC_LINEAR_GMSH_MESH = "resources/unit-cube-fully-periodic-linear-gmsh.msh";

  bool isEquivalentAngleDegrees(double a, double b, double epsilon = 1e-6);

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
