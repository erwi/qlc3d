#ifndef PROJECT_QLC3D_ELEMENT_SPLIT_H
#define PROJECT_QLC3D_ELEMENT_SPLIT_H
#include <vector>
#include <stdexcept>

/**
 * Functions for splitting and recombining quadratic tetrahedra/triangles to linear tetrahedra/triangles and vice versa. Node numbering
 * is as in fig 6.1 (b) in Eero's thesis.
 */

class ElementSplitCombineException : public std::runtime_error {
public:
  ElementSplitCombineException(const char *what) : std::runtime_error(what) {}
};


std::vector<std::vector<unsigned int>> splitQuadraticTetrahedronToLinear(const std::vector<unsigned int> &quadraticTetrahedron);
std::vector<unsigned int> recombineLinearTetsToQuadratic(const std::vector<std::vector<unsigned int>> &linearTets);

std::vector<std::vector<unsigned int>> splitQuadraticTriangleToLinear(const std::vector<unsigned int> &quadraticTriangle);
std::vector<unsigned int> recombineLinearTrianglesToQuadratic(const std::vector<std::vector<unsigned int>> &linearTriangles);

#endif //PROJECT_QLC3D_ELEMENT_SPLIT_H
