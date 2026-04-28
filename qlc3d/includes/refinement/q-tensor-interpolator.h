#ifndef Q_TENSOR_INTERPOLATOR_H
#define Q_TENSOR_INTERPOLATOR_H

#include <geometry.h>
#include <solutionvector.h>

/**
 * @brief Linearly interpolates a scalar field using four barycentric weights.
 *
 * @param loc  Four barycentric coordinates summing to 1.
 * @param S    Scalar field values at the four tet corners.
 * @return     Interpolated scalar value.
 */
double interpolateScalar(const double loc[4], const double S[4]);

/**
 * @brief Evaluates all 10 TET10 quadratic shape functions at a barycentric point.
 *
 * The TET10 element uses volume (barycentric) coordinates @c loc[0..3] where
 * @c loc[0] = 1 - loc[1] - loc[2] - loc[3].  Node ordering follows the Gmsh
 * TET10 convention: nodes 0–3 are corner nodes and nodes 4–9 are mid-edge nodes
 * in the order AB, BC, AC, AD, CD, BD.
 *
 * This function is exposed in the header to allow direct unit testing of the
 * algebraic correctness of the shape functions.
 *
 * @param loc  Four barycentric (volume) coordinates summing to 1.
 * @param N    Output array of length 10; N[i] receives the i-th shape function value.
 */
void evaluateTet10ShapeFunctions(const double loc[4], double N[10]);

/**
 * @brief Interpolates Q-tensor values from an old mesh onto a new mesh.
 *
 * For each LC node in @p geomNew, the enclosing tetrahedron in @p geomOld
 * is found and barycentric coordinates are computed.
 *
 * - If the new node coincides with an existing corner node of the enclosing
 *   tet (largest barycentric coordinate >= 0.99999) the value is copied
 *   directly without interpolation.
 * - If the source geometry is a TET10 (QUADRATIC_TETRAHEDRON) mesh all 10
 *   nodes of the enclosing element contribute via the quadratic shape
 *   functions.  This reproduces degree-2 polynomial fields exactly.
 * - If the source geometry is a TET4 (LINEAR_TETRAHEDRON) mesh the four
 *   corner nodes are interpolated linearly using barycentric weights.
 *
 * @param geomNew  Target geometry whose LC nodes need Q values.
 * @param geomOld  Source geometry from which Q values are read.
 * @param qOld     Q-tensor solution on @p geomOld; must have 5 degrees of
 *                 freedom per node.
 * @return         A new SolutionVector sized to geomNew.getnpLC() x 5
 *                 containing the interpolated Q-tensor values.
 */
SolutionVector interpolateQTensor(Geometry &geomNew,
                                  Geometry &geomOld,
                                  const SolutionVector &qOld);

#endif // Q_TENSOR_INTERPOLATOR_H

