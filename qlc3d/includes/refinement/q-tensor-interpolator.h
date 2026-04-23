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
 * @brief Interpolates Q-tensor values from an old mesh onto a new mesh.
 *
 * For each LC node in @p geomNew, the enclosing tetrahedron in @p geomOld
 * is found, barycentric coordinates are computed, and the five Q-tensor
 * components are interpolated linearly.  If the new node coincides with an
 * existing corner node of the enclosing tet (largest barycentric coordinate
 * >= 0.99999) the value is copied directly without interpolation.
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
