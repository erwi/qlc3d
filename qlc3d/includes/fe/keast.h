#ifndef PROJECT_QLC3D_KEAST_H
#define PROJECT_QLC3D_KEAST_H
#include <vector>
#include <cassert>
/*
 * "Moderate degree tetrahedral quadrature formulas"
 * Keast, P.
 * Computer Methods in Applied Mechanics and Engineering, 55(3), 339-348, 1986
 * https://doi.org/10.1016/0045-7825(86)90059-2
 *
 *
 * NOTE: An additional division by 6 of the weights to account for the volume of the tetrahedron.
 */


struct IntegrationPoints {
  /** max polynomial order */
  const int p;
  const std::vector<double> points;
  const std::vector<double> weights;

  [[nodiscard]] unsigned int numGaussPoints() const { return weights.size(); }

  /**
   * Number of dimensions of the integration points. For example, 3 for 3D tetrahedron points, 2 for 2D triangle points.
   */
  [[nodiscard]] unsigned int numDimensions() const {
    assert(points.size() % numGaussPoints() == 0);
    return points.size() / numGaussPoints();
  }
};

/** Keast integration points for tetrahedral element with 1st order polynomial, with 1 gauss point */
const IntegrationPoints Keast0 = {
        1,
        {0.25, 0.25, 0.25},
        {1. / 6}
};

/** Keast integration points for tetrahedral element with 4th order polynomial, with 11 gauss points */
const IntegrationPoints Keast4 = {
        4,
        {0.2500000000000000,  0.2500000000000000,  0.2500000000000000,
         0.7857142857142857,  0.0714285714285714,  0.0714285714285714,
         0.0714285714285714,  0.0714285714285714,  0.0714285714285714,
         0.0714285714285714,  0.0714285714285714,  0.7857142857142857,
         0.0714285714285714,  0.7857142857142857,  0.0714285714285714,
         0.1005964238332008,  0.3994035761667992,  0.3994035761667992,
         0.3994035761667992,  0.1005964238332008,  0.3994035761667992,
         0.3994035761667992,  0.3994035761667992,  0.1005964238332008,
         0.3994035761667992,  0.1005964238332008,  0.1005964238332008,
         0.1005964238332008,  0.3994035761667992,  0.1005964238332008,
         0.1005964238332008,  0.1005964238332008,  0.3994035761667992},

        {
                -0.0789333333333333 / 6,
                0.0457333333333333 / 6,
                0.0457333333333333 / 6,
                0.0457333333333333 / 6,
                0.0457333333333333 / 6,
                0.1493333333333333 / 6,
                0.1493333333333333 / 6,
                0.1493333333333333 / 6,
                0.1493333333333333 / 6,
                0.1493333333333333 / 6,
                0.1493333333333333 / 6,
        }
};


/*
 * Weights and Gauss points from table 10.3 in J. E. Akin
 */
const IntegrationPoints Tri4thOrder {
  4,
  {0    , 0,
   0.5  , 0,
   1.0  , 0,
   0.5  , 0.5,
   0    , 1.0,
   0    , 0.5,
   1./3., 1./3.},
  {1. / 40, 1./ 15, 1. / 40, 1. / 15, 1. / 40, 1. / 15, 9. / 40}
};

#endif //PROJECT_QLC3D_KEAST_H
