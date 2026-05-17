#ifndef PROJECT_QLC3D_GAUSSIAN_QUADRATURE_H
#define PROJECT_QLC3D_GAUSSIAN_QUADRATURE_H
#include <geom/vec3.h>
#include <lc-representation.h>
#include <fe/keast.h>
#include <util/exception.h>

#include <vector>
#include <cassert>

/**
 * Common base class for both triangle and tetrahedron shape functions.
 */
class ShapeFunction {
protected:
  unsigned int elementOrder;
  unsigned int nodesPerElement = 0;
  unsigned int numGaussPoints = 0;
  unsigned int currentPoint = 0;

  const IntegrationPoints *integrationPoints = nullptr;

  /* Local coordinates sampled at integration points of the element */
  std::vector<double> sh;

  explicit ShapeFunction(unsigned int elementOrder) : elementOrder(elementOrder) {}

  [[nodiscard]] const double& get(const std::vector<double> &vec, unsigned int i) const { return vec[currentPoint * nodesPerElement + i]; }
public:
  virtual void setIntegrationPoints(const IntegrationPoints &integrationPoints) = 0;
  [[nodiscard]] unsigned int getNumGaussPoints() const { return numGaussPoints; }
  [[nodiscard]] double getWeight() const { return integrationPoints->weights[currentPoint]; }
  [[nodiscard]] unsigned int getNumPointsPerElement() const { return nodesPerElement; }
  [[nodiscard]] bool hasNextPoint() {
    bool hasNext = currentPoint < numGaussPoints;
    if (!hasNext) {
      currentPoint = 0; // side effect: reset to start so that next shape function can be calculated. TODO: should we reset it explicitly?
    }
    return hasNext;

  }
  void nextPoint() { currentPoint++; }

  /**
   * Get the i'th shape function value for the current integration point
   * @param i = 0..3 for linear tetrahedron, 0..9 for quadratic tetrahedron
   */
  [[nodiscard]] const double& N(unsigned int i) const { return get(sh, i); }
  [[nodiscard]] double sample(const double *value) const {
    double sum = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      sum += value[i] * N(i);
    }
    return sum;
  }


  template<typename Src>
  void sampleQ(const Src &source, double &q1, double &q2, double &q3, double &q4, double &q5) const {
    q1 = q2 = q3 = q4 = q5 = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      q1 += source[i][0] * N(i);
      q2 += source[i][1] * N(i);
      q3 += source[i][2] * N(i);
      q4 += source[i][3] * N(i);
      q5 += source[i][4] * N(i);
    }
  }

  template<typename Src>
  void sampleAll(const Src* source, double &v1, double &v2, double &v3, double &v4, double &v5, double &v6) const {
    v1 = v2 = v3 = v4 = v5 = v6 = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      v1 += source[i][0] * N(i);
      v2 += source[i][1] * N(i);
      v3 += source[i][2] * N(i);
      v4 += source[i][3] * N(i);
      v5 += source[i][4] * N(i);
      v6 += source[i][5] * N(i);
    }
  }

  template<typename Src>
  void sample(const Src source, Vec3 &destination) const {
    destination.set(0, 0, 0);
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      destination.add(source[i].x() * N(i),
                      source[i].y() * N(i),
                      source[i].z() * N(i));
    }
  }
};

/**
 * Shape function for Triangle element. This is pretty simple as no gradients are currently needed.
 */
class TriShapeFunction : public ShapeFunction {
  std::vector<double> shR;
  std::vector<double> shS;

  void setLinearTrianglePoints(const IntegrationPoints &integrationPoints) {
    nodesPerElement = 3;

    sh.resize(numGaussPoints * nodesPerElement, 0);
    shR.resize(numGaussPoints * nodesPerElement, 0);
    shS.resize(numGaussPoints * nodesPerElement, 0);

    for(unsigned int i = 0; i < numGaussPoints; i++) {
      double r = integrationPoints.points[i * 2 + 0];
      double s = integrationPoints.points[i * 2 + 1];

      sh[i * nodesPerElement + 0] = 1 - r - s;
      sh[i * nodesPerElement + 1] = r;
      sh[i * nodesPerElement + 2] = s;

      shR[i * nodesPerElement + 0] = -1.0;
      shR[i * nodesPerElement + 1] = 1.0;
      shR[i * nodesPerElement + 2] = 0.0;

      shS[i * nodesPerElement + 0] = -1.0;
      shS[i * nodesPerElement + 1] = 0.0;
      shS[i * nodesPerElement + 2] = 1.0;
    }
  }

  void setQuadraticTrianglePoints(const IntegrationPoints &integrationPoints) {
    nodesPerElement = 6;
    sh.resize(numGaussPoints * nodesPerElement, 0);
    shR.resize(numGaussPoints * nodesPerElement, 0);
    shS.resize(numGaussPoints * nodesPerElement, 0);

    for (unsigned int i = 0; i < numGaussPoints; i++) {
      double r = integrationPoints.points[i * 2 + 0];
      double s = integrationPoints.points[i * 2 + 1];

      sh[i * nodesPerElement + 0] = (1 - r - s) * (1 - 2 * r - 2 * s);
      sh[i * nodesPerElement + 1] = r * (2 * r - 1);
      sh[i * nodesPerElement + 2] = s * (2 * s - 1);
      sh[i * nodesPerElement + 3] = 4 * r * (1 - r - s);
      sh[i * nodesPerElement + 4] = 4 * r * s;
      sh[i * nodesPerElement + 5] = 4 * s * (1 - r - s);

      shR[i * nodesPerElement + 0] = -3 + 4 * r + 4 * s;
      shR[i * nodesPerElement + 1] = 4 * r - 1;
      shR[i * nodesPerElement + 2] = 0;
      shR[i * nodesPerElement + 3] = 4 - 8 * r - 4 * s;
      shR[i * nodesPerElement + 4] = 4 * s;
      shR[i * nodesPerElement + 5] = -4 * s;

      shS[i * nodesPerElement + 0] = -3 + 4 * r + 4 * s;
      shS[i * nodesPerElement + 1] = 0;
      shS[i * nodesPerElement + 2] = 4 * s - 1;
      shS[i * nodesPerElement + 3] = -4 * r;
      shS[i * nodesPerElement + 4] = 4 * r;
      shS[i * nodesPerElement + 5] = 4 - 4 * r - 8 * s;
    }
  }

public:
  TriShapeFunction(unsigned int elementOrder) : ShapeFunction(elementOrder) {}

  void setIntegrationPoints(const IntegrationPoints &integrationPoints) override {
    if (this->integrationPoints != nullptr) {
      return; // already initialised
    }
    this->integrationPoints = &integrationPoints;
    numGaussPoints = integrationPoints.weights.size();
    assert(numGaussPoints == integrationPoints.points.size() / 2);

    if (elementOrder == 1) {
      setLinearTrianglePoints(integrationPoints);
    } else if (elementOrder == 2) {
      setQuadraticTrianglePoints(integrationPoints);
    } else {
      RUNTIME_ERROR("Unsupported element order " + std::to_string(elementOrder));
    }
  }
};

/**
 * Shape function for quadratic or linear tetrahedral elements (TET10 or TET4).
 *
 * Stores pre-computed shape function values and their reference-space derivatives at all Gauss
 * points.  The global (x,y,z) derivatives are recomputed each time @c initialiseElement is
 * called by evaluating the isoparametric Jacobian at the current Gauss point.
 */
class TetShapeFunction : public ShapeFunction {
protected:
  std::vector<double> shR;
  std::vector<double> shS;
  std::vector<double> shT;

  std::vector<double> shX;
  std::vector<double> shY;
  std::vector<double> shZ;

  [[nodiscard]] const double& getShR(unsigned int i) const { return get(shR, i); }
  [[nodiscard]] const double& getShS(unsigned int i) const { return get(shS, i); }
  [[nodiscard]] const double& getShT(unsigned int i) const { return get(shT, i); }

private:
  void initialiseLinearTet() {
    assert(integrationPoints != nullptr);
    nodesPerElement = 4;
    unsigned int numDimensions = integrationPoints->numDimensions();
    assert(numDimensions == 2 || numDimensions == 3);
    sh.resize(numGaussPoints * nodesPerElement, 0);
    shR.resize(numGaussPoints * nodesPerElement, 0);
    shS.resize(numGaussPoints * nodesPerElement, 0);
    shT.resize(numGaussPoints * nodesPerElement, 0);

    for (unsigned int i = 0; i < numGaussPoints; ++i) {
      double r = integrationPoints->points[i * numDimensions + 0];
      double s = integrationPoints->points[i * numDimensions + 1];
      double t = numDimensions == 2 ? 0 : integrationPoints->points[i * numDimensions + 2];

      sh[i * nodesPerElement + 0] = 1 - r - s - t;
      sh[i * nodesPerElement + 1] = r;
      sh[i * nodesPerElement + 2] = s;
      sh[i * nodesPerElement + 3] = t;

      shR[i * nodesPerElement + 0] = -1.0;
      shR[i * nodesPerElement + 1] = 1.0;
      shR[i * nodesPerElement + 2] = 0.0;
      shR[i * nodesPerElement + 3] = 0.0;

      shS[i * nodesPerElement + 0] = -1.0;
      shS[i * nodesPerElement + 1] = 0.0;
      shS[i * nodesPerElement + 2] = 1.0;
      shS[i * nodesPerElement + 3] = 0.0;

      shT[i * nodesPerElement + 0] = -1.0;
      shT[i * nodesPerElement + 1] = 0.0;
      shT[i * nodesPerElement + 2] = 0.0;
      shT[i * nodesPerElement + 3] = 1.0;
    }

    shX.resize(nodesPerElement, 0);
    shY.resize(nodesPerElement, 0);
    shZ.resize(nodesPerElement, 0);

    for (unsigned int i = 0; i < nodesPerElement; i++) {
      shX[i] = 0.;
      shY[i] = 0.;
      shZ[i] = 0.;
    }
  }

  void initialiseQuadraticTet() {
    assert(integrationPoints != nullptr);

    unsigned int dimensions = integrationPoints->numDimensions();
    assert(dimensions == 3 || dimensions == 2);

    nodesPerElement = 10;

    sh.resize(numGaussPoints * nodesPerElement, 0);
    shR.resize(numGaussPoints * nodesPerElement, 0);
    shS.resize(numGaussPoints * nodesPerElement, 0);
    shT.resize(numGaussPoints * nodesPerElement, 0);

    for (unsigned int i = 0; i < numGaussPoints; ++i) {
      double r = integrationPoints->points[i * dimensions + 0];
      double s = integrationPoints->points[i * dimensions + 1];
      // for triangle weights, local coordinate t is always 0. This case is used for surface integrals
      double t = dimensions == 2 ? 0 : integrationPoints->points[i * dimensions + 2];

      // corner nodes expressed in natural coordinates
      double N1 = 1 - r - s - t;
      double N2 = r;
      double N3 = s;
      double N4 = t;

      sh[i * nodesPerElement + 0] = N1 * (2 * N1 - 1);
      sh[i * nodesPerElement + 1] = N2 * (2 * N2 - 1);
      sh[i * nodesPerElement + 2] = N3 * (2 * N3 - 1);
      sh[i * nodesPerElement + 3] = N4 * (2 * N4 - 1);

      // mid-edge nodes — Gmsh TET10 ordering: AB, BC, AC, AD, CD, BD
      sh[i * nodesPerElement + 4] = 4 * N1 * N2;  // AB
      sh[i * nodesPerElement + 5] = 4 * N2 * N3;  // BC
      sh[i * nodesPerElement + 6] = 4 * N3 * N1;  // AC
      sh[i * nodesPerElement + 7] = 4 * N1 * N4;  // AD
      sh[i * nodesPerElement + 8] = 4 * N3 * N4;  // CD  (Gmsh [8])
      sh[i * nodesPerElement + 9] = 4 * N2 * N4;  // BD  (Gmsh [9])


      shR[i * nodesPerElement + 0] = 4 * r + 4 * s + 4 * t - 3;
      shS[i * nodesPerElement + 0] = 4 * r + 4 * s + 4 * t - 3;
      shT[i * nodesPerElement + 0] = 4 * r + 4 * s + 4 * t - 3;

      shR[i * nodesPerElement + 1] = 4 * r - 1;
      shS[i * nodesPerElement + 1] = 0;
      shT[i * nodesPerElement + 1] = 0;

      shR[i * nodesPerElement + 2] = 0;
      shS[i * nodesPerElement + 2] = 4 * s - 1;
      shT[i * nodesPerElement + 2] = 0;

      shR[i * nodesPerElement + 3] = 0;
      shS[i * nodesPerElement + 3] = 0;
      shT[i * nodesPerElement + 3] = 4 * t - 1;

      // [4] = AB: 4*N1*N2 = 4*(1-r-s-t)*r
      shR[i * nodesPerElement + 4] = -8 * r - 4 * s - 4 * t + 4;
      shS[i * nodesPerElement + 4] = -4 * r;
      shT[i * nodesPerElement + 4] = -4 * r;

      // [5] = BC: 4*N2*N3 = 4*r*s
      shR[i * nodesPerElement + 5] = 4 * s;
      shS[i * nodesPerElement + 5] = 4 * r;
      shT[i * nodesPerElement + 5] = 0;

      // [6] = AC: 4*N3*N1 = 4*s*(1-r-s-t)
      shR[i * nodesPerElement + 6] = -4 * s;
      shS[i * nodesPerElement + 6] = -4 * r - 8 * s - 4 * t + 4;
      shT[i * nodesPerElement + 6] = -4 * s;

      // [7] = AD: 4*N1*N4 = 4*(1-r-s-t)*t
      shR[i * nodesPerElement + 7] = -4 * t;
      shS[i * nodesPerElement + 7] = -4 * t;
      shT[i * nodesPerElement + 7] = -4 * r - 4 * s - 8 * t + 4;

      // [8] = CD: 4*N3*N4 = 4*s*t  (Gmsh [8])
      shR[i * nodesPerElement + 8] = 0;
      shS[i * nodesPerElement + 8] = 4 * t;
      shT[i * nodesPerElement + 8] = 4 * s;

      // [9] = BD: 4*N2*N4 = 4*r*t  (Gmsh [9])
      shR[i * nodesPerElement + 9] = 4 * t;
      shS[i * nodesPerElement + 9] = 0;
      shT[i * nodesPerElement + 9] = 4 * r;
    }

    shX.resize(nodesPerElement, 0);
    shY.resize(nodesPerElement, 0);
    shZ.resize(nodesPerElement, 0);

    for (unsigned int i = 0; i < nodesPerElement; i++) {
      shX[i] = 0.;
      shY[i] = 0.;
      shZ[i] = 0.;
    }

  }

public:
  explicit TetShapeFunction(unsigned int elementOrder) : ShapeFunction(elementOrder)
  {
    // do separate initialisation to simplify use with openmp parallel for loops where separate instances of this
    // class are used
  }

  void setIntegrationPoints(const IntegrationPoints &integrationPoints) override {
    if (this->integrationPoints != nullptr ) {
      return; // already initialised
    }
    this->integrationPoints = &integrationPoints;

    //assert(integrationPoints.weights.size() == integrationPoints.points.size() / 3);
    numGaussPoints = integrationPoints.numGaussPoints();

    switch (elementOrder) {
      case 1:
        initialiseLinearTet();
        break;
      case 2:
        initialiseQuadraticTet();
        break;
      default:
        RUNTIME_ERROR("Unsupported element order " + std::to_string(elementOrder));
    }
  }


  [[nodiscard]] unsigned int getNumGaussPoints() const { return numGaussPoints; }

  /**
   * Prepare shape-function state for the current Gauss point within an element.
   *
   * Computes the element inverse Jacobian via the standard isoparametric summation over all
   * @c nodesPerElement nodes at the current quadrature location, then transforms reference-space
   * derivatives into global (x,y,z) derivatives.
   *
   * Although the Jacobian is mathematically constant for straight-sided elements, the quadratic
   * shape-function contributions produce slightly different floating-point sums at different Gauss
   * points (difference ≈ 2×10⁻¹⁵ relative).  Recomputing at every point avoids introducing a
   * systematic bias into the assembled stiffness matrix that would otherwise accumulate across
   * the mesh and degrade solution accuracy.
   *
   * @param nodes        Pointer to the nodal coordinates for the active element.
   *                     Must have at least @c nodesPerElement valid entries.
   * @param determinant  Pre-computed element determinant (scalar triple product of edge vectors,
   *                     equal to 6 × signed volume) in the same coordinate units as @p nodes.
   */
  void initialiseElement(Vec3 *nodes, double determinant) {
    double xr, xs, xt, yr, ys, yt, zr, zs, zt;
    xr = xs = xt = yr = ys = yt = zr = zs = zt = 0.0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      double x = nodes[i].x();
      double y = nodes[i].y();
      double z = nodes[i].z();

      double r = getShR(i);
      double s = getShS(i);
      double t = getShT(i);

      xr += x * r;
      xs += x * s;
      xt += x * t;
      yr += y * r;
      ys += y * s;
      yt += y * t;
      zr += z * r;
      zs += z * s;
      zt += z * t;
    }

    double Jinv[3][3] = {
              {(zt * ys - yt * zs) / determinant, (xt * zs - zt * xs) / determinant, (xs * yt - ys * xt) / determinant}
            , {(yt * zr - zt * yr) / determinant, (zt * xr - xt * zr) / determinant, (xt * yr - yt * xr) / determinant}
            , {(yr * zs - ys * zr) / determinant, (xs * zr - xr * zs) / determinant, (ys * xr - xs * yr) / determinant}
    };

    // x,y,z derivatives of shape functions
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      double r = getShR(i);
      double s = getShS(i);
      double t = getShT(i);
      shX[i] = r * Jinv[0][0] +
               s * Jinv[1][0] +
               t * Jinv[2][0];
      shY[i] = r * Jinv[0][1] +
               s * Jinv[1][1] +
               t * Jinv[2][1];
      shZ[i] = r * Jinv[0][2] +
               s * Jinv[1][2] +
               t * Jinv[2][2];
    }
  }

  [[nodiscard]] const double& Nx(int i) const { return shX[i]; }
  [[nodiscard]] const double& Ny(int i) const { return shY[i]; }
  [[nodiscard]] const double& Nz(int i) const { return shZ[i]; }

  [[nodiscard]] double sampleX(const double *values) const {
    double sum = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      sum += values[i] * Nx(i);
    }
    return sum;
  }

  [[nodiscard]] double sampleY(const double *values) const {
    double sum = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      sum += values[i] * Ny(i);
    }
    return sum;
  }

  [[nodiscard]] double sampleZ(const double *values) const {
    double sum = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      sum += values[i] * Nz(i);
    }
    return sum;
  }

  template<typename Src>
  void sampleQX(const Src &source, double &v1x, double &v2x, double &v3x, double &v4x, double &v5x) const {
    // gradient along x
    v1x = v2x = v3x = v4x = v5x = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      v1x += source[i][0] * Nx(i);
      v2x += source[i][1] * Nx(i);
      v3x += source[i][2] * Nx(i);
      v4x += source[i][3] * Nx(i);
      v5x += source[i][4] * Nx(i);
    }
  }

  template<typename Src>
  void sampleQY(const Src &source, double &v1y, double &v2y, double &v3y, double &v4y, double &v5y) const {
    // gradient along y
    v1y = v2y = v3y = v4y = v5y = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      v1y += source[i][0] * Ny(i);
      v2y += source[i][1] * Ny(i);
      v3y += source[i][2] * Ny(i);
      v4y += source[i][3] * Ny(i);
      v5y += source[i][4] * Ny(i);
    }
  }

  template<typename Src>
  void sampleQZ(const Src &source, double &v1z, double &v2z, double &v3z, double &v4z, double &v5z) const {
    // gradient along z
    v1z = v2z = v3z = v4z = v5z = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      v1z += source[i][0] * Nz(i);
      v2z += source[i][1] * Nz(i);
      v3z += source[i][2] * Nz(i);
      v4z += source[i][3] * Nz(i);
      v5z += source[i][4] * Nz(i);
    }
  }



};

#endif //PROJECT_QLC3D_GAUSSIAN_QUADRATURE_H
