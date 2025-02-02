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
  [[nodiscard]] const double& N(int i) const { return get(sh, i); }

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

    sh.resize(numGaussPoints * nodesPerElement, 0);
    shR.resize(numGaussPoints * nodesPerElement, 0);
    shS.resize(numGaussPoints * nodesPerElement, 0);
    shT.resize(numGaussPoints * nodesPerElement, 0);

    for (unsigned int i = 0; i < numGaussPoints; ++i) {
      double r = integrationPoints->points[i * 3 + 0];
      double s = integrationPoints->points[i * 3 + 1];
      double t = integrationPoints->points[i * 3 + 2];

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
    nodesPerElement = 10;

    sh.resize(numGaussPoints * nodesPerElement, 0);
    shR.resize(numGaussPoints * nodesPerElement, 0);
    shS.resize(numGaussPoints * nodesPerElement, 0);
    shT.resize(numGaussPoints * nodesPerElement, 0);

    for (unsigned int i = 0; i < numGaussPoints; ++i) {
      double r = integrationPoints->points[i * 3 + 0];
      double s = integrationPoints->points[i * 3 + 1];
      double t = integrationPoints->points[i * 3 + 2];

      // corner nodes expressed in natural coordinates
      double N1 = 1 - r - s - t;
      double N2 = r;
      double N3 = s;
      double N4 = t;

      sh[i * nodesPerElement + 0] = N1 * (2 * N1 - 1);
      sh[i * nodesPerElement + 1] = N2 * (2 * N2 - 1);
      sh[i * nodesPerElement + 2] = N3 * (2 * N3 - 1);
      sh[i * nodesPerElement + 3] = N4 * (2 * N4 - 1);

      // mid-edge nodes
      // TODO: probably we must swap some to match node ordering of GMSH
      sh[i * nodesPerElement + 4] = 4 * N1 * N2;
      sh[i * nodesPerElement + 5] = 4 * N2 * N3;
      sh[i * nodesPerElement + 6] = 4 * N3 * N1;
      sh[i * nodesPerElement + 7] = 4 * N1 * N4;
      sh[i * nodesPerElement + 8] = 4 * N2 * N4;
      sh[i * nodesPerElement + 9] = 4 * N3 * N4;


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

      shR[i * nodesPerElement + 4] = -8 * r - 4 * s - 4 * t + 4;
      shS[i * nodesPerElement + 4] = -4 * r;
      shT[i * nodesPerElement + 4] = -4 * r;

      shR[i * nodesPerElement + 5] = 4 * s;
      shS[i * nodesPerElement + 5] = 4 * r;
      shT[i * nodesPerElement + 5] = 0;

      shR[i * nodesPerElement + 6] = -4 * s;
      shS[i * nodesPerElement + 6] = -4 * r - 8 * s - 4 * t + 4;
      shT[i * nodesPerElement + 6] = -4 * s;

      shR[i * nodesPerElement + 7] = -4 * t;
      shS[i * nodesPerElement + 7] = -4 * t;
      shT[i * nodesPerElement + 7] = -4 * r - 4 * s - 8 * t + 4;

      shR[i * nodesPerElement + 8] = 4 * t;
      shS[i * nodesPerElement + 8] = 0;
      shT[i * nodesPerElement + 8] = 4 * r;

      shR[i * nodesPerElement + 9] = 0;
      shS[i * nodesPerElement + 9] = 4 * t;
      shT[i * nodesPerElement + 9] = 4 * s;
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

    assert(integrationPoints.weights.size() == integrationPoints.points.size() / 3);
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

  [[nodiscard]] double sample(const double *value) const {
    double sum = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      sum += value[i] * N(i);
    }
    return sum;
  }

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

  template<typename Src>
  void sampleQ(const Src* source, double &q1, double &q2, double &q3, double &q4, double &q5) const {
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
};

class BoundaryIntegralShapeFunction : public TetShapeFunction {

public:
  BoundaryIntegralShapeFunction(unsigned int elementOrder) : TetShapeFunction(elementOrder) {
    currentPoint = 0;
  }

  void setIntegrationPoints(const IntegrationPoints &integrationPoints) override {

    this->numGaussPoints = integrationPoints.numGaussPoints();
    this->integrationPoints = &integrationPoints;
    // expect 2D integration points since this is a boundary integral
    assert(integrationPoints.points.size() / numGaussPoints == 2);

    nodesPerElement = 4;

    sh.resize(numGaussPoints * nodesPerElement, 0);
    shR.resize(numGaussPoints * nodesPerElement, 0);
    shS.resize(numGaussPoints * nodesPerElement, 0);
    shT.resize(numGaussPoints * nodesPerElement, 0);

    for (unsigned int i = 0; i < numGaussPoints; ++i) {
      double r = integrationPoints.points[i * 3 + 0];
      double s = integrationPoints.points[i * 3 + 1];
      double t = 0; // the t-node is the internal node in the element, which is always 0 since integration is on the boundary

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

      shT[i * nodesPerElement + 0] = 0.0; //-1.0;
      shT[i * nodesPerElement + 1] = 0.0;
      shT[i * nodesPerElement + 2] = 0.0;
      shT[i * nodesPerElement + 3] = 0.0; //1.0;
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

};

/*
inline TetShapeFunction createLinearTetShapeFunction() {
  // Weights and Gauss points from table 10.4 in
  // J. E. Akin Finite Element Analysis with Error Estimators
  const double w1 = -74.0/5625.0;
  const double w2 = 343.0/45000.0;
  const double w3 = 56.0/2250.0;
  std::vector<double> weights = {w1, w2, w2, w2, w2, w3, w3, w3, w3, w3, w3};
  const double a = (1 + sqrt(5.0 / 14.0)) / 4.0;
  const double b = (1 - sqrt(5.0 / 14.0)) / 4.0;

  std::vector<double> r = {
          0.25,
          11.0 / 14.0, 1.0 / 14.0, 1.0 / 14.0, 1.0 / 14.0,
          a, a, a, b, b, b};

  std::vector<double> s = {
          0.25,
          1.0 / 14.0, 11.0 / 14.0, 1.0 / 14.0, 1.0 / 14.0,
          a, b, b, a, a, b};

  std::vector <double> t = {
          0.25,
          1.0 / 14.0, 1.0 / 14.0, 11.0 / 14.0, 1.0 / 14.0,
          b, a, b, a, b, a};

  return TetShapeFunction(4, weights, r, s, t);
}
 */

template<unsigned int NGP>
class GaussianQuadratureTet {
  static const unsigned int NPE = 4;
  double w[NGP];

  double shX[NGP];
  double shY[NGP];
  double shZ[NGP];

  double sh[NGP][NPE];
  double shR[NGP][NPE];
  double shS[NGP][NPE];
  double shT[NGP][NPE];

  unsigned int currentPoint = 0;
public:
  GaussianQuadratureTet(const double* weights,
                        const double* r,
                        const double* s,
                        const double* t) {

    for (unsigned int i = 0; i < NGP; ++i) {
      w[i] = weights[i];

      sh[i][0] = 1 - r[i] - s[i] - t[i];
      sh[i][1] = r[i];
      sh[i][2] = s[i];
      sh[i][3] = t[i];

      shR[i][0] = -1.0;
      shR[i][1] = 1.0;
      shR[i][2] = 0.0;
      shR[i][3] = 0.0;

      shS[i][0] = -1.0;
      shS[i][1] = 0.0;
      shS[i][2] = 1.0;
      shS[i][3] = 0.0;

      shT[i][0] = -1.0;
      shT[i][1] = 0.0;
      shT[i][2] = 0.0;
      shT[i][3] = 1.0;
    }

    for (unsigned int i = 0; i < NPE; i++) {
      shX[i] = 0.;
      shY[i] = 0.;
      shZ[i] = 0.;
    }
  }

  [[nodiscard]] unsigned int numGaussPoints() const { return NGP; }

  void initialiseElement(Vec3 *nodes, double determinant) {
    currentPoint = 0;

    double xr, xs, xt, yr, ys, yt, zr, zs, zt;
    xr = xs = xt = yr = ys = yt = zr = zs = zt = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      double x = nodes[i].x();
      double y = nodes[i].y();
      double z = nodes[i].z();
      xr += x * shR[0][i];
      xs += x * shS[0][i];
      xt += x * shT[0][i];
      yr += y * shR[0][i];
      ys += y * shS[0][i];
      yt += y * shT[0][i];
      zr += z * shR[0][i];
      zs += z * shS[0][i];
      zt += z * shT[0][i];
    }

    double Jinv[3][3] = {
              {(zt * ys - yt * zs) / determinant, (xt * zs - zt * xs) / determinant, (xs * yt - ys * xt) / determinant}
            , {(yt * zr - zt * yr) / determinant, (zt * xr - xt * zr) / determinant, (xt * yr - yt * xr) / determinant}
            , {(yr * zs - ys * zr) / determinant, (xs * zr - xr * zs) / determinant, (ys * xr - xs * yr) / determinant}
    };

    // x,y,z derivatives of shape functions
    for (int i = 0; i < 4; i++) {
      shX[i] = shR[0][i] * Jinv[0][0] +
               shS[0][i] * Jinv[1][0] +
               shT[0][i] * Jinv[2][0];
      shY[i] = shR[0][i] * Jinv[0][1] +
               shS[0][i] * Jinv[1][1] +
               shT[0][i] * Jinv[2][1];
      shZ[i] = shR[0][i] * Jinv[0][2] +
               shS[0][i] * Jinv[1][2] +
               shT[0][i] * Jinv[2][2];
    }
  }

  [[nodiscard]] double weight() const { return w[currentPoint]; }
  [[nodiscard]] double N(int i) const { return sh[currentPoint][i]; }
  [[nodiscard]] double Nx(int i) const { return shX[i]; }
  [[nodiscard]] double Ny(int i) const { return shY[i]; }
  [[nodiscard]] double Nz(int i) const { return shZ[i]; }
  [[nodiscard]] bool hasNextPoint() const { return currentPoint < NGP; }

  void nextPoint() { currentPoint++; }

  [[nodiscard]] double sample(const double *values) const {
    return values[0] * N(0) +
           values[1] * N(1) +
           values[2] * N(2) +
           values[3] * N(3);
  }

  [[nodiscard]] double sampleX(const double *values) const {
    return values[0] * Nx(0) + values[1] * Nx(1) + values[2] * Nx(2) + values[3] * Nx(3);
  }

  [[nodiscard]] double sampleY(const double *values) const {
    return values[0] * Ny(0) + values[1] * Ny(1) + values[2] * Ny(2) + values[3] * Ny(3);
  }

  [[nodiscard]] double sampleZ(const double *values) const {
    return values[0] * Nz(0) + values[1] * Nz(1) + values[2] * Nz(2) + values[3] * Nz(3);
  }

  template<typename Src>
  void sampleAll(const Src* source, double &v1, double &v2, double &v3, double &v4, double &v5) const {
    v1 = v2 = v3 = v4 = v5 = 0;
    for (unsigned int i = 0; i < NPE; i++) {
      v1 += source[i][0] * N(i);
      v2 += source[i][1] * N(i);
      v3 += source[i][2] * N(i);
      v4 += source[i][3] * N(i);
      v5 += source[i][4] * N(i);
    }
  }

  template<typename Src>
  void sampleQ(const Src &source, double &v1, double &v2, double &v3, double &v4, double &v5) const {
    v1 = v2 = v3 = v4 = v5 = 0;
    const unsigned int nodesPerElement = source.size();
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      v1 += source[i][0] * N(i);
      v2 += source[i][1] * N(i);
      v3 += source[i][2] * N(i);
      v4 += source[i][3] * N(i);
      v5 += source[i][4] * N(i);
    }
  }

  /** Deprecated */
  template<typename Src>
  void sampleAllX(const Src* source, double &v1x, double &v2x, double &v3x, double &v4x, double &v5x) const {
    // gradient along x
    v1x = v2x = v3x = v4x = v5x = 0;
    for (unsigned int i = 0; i < NPE; i++) {
      v1x += source[i][0] * Nx(i);
      v2x += source[i][1] * Nx(i);
      v3x += source[i][2] * Nx(i);
      v4x += source[i][3] * Nx(i);
      v5x += source[i][4] * Nx(i);
    }
  }

  /** Deprecated */
  template<typename Src>
  void sampleAllY(const Src* source, double &v1y, double &v2y, double &v3y, double &v4y, double &v5y) const {
    // gradient along y
    v1y = v2y = v3y = v4y = v5y = 0;
    for (unsigned int i = 0; i < NPE; i++) {
      v1y += source[i][0] * Ny(i);
      v2y += source[i][1] * Ny(i);
      v3y += source[i][2] * Ny(i);
      v4y += source[i][3] * Ny(i);
      v5y += source[i][4] * Ny(i);
    }
  }

  /** Deprecated */
  template<typename Src>
  void sampleAllZ(const Src* source, double &v1z, double &v2z, double &v3z, double &v4z, double &v5z) const {
    // gradient along z
    v1z = v2z = v3z = v4z = v5z = 0;
    for (unsigned int i = 0; i < NPE; i++) {
      v1z += source[i][0] * Nz(i);
      v2z += source[i][1] * Nz(i);
      v3z += source[i][2] * Nz(i);
      v4z += source[i][3] * Nz(i);
      v5z += source[i][4] * Nz(i);
    }
  }

  template<typename Src>
  void sampleQX(const Src &source, double &v1x, double &v2x, double &v3x, double &v4x, double &v5x) const {
    // gradient along x
    const unsigned int nodesPerElement = source.size();
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
    const unsigned int nodesPerElement = source.size();
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
    const unsigned int nodesPerElement = source.size();
    v1z = v2z = v3z = v4z = v5z = 0;
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      v1z += source[i][0] * Nz(i);
      v2z += source[i][1] * Nz(i);
      v3z += source[i][2] * Nz(i);
      v4z += source[i][3] * Nz(i);
      v5z += source[i][4] * Nz(i);
    }
  }

  template<typename Src>
  void sampleAll(const Src* source, double &v1, double &v2, double &v3, double &v4, double &v5, double &v6) const {
    v1 = v2 = v3 = v4 = v5 = v6 = 0;
    for (unsigned int i = 0; i < NPE; i++) {
      v1 += source[i][0] * N(i);
      v2 += source[i][1] * N(i);
      v3 += source[i][2] * N(i);
      v4 += source[i][3] * N(i);
      v5 += source[i][4] * N(i);
      v6 += source[i][5] * N(i);
    }
  }

  template<typename Src>
  void samplePermittivity(const Src &source, double &v1, double &v2, double &v3, double &v4, double &v5, double &v6) const {
    v1 = v2 = v3 = v4 = v5 = v6 = 0;
    const unsigned int nodesPerElement = source.size();
    for (unsigned int i = 0; i < nodesPerElement; i++) {
      v1 += source[i][0] * N(i);
      v2 += source[i][1] * N(i);
      v3 += source[i][2] * N(i);
      v4 += source[i][3] * N(i);
      v5 += source[i][4] * N(i);
      v6 += source[i][5] * N(i);
    }
  }

};

template<unsigned int NGP>
class GaussianQuadratureTri {
  static const unsigned int NPE = 3;
  double w[NGP];
  double gp[NGP][NPE];

  double shX[NGP];
  double shY[NGP];

  double sh[NGP][NPE];
  double shR[NGP][NPE];
  double shS[NGP][NPE];

  unsigned int currentPoint = 0;
public:
  GaussianQuadratureTri(const double* weights,
                        const double* r,
                        const double* s) {

    for (unsigned int i = 0; i < NGP; ++i) {
      w[i] = weights[i];

      sh[i][0] = 1 - r[i] - s[i];
      sh[i][1] = r[i];
      sh[i][2] = s[i];

      shR[i][0] = -1.0;
      shR[i][1] = 1.0;
      shR[i][2] = 0.0;

      shS[i][0] = -1.0;
      shS[i][1] = 0.0;
      shS[i][2] = 1.0;
    }

    for (unsigned int i = 0; i < NPE; i++) {
      shX[i] = 0.;
      shY[i] = 0.;
    }
  }

  [[nodiscard]] double weight(unsigned int i) const { return w[i]; }
  [[nodiscard]] double gaussPoint(unsigned int i, unsigned int j) const { return gp[i][j]; }
  [[nodiscard]] unsigned int numGaussPoints() const { return NGP; }
  [[nodiscard]] unsigned int numPointsPerElement() const { return NPE; }

  void initialiseElement(Vec3 *nodes, double determinant) {
    currentPoint = 0;
  }

  [[nodiscard]] double weight() const { return w[currentPoint]; }
  [[nodiscard]] double gaussPoint(int i) const { return gp[currentPoint][i]; }
  [[nodiscard]] double N(int i) const { return sh[currentPoint][i]; }
  //[[nodiscard]] double Nx(int i) const { return shX[i]; }
  //[[nodiscard]] double Ny(int i) const { return shY[i]; }
  //[[nodiscard]] double Nz(int i) const { return shZ[i]; }
  [[nodiscard]] bool hasNextPoint() const { return currentPoint < NGP; }

  void nextPoint() { currentPoint++; }

  [[nodiscard]] double sample(const double *values) const {
    return values[0] * N(0) +
           values[1] * N(1) +
           values[2] * N(2) +
           values[3] * N(3);
  }

  template<typename Src>
  void sampleAll(const Src* source, double &v1, double &v2, double &v3, double &v4, double &v5) const {
    v1 = v2 = v3 = v4 = v5 = 0;
    for (unsigned int i = 0; i < NPE; i++) {
      v1 += source[i][0] * N(i);
      v2 += source[i][1] * N(i);
      v3 += source[i][2] * N(i);
      v4 += source[i][3] * N(i);
      v5 += source[i][4] * N(i);
    }
  }

  /*
  void sample(const Vec3 &source, Vec3 &destination) const {
    destination.set(0, 0, 0);
    for (unsigned int i = 0; i < NPE; i++) {
      destination.add(source.x() * N(i),
                      source.y() * N(i),
                      source.z() * N(i));
    }
  }
   */

  void sample(const Vec3 source[3], Vec3 &destination) const {
    destination.set(0, 0, 0);
    for (unsigned int i = 0; i < NPE; i++) {
      destination.add(source[i].x() * N(i),
                      source[i].y() * N(i),
                      source[i].z() * N(i));
    }
  }
};


/** Create 4th order linear tet gaussian quadrature integrator */
inline GaussianQuadratureTet<11> gaussQuadratureTet4thOrder() {
  // Weights and Gauss points from table 10.4 in
  // J. E. Akin Finite Element Analysis with Error Estimators
  const double w1 = -74.0/5625.0;
  const double w2 = 343.0/45000.0;
  const double w3 = 56.0/2250.0;
  const double weights[11] = {w1,w2, w2, w2, w2,w3, w3, w3, w3, w3, w3 };
  const double a = (1 + sqrt(5.0 / 14.0)) / 4.0;
  const double b = (1 - sqrt(5.0 / 14.0)) / 4.0;

  const double r[11] = {
          0.25,
          11.0 / 14.0, 1.0 / 14.0, 1.0 / 14.0, 1.0 / 14.0,
          a, a, a, b, b, b};

  const double s[11] = {
          0.25,
          1.0 / 14.0, 11.0 / 14.0, 1.0 / 14.0, 1.0 / 14.0,
          a, b, b, a, a, b};

  const double t[11] = {
          0.25,
          1.0 / 14.0, 1.0 / 14.0, 11.0 / 14.0, 1.0 / 14.0,
          b, a, b, a, b, a};

  return GaussianQuadratureTet<11>(weights, r, s, t);
}

inline GaussianQuadratureTet<7> gaussQuadratureTetBoundaryIntegral4thOrder() {
  // Weights and Gauss points from table 10.3 in J. E. Akin
  // Note that these are triangle element quadrature points, used on the boundary of the tet element
  // so the t coordinate is always 0
  const double w[7] = {1. / 40, 1./ 15, 1. / 40, 1. / 15, 1. / 40, 1. / 15, 9. / 40};
  const double r[7] = {0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 1./3.};
  const double s[7] = {0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 1./3.};
  const double t[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // all-0 because this is integral along the t = 0 surface

  return GaussianQuadratureTet<7>(w, r, s, t);
}

inline GaussianQuadratureTri<7> gaussianQuadratureTri4thOrder() {
  // Weights and Gauss points from table 10.3 in J. E. Akin
  // TODO: do we need 4th order for the triangles? The table also has
  // 3rd, 2nd and 1st order quadrature points
  const double w[7] = {1. / 40, 1./ 15, 1. / 40, 1. / 15, 1. / 40, 1. / 15, 9. / 40};
  const double r[7] = {0.0, 0.5, 1.0, 0.5, 0.0, 0.0, 1./3.};
  const double s[7] = {0.0, 0.0, 0.0, 0.5, 1.0, 0.5, 1./3.};
  return GaussianQuadratureTri<7>(w, r, s);
}



#endif //PROJECT_QLC3D_GAUSSIAN_QUADRATURE_H
