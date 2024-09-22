#ifndef PROJECT_QLC3D_GAUSSIAN_QUADRATURE_H
#define PROJECT_QLC3D_GAUSSIAN_QUADRATURE_H
#include <geom/vec3.h>
#include <lc-representation.h>

template<unsigned int NGP>
class GaussianQuadratureTet {
  static const unsigned int NPE = 4;
  double w[NGP];
  double gp[NGP][NPE];

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

  [[nodiscard]] double weight(unsigned int i) const { return w[i]; }
  [[nodiscard]] double gaussPoint(unsigned int i, unsigned int j) const { return gp[i][j]; }
  [[nodiscard]] unsigned int numGaussPoints() const { return NGP; }
  [[nodiscard]] unsigned int numPointsPerElement() const { return NPE; }

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
    }//end for i
  }

  [[nodiscard]] double weight() const { return w[currentPoint]; }
  [[nodiscard]] double gaussPoint(int i) const { return gp[currentPoint][i]; }
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

  /*
  [[nodiscard]] double sampleX(const double *values) const {
    return values[0] * Nx(0) + values[1] * Nx(1) + values[2] * Nx(2) + values[3] * Nx(3);
  }

  //[[nodiscard]] double sampleY(const double *values) const {
    return values[0] * Ny(0) + values[1] * Ny(1) + values[2] * Ny(2) + values[3] * Ny(3);
  }

  //[[nodiscard]] double sampleZ(const double *values) const {
    return values[0] * Nz(0) + values[1] * Nz(1) + values[2] * Nz(2) + values[3] * Nz(3);
  }
  */

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
   */
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
