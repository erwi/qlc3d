
//
//  STACK ALLOCATED 3D TETRAHEDRAL SHAPE FUNCTIONS
//


#include <vector>
#include "geom/vec3.h"

#ifndef SHAPEFUNCTION3D_H
#define SHAPEFUNCTION3D_H

// Table 10.4 Quadrature for unit tetrahedra in Zienkiewicz
#define NGP4 11
#define W4_11 -74.0/5625.0
#define W4_12 343.0/45000.0
#define W4_13 56.0/2250.0

#define A4  (1+sqrt(5.0/14.0))/4.0
#define B4  (1-sqrt(5.0/14.0))/4.0
    const double W4[NGP4] = {W4_11, W4_12, W4_12, W4_12, W4_12, W4_13, W4_13, W4_13, W4_13, W4_13, W4_13};
    const double gp[NGP4][4]={	{0.25	  , 0.25	,	0.25	,0.25}
                            ,{11.0/14.0     ,	1.0/14.0	,	1.0/14.0	,1.0/14.0}
                            ,{1.0/14.0      ,	11.0/14.0	,	1.0/14.0	,1.0/14.0}
                            ,{1.0/14.0	  ,	1.0/14.0	,	11.0/14.0   ,1.0/14.0},
                            {1.0/14.0  , 1.0/14.0	,	1.0/14.0	,11.0/14.0},
                            {A4	  , A4	,  B4 , B4},
                            {A4	  , B4	,  A4 , B4},
                            {A4   , B4  ,  B4 , A4},
                            {B4	  , A4  ,  A4 , B4},
                            {B4	  , A4  ,  B4 , A4},
                            {B4	  , B4  ,  A4 , A4}};

class Shape4thOrder {
    public:
    const unsigned int ngp;
    double w[NGP4];
    double sh1[NGP4][4];
    double sh1r[NGP4][4];
    double sh1s[NGP4][4];
    double sh1t[NGP4][4];



    Shape4thOrder(): ngp(NGP4) {
        for (unsigned int i = 0 ; i < ngp ; i++) {
            w[i] = W4[i];
            // P1 Shape functions
            sh1[i][0]=1-gp[i][0]-gp[i][1]-gp[i][2];
            sh1[i][1]=gp[i][0];
            sh1[i][2]=gp[i][1];
            sh1[i][3]=gp[i][2];
            // P1 Shape functions r-derivatives
            sh1r[i][0]=-1.0;
            sh1r[i][1]=1.0;
            sh1r[i][2]=0.0;
            sh1r[i][3]=0.0;
            // P1 Shape functions s-derivatives
            sh1s[i][0]=-1.0;
            sh1s[i][1]=0.0;
            sh1s[i][2]=1.0;
            sh1s[i][3]=0.0;
            // P1 Shape functions t-derivatives
            sh1t[i][0]=-1.0;
            sh1t[i][1]=0.0;
            sh1t[i][2]=0.0;
            sh1t[i][3]=1.0;
        }
    }
};

class ShapeFunction {
  unsigned int currentPoint = 0;
  std::vector<double*> sh1;
  std::vector<double*> sh1r;
  std::vector<double*> sh1s;
  std::vector<double*> sh1t;

  std::vector<double> shape;
  std::vector<double> shapeX;
  std::vector<double> shapeY;
  std::vector<double> shapeZ;

public:
ShapeFunction() {
    for (unsigned int i = 0 ; i < NGP4 ; i++) {
      double r = gp[i][0];
      double s = gp[i][1];
      double t = gp[i][2];

      sh1.push_back(new double[4]);
      sh1[i][0] = 1 - r - s - t;
      sh1[i][1] = r;
      sh1[i][2] = s;
      sh1[i][3] = t;

      sh1r.push_back(new double[4]);
      sh1r[i][0] = -1.0;
      sh1r[i][1] = 1.0;
      sh1r[i][2] = 0.0;
      sh1r[i][3] = 0.0;

      sh1s.push_back(new double[4]);
      sh1s[i][0] = -1.0;
      sh1s[i][1] = 0.0;
      sh1s[i][2] = 1.0;
      sh1s[i][3] = 0.0;

      sh1t.push_back(new double[4]);
      sh1t[i][0] = -1.0;
      sh1t[i][1] = 0.0;
      sh1t[i][2] = 0.0;
      sh1t[i][3] = 1.0;
    }

    shape = {0, 0, 0, 0};
    shapeX = {0, 0, 0, 0};
    shapeY = {0, 0, 0, 0};
    shapeZ = {0, 0, 0, 0};
  }

  [[nodiscard]] unsigned int numGaussPoints() const { return NGP4; }

  void initialiseElement(Vec3 *nodes, double determinant) {
    currentPoint = 0;

    double xr, xs, xt, yr, ys, yt, zr, zs, zt;
    xr = xs = xt = yr = ys = yt = zr = zs = zt = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      double x = nodes[i].x();
      double y = nodes[i].y();
      double z = nodes[i].z();
      xr += x * sh1r[0][i];
      xs += x * sh1s[0][i];
      xt += x * sh1t[0][i];
      yr += y * sh1r[0][i];
      ys += y * sh1s[0][i];
      yt += y * sh1t[0][i];
      zr += z * sh1r[0][i];
      zs += z * sh1s[0][i];
      zt += z * sh1t[0][i];
    }

    double Jinv[3][3] = {
            {(zt * ys - yt * zs) / determinant, (xt * zs - zt * xs) / determinant, (xs * yt - ys * xt) / determinant}
            , {(yt * zr - zt * yr) / determinant, (zt * xr - xt * zr) / determinant, (xt * yr - yt * xr) / determinant}
            , {(yr * zs - ys * zr) / determinant, (xs * zr - xr * zs) / determinant, (ys * xr - xs * yr) / determinant}
    };

    // x,y,z derivatives of shape functions
    for (int i = 0; i < 4; i++) {
      shapeX[i] = sh1r[0][i] * Jinv[0][0] +
                  sh1s[0][i] * Jinv[1][0] +
                  sh1t[0][i] * Jinv[2][0];
      shapeY[i] = sh1r[0][i] * Jinv[0][1] +
                  sh1s[0][i] * Jinv[1][1] +
                  sh1t[0][i] * Jinv[2][1];
      shapeZ[i] = sh1r[0][i] * Jinv[0][2] +
                  sh1s[0][i] * Jinv[1][2] +
                  sh1t[0][i] * Jinv[2][2];
    }//end for i
  }

  double N(int i) const {
    return sh1[currentPoint][i];
  }

  double Nx(int i) const {
    return shapeX[i];
  }

  double Ny(int i) const {
    return shapeY[i];
  }

  double Nz(int i) const {
    return shapeZ[i];
  }

  double getWeight() const {
    return W4[currentPoint];
  }

  void nextPoint() {
    currentPoint++;
  }

  bool hasNextPoint() const {
    return currentPoint < NGP4;
  }

  double calculate(const double *values) const {
    return values[0] * sh1[currentPoint][0] +
           values[1] * sh1[currentPoint][1] +
           values[2] * sh1[currentPoint][2] +
           values[3] * sh1[currentPoint][3];
  }

  double gradientX(const double *values) const {
    return values[0] * shapeX[0] +
           values[1] * shapeX[1] +
           values[2] * shapeX[2] +
           values[3] * shapeX[3];
  }

  double gradientY(const double *values) const {
    return values[0] * shapeY[0] +
           values[1] * shapeY[1] +
           values[2] * shapeY[2] +
           values[3] * shapeY[3];
  }

  double gradientZ(const double *values) const {
    return values[0] * shapeZ[0] +
           values[1] * shapeZ[1] +
           values[2] * shapeZ[2] +
           values[3] * shapeZ[3];
  }
};


#endif // SHAPEFUNCTION3D_H
#undef NGP4
#undef W4_11
#undef W4_12
#undef W4_13
#undef A4
#undef B4
