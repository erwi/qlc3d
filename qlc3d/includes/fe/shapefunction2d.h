#ifndef SHAPEFUNCTION2D_H
#define SHAPEFUNCTION2D_H


#define NGPS4 7
namespace SurfaceShapes
{
  /*
    const double gps4[NGPS4][2]= {
        {0.8168476, 0.09157621},
        {0.09157621,0.8168476},
        {0.09157621,0.09157621},
        {0.1081030, 0.4459485},
        {0.4459485, 0.1081030},
        {0.4459485, 0.4459485}};
        */

  // From table 10.3 in Zienkiewicz
  const double gps4[NGPS4][2] = {
          {0.0, 0.0},
          {0.5, 0.0},
          {1.0, 0.0},
          {0.5, 0.5},
          {0.0, 1.0},
          {0.0, 0.5},
          {1./3., 1./3.}
  };

}
class ShapeSurf4thOrder{
  const unsigned int ngps;
  double w[NGPS4] = {1. / 40, 1./ 15, 1. / 40, 1. / 15, 1. / 40, 1. / 15, 9. / 40};
  double sh1[NGPS4][4];
  double sh1r[NGPS4][4];
  double sh1s[NGPS4][4];
  double sh1t[NGPS4][4];

  double shape[4] = {0, 0, 0, 0};
  double shapeX[4] = {0, 0, 0, 0};
  double shapeY[4] = {0, 0, 0, 0};
  double shapeZ[4] = {0, 0, 0, 0};

  unsigned int currentPoint = 0;

public:

  ShapeSurf4thOrder() : ngps(NGPS4) {
    /*
    w[0] = 0.05497587;
    w[1] = 0.05497587;
    w[2] = 0.05497587;
    w[3] = 0.1116908;
    w[4] = 0.1116908;
    w[5] = 0.1116908;
     */
    for(unsigned int i = 0 ; i < NGPS4 ; i++) {
      double r = SurfaceShapes::gps4[i][0];
      double s = SurfaceShapes::gps4[i][1];
      double t = 0; // always 0 because this is integral along the t = 0 surface
      // P1 Shape functions
      sh1[i][0] = 1 - r - s - t;
      sh1[i][1] = r;
      sh1[i][2] = s;
      sh1[i][3] = t;
      // P1 Shape functions r-derivatives
      sh1r[i][0]=-1;
      sh1r[i][1]=1;
      sh1r[i][2]=0;
      sh1r[i][3]=0;
      // P1 Shape functions s-derivatives
      sh1s[i][0]=-1;
      sh1s[i][1]=0;
      sh1s[i][2]=1;
      sh1s[i][3]=0;
      //P1 Shape functions t-derivatives
      sh1t[i][0]=-1;
      sh1t[i][1]=0;
      sh1t[i][2]=0;
      sh1t[i][3]=1;
    }
  }

  [[ndiscard]] unsigned int numGaussPoints() const { return NGPS4; }

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
    return w[currentPoint];
  }

  void nextPoint() {
    currentPoint++;
  }

  bool hasNextPoint() const {
    return currentPoint < NGPS4;
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


#endif // SHAPEFUNCTION2D_H
