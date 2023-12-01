#ifndef PROJECT_QLC3D_LC_ENERGY_TERMS_H
#define PROJECT_QLC3D_LC_ENERGY_TERMS_H
#include <fe/gaussian-quadrature.h>

namespace LcEnergyTerms {
  const double D2 = 1. / 2;
  const double D3 = 1. / 3;
  const double D4 = 1. / 4;
  const double D6 = 1. / 6;
  const double rt2 = std::sqrt(2.0);
  const double rt3 = std::sqrt(3.0);
  const double rt6 = std::sqrt(6.0);
  const double eps0 = 8.854187817e-12;

  inline void assembleThermotropic(double lK[20][20],
                                   double lL[20],
                                   const GaussianQuadratureTet<11> &shapes,
                                   const double &tetDeterminant,
                                   const double &q1,
                                   const double &q2,
                                   const double &q3,
                                   const double &q4,
                                   const double &q5,
                                   const double &A,
                                   const double &B,
                                   const double &C) {

    const double R = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4 + q5 * q5;
    const double mul = shapes.weight() * tetDeterminant;
    for (int i = 0; i < 4; i++) {
      const double ShR = shapes.N(i) * mul;
      // RHS terms
      lL[i + 0] += (A * q1 + D3 * 0.5 * B * (q5 * q5 * rt6 * 0.5 - q3 * q3 * rt6 - rt6 * q2 * q2 + q1 * q1 * rt6 + q4 * q4 * rt6 * 0.5) + C * R * q1) * ShR;
      lL[i + 4] += (A * q2 + D3 * B * (0.75 * q5 * q5 * rt2 - q1 * rt6 * q2 - 0.75 * q4 * q4 * rt2) + C * R * q2) * ShR;
      lL[i + 8] += (A * q3 + D3 * B * (-q3 * q1 * rt6 + 1.5 * rt2 * q5 * q4) + C * R * q3) * mul * shapes.N(i);
      lL[i + 12] += (A * q4 + D3 * 0.5 * B * (3.0 * q3 * rt2 * q5 + q4 * q1 * rt6 - 3.0 * q4 * q2 * rt2) + C * R * q4) * ShR;
      lL[1 + 16] += (A * q5 + D3 * 0.5 * B * (q5 * q1 * rt6 + 3.0 * q5 * q2 * rt2 + 3.0 * q3 * rt2 * q4) + C * R * q5) * ShR;

      for (int j = 0; j < 4; j++) {
        const double ShRC = shapes.N(i) * shapes.N(j) * mul;
        // Matrix diagonal terms 11, 22, 33, 44, 55
        lK[i+ 0][i+ 0] += (A + D3 * B * q1 * rt6 + 2.0 * C * q1 * q1 + C * R) * ShRC;
        lK[i+ 4][i+ 4] += (A - D3 * B * q1 * rt6 + 2.0 * C * q2 * q2 + C * R) * ShRC;
        lK[i+ 8][i+ 8] += (A - B * q1 * rt6 * D3 + 2.0 * C * q3 * q3 + C * R) * ShRC;
        lK[i+12][i+12] += (A + B * (q1 * rt6 - 3.0 * q2 * rt2) * D6 + 2.0 * C * q4 * q4 + C * R) * ShRC;
        lK[i+16][i+16] += (A + B * (q1 * rt6 + 3.0 * q2 * rt2) * D6 + 2.0 * C * q5 * q5 + C * R) * ShRC;

        // Matrix off-diagonal terms
        // 12
        lK[i+0][j+4] += (-D3 * B * rt6 * q2 + 2.0 * C * q2 * q1) * ShRC;
        lK[i+4][j+0] += (-D3 * B * rt6 * q2 + 2.0 * C * q2 * q1) * ShRC;
        // 13
        lK[i+0][j+8] += (-B * q3 * rt6 * D3 + 2.0 * C * q3 * q1) * ShRC;
        lK[i+8][j+0] += (-B * q3 * rt6 * D3 + 2.0 * C * q3 * q1) * ShRC;
        // 14
        lK[i+0 ][j+12] += (B * q4 * rt6 * D6 + 2.0 * C * q4 * q1) * ShRC;
        lK[i+12][j+0 ] += (B * q4 * rt6 * D6 + 2.0 * C * q4 * q1) * ShRC;
        // 15
        lK[i+0 ][j+16] += (B * q5 * rt6 * D6 + 2.0 * C * q5 * q1) * ShRC;
        lK[i+16][j+ 0] += (B * q5 * rt6 * D6 + 2.0 * C * q5 * q1) * ShRC;

        // 23
        lK[i +4][j + 8] += (2.0 * C * q3 * q2) * ShRC;
        lK[i +8][j + 4] += (2.0 * C * q3 * q2) * ShRC;
        // 24
        lK[i +4][j +12] += (-B * q4 * rt2 * D2 + 2.0 * C * q4 * q2) * ShRC;
        lK[i+12][j + 4] += (-B * q4 * rt2 * D2 + 2.0 * C * q4 * q2) * ShRC;
        // 25
        lK[i+4 ][j +16] += (B * q5 * rt2 * D2  + 2.0 * C * q5 * q2) * ShRC;
        lK[i+16][j + 4] += (B * q5 * rt2 * D2  + 2.0 * C * q5 * q2) * ShRC;

        // 34
        lK[i+ 8][j +12] += (B * q5 * rt2 * D2  + 2.0 * C * q4 * q3) * ShRC;
        lK[i+12][j + 8] += (B * q5 * rt2 * D2  + 2.0 * C * q4 * q3) * ShRC;
        // 35
        lK[i+ 8][j +16] += (B * q4 * rt2 * D2  + 2.0 * C * q5 * q3) * ShRC;
        lK[i+16][j + 8] += (B * q4 * rt2 * D2  + 2.0 * C * q5 * q3) * ShRC;

        // 45
        lK[i+12][j+16] += (B * q3 * rt2 * D2 + 2.0 * C * q5 * q4) * ShRC;
        lK[i+16][j+12] += (B * q3 * rt2 * D2 + 2.0 * C * q5 * q4) * ShRC;
      }
    }
  }

  inline void assembleElasticL1(double lK[20][20],
                                double lL[20],
                                const GaussianQuadratureTet<11> &shapes,
                                const double &tetDeterminant,
                                const double &q1x, const double &q1y, const double &q1z,
                                const double &q2x, const double &q2y, const double &q2z,
                                const double &q3x, const double &q3y, const double &q3z,
                                const double &q4x, const double &q4y, const double &q4z,
                                const double &q5x, const double &q5y, const double &q5z,
                                const double &L1) {
    const double mul = shapes.weight() * tetDeterminant;
    for (int i = 0; i < 4; i++) {
      const double ShRx = shapes.Nx(i);
      const double ShRy = shapes.Ny(i);
      const double ShRz = shapes.Nz(i);

      lL[i + 0] += (ShRx * q1x + ShRy * q1y + ShRz * q1z) * L1 * mul;
      lL[i + 4] += (ShRx * q2x + ShRy * q2y + ShRz * q2z) * L1 * mul;
      lL[i + 8] += (ShRx * q3x + ShRy * q3y + ShRz * q3z) * L1 * mul;
      lL[i + 12] += (ShRx * q4x + ShRy * q4y + ShRz * q4z) * L1 * mul;
      lL[i + 16] += (ShRx * q5x + ShRy * q5y + ShRz * q5z) * L1 * mul;

      for (int j = 0; j < 4; j++) {
        const double ShCx = shapes.Nx(j);
        const double ShCy = shapes.Ny(j);
        const double ShCz = shapes.Nz(j);

        // L1- matrix term 'dot' only appears on diagonal
        const double dot = (ShRx * ShCx + ShRy * ShCy + ShRz * ShCz) * L1 * mul;
        lK[i   ][j   ] += dot;
        lK[i + 4 ][j + 4 ] += dot;
        lK[i + 8 ][j + 8 ] += dot;
        lK[i + 12][j + 12] += dot;
        lK[i + 16][j + 16] += dot;
      }
    }
  }

  inline void assembleThreeElasticConstants(double lK[20][20],
                                            double lL[20],
                                            const GaussianQuadratureTet<11> &shapes,
                                            const double &tetDeterminant,
                                            const double &q1, const double &q2, const double &q3, const double &q4, const double &q5,
                                            const double &q1x, const double &q1y, const double &q1z,
                                            const double &q2x, const double &q2y, const double &q2z,
                                            const double &q3x, const double &q3y, const double &q3z,
                                            const double &q4x, const double &q4y, const double &q4z,
                                            const double &q5x, const double &q5y, const double &q5z,
                                            const double &L2, const double &L3, const double &L6
                                            ) {
    const double rt23 = rt2 * rt3;
    const double mul = shapes.weight() * tetDeterminant;
    for (int i = 0; i < 4; i++) {
      const double ShR = shapes.N(i);
      const double ShRx = shapes.Nx(i);
      const double ShRy = shapes.Ny(i);
      const double ShRz = shapes.Nz(i);

      //L2 and L6 q1 terms
      double temp = (ShRx * q1x * D6 - ShRx * rt3 * q2x * D6 - ShRx * rt3 * q3y * D6 - ShRx * rt3 * q5z * D6 -
                     ShRy * q3x * rt3 * D6 + ShRy * q1y * D6 + ShRy * rt3 * q2y * D6 - ShRy * rt3 * q4z * D6 +
                     ShRz * q5x * rt3 * D3 + ShRz * q4y * rt3 * D3 + 2.0 * D3 * ShRz * q1z) * L2;
      temp += L3 / 6 * ((ShRx * (2 * q5z - q2x - q3y) + ShRy * (q2y - q3x + 2 * q4z) - ShRz * (q4y + q5x)) * rt3 +
                        (ShRx * q1x + ShRy * q1y + 4 * ShRz * q1z));
      // L6 q1 term
      temp += (-ShRx * q1x * q1 * rt23 * D6 - ShRy * q1y * q1 * rt23 * D6 + ShRz * q1 * rt23 * q1z * D3 -
               ShR * rt23 * q5y * q5y * D2 * D6 + ShRx * q5 * rt2 * q1z * D2 + ShRy * q3 * rt2 * q1x * D2 +
               ShRy * q4 * rt2 * q1z * D2 + ShRz * q5 * rt2 * q1x * D2 - ShR * rt23 * q4x * q4x * D6 * D2 -
               ShR * rt23 * q4y * q4y * D6 * D2 - ShR * rt23 * q3x * q3x * D2 * D6 - ShR * rt23 * q2y * q2y * D2 * D6 +
               ShR * rt23 * q2z * q2z * D6 + ShR * rt23 * q3z * q3z * D6 + ShR * rt23 * q5z * q5z * D6 +
               ShR * rt23 * q4z * q4z * D6 + ShR * rt23 * q1z * q1z * D6 + ShRz * q4 * rt2 * q1y * D2 -
               ShR * rt23 * q1x * q1x * D2 * D6 - ShR * rt23 * q3y * q3y * D2 * D6 - ShR * rt23 * q2x * q2x * D2 * D6 -
               ShR * rt23 * q1y * q1y * D2 * D6 + ShRx * q1x * q2 * rt2 * D2 + ShRx * q3 * rt2 * q1y * D2 -
               ShRy * q1y * q2 * rt2 * D2 - ShR * rt23 * q5x * q5x * D2 * D6) * L6;
      lL[i + 0] += temp * mul;
      // L2 and L6 q2 term
      temp = (-ShRx * q1x * rt3 * D6 + q2x * ShRx * D2 + q3y * ShRx * D2 + q5z * ShRx * D2 - q3x * ShRy * D2 +
              ShRy * q1y * rt3 * D6 + q2y * ShRy * D2 - q4z * ShRy * D2) * L2;
      temp += L3 / 2 * (ShRx * (q2x - q3y - q1x / rt3) + ShRy * (q2y + q3x + q1y / rt3) + ShRz * (q5x - q4y));
      temp += (-ShR * rt2 * q5y * q5y * D4 - ShR * rt2 * q1y * q1y * D4 + ShR * rt2 * q4x * q4x * D4 +
               ShR * rt2 * q1x * q1x * D4 - ShR * rt2 * q3y * q3y * D4 + ShR * rt2 * q2x * q2x * D4 +
               ShR * rt2 * q5x * q5x * D4 + ShR * rt2 * q3x * q3x * D4 - ShR * rt2 * q4y * q4y * D4 -
               ShR * rt2 * q2y * q2y * D4 - ShRx * q1 * rt23 * q2x * D6 + ShRx * rt2 * q2 * q2x * D2 +
               ShRx * q3 * q2y * rt2 * D2 + ShRx * q5 * q2z * rt2 * D2 - ShRy * q1 * rt23 * q2y * D6 -
               ShRy * rt2 * q2 * q2y * D2 + ShRy * q3 * q2x * rt2 * D2 + ShRy * q4 * q2z * rt2 * D2 +
               ShRz * q5 * q2x * rt2 * D2 + ShRz * q4 * q2y * rt2 * D2 + ShRz * q1 * rt23 * q2z * D3) * L6;
      lL[i + 4] += temp * mul;
      // L2 and L6 q3 terms
      temp = (ShRx * q3x * D2 - ShRx * q1y * rt3 * D6 - ShRx * q2y * D2 + ShRx * q4z * D2 - ShRy * q1x * rt3 * D6 +
              ShRy * q2x * D2 + ShRy * q3y * D2 + ShRy * q5z * D2) * L2;
      temp += L3 / 2 * (ShRx * (q2y + q3x - q1y / rt3) + ShRy * (q3y - q2x - q1x / rt3) + ShRz * (q4x + q5y));
      temp += (ShR * rt2 * q1x * q1y * D2 + ShR * rt2 * q2x * q2y * D2 + ShR * rt2 * q3x * q3y * D2 +
               ShR * rt2 * q5x * q5y * D2 + ShR * rt2 * q4x * q4y * D2 - ShRx * q3x * q1 * rt23 * D6 +
               ShRx * q3x * q2 * rt2 * D2 + ShRx * q3 * rt2 * q3y * D2 + ShRx * q5 * rt2 * q3z * D2 -
               ShRy * q3y * q1 * rt23 * D6 - ShRy * q3y * q2 * rt2 * D2 + ShRy * q3 * rt2 * q3x * D2 +
               ShRy * q4 * rt2 * q3z * D2 + ShRz * q5 * rt2 * q3x * D2 + ShRz * q4 * rt2 * q3y * D2 +
               ShRz * q1 * rt23 * q3z * D3) * L6;
      lL[i + 8] += temp * mul;
      // L2 and L6 q4 terms
      temp = (ShRy * q5x * D2 + ShRy * q4y * D2 + ShRy * q1z * rt3 * D3 + ShRz * q3x * D2 - ShRz * q1y * rt3 * D6 -
              ShRz * q2y * D2 + ShRz * q4z * D2) * L2;
      temp += L3 / 2 * (ShRx * (q3z + q5y) + ShRy * (q4y - q2z - q1z / rt3) + ShRz * (q4z + 2 * q1y / rt3));
      temp += (ShR * rt2 * q1y * q1z * D2 + ShR * rt2 * q2y * q2z * D2 + ShR * rt2 * q3y * q3z * D2 +
               ShR * rt2 * q5y * q5z * D2 + ShR * rt2 * q4y * q4z * D2 - ShRx * q4x * q1 * rt23 * D6 +
               ShRx * q4x * q2 * rt2 * D2 + ShRx * q3 * rt2 * q4y * D2 + ShRx * q5 * rt2 * q4z * D2 -
               ShRy * q4y * q1 * rt23 * D6 - ShRy * q4y * q2 * rt2 * D2 + ShRy * q3 * rt2 * q4x * D2 +
               ShRy * q4 * rt2 * q4z * D2 + ShRz * q5 * rt2 * q4x * D2 + ShRz * q4 * rt2 * q4y * D2 +
               ShRz * q1 * rt23 * q4z * D3) * L6;
      lL[i + 12] += temp * mul;
      // L2 and L6 q5 terms
      temp = (ShRx * q5x * D2 + ShRx * q4y * D2 + ShRx * q1z * rt3 * D3 - ShRz * q1x * rt3 * D6 + ShRz * q2x * D2 +
              ShRz * q3y * D2 + ShRz * q5z * D2) * L2;
      temp += L3 / 2 * (ShRx * (q2z + q5x - q1z / rt3) + ShRy * (q3z + q4x) + ShRz * (q5z + 2 * q1x / rt3));
      temp += (ShR * rt2 * q1x * q1z * D2 + ShR * rt2 * q2x * q2z * D2 + ShR * rt2 * q3x * q3z * D2 +
               ShR * rt2 * q5x * q5z * D2 + ShR * rt2 * q4x * q4z * D2 - ShRx * q5x * q1 * rt23 * D6 +
               ShRx * q5x * q2 * rt2 * D2 + ShRx * q3 * rt2 * q5y * D2 + ShRx * q5 * rt2 * q5z * D2 -
               ShRy * q5y * q1 * rt23 * D6 - ShRy * q5y * q2 * rt2 * D2 + ShRy * q3 * rt2 * q5x * D2 +
               ShRy * q4 * rt2 * q5z * D2 + ShRz * q5 * rt2 * q5x * D2 + ShRz * q4 * rt2 * q5y * D2 +
               ShRz * q1 * rt23 * q5z * D3) * L6;
      lL[i + 16] += temp * mul;

      for (int j = 0; j < 4; j++) {
        const double ShC = shapes.N(j);
        const double ShCx = shapes.Nx(j);
        const double ShCy = shapes.Ny(j);
        const double ShCz = shapes.Nz(j);

        // dlL[0]/dq1 -> dlL[0]/dq5 L2 and L6 terms
        temp = (ShRx * ShCx * D6 + ShRy * ShCy * D6 + 2.0 * D3 * ShRz * ShCz) * L2;
        temp += L3/6*(ShCx*ShRx + ShCy*ShRy + 4*ShCz*ShRz);
        temp += (-ShC * ShRx * q1x * rt23 * D6 - ShC * ShRy * q1y * rt23 * D6 + ShC * ShRz * rt23 * q1z * D3 - ShCx * ShRx * q1 * rt23 * D6 + ShCx * ShRy * q3 * rt2 / 2.0 + ShCx * ShRz * q5 * rt2 * D2 - ShCx * ShR * rt23 * q1x * D6 + ShCx * ShRx * q2 * rt2 * D2 - ShCy * ShRy * q1 * rt23 * D6 + ShCy * ShRz * q4 * rt2 * D2 - ShCy * ShR * rt23 * q1y * D6 + ShCy * ShRx * q3 * rt2 * D2 - ShCy * ShRy * q2 * rt2 * D2 + ShCz * ShRz * q1 * rt23 * D3 + ShCz * ShRx * q5 * rt2 * D2 + ShCz * ShRy * q4 * rt2 * D2 + ShCz * ShR * rt23 * q1z * D3) * L6;
        lK[i][j   ] += temp * mul;

        temp = (-ShRx * rt3 * ShCx * D6 + ShRy * rt3 * ShCy * D6) * L2;
        temp += L3*rt3/6*(ShCy*ShRy - ShCx*ShRx);
        temp += (ShC * ShRx * q1x * rt2 * D2 - ShC * ShRy * q1y * rt2 * D2 - ShR * rt23 * q2x * ShCx * D6 - ShR * rt23 * q2y * ShCy * D6 + ShR * rt23 * q2z * ShCz * D3) * L6;;
        lK[i  ][j + 4] += temp * mul;
        lK[j + 4][i  ] += temp * mul;

        temp = (-ShRy * rt3 * ShCx * D6 - ShRx * rt3 * ShCy * D6) * L2;
        temp += -L3*rt3/6*(ShCx*ShRy + ShCy*ShRx);
        temp += (ShC * ShRy * rt2 * q1x * D2 + ShC * ShRx * rt2 * q1y * D2 - ShR * rt23 * q3x * ShCx * D6 - ShR * rt23 * q3y * ShCy * D6 + ShR * rt23 * q3z * ShCz * D3) * L6;
        lK[i  ][j + 8] += temp * mul;
        lK[j + 8][i  ] += temp * mul;

        temp = (ShRz * rt3 * ShCy * D3 - ShRy * rt3 * ShCz * D6) * L2;
        temp += L3*rt3/6*(2*ShCz*ShRy - ShCy*ShRz);
        temp += (ShC * ShRy * rt2 * q1z * D2 + ShC * ShRz * rt2 * q1y * D2 - ShR * rt23 * q4x * ShCx * D6 - ShR * rt23 * q4y * ShCy * D6 + ShR * rt23 * q4z * ShCz * D3) * L6;
        lK[i   ][j + 12] += temp * mul;
        lK[j + 12][i   ] += temp * mul;

        temp = (ShRz * rt3 * ShCx * D3 - ShRx * rt3 * ShCz * D6) * L2;
        temp += L3*rt3/6*(2*ShCz*ShRx - ShCx*ShRz);
        temp += (ShC * ShRx * rt2 * q1z * D2 + ShC * ShRz * rt2 * q1x * D2 - ShR * rt23 * q5x * ShCx * D6 - ShR * rt23 * q5y * ShCy * D6 + ShR * rt23 * q5z * ShCz * D3) * L6;
        lK[i   ][j + 16] += temp * mul;
        lK[j + 16][i   ] += temp * mul;

        // dlL[4]/dq2 -> dlL[4]/dq5 L2 and L6 terms
        temp = 0.5 * (ShRx * ShCx + ShRy * ShCy) * L2;
        temp += L3/2*(ShCx*ShRx + ShCy*ShRy);
        temp += (ShC * ShRx * rt2 * q2x * D2 - ShC * ShRy * rt2 * q2y * D2 + ShCx * ShR * rt2 * q2x * D2 - ShCx * ShRx * q1 * rt23 * D6 + ShCx * ShRx * q2 * rt2 * D2 + ShCx * ShRy * q3 * rt2 * D2 + ShCx * ShRz * q5 * rt2 * D2 - ShCy * ShR * rt2 * q2y * D2 + ShCy * ShRx * q3 * rt2 * D2 - ShCy * ShRy * q1 * rt23 * D6 - ShCy * ShRy * q2 * rt2 * D2 + ShCy * ShRz * q4 * rt2 * D2 + ShCz * ShRx * q5 * rt2 * D2 + ShCz * ShRy * q4 * rt2 * D2 + ShCz * ShRz * q1 * rt23 * D3) * L6;
        lK[i + 4][j + 4] += temp * mul;

        temp = 0.5 * (-ShRy * ShCx + ShRx * ShCy) * L2;
        temp += L3/2*(ShCx*ShRy - ShCy*ShRx);
        temp += (ShC * ShRx * q2y * rt2 * D2 + ShC * ShRy * q2x * rt2 * D2 + ShR * rt2 * q3x * ShCx * D2 - ShR * rt2 * q3y * ShCy * D2) * L6;
        lK[i + 4][j + 8] += temp * mul;
        lK[j + 8][i + 4] += temp * mul;

        temp = -0.5 * ShCz * L2 * ShRy;
        temp += -L3/2*ShCy*ShRz;
        temp += (ShC * ShRy * q2z * rt2 * D2 + ShC * ShRz * q2y * rt2 * D2 + ShR * rt2 * q4x * ShCx * D2 - ShR * rt2 * q4y * ShCy * D2) * L6;
        lK[i + 4 ][j + 12] += temp * mul;
        lK[j + 12][i + 4 ] += temp * mul;

        temp = 0.5 * ShCz * L2 * ShRx;
        temp += L3/2*ShCx*ShRz;
        temp += (ShC * ShRx * q2z * rt2 * D2 + ShC * ShRz * q2x * rt2 * D2 + ShR * rt2 * q5x * ShCx * D2 - ShR * rt2 * q5y * ShCy * D2) * L6;
        lK[i + 4 ][j + 16] += temp * mul;
        lK[j + 16][i + 4 ] += temp * mul;

        // dlL[8] / dq3 -> dlL[8]/dq5 L2 and L6 terms
        temp = 0.5 * (ShRx * ShCx + ShRy * ShCy) * L2;
        temp += L3/2*(ShCx*ShRx + ShCy*ShRy);
        temp += (ShC * ShRx * rt2 * q3y * D2 + ShC * ShRy * rt2 * q3x * D2 + ShCx * ShR * rt2 * q3y * D2 - ShCx * ShRx * q1 * rt23 * D6 + ShCx * ShRx * q2 * rt2 * D2 + ShCx * ShRy * q3 * rt2 * D2 + ShCx * ShRz * q5 * rt2 * D2 + ShCy * ShR * rt2 * q3x * D2 + ShCy * ShRx * q3 * rt2 * D2 - ShCy * ShRy * q1 * rt23 * D6 - ShCy * ShRy * q2 * rt2 * D2 + ShCy * ShRz * q4 * rt2 * D2 + ShCz * ShRx * q5 * rt2 * D2 + ShCz * ShRy * q4 * rt2 * D2 + ShCz * ShRz * q1 * rt23 * D3) * L6;
        lK[i + 8][j + 8] += temp  * mul;

        temp = ShCz * L2 * ShRx * 0.5;
        temp += L3/2*ShCx*ShRz;
        temp += (ShC * ShRy * rt2 * q3z * D2 + ShC * ShRz * rt2 * q3y * D2 + ShR * rt2 * q4y * ShCx * D2 + ShR * rt2 * q4x * ShCy * D2) * L6;
        lK[i + 8 ][j + 12] += temp * mul;
        lK[j + 12][i + 8 ] += temp * mul;

        temp = 0.5 * ShCz * L2 * ShRy;
        temp += L3/2*ShCy*ShRz;
        temp += (ShC * ShRx * rt2 * q3z * D2 + ShC * ShRz * rt2 * q3x * D2 + ShR * rt2 * q5y * ShCx * D2 + ShR * rt2 * q5x * ShCy * D2) * L6;
        lK[i + 8 ][j + 16] += temp * mul;
        lK[j + 16][i + 8 ] += temp * mul;

        // dlL[12] / dq4 -> dlL[8]/dq5 L2 and L6 terms
        temp = 0.5 * (ShRy * ShCy + ShRz * ShCz) * L2;
        temp += L3/2*(ShCy*ShRy + ShCz*ShRz);
        temp += (ShC * ShRy * rt2 * q4z * D2 + ShC * ShRz * rt2 * q4y * D2 - ShCx * ShRx * q1 * rt23 * D6 + ShCx * ShRx * q2 * rt2 * D2 + ShCx * ShRy * q3 * rt2 * D2 + ShCx * ShRz * q5 * rt2 * D2 + ShCy * ShR * rt2 * q4z * D2 + ShCy * ShRx * q3 * rt2 * D2 - ShCy * ShRy * q1 * rt23 * D6 - ShCy * ShRy * q2 * rt2 * D2 + ShCy * ShRz * q4 * rt2 * D2 + ShCz * ShR * rt2 * q4y * D2 + ShCz * ShRx * q5 * rt2 * D2 + ShCz * ShRy * q4 * rt2 * D2 + ShCz * ShRz * q1 * rt23 * D3) * L6;
        lK[i + 12][j + 12] += temp * mul;

        temp = 0.5 * ShCx * L2 * ShRy;
        temp += L3/2*ShCy*ShRx;
        temp += (ShC * ShRx * rt2 * q4z * D2 + ShC * ShRz * rt2 * q4x * D2 + ShR * rt2 * q5z * ShCy * D2 + ShR * rt2 * q5y * ShCz * D2) * L6;
        lK[i + 12][j + 16] += temp * mul;
        lK[j + 16][i + 12] += temp * mul;
        // dlL[12] / dq5 L2 and L6 terms
        temp = 0.5 * (ShRx * ShCx + ShRz * ShCz) * L2;
        temp += L3/2*(ShCx*ShRx + ShCz*ShRz);
        temp += (ShC * ShRx * rt2 * q5z * D2 + ShC * ShRz * rt2 * q5x * D2 + ShCx * ShR * rt2 * q5z * D2 - ShCx * ShRx * q1 * rt23 * D6 + ShCx * ShRx * q2 * rt2 * D2 + ShCx * ShRy * q3 * rt2 * D2 + ShCx * ShRz * q5 * rt2 * D2 + ShCy * ShRx * q3 * rt2 * D2 - ShCy * ShRy * q1 * rt23 * D6 - ShCy * ShRy * q2 * rt2 * D2 + ShCy * ShRz * q4 * rt2 * D2 + ShCz * ShR * rt2 * q5x * D2 + ShCz * ShRx * q5 * rt2 * D2 + ShCz * ShRy * q4 * rt2 * D2 + ShCz * ShRz * q1 * rt23 * D3) * L6;
        lK[i + 16][j + 16] += temp * mul;
      }
    }
  }


  inline void assembleDielectric(double lL[20],
                                 const GaussianQuadratureTet<11> &shapes,
                                 const double &tetDeterminant,
                                 const double &Vx,
                                 const double &Vy,
                                 const double &Vz,
                                 const double &deleps) {
    const double mul = shapes.weight() * tetDeterminant;
    for (int i = 0; i < 4; i++) {
      const double ShR = shapes.N(i) * mul;
      lL[i + 0] += (rt6 * (Vx * Vx + Vy * Vy - 2.0 * Vz * Vz) * deleps * D3 * D6 * eps0) * ShR;
      lL[i + 4] += (-rt2 * (Vx * Vx - Vy * Vy) * deleps * D6 * eps0) * ShR;
      lL[i + 8] += (-rt2 * Vx * Vy * deleps * D3 * eps0) * ShR;
      lL[i + 12] += (-rt2 * Vz * Vy * deleps * D3 * eps0) * ShR;
      lL[i + 16] += (-rt2 * Vx * Vz * deleps * D3 * eps0) * ShR;
    }
  }

}

#endif //PROJECT_QLC3D_LC_ENERGY_TERMS_H
