#include <math.h>
#include <solutionvector.h>
#include <lc.h>
#include <solver-settings.h>
#include <geometry.h>
#include <shapefunction3d.h>
#include <shapefunction2d.h>
#include <util/logging.h>
#include <util/hash.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <potential/potential-solver.h>

#include <random>

// SPAMTRIX INCLUDES
#include "spamtrix_ircmatrix.hpp"
#include "spamtrix_vector.hpp"
#include "spamtrix_iterativesolvers.hpp"
#include "spamtrix_luincpreconditioner.hpp"

const idx npt = 4; //Number of Points per Tetrahedra

double rt2 = sqrt(2.0);
double rt3 = sqrt(3.0);
double rt6 = sqrt(6.0);

void Pot_GMRES(const SpaMtrix::IRCMatrix &K,
               const SpaMtrix::Vector &B,
               SpaMtrix::Vector &X,
               const SolverSettings &settings);

// DECLARATION ONLY
void setUniformEField(const Electrodes &electrodes, SolutionVector &v, const Coordinates &coordinates);

inline void localKL(
        const Geometry &geom,
        double lK[npt][npt],
        double lL[npt],
        unsigned int elementIndex,
        const SolutionVector &q,
        const LC &lc,
        const Electrodes &electrodes,
        const Shape4thOrder &shapes)

{
  idx i, j;
  double eper, deleps;
  double S0 = lc.S0();
  eper    = 0;
  deleps  = 0;
  double efe = 0.0, efe2 = 0.0;
  const Mesh &mesh = geom.getTetrahedra();
  if (mesh.getMaterialNumber(elementIndex) == MAT_DOMAIN1) { // if LC element
    eper = lc.eps_per() / S0;
    deleps = (lc.eps_par() - lc.eps_per()) / S0;
    efe  = 2.0 / (9 * S0) * (lc.e1() + 2 * lc.e3());
    efe2 = 4.0 / (9 * S0 * S0) * (lc.e1() - lc.e3());
  } else { // otherwise dielectric
    idx ind_de = mesh.getDielectricNumber(elementIndex) - 1; // -1 for 0 indexing
    eper = electrodes.getDielectricPermittivity(ind_de);
  }
  idx tt[4];
  mesh.loadNodes(elementIndex, tt);
  memset(lK, 0, 4 * 4 * sizeof(double));
  memset(lL, 0, 4 * sizeof(double));

  double Jdet = mesh.getDeterminant(elementIndex);
  const Coordinates &coordinates = geom.getCoordinates();
  Vec3 p[4] = {coordinates.getPoint(tt[0]),
               coordinates.getPoint(tt[1]),
               coordinates.getPoint(tt[2]),
               coordinates.getPoint(tt[3])};
  // Jacobian
  double xr, xs, xt, yr, ys, yt, zr, zs, zt;
  xr = xs = xt = yr = ys = yt = zr = zs = zt = 0.0;
  for (i = 0; i < npt; i++) {
    double x = p[i].x() * 1e-6;
    double y = p[i].y() * 1e-6;
    double z = p[i].z() * 1e-6;
    xr += x * shapes.sh1r[0][i];
    xs += x * shapes.sh1s[0][i];
    xt += x * shapes.sh1t[0][i];
    yr += y * shapes.sh1r[0][i];
    ys += y * shapes.sh1s[0][i];
    yt += y * shapes.sh1t[0][i];
    zr += z * shapes.sh1r[0][i];
    zs += z * shapes.sh1s[0][i];
    zt += z * shapes.sh1t[0][i];
  }//end for i
  if (Jdet < 0) Jdet = -Jdet;
  double Jinv[3][3] = {{(zt * ys - yt * zs) / Jdet, (xt * zs - zt * xs) / Jdet, (xs * yt - ys * xt) / Jdet}
          , {(yt * zr - zt * yr) / Jdet, (zt * xr - xt * zr) / Jdet, (xt * yr - yt * xr) / Jdet}
          , {(yr * zs - ys * zr) / Jdet, (xs * zr - xr * zs) / Jdet, (ys * xr - xs * yr) / Jdet}
  };
  double Sh[4], dSh[4][3];
  // x,y,z derivatives of shape functions
  for (i = 0; i < 4; i++) {
    dSh[i][0] = shapes.sh1r[0][i] * Jinv[0][0] +
                shapes.sh1s[0][i] * Jinv[1][0] +
                shapes.sh1t[0][i] * Jinv[2][0];
    dSh[i][1] = shapes.sh1r[0][i] * Jinv[0][1] +
                shapes.sh1s[0][i] * Jinv[1][1] +
                shapes.sh1t[0][i] * Jinv[2][1];
    dSh[i][2] = shapes.sh1r[0][i] * Jinv[0][2] +
                shapes.sh1s[0][i] * Jinv[1][2] +
                shapes.sh1t[0][i] * Jinv[2][2];
  }//end for i
  for (unsigned int igp = 0; igp < shapes.ngp; igp++) {
    Sh[0] = shapes.sh1[igp][0];
    Sh[1] = shapes.sh1[igp][1];
    Sh[2] = shapes.sh1[igp][2];
    Sh[3] = shapes.sh1[igp][3];
    double e11, e22, e33, e12, e13, e23;
    e11 = 0; e22 = 0; e33 = 0; e12 = 0; e13 = 0; e23 = 0;
    // Function variables and derivatives
    double q1 = 0, q2 = 0, q3 = 0, q4 = 0, q5 = 0;
    double q1x = 0, q2x = 0, q3x = 0, q4x = 0, q5x = 0;
    double q1y = 0, q2y = 0, q3y = 0, q4y = 0, q5y = 0;
    double q1z = 0, q2z = 0, q3z = 0, q4z = 0, q5z = 0;
    //if (1){//
    if (mesh.getMaterialNumber(elementIndex) != MAT_DOMAIN1) { // if this element is not LC
      // only set diagonal permittivities to non-zero
      for (i = 0 ; i < npt ; ++i) {
        e11 += Sh[i] * eper;
        e22 += Sh[i] * eper;
        e33 += Sh[i] * eper;
      }
    } else { // otherwise LC ->
      for (i = 0; i < npt; i++) {
        q1 += Sh[i] * q.getValue(tt[i], 0);
        q2 += Sh[i] * q.getValue(tt[i], 1);
        q3 += Sh[i] * q.getValue(tt[i], 2);
        q4 += Sh[i] * q.getValue(tt[i], 3);
        q5 += Sh[i] * q.getValue(tt[i], 4);
        q1x += dSh[i][0] * q.getValue(tt[i], 0);
        q2x += dSh[i][0] * q.getValue(tt[i], 1);
        q3x += dSh[i][0] * q.getValue(tt[i], 2);
        q4x += dSh[i][0] * q.getValue(tt[i], 3);
        q5x += dSh[i][0] * q.getValue(tt[i], 4);
        q1y += dSh[i][1] * q.getValue(tt[i], 0);
        q2y += dSh[i][1] * q.getValue(tt[i], 1);
        q3y += dSh[i][1] * q.getValue(tt[i], 2);
        q4y += dSh[i][1] * q.getValue(tt[i], 3);
        q5y += dSh[i][1] * q.getValue(tt[i], 4);
        q1z += dSh[i][2] * q.getValue(tt[i], 0);
        q2z += dSh[i][2] * q.getValue(tt[i], 1);
        q3z += dSh[i][2] * q.getValue(tt[i], 2);
        q4z += dSh[i][2] * q.getValue(tt[i], 3);
        q5z += dSh[i][2] * q.getValue(tt[i], 4);

        e11 += Sh[i] * (((2.0 / 3.0 / S0) * (-q1 / rt6 + q2 / rt2) + (1.0 / 3.0)) * deleps + eper); //~nx*nx
        e22 += Sh[i] * (((2.0 / 3.0 / S0) * (-q1 / rt6 - q2 / rt2) + (1.0 / 3.0)) * deleps + eper); //~ny*ny
        e33 += Sh[i] * (((2.0 / 3.0 / S0) * (2.0 * q1 / rt6)        + (1.0 / 3.0)) * deleps + eper); //~nz*nz
        e12 += Sh[i] * (2.0 / 3.0 / S0) * (q3 / rt2) * deleps;           //~nx*ny
        e13 += Sh[i] * (2.0 / 3.0 / S0) * (q5 / rt2) * deleps;           //~nx*nz
        e23 += Sh[i] * (2.0 / 3.0 / S0) * (q4 / rt2) * deleps;
      }
    }
    // Local K and L
    double mul = shapes.w[igp] * Jdet;
    for (i = 0; i < 4; i++) {
      const double ShRx = mul * dSh[i][0];
      const double ShRy = mul * dSh[i][1];
      const double ShRz = mul * dSh[i][2];

      // Flexoelectric polarisation terms, minus, since minus residual formed
      lL[i] -= -efe*( (ShRx*(q2x + q3y + q5z) + ShRy*(q3x + q4z - q2y)
                       + ShRz*(q4y + q5x))/rt2 + (2*ShRz*q1z - ShRx*q1x - ShRy*q1y)/rt6 );

      lL[i] -= -efe2*((ShRx*q1*q1x + ShRy*q1*q1y + 4*ShRz*q1*q1z)/6
                      + ( ShRx*(q2*(q2x+q3y+q5z) + q3*(q3x-q2y+q4z) + q5*(q4y+q5x))
                          + ShRy*(q2*(q2y-q3x-q4z) + q3*(q2x+q3y+q5z) + q4*(q4y+q5x))
                          + ShRz*(q4*(q4z-q2y+q3x) + q5*(q3y+q2x+q5z)))/2
                      + ( ShRx*(-q1*(q2x+q3y+q5z) - q2*q1x - q3*q1y + 2*q5*q1z)
                          + ShRy*( q1*(q2y-q3x-q4z) + q2*q1y - q3*q1x + 2*q4*q1z)
                          + ShRz*(2*q1*(q4y+q5x) - q4*q1y - q5*q1x))*rt3/6);

      for (j = 0; j < 4; j++) {
        lK[i][j] += mul * (
                dSh[i][0] * dSh[j][0] * e11 +
                dSh[i][1] * dSh[j][1] * e22 +
                dSh[i][2] * dSh[j][2] * e33 +
                dSh[i][0] * dSh[j][1] * (e12) +
                dSh[i][1] * dSh[j][0] * (e12) +
                dSh[i][1] * dSh[j][2] * (e23) +
                dSh[i][2] * dSh[j][1] * (e23) +
                dSh[i][0] * dSh[j][2] * (e13) +
                dSh[i][2] * dSh[j][0] * (e13)
        );
      }//end for j
    }//end for i
  }//end for igp
}// end void localKL

void localKL_N(
        const Coordinates &coordinates,
        idx *tt,
        double lK[npt][npt],
        double lL[npt],
        int it,
        int index_to_Neumann,
        const Mesh &mesh,
        const Mesh &surf_mesh,
        const SolutionVector &q,
        const LC &lc,
        const ShapeSurf4thOrder &shapes) {
  int i, j;
  double S0 = lc.S0();
  double eper   = 4.5;
  double deleps = 0.0;    // default for dielectric material
  double efe = 0.0, efe2 = 0.0;
  if (mesh.getMaterialNumber(index_to_Neumann) == MAT_DOMAIN1) {
    eper   = lc.eps_per() / S0;
    deleps = (lc.eps_par() - lc.eps_per()) / S0;
    efe  = 2.0 / (9 * S0) * (lc.e1() + 2 * lc.e3());
    efe2 = 4.0 / (9 * S0 * S0) * (lc.e1() - lc.e3());
  } else {
    throw std::runtime_error(fmt::format("Expected DOMAIN1 material in {}, {}.", __FILE__, __func__));
  }
  memset(lK, 0, npt * npt * sizeof(double));
  memset(lL, 0, 4 * sizeof(double));
  Vec3 n = surf_mesh.getSurfaceNormal(it);
  double eDet = surf_mesh.getDeterminant(it);
  double Jdet = mesh.getDeterminant(index_to_Neumann);
#ifndef NDEBUG
  assert(eDet > 0);
  assert(Jdet > 0);
  if (abs(n.norm2() - 1.0) > 0.01) {
    throw std::runtime_error(fmt::format("Surface normal vector for triangle {} is not of unit length.", it));
  }
#endif

  double xr, xs, xt, yr, ys, yt, zr, zs, zt;
  xr = xs = xt = yr = ys = yt = zr = zs = zt = 0.0;
  for (i = 0; i < 4; i++) {
    Vec3 p = coordinates.getPoint(tt[i]); // tt is reordered volume element
    double x = p.x() * 1e-6;
    double y = p.y() * 1e-6;
    double z = p.z() * 1e-6;

    xr += x * shapes.sh1r[0][i];
    xs += x * shapes.sh1s[0][i];
    xt += x * shapes.sh1t[0][i];
    yr += y * shapes.sh1r[0][i];
    ys += y * shapes.sh1s[0][i];
    yt += y * shapes.sh1t[0][i];
    zr += z * shapes.sh1r[0][i];
    zs += z * shapes.sh1s[0][i];
    zt += z * shapes.sh1t[0][i];
  }//end for i
  double Jinv[3][3] = {{ (zt * ys - yt * zs) / Jdet , (xt * zs - zt * xs) / Jdet , (xs * yt - ys * xt) / Jdet}
          , { (yt * zr - zt * yr) / Jdet , (zt * xr - xt * zr) / Jdet , (xt * yr - yt * xr) / Jdet}
          , { (yr * zs - ys * zr) / Jdet , (xs * zr - xr * zs) / Jdet , (ys * xr - xs * yr) / Jdet}
  };
  double Sh[4], dSh[4][3];
  for (i = 0; i < 4; i++) {
    dSh[i][0] = shapes.sh1r[0][i] * Jinv[0][0] +
                shapes.sh1s[0][i] * Jinv[1][0] +
                shapes.sh1t[0][i] * Jinv[2][0];
    dSh[i][1] = shapes.sh1r[0][i] * Jinv[0][1] +
                shapes.sh1s[0][i] * Jinv[1][1] +
                shapes.sh1t[0][i] * Jinv[2][1];
    dSh[i][2] = shapes.sh1r[0][i] * Jinv[0][2] +
                shapes.sh1s[0][i] * Jinv[1][2] +
                shapes.sh1t[0][i] * Jinv[2][2];
  }//end for i
  // Jacobian
  for (unsigned int igp = 0; igp < shapes.ngps; igp++) {
    Sh[0] = shapes.sh1[igp][0];
    Sh[1] = shapes.sh1[igp][1];
    Sh[2] = shapes.sh1[igp][2];
    Sh[3] = shapes.sh1[igp][3];
    double e11 = 0, e22 = 0, e33 = 0, e12 = 0, e23 = 0, e13 = 0;
    // Function variables and derivatives
    double q1 = 0, q2 = 0, q3 = 0, q4 = 0, q5 = 0;
    double q1x = 0, q2x = 0, q3x = 0, q4x = 0, q5x = 0;
    double q1y = 0, q2y = 0, q3y = 0, q4y = 0, q5y = 0;
    double q1z = 0, q2z = 0, q3z = 0, q4z = 0, q5z = 0;

    for (i = 0; i < 4; i++) {
      q1 += Sh[i] * q.getValue(tt[i], 0);
      q2 += Sh[i] * q.getValue(tt[i], 1);
      q3 += Sh[i] * q.getValue(tt[i], 2);
      q4 += Sh[i] * q.getValue(tt[i], 3);
      q5 += Sh[i] * q.getValue(tt[i], 4);
      q1x += dSh[i][0] * q.getValue(tt[i], 0);
      q2x += dSh[i][0] * q.getValue(tt[i], 1);
      q3x += dSh[i][0] * q.getValue(tt[i], 2);
      q4x += dSh[i][0] * q.getValue(tt[i], 3);
      q5x += dSh[i][0] * q.getValue(tt[i], 4);
      q1y += dSh[i][1] * q.getValue(tt[i], 0);
      q2y += dSh[i][1] * q.getValue(tt[i], 1);
      q3y += dSh[i][1] * q.getValue(tt[i], 2);
      q4y += dSh[i][1] * q.getValue(tt[i], 3);
      q5y += dSh[i][1] * q.getValue(tt[i], 4);
      q1z += dSh[i][2] * q.getValue(tt[i], 0);
      q2z += dSh[i][2] * q.getValue(tt[i], 1);
      q3z += dSh[i][2] * q.getValue(tt[i], 2);
      q4z += dSh[i][2] * q.getValue(tt[i], 3);
      q5z += dSh[i][2] * q.getValue(tt[i], 4);

      e11 += shapes.sh1[igp][i] * (((2.0 / 3.0 / S0) * (-q1 / rt6 + q2 / rt2) + (1.0 / 3.0)) * deleps + eper); //~nx*nx
      e22 += shapes.sh1[igp][i] * (((2.0 / 3.0 / S0) * (-q1 / rt6 - q2 / rt2) + (1.0 / 3.0)) * deleps + eper); //~ny*ny
      e33 += shapes.sh1[igp][i] * (((2.0 / 3.0 / S0) * (2.0 * q1 / rt6)       + (1.0 / 3.0)) * deleps + eper); //~nz*nz
      e12 += shapes.sh1[igp][i] * (2.0 / 3.0 / S0) * (q3 / rt2) * deleps;           //~nx*ny
      e13 += shapes.sh1[igp][i] * (2.0 / 3.0 / S0) * (q5 / rt2) * deleps;           //~nx*nz
      e23 += shapes.sh1[igp][i] * (2.0 / 3.0 / S0) * (q4 / rt2) * deleps;       //~ny*nz
    }//end for i
    double mul = shapes.w[igp] * eDet;
    double nx = n.x(), ny = n.y(), nz = n.z(); // interior normal?
    for (i = 0; i < 4; i++) {
      const double ShR  = mul * Sh[i];

      // Flexoelectric polarisation terms, minus, since minus residual formed
      lL[i] -= -efe*ShR*( (nx*(q2x + q3y + q5z) + ny*(q3x + q4z - q2y)
                           + nz*(q4y + q5x))/rt2 + (2*nz*q1z - nx*q1x - ny*q1y)/rt6 );

      lL[i] -= -efe2*ShR*((nx*q1*q1x + ny*q1*q1y + 4*nz*q1*q1z)/6
                          + ( nx*(q2*(q2x+q3y+q5z) + q3*(q3x-q2y+q4z) + q5*(q4y+q5x))
                              + ny*(q2*(q2y-q3x-q4z) + q3*(q2x+q3y+q5z) + q4*(q4y+q5x))
                              + nz*(q4*(q4z-q2y+q3x) + q5*(q3y+q2x+q5z)))/2
                          + ( nx*(-q1*(q2x+q3y+q5z) - q2*q1x - q3*q1y + 2*q5*q1z)
                              + ny*( q1*(q2y-q3x-q4z) + q2*q1y - q3*q1x + 2*q4*q1z)
                              + nz*(2*q1*(q4y+q5x) - q4*q1y - q5*q1x))*rt3/6);

      for (j = 0; j < 4; j++) {
        lK[i][j] += mul * Sh[i] * (((e11 - 1) * dSh[j][0] + e12 * dSh[j][1] + e13 * dSh[j][2]) * nx
                                   + (e12 * dSh[j][0] + (e22 - 1) * dSh[j][1] + e23 * dSh[j][2]) * ny
                                   + (e13 * dSh[j][0] + e23 * dSh[j][1] + (e33 - 1) * dSh[j][2]) * nz);
      }//end for j
    }//end for i
  }//end for igp
}
// end void localKL


bool isFixedNode(idx i) {
  return i == NOT_AN_INDEX;
}

bool isFreeNode(idx i) {
  return i < NOT_AN_INDEX;
}

void assemble_volume(
        const Geometry &geometry,
        const SolutionVector &v,
        const SolutionVector &q,
        const LC &lc,
        SpaMtrix::IRCMatrix &K,
        SpaMtrix::Vector &L,
        const Electrodes &electrodes) {
  Shape4thOrder shapes;
  const unsigned int elementCount = geometry.getTetrahedra().getnElements();

  double lK[npt][npt];
  double lL[npt];
  idx t[npt];
  idx mapped[npt];

#pragma omp parallel for default(none) shared(geometry, v, q, lc, K, L, shapes, electrodes, elementCount) private(lK, lL, t, mapped)
  for (idx elementIndex = 0; elementIndex < elementCount; elementIndex++) {
    geometry.getTetrahedra().loadNodes(elementIndex, t);
    v.loadEquNodes(&t[0], &t[npt], mapped);

    localKL(geometry, lK, lL, elementIndex, q, lc, electrodes, shapes);

    ///*
    for (idx rowCounter = 0; rowCounter < npt; rowCounter++) {
      idx rowDof = mapped[rowCounter];

      for (idx colCounter = 0; colCounter < npt; colCounter++) {
        idx colDof = mapped[colCounter];

        if (isFreeNode(rowDof) && isFreeNode(colDof)) {
          // Free node at row contributes to the LHS matrix only
          //K.sparse_add(colDof, rowDof, lK[colCounter][rowCounter]); // OK periodic
          K.sparse_add(rowDof, colDof, lK[rowCounter][colCounter]); // OK periodic
        } else if (isFixedNode(rowDof) && isFreeNode(colDof)) {
          // Fixed note at row contributes to the RHS vector only for free col values: L = -K * v
          double fixedValue = v.getValue(t[rowCounter]);
          const double update = lL[rowCounter] - lK[rowCounter][colCounter] * fixedValue;
#pragma omp atomic
          L[colDof] += update;
        }
      }
    }

  }
}// end void assemble_volume

////////////////////////////////////////////////////////////////////////////////
// assemble Neumann boundaries
void assemble_Neumann(
        const Coordinates &coordinates,
        SolutionVector &v,
        const SolutionVector &q,
        const LC &lc,
        const Mesh &mesh,
        const Mesh &surf_mesh,
        SpaMtrix::IRCMatrix &K,
        SpaMtrix::Vector &L) {
  ShapeSurf4thOrder shapes;

  for (idx it = 0; it < surf_mesh.getnElements(); it++) {
    bool isNeumann = surf_mesh.getMaterialNumber(it) == MAT_NEUMANN;
    if (!isNeumann) {
      continue;
    }
    int index_to_Neumann = surf_mesh.getConnectedVolume(it);
    if (index_to_Neumann == -1) {
      throw std::runtime_error(fmt::format("Surface element {} mit material Neumann is not connected to any volume element.", it));
    }
      double lK[4][4];
      double lL[4];
      idx ee[3] = {   surf_mesh.getNode(it, 0) ,
                      surf_mesh.getNode(it, 1) ,
                      surf_mesh.getNode(it, 2)
      } ;
      idx tt[4] = { tt[0] =   mesh.getNode(index_to_Neumann, 0),
              mesh.getNode(index_to_Neumann, 1),
              mesh.getNode(index_to_Neumann, 2),
              mesh.getNode(index_to_Neumann, 3)
      };
      int intr = 0;//find  index to internal node
      for (int i = 0; i < 4; i++) {
        if ((tt[i] != ee[0]) && (tt[i] != ee[1]) && (tt[i] != ee[2])) {
          intr = i;
          break;
        }
      }
      idx ti[4] = { ee[0], ee[1], ee[2], tt[intr] }; // reordered local element, internal node is always last
      localKL_N(coordinates, &ti[0], lK , lL, it, index_to_Neumann, mesh, surf_mesh, q, lc, shapes);
      for (int rowCounter = 0; rowCounter < 4; rowCounter++) {
        idx rowDof = v.getEquNode(ti[rowCounter]);
        if (rowDof == NOT_AN_INDEX) { // HANDLE FIXED NODE
          double fixedValue = v.getValue(ti[rowCounter]);

          for (int colCounter = 0; colCounter < 4 ; colCounter++) {
            idx colDof = v.getEquNode(ti[colCounter]);   // CONNECTED NODE DOF ORDER
            if (colDof != NOT_AN_INDEX) {
              //L[colDof ] += lL[rowCounter] - lK[colCounter][rowCounter] * v.getValue(ti[rowCounter]);
              // Fixed note at row contributes to the RHS vector only for free col values: L = -K * v

              const double update = lL[rowCounter] - lK[colCounter][rowCounter] * fixedValue;
              L[colDof] += update;
            }
          }
        }// END HANDLE FIXED NODE
        else { // HANDLE FREE NODE
          for (int colCounter = 0; colCounter < 4; colCounter++) { // FOR COLUMNS
            idx colDof = v.getEquNode(ti[colCounter]);
            if (colDof != NOT_AN_INDEX) {// NON-FIXED NODE
              K.sparse_add(rowDof, colDof, lK[rowCounter][colCounter]);
            }
          }//end for j
        }// END HANDLE FREE NODES
      }//end for i
  }//end for it
}//end void assemble_Neumann

void calcpot3d(
    SpaMtrix::IRCMatrix &Kpot,
    SolutionVector &v,
    const SolutionVector &q,
    const LC &lc,
    const Geometry &geom,
    const SolverSettings &settings,
    const Electrodes &electrodes) {
    // First check whether potential calculation is actually needed...
    // NO NEED TO CALCULATE POTENTIAL IF...
    if ((v.getnFixed() == 0) ||        // no fixed potential nodes OR
            (!electrodes.getCalcPot())) {  // no need to calculate potential
        v.setValuesTo(0.0); // if no potential calculation, set all values to zero
        if (electrodes.hasElectricField()) {
            setUniformEField(electrodes, v, geom.getCoordinates());
        }
        return;
    }
    Kpot = 0.0; // clears values but keeps sparsity structure
    SpaMtrix::Vector L(v.getnFreeNodes());
    SpaMtrix::Vector V(v.getnFreeNodes());
    // PROVIDE POTENTIAL FROM PREVIOUS STEP AS INITIAL GUESS
    for (idx i = 0 ; i < v.getnDoF(); i++) {
        const idx ind = v.getEquNode(i);
        if (ind != NOT_AN_INDEX) {
            V[ind] = v.getValue(i);
        }
    }
    // Assemble system
    omp_set_num_threads(1);
    assemble_volume(geom, v, q, lc, Kpot , L, electrodes);
    assemble_Neumann(geom.getCoordinates() , v , q , lc , geom.getTetrahedra() , geom.getTriangles(), Kpot , L);

#ifdef LOG_DEBUG_HASH
    int64_t lHash = hashCode64(&L[0], &L[L.getLength() - 1]);
    Log::info("Calcpot. lHash={:X}", lHash);
#endif
    // GMRES SOLVER
    Pot_GMRES(Kpot, L, V, settings);
    // COPY NON-FIXED VALUES BACK TO SOLUTIONVECTOR
    for (idx i = 0 ; i < v.getnDoF() ; i++) {
        idx ind = v.getEquNode(i);
        // EQU NODES OF FIXED DOFS ARE ALL "NOT_AN_INDEX"
        if (ind != NOT_AN_INDEX) {
            v.setValue(i, 0, V[ind]);
        }
    }
}

void setUniformEField(const Electrodes &electrodes, SolutionVector &v, const Coordinates &coordinates) {
  // SETS POTENTIAL VALUES IN V SUCH THA A UNIFORM E-FIELD IS CREATED
  // CENTRE OF STRUCTURE IS ASSUMED TO BE AT 0V, VOLTAGE VALUES FOR NODES
  // IN THE DIRECTION OF THE EFIELD ARE SET ACCORDING TO THEIR DISTANCE
  // TO THE CENTRE
  // GET E-FIELD DIRECTION VECTOR AND MAGNITUDE
  Vec3 E = electrodes.getElectricField();
  double Emag = E.norm();
  Vec3 Ehat = E.normalized();
  //1. CALCULATE CENTRE OF STRUCTURE
  int np = v.getnDoF();
  double xmax = 0.0;
  double ymax = 0.0;
  double zmax = 0.0;
  for (int i = 0 ; i < np ; i++) {
    Vec3 p = coordinates.getPoint(i);
    xmax = p.x() > xmax ? p.x() : xmax;
    ymax = p.y() > ymax ? p.y() : ymax;
    zmax = p.z() > zmax ? p.z() : zmax;
  }
  Vec3 centre = {xmax / 2.0 , ymax / 2.0 , zmax / 2.0}; // centre of structure, assuming min is (0, 0, 0)

  //2. LOOP OVER EACH NODE AND CACLCULATE ITS DISTANCE TO CENTRE
  for (int i = 0 ; i < np ; i++) {
    Vec3 pos = coordinates.getPoint(i);
    Vec3 vec = centre - pos;

    // want distance along EField, i.e. dot product
    double dist = vec.dot(Ehat);
    // set potential value as distance*magnitude
    v.setValue(i, 0, dist * Emag + v.getValue(i));
  }
}

/**
 * Solves the Linear simulatenous equation Ax=b using the GMRES method.
 */
void Pot_GMRES(const SpaMtrix::IRCMatrix &K,
               const SpaMtrix::Vector &B,
               SpaMtrix::Vector &X,
               const SolverSettings &settings) {
    idx size = K.getNumRows();
    idx maxiter     = settings.getV_GMRES_Maxiter();
    idx restart     = settings.getV_GMRES_Restart();
    maxiter = maxiter < size ? maxiter : size;
    restart = restart < maxiter ? restart : maxiter;
    double toler    = settings.getV_GMRES_Toler();
    SpaMtrix::LUIncPreconditioner LU(K); // DOES THIS HAVE TO BE RECOMPUTED EACH TIME??
    SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);
    if (!solver.gmres(K, X, B, LU)) {
        Log::warn("GMRES did not converge in {} iterations when solving for potential. Tolerance achieved is {}.",
                  solver.maxIter, solver.toler);
    }
}
