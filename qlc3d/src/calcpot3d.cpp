#include <calcpot3d.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#include <solutionvector.h>
#include <lc.h>
#include <settings.h>
#include <geometry.h>
#include <shapefunction3d.h>
#include <shapefunction2d.h>
#include <util/logging.h>

// SPAMTRIX INCLUDES
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <spamtrix_luincpreconditioner.hpp>

const idx npt = 4; //Number of Points per Tetrahedra

double rt2 = sqrt(2.0);
double rt3 = sqrt(3.0);
double rt6 = sqrt(6.0);

void assemble_volume(double *p,
                     SolutionVector *v,
                     SolutionVector *q,
                     const LC &lc,
                     Mesh *mesh,
                     SpaMtrix::IRCMatrix &K,
                     SpaMtrix::Vector &L,
                     Electrodes *electrodes);

void assemble_Neumann(double *p,
                      SolutionVector *v,
                      SolutionVector *q,
                      const LC &lc,
                      Mesh *mesh,
                      Mesh *surf_mesh,
                      SpaMtrix::IRCMatrix &K,
                      SpaMtrix::Vector &L);

void Pot_GMRES(SpaMtrix::IRCMatrix &K,
               SpaMtrix::Vector &B,
               SpaMtrix::Vector &X,
               Settings *settings);

// DECLARATION ONLY
void setUniformEField(Electrodes &electrodes, SolutionVector &v, double *p);

void calcpot3d(
    SpaMtrix::IRCMatrix &K,
    SolutionVector *v,
    SolutionVector *q,
    const LC &lc,
    Geometry &geom,
    Settings *settings,
    Electrodes *electrodes) {
    // First check whether potential calculation is actually needed...
    // NO NEED TO CALCULATE POTENTIAL IF...
    if ((v->getnFixed() == 0) ||        // no fixed potential nodes OR
            (!electrodes->getCalcPot())) {  // no need to calculate potential
        v->setValuesTo(0.0); // if no potential calculation, set all values to zero
        if (electrodes->isEField()) {
            setUniformEField(*electrodes, *v, geom.getPtrTop());
        }
        return;
    }
    K = 0.0; // clears values but keeps sparsity structure
    SpaMtrix::Vector L(v->getnFreeNodes());
    SpaMtrix::Vector V(v->getnFreeNodes());
    // PROVIDE POTENTIAL FROM PREVIOUS STEP AS INITIAL GUESS
    for (idx i = 0 ; i < v->getnDoF(); i++) {
        const idx ind = v->getEquNode(i);
        if (ind != NOT_AN_INDEX) {
            V[ind] = v->getValue(i);
        }
    }
    // Assemble system
    assemble_volume(geom.getPtrTop(), v, q, lc, geom.t, K , L, electrodes);
    assemble_Neumann(geom.getPtrTop() , v , q , lc , geom.t , geom.e , K , L);
//#ifdef DEBUG
//    K->DetectZeroDiagonals();
//#endif
    // GMRES SOLVER
    Pot_GMRES(K, L, V, settings);
    // COPY NON-FIXED VALUES BACK TO SOLUTIONVECTOR
    for (idx i = 0 ; i < v->getnDoF() ; i++) {
        idx ind = v->getEquNode(i);
        // EQU NODES OF FIXED DOFS ARE ALL "NOT_AN_INDEX"
        if (ind != NOT_AN_INDEX) {
            v->setValue(i, 0, V[ind]);
        }
    }// end for i
}//end calcpot3d

inline void localKL(
    double *p,
    int *tt,
    double lK[npt][npt],
    double lL[npt],
    int it,
    Mesh *mesh,
    SolutionVector *q,
    const LC &lc,
    Electrodes *electrodes,
    const Shape4thOrder &shapes)

{
    idx i, j;
    double eper, deleps;
    double S0 = lc.S0();
    eper    = 0;
    deleps  = 0;
    double efe = 0.0, efe2 = 0.0;
    if (mesh->getMaterialNumber(it) == MAT_DOMAIN1) { // if LC element
        eper = lc.eps_per() / S0;
        deleps = (lc.eps_par() - lc.eps_per()) / S0;
        efe  = 2.0 / (9 * S0) * (lc.e1() + 2 * lc.e3());
        efe2 = 4.0 / (9 * S0 * S0) * (lc.e1() - lc.e3());
    } else { // otherwise dielectric
        idx ind_de = mesh->getDielectricNumber(it) - 1; // -1 for 0 indexing
        eper = electrodes->getDielectricPermittivity(ind_de);
    }

    memset(lK, 0, 4 * 4 * sizeof(double));
    memset(lL, 0, 4 * sizeof(double));

    double Jdet = mesh->getDeterminant(it);
    // Jacobian
    double xr, xs, xt, yr, ys, yt, zr, zs, zt;
    xr = xs = xt = yr = ys = yt = zr = zs = zt = 0.0;
    for (i = 0; i < npt; i++) {
        xr += shapes.sh1r[0][i] * p[ tt[i] * 3 + 0 ] * 1e-6;
        xs += shapes.sh1s[0][i] * p[ tt[i] * 3 + 0 ] * 1e-6;
        xt += shapes.sh1t[0][i] * p[ tt[i] * 3 + 0 ] * 1e-6;
        yr += shapes.sh1r[0][i] * p[ tt[i] * 3 + 1 ] * 1e-6;
        ys += shapes.sh1s[0][i] * p[ tt[i] * 3 + 1 ] * 1e-6;
        yt += shapes.sh1t[0][i] * p[ tt[i] * 3 + 1 ] * 1e-6;
        zr += shapes.sh1r[0][i] * p[ tt[i] * 3 + 2 ] * 1e-6;
        zs += shapes.sh1s[0][i] * p[ tt[i] * 3 + 2 ] * 1e-6;
        zt += shapes.sh1t[0][i] * p[ tt[i] * 3 + 2 ] * 1e-6;
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
        if (mesh->getMaterialNumber(it) != MAT_DOMAIN1) { // if this element is not LC
            // only set diagonal permittivities to non-zero
            for (i = 0 ; i < npt ; ++i) {
                e11 += Sh[i] * eper;
                e22 += Sh[i] * eper;
                e33 += Sh[i] * eper;
            }
        } else { // otherwise LC ->
            for (i = 0; i < npt; i++) {
                q1 += Sh[i] * q->getValue(tt[i], 0);
                q2 += Sh[i] * q->getValue(tt[i], 1);
                q3 += Sh[i] * q->getValue(tt[i], 2);
                q4 += Sh[i] * q->getValue(tt[i], 3);
                q5 += Sh[i] * q->getValue(tt[i], 4);
                q1x += dSh[i][0] * q->getValue(tt[i], 0);
                q2x += dSh[i][0] * q->getValue(tt[i], 1);
                q3x += dSh[i][0] * q->getValue(tt[i], 2);
                q4x += dSh[i][0] * q->getValue(tt[i], 3);
                q5x += dSh[i][0] * q->getValue(tt[i], 4);
                q1y += dSh[i][1] * q->getValue(tt[i], 0);
                q2y += dSh[i][1] * q->getValue(tt[i], 1);
                q3y += dSh[i][1] * q->getValue(tt[i], 2);
                q4y += dSh[i][1] * q->getValue(tt[i], 3);
                q5y += dSh[i][1] * q->getValue(tt[i], 4);
                q1z += dSh[i][2] * q->getValue(tt[i], 0);
                q2z += dSh[i][2] * q->getValue(tt[i], 1);
                q3z += dSh[i][2] * q->getValue(tt[i], 2);
                q4z += dSh[i][2] * q->getValue(tt[i], 3);
                q5z += dSh[i][2] * q->getValue(tt[i], 4);
                
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
    double *p,
    idx *tt,
    double lK[npt][npt],
    double lL[npt],
    int it,
    int index_to_Neumann,
    Mesh  *mesh,
    Mesh *surf_mesh,
    SolutionVector *q,
    const LC &lc,
    const ShapeSurf4thOrder &shapes) {
    int i, j;
    double S0 = lc.S0();
    double eper   = 4.5;
    double deleps = 0.0;    // default for dielectric material
    double efe = 0.0, efe2 = 0.0;
    if (mesh->getMaterialNumber(index_to_Neumann) == MAT_DOMAIN1) {
        eper   = lc.eps_per() / S0;
        deleps = (lc.eps_par() - lc.eps_per()) / S0;
        efe  = 2.0 / (9 * S0) * (lc.e1() + 2 * lc.e3());
        efe2 = 4.0 / (9 * S0 * S0) * (lc.e1() - lc.e3());
    } else {
        throw std::runtime_error(fmt::format("Expected DOMAIN1 material in {}, {}.", __FILE__, __func__));
    }
    memset(lK, 0, npt * npt * sizeof(double));
    memset(lL, 0, 4 * sizeof(double));
    double n[3];
    surf_mesh->CopySurfaceNormal(it, n);
    double eDet = surf_mesh->getDeterminant(it);
    double Jdet = mesh->getDeterminant(index_to_Neumann);
#ifdef DEBUG
    assert(eDet > 0);
    assert(Jdet > 0);
    if (abs(n[0]*n[0] + n[1]*n[1] + n[2]*n[2] - 1.0) > 0.01) {
        throw std::runtime_error(fmt::format("Surface normal vector for triangle {} is not of unit length.", it));
    }
#endif

    double xr, xs, xt, yr, ys, yt, zr, zs, zt;
    xr = xs = xt = yr = ys = yt = zr = zs = zt = 0.0;
    for (i = 0; i < 4; i++) {
        xr += shapes.sh1r[0][i] * p[ tt[i] * 3 + 0] * 1e-6; //  <- tt is reordered volume element
        xs += shapes.sh1s[0][i] * p[ tt[i] * 3 + 0] * 1e-6;
        xt += shapes.sh1t[0][i] * p[ tt[i] * 3 + 0] * 1e-6;
        yr += shapes.sh1r[0][i] * p[ tt[i] * 3 + 1] * 1e-6;
        ys += shapes.sh1s[0][i] * p[ tt[i] * 3 + 1] * 1e-6;
        yt += shapes.sh1t[0][i] * p[ tt[i] * 3 + 1] * 1e-6;
        zr += shapes.sh1r[0][i] * p[ tt[i] * 3 + 2] * 1e-6;
        zs += shapes.sh1s[0][i] * p[ tt[i] * 3 + 2] * 1e-6;
        zt += shapes.sh1t[0][i] * p[ tt[i] * 3 + 2] * 1e-6;
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
            q1 += Sh[i] * q->getValue(tt[i], 0);
            q2 += Sh[i] * q->getValue(tt[i], 1);
            q3 += Sh[i] * q->getValue(tt[i], 2);
            q4 += Sh[i] * q->getValue(tt[i], 3);
            q5 += Sh[i] * q->getValue(tt[i], 4);
            q1x += dSh[i][0] * q->getValue(tt[i], 0);
            q2x += dSh[i][0] * q->getValue(tt[i], 1);
            q3x += dSh[i][0] * q->getValue(tt[i], 2);
            q4x += dSh[i][0] * q->getValue(tt[i], 3);
            q5x += dSh[i][0] * q->getValue(tt[i], 4);
            q1y += dSh[i][1] * q->getValue(tt[i], 0);
            q2y += dSh[i][1] * q->getValue(tt[i], 1);
            q3y += dSh[i][1] * q->getValue(tt[i], 2);
            q4y += dSh[i][1] * q->getValue(tt[i], 3);
            q5y += dSh[i][1] * q->getValue(tt[i], 4);
            q1z += dSh[i][2] * q->getValue(tt[i], 0);
            q2z += dSh[i][2] * q->getValue(tt[i], 1);
            q3z += dSh[i][2] * q->getValue(tt[i], 2);
            q4z += dSh[i][2] * q->getValue(tt[i], 3);
            q5z += dSh[i][2] * q->getValue(tt[i], 4);

            e11 += shapes.sh1[igp][i] * (((2.0 / 3.0 / S0) * (-q1 / rt6 + q2 / rt2) + (1.0 / 3.0)) * deleps + eper); //~nx*nx
            e22 += shapes.sh1[igp][i] * (((2.0 / 3.0 / S0) * (-q1 / rt6 - q2 / rt2) + (1.0 / 3.0)) * deleps + eper); //~ny*ny
            e33 += shapes.sh1[igp][i] * (((2.0 / 3.0 / S0) * (2.0 * q1 / rt6)       + (1.0 / 3.0)) * deleps + eper); //~nz*nz
            e12 += shapes.sh1[igp][i] * (2.0 / 3.0 / S0) * (q3 / rt2) * deleps;           //~nx*ny
            e13 += shapes.sh1[igp][i] * (2.0 / 3.0 / S0) * (q5 / rt2) * deleps;           //~nx*nz
            e23 += shapes.sh1[igp][i] * (2.0 / 3.0 / S0) * (q4 / rt2) * deleps;       //~ny*nz
        }//end for i
        double mul = shapes.w[igp] * eDet;
        double nx = n[0], ny = n[1], nz = n[2]; // interior normal?
        for (i = 0; i < 4; i++) {
            const double ShR  = mul * Sh[i];
            const double ShRx = mul * dSh[i][0];
            const double ShRy = mul * dSh[i][1];
            const double ShRz = mul * dSh[i][2];
    
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
                lK[i][j] += mul * Sh[i] * (((e11 - 1) * dSh[j][0] + e12 * dSh[j][1] + e13 * dSh[j][2]) * n[0]
                                           + (e12 * dSh[j][0] + (e22 - 1) * dSh[j][1] + e23 * dSh[j][2]) * n[1]
                                           + (e13 * dSh[j][0] + e23 * dSh[j][1] + (e33 - 1) * dSh[j][2]) * n[2]);
            }//end for j
        }//end for i
    }//end for igp
}
// end void localKL




void assemble_volume(
    double *p,
    SolutionVector *v,
    SolutionVector *q,
    const LC &lc,
    Mesh *mesh,
    SpaMtrix::IRCMatrix &K,
    SpaMtrix::Vector &L,
    Electrodes *electrodes) {
    Shape4thOrder shapes;
#ifndef DEBUG
    #pragma omp parallel for
#endif
    for (idx it = 0; it < mesh->getnElements(); it++) {
        double lK[npt][npt];
        double lL[npt];
        int t[4] = {0, 0, 0, 0};
        t[0] =  mesh->getNode(it, 0);
        t[1] =  mesh->getNode(it, 1);
        t[2] =  mesh->getNode(it, 2);
        t[3] =  mesh->getNode(it, 3);
        localKL(p, t, lK, lL, it, mesh, q, lc, electrodes, shapes);
        //printlK( &lK[0][0] , npt);
        for (idx i = 0; i < npt; i++) { // FOR ROWS
            idx ri = v->getEquNode(t[i]);
            // RHS FIXED NODE HANDLING
            if (ri == NOT_AN_INDEX) {  // IF THIS NODE IS FIXED
                for (int j = 0; j < 4 ; j++) { // SET CONTRIBUTION TO CONNECTED *FREE* NODES
                    idx nc = v->getEquNode(t[j]);   // INDEX TO CONNECTED NODE DEGREE OF FREEDOM POSITION
                    if (nc != NOT_AN_INDEX) {
                        #pragma omp atomic
                        L[ nc ] += lL[i] - lK[i][j] * v->getValue(t[i]); // L = -K*v
                    }
                }
            }// END IF ROW NODE IS FIXED
            if (ri != NOT_AN_INDEX)
                for (idx j = 0; j < npt  ; j++) { // FOR COLUMNS
                    idx rj = v->getEquNode(t[j]);
                    if (rj != NOT_AN_INDEX) {
                        K.sparse_add(rj, ri, lK[j][i]);
                    }
                }//end for j
        }//end for i
    }//end for it
}// end void assemble_volume

////////////////////////////////////////////////////////////////////////////////
// assemble Neumann boundaries
void assemble_Neumann(
    double *p,
    SolutionVector *v,
    SolutionVector *q,
    const LC &lc,
    Mesh *mesh,
    Mesh *surf_mesh,
    SpaMtrix::IRCMatrix &K,
    SpaMtrix::Vector &L) {
    ShapeSurf4thOrder shapes;
    //#pragma omp parallel for
    for (idx it = 0; it < surf_mesh->getnElements(); it++) {
        int index_to_Neumann = surf_mesh->getConnectedVolume(it);
        if ((index_to_Neumann > -1) && (surf_mesh->getMaterialNumber(it) == MAT_NEUMANN)) { // if connected to LC tet
            double lK[4][4];
            double lL[4];
            idx ee[3] = {   surf_mesh->getNode(it, 0) ,
                            surf_mesh->getNode(it, 1) ,
                            surf_mesh->getNode(it, 2)
                        } ;
            idx tt[4] = { tt[0] =   mesh->getNode(index_to_Neumann, 0),
                          mesh->getNode(index_to_Neumann, 1),
                          mesh->getNode(index_to_Neumann, 2),
                          mesh->getNode(index_to_Neumann, 3)
                        };
            int intr = 0;//find  index to internal node
            for (int i = 0; i < 4; i++) {
                if ((tt[i] != ee[0]) && (tt[i] != ee[1]) && (tt[i] != ee[2])) {
                    intr = i;
                    break;
                }
            }
            idx ti[4] = { ee[0], ee[1], ee[2], tt[intr] }; // reordered local element, internal node is always last
            localKL_N(p, &ti[0], lK , lL, it, index_to_Neumann, mesh, surf_mesh, q, lc, shapes);
            for (int i = 0; i < 4; i++) {
                idx ri = v->getEquNode(ti[i]);
                if (ri == NOT_AN_INDEX) { // HANDLE FIXED NODE
                    for (int j = 0; j < 4 ; j++) {
                        idx cr = v->getEquNode(ti[j]);   // CONNECTED NODE DOF ORDER
                        if (cr != NOT_AN_INDEX) {
                            #pragma omp atomic
                            L[cr ] += lL[i] - lK[j][i] * v->getValue(ti[i]);
                        }
                    }
                }// END HANDLE FIXED NODE
                else { // HANDLE FREE NODE
                    for (int j = 0; j < 4; j++) { // FOR COLUMNS
                        idx rj = v->getEquNode(ti[j]);
                        if (rj != NOT_AN_INDEX) {// NON-FIXED NODE
                            K.sparse_add(rj, ri, lK[j][i]);
                        }
                    }//end for j
                }// END HANDLE FREE NODES
            }//end for i
        }//end if LC
    }//end for it
}//end void assemble_Neumann


void setUniformEField(Electrodes &electrodes, SolutionVector &v, double *p) {
    // SETS POTENTIAL VALUES IN V SUCH THA A UNIFORM E-FIELD IS CREATED
    // CENTRE OF STRUCTURE IS ASSUMED TO BE AT 0V, VOLTAGE VALUES FOR NODES
    // IN THE DIRECTION OF THE EFIELD ARE SET ACCORDING TO THEIR DISTANCE
    // TO THE CENTRE
    // GET E-FIELD DIRECTION VECTOR AND MAGNITUDE
    double *E = &(electrodes.EField[0]);
    double Emag = sqrt(E[0] * E[0] + E[1] * E[1] + E[2] * E[2]);
    double Ehat[3] = { E[0] / Emag,
                       E[1] / Emag,
                       E[2] / Emag
                     };
    //1. CALCULATE CENTRE OF STRUCTURE
    int np = v.getnDoF();
    double xmax = 0.0;
    double ymax = 0.0;
    double zmax = 0.0;
    for (int i = 0 ; i < np ; i++) {
        xmax = p[3 * i + 0 ] > xmax ? p[3 * i + 0] : xmax;
        ymax = p[3 * i + 1 ] > ymax ? p[3 * i + 1] : ymax;
        zmax = p[3 * i + 2 ] > zmax ? p[3 * i + 2] : zmax;
    }
    double centre[3] = {xmax / 2.0 , ymax / 2.0 , zmax / 2.0}; // centre of structure
    //2. LOOP OVER EACH NODE AND CACLCULATE ITS DISTANCE TO CENTRE
    for (int i = 0 ; i < np ; i++) {
        double *pos = &p[i * 3]; // shortcut to this node
        double vec[3] = { centre[0] - pos[0],
                          centre[1] - pos[1],
                          centre[2] - pos[2]
                        }; // vector from centre to node i
        // want distance along EField, i.e. dot product
        double dist = vec[0] * Ehat[0] + vec[1] * Ehat[1] + vec[2] * Ehat[2];
        // set potential value as distance*magnitude
        v.setValue(i, 0, dist * Emag + v.getValue(i));
    }
}

/**
 * Solves the Linear simulatenous equation Ax=b using the GMRES method.
 */
void Pot_GMRES(SpaMtrix::IRCMatrix &K,
               SpaMtrix::Vector &B,
               SpaMtrix::Vector &X,
               Settings *settings) {
    idx size = K.getNumRows();
    idx maxiter     = settings->getV_GMRES_Maxiter();
    idx restart     = settings->getV_GMRES_Restart();
    maxiter = maxiter < size ? maxiter : size;
    restart = restart < maxiter ? restart : maxiter;
    double toler    = settings->getV_GMRES_Toler();
    SpaMtrix::LUIncPreconditioner LU(K); // DOES THIS HAVE TO BE RECOMPUTED EACH TIME??
    //SpaMtrix::DiagPreconditioner LU(K);
    SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);
    if (!solver.gmres(K, X, B, LU)) {
        Log::warn("GMRES did not converge in {} iterations when solving for potential. Tolerance achieved is {}.",
                  solver.maxIter, solver.toler);
    }
}
