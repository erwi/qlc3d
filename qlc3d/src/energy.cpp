#include <energy.h>
#include <cstdio>
#include <iostream>
#include <math.h>
#include <lc.h>
#include <geometry.h>
#include <solutionvector.h>
#include <shapefunction3d.h>
namespace Energy {
#define SIGN(x)     (x)>=0 ? 1:-1


// Gauss integration
// ---------------------------------------------------------
//     3D Gauss-Legendre weights for N = 11, D = 4
// ---------------------------------------------------------
const int ngp = 11;
const double a = (1 + sqrt(5.0 / 14.0)) / 4.0;
const double b = (1 - sqrt(5.0 / 14.0)) / 4.0;

static double gp[ngp][4] = {
    {0.25           , 0.25      ,   0.25        , 0.25},
    {11.0 / 14.0      , 1.0 / 14.0  ,   1.0 / 14.0    , 1.0 / 14.0},
    {1.0 / 14.0      ,    11.0 / 14.0   ,   1.0 / 14.0    , 1.0 / 14.0},
    {1.0 / 14.0     , 1.0 / 14.0    ,   11.0 / 14.0   , 1.0 / 14.0},
    {1.0 / 14.0     , 1.0 / 14.0    ,   1.0 / 14.0    , 11.0 / 14.0},
    {a        , a       ,   b       , b},
    {a        , b       ,   a       , b},
    {a         , b       ,   b      , a},
    {b        , a       ,   a       , b},
    {b        , a       ,   b       , a},
    {b        , b       ,   a       , a}
};

const double w11 = -74.0 / 5625.0;
const double w12 = 343.0 / 45000.0;
const double w13 = 56.0 / 2250.0;
static double w[ngp] = {w11, w12, w12, w12, w12, w13, w13, w13, w13, w13, w13};

static double sh1[ngp][4]; // P1 Shape functions
static double sh1r[ngp][4]; // P1 Shape functions r-derivatives
static double sh1s[ngp][4]; // P1 Shape functions s-derivatives
static double sh1t[ngp][4]; //P1 shape functions t-derivative

static double rt2 = sqrt(2.0);
static double rt3 = sqrt(3.0);
static double rt6 = sqrt(6.0);

//static double rt3 = sqrt(3.0);

void init_shape() {
    for (int i = 0; i < ngp; i++) {
        // P1 Shape functions
        sh1[i][0] = 1 - gp[i][0] - gp[i][1] - gp[i][2];
        sh1[i][1] = gp[i][0];
        sh1[i][2] = gp[i][1];
        sh1[i][3] = gp[i][2];
        // P1 Shape functions r-derivatives
        sh1r[i][0] = -1.0;
        sh1r[i][1] = 1.0;
        sh1r[i][2] = 0.0;
        sh1r[i][3] = 0.0;
        // P1 Shape functions s-derivatives
        sh1s[i][0] = -1.0;
        sh1s[i][1] = 0.0;
        sh1s[i][2] = 1.0;
        sh1s[i][3] = 0.0;
        // P1 Shape functions t-derivatives
        sh1t[i][0] = -1.0;
        sh1t[i][1] = 0.0;
        sh1t[i][2] = 0.0;
        sh1t[i][3] = 1.0;
    }
}// enf void init_shape

}// end namespace Energy


void CalculateFreeEnergy(FILE *fid,
                         int currentIteration,
                         double currentTime,
                         const LC &lc,
                         Geometry *geom,
                         SolutionVector *v,
                         SolutionVector *q) {
    using std::cout;
    using std::endl;
    cout << " Calculating free energy..."; fflush(stdout);
    // IF FIRST ITERATION, PRINT HEADERS
    if (currentIteration == 1) {
        fprintf(fid, "%% columns are:\n");
        fprintf(fid, "%% time[s],splay,twist,bend,thermotropic,dielectric\n");
        fprintf(fid, "F = [ ...\n");
    }
    //
    //
    //  BULK ENERGY
    //
    //
    using namespace Energy;
    init_shape();
    double e0 = 8.8541878176 * 1e-12;
    double S0 = lc.S0();
    double epsav = lc.eps_per() / S0;
    double deleps = (lc.eps_par() - lc.eps_per()) / S0 ;
    double A = lc.A();
    double B = lc.B();
    double C = lc.C();
    double efe  = 2.0 / (9 * S0) * (lc.e1() + 2 * lc.e3());
    double efe2 = 4.0 / (9 * S0 * S0) * (lc.e1() - lc.e3());
    double f0 = (3.0 * A / 4.0) * (S0 * S0) + (B / 4.0) * (S0 * S0 * S0) + (9.0 * C / 16.0) * (S0 * S0 * S0 * S0);
    double pi = 3.141569;
    double q0(0);           // CHIRALITY
    if (lc.p0() > 0.0) {
        q0 = 2 * pi / lc.p0() ;
    }
    // IF SINGLE ELASTIC COEFFICIENT
    //double L1(0);
    //if ( (lc->K11 == lc->K22) && (lc->K11 == lc->K33 ) )
    //    L1 = 2.0*(lc->K33-lc->K11+3.0*lc->K22)/(S0*S0*27.0);
    double *p = geom->getPtrTop();
    // energy variables
    double Fe = 0;  // electric energy
    double Fflx = 0;    // flexoelectric enery
    double Fth  = 0;    // thermotropic energy
    double F11 = 0; // SPLAY, TWIST AND BEND ENERGIES
    double F22 = 0;
    double F33 = 0;
    //loop over each element and calculate elastic energy contribution
    for (idx x = 0 ; x < geom->t->getnElements() ; x++) {
        if (geom->t->getMaterialNumber(x) <= MAT_DOMAIN7) { //IF LC ELEMENT
            idx tt[4] = {geom->t->getNode(x, 0),
                         geom->t->getNode(x, 1),
                         geom->t->getNode(x, 2),
                         geom->t->getNode(x, 3)
                        };
            for (int igp = 0; igp < ngp; igp++) { //loop over each gauss point
                double q1(0), q2(0), q3(0), q4(0), q5(0);
                double q1x = 0, q2x = 0, q3x = 0, q4x = 0, q5x = 0;
                double q1y = 0, q2y = 0, q3y = 0, q4y = 0, q5y = 0;
                double q1z = 0, q2z = 0, q3z = 0, q4z = 0, q5z = 0;
                double Vx = 0, Vy = 0, Vz = 0;
                // Inverse Jacobian
                double xr, xs, xt, yr, ys, yt, zr, zs, zt;
                xr = xs = xt = yr = ys = yt = zr = zs = zt = 0.0;
                for (int i = 0; i < 4 ; i++) {
                    xr += sh1r[igp][i] * p[tt[i] * 3 + 0] * 1e-6;
                    xs += sh1s[igp][i] * p[tt[i] * 3 + 0] * 1e-6;
                    xt += sh1t[igp][i] * p[tt[i] * 3 + 0] * 1e-6;
                    yr += sh1r[igp][i] * p[tt[i] * 3 + 1] * 1e-6;
                    ys += sh1s[igp][i] * p[tt[i] * 3 + 1] * 1e-6;
                    yt += sh1t[igp][i] * p[tt[i] * 3 + 1] * 1e-6;
                    zr += sh1r[igp][i] * p[tt[i] * 3 + 2] * 1e-6;
                    zs += sh1s[igp][i] * p[tt[i] * 3 + 2] * 1e-6;
                    zt += sh1t[igp][i] * p[tt[i] * 3 + 2] * 1e-6;
                }
                //----------------
                // Jacobian
                double Jdet = geom->t->getDeterminant(x);
                double Jinv[3][3] = {
                    {(zt * ys - yt * zs) / Jdet, (xt * zs - zt * xs) / Jdet, (xs * yt - ys * xt) / Jdet},
                    {(yt * zr - zt * yr) / Jdet, (zt * xr - xt * zr) / Jdet, (xt * yr - yt * xr) / Jdet},
                    {(yr * zs - ys * zr) / Jdet, (xs * zr - xr * zs) / Jdet, (ys * xr - xs * yr) / Jdet}
                };
                // Shape function derivatives
                double dSh[4][3];
                for (int i = 0 ; i < 4 ; i++) {
                    dSh[i][0] = sh1r[igp][i] * Jinv[0][0] + sh1s[igp][i] * Jinv[1][0] + sh1t[igp][i] * Jinv[2][0];
                    dSh[i][1] = sh1r[igp][i] * Jinv[0][1] + sh1s[igp][i] * Jinv[1][1] + sh1t[igp][i] * Jinv[2][1];
                    dSh[i][2] = sh1r[igp][i] * Jinv[0][2] + sh1s[igp][i] * Jinv[1][2] + sh1t[igp][i] * Jinv[2][2];
                }
                for (int i = 0; i < 4 ; i++) {
                    q1 += sh1[igp][i] * q->getValue(geom->t->getNode(x, i) , 0) ;            // q1i * Ni  = A1
                    q2 += sh1[igp][i] * q->getValue(geom->t->getNode(x, i) , 1) ;
                    q3 += sh1[igp][i] * q->getValue(geom->t->getNode(x, i) , 2) ;
                    q4 += sh1[igp][i] * q->getValue(geom->t->getNode(x, i) , 3) ;
                    q5 += sh1[igp][i] * q->getValue(geom->t->getNode(x, i) , 4) ;
                    q1x += dSh[i][0] * q->getValue(geom->t->getNode(x, i) , 0);
                    q2x += dSh[i][0] * q->getValue(geom->t->getNode(x, i) , 1);
                    q3x += dSh[i][0] * q->getValue(geom->t->getNode(x, i) , 2);
                    q4x += dSh[i][0] * q->getValue(geom->t->getNode(x, i) , 3);
                    q5x += dSh[i][0] * q->getValue(geom->t->getNode(x, i) , 4);
                    q1y += dSh[i][1] * q->getValue(geom->t->getNode(x, i) , 0);
                    q2y += dSh[i][1] * q->getValue(geom->t->getNode(x, i) , 1);
                    q3y += dSh[i][1] * q->getValue(geom->t->getNode(x, i) , 2);
                    q4y += dSh[i][1] * q->getValue(geom->t->getNode(x, i) , 3);
                    q5y += dSh[i][1] * q->getValue(geom->t->getNode(x, i) , 4);
                    q1z += dSh[i][2] * q->getValue(geom->t->getNode(x, i) , 0);
                    q2z += dSh[i][2] * q->getValue(geom->t->getNode(x, i) , 1);
                    q3z += dSh[i][2] * q->getValue(geom->t->getNode(x, i) , 2);
                    q4z += dSh[i][2] * q->getValue(geom->t->getNode(x, i) , 3);
                    q5z += dSh[i][2] * q->getValue(geom->t->getNode(x, i) , 4);
                    // electric fields
                    Vx += dSh[i][0] * v->getValue(geom->t->getNode(x, i));
                    Vy += dSh[i][1] * v->getValue(geom->t->getNode(x, i));
                    Vz += dSh[i][2] * v->getValue(geom->t->getNode(x, i));
                }//end for i
                double R = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4 + q5 * q5;
                double mul = w[igp] * Jdet;
                // IF 3 ELASTIC COEFFICIENTS
                // CALCULATE TWIST ENERGY: 1/2*K22(n.curl(n) - q0 )^2
                // G4 = (9*S^2 / 4 )* n.curl(n)
                double G4 = (q2 * q4x - q4 * q2x - q3 * q5x + q5 * q3x + q2 * q5y + q3 * q4y - q4 * q3y - q5 * q2y - 2 * q2 * q3z + 2 * q3 * q2z + q4 * q5z - q5 * q4z) * 0.5
                            + (3 * q1 * q4x - 3 * q4 * q1x - 3 * q1 * q5y + 3 * q5 * q1y) / (rt2 * rt6);
                // chiral twist offset, scaled by S^2
                double aa = 3.0 / 2.0 * R;
                double F_twist = (G4 - aa * q0) / (9.0 * S0 * S0) * 4; // DIVIDE OUT (9S^2/4) term
                F_twist *= F_twist; // (twist^2)
                double G1 = 1.0 / (rt6 * rt6) * (q1x * q1x + q1y * q1y + q1z * q1z) * 6.0 + 1.0 / (rt2 * rt2) * ((q2x * q2x) * 2.0 + (q3x * q3x) * 2.0 + (q4x * q4x) * 2.0 + (q5x * q5x) * 2.0 + (q2y * q2y) * 2.0 + (q3y * q3y) * 2.0 + (q4y * q4y) * 2.0 + (q5y * q5y) * 2.0 + (q2z * q2z) * 2.0 + (q3z * q3z) * 2.0 + (q4z * q4z) * 2.0 + (q5z * q5z) * 2.0);
                double G2 = 1.0 / (rt6 * rt6) * (q1x * q1x + q1y * q1y + (q1z * q1z) * 4.0) + 1.0 / (rt2 * rt2) * (q2x * q3y * 2.0 - q3x * q2y * 2.0 + q5x * q4y * 2.0 + q2x * q5z * 2.0 + q3x * q4z * 2.0 - q2y * q4z * 2.0 + q3y * q5z * 2.0 + q2x * q2x + q3x * q3x + q5x * q5x + q2y * q2y + q3y * q3y + q4y * q4y + q4z * q4z + q5z * q5z) - (q1x * q2x * 2.0 + q1x * q3y * 2.0 + q3x * q1y * 2.0 - q1y * q2y * 2.0 + q1x * q5z * 2.0 - q5x * q1z * 4.0 + q1y * q4z * 2.0 - q4y * q1z * 4.0) / (rt2 * rt6);
                double G6 = 1.0 / (rt2 * rt2 * rt2) * (q2 * (q2x * q2x) * 2.0 + q2 * (q3x * q3x) * 2.0 + q2 * (q4x * q4x) * 2.0 + q2 * (q5x * q5x) * 2.0 - q2 * (q2y * q2y) * 2.0 - q2 * (q3y * q3y) * 2.0 - q2 * (q4y * q4y) * 2.0 - q2 * (q5y * q5y) * 2.0 + q3 * q2x * q2y * 4.0 + q3 * q3x * q3y * 4.0 + q3 * q4x * q4y * 4.0 + q3 * q5x * q5y * 4.0 + q5 * q2x * q2z * 4.0 + q5 * q3x * q3z * 4.0 + q5 * q4x * q4z * 4.0 + q5 * q5x * q5z * 4.0 + q4 * q2y * q2z * 4.0 + q4 * q3y * q3z * 4.0 + q4 * q4y * q4z * 4.0 + q4 * q5y * q5z * 4.0) + (1.0 / (rt6 * rt6) * (q2 * (q1x * q1x) - q2 * (q1y * q1y) + q3 * q1x * q1y * 2.0 + q5 * q1x * q1z * 2.0 + q4 * q1y * q1z * 2.0) * 6.0) / rt2 - q1 * 1.0 / (rt6 * rt6 * rt6) * (q1x * q1x + q1y * q1y - (q1z * q1z) * 2.0) * 6.0 - (q1 * 1.0 / (rt2 * rt2) * (q2x * q2x + q3x * q3x + q4x * q4x + q5x * q5x + q2y * q2y + q3y * q3y + q4y * q4y + q5y * q5y - (q2z * q2z) * 2.0 - (q3z * q3z) * 2.0 - (q4z * q4z) * 2.0 - (q5z * q5z) * 2.0) * 2.0) / rt6;
                double F_splay =  4 * G2 / (9 * S0 * S0) - 2 * G1 / (27 * S0 * S0) - 4 * G6 / (27 * S0 * S0 * S0);
                double F_bend  =   2 * G1 / (27.0 * S0 * S0) + 4 * G6 / (27 * S0 * S0 * S0);
                F11 += 0.5 * lc.K11() * mul * F_splay;
                F22 += 0.5 * lc.K22() * mul * F_twist;
                F33 += 0.5 * lc.K33() * mul * F_bend;
                double Fel_elem = e0 * (-Vx * Vx - Vy * Vy - Vz * Vz) * epsav * 0.5 +
                                  e0 * deleps * (Vx * Vx * q1 * rt6 / 12.0  - Vx * Vx * q2 * rt2 / 4.0 - Vx * Vy * q3 * rt2 / 2.0  -  Vx * Vz * q5 * rt2 / 2.0
                                                 - Vy * Vz * q4 * rt2 / 2.0   - Vz * Vz * q1 * rt6 / 6.0 + Vy * Vy * q2 / 4.0);
                Fe += mul * Fel_elem;
                if (efe != 0) { // if flexoelectric terms (actually from -f_E)
                    // N.B. this is identical to old version, just simplified
                    double Fflexo = efe*( (2*Vz*q1z - Vx*q1x - Vy*q1y)/rt6
                                  + (Vx*(q2x+q3y+q5z) + Vy*(-q2y+q3x+q4z) + Vz*(q4y+q5x))/rt2 );
                    Fflx += mul * Fflexo;
                }// end if flexoelectric terms
                if (efe2 != 0) { // if flexoelectric terms (actually from -f_E)
                    double Fflexo = efe2*( ( Vx*(-q1*(q2x+q3y+q5z) - q2*q1x - q3*q1y + 2*q5*q1z)
                        + Vy*(q1*(q2y-q3x-q4z) + q2*q1y - q3*q1x + 2*q4*q1z)
                        + Vz*(2*q1*(q4y+q5x) - q4*q1y - q5*q1x) )*rt3/6
                        + Vx*(q2*(q2x+q3y+q5z) + q3*(-q2y+q3x+q4z) + q5*(q4y+q5x))/2
                        + (Vy*q2-Vz*q4)*(q2y-q3x-q4z)/2 + (Vy*q3+Vz*q5)*(q2x+q3y+q5z)/2
                        + Vy*q4*(q4y+q5x)/2 + q1*(Vx*q1x + 4*Vz*q1z + Vy*q1y)/6 );
                    Fflx += mul * Fflexo;
                }// end if flexoelectric terms
                double Fth_elem = A * (R) / 2.0 +
                                  B * (q5 * q5 * q1 * rt6 / 4.0 - q1 * rt6 * q2 * q2 / 2.0
                                       - q3 * q3 * q1 * rt6 / 2.0 + 3.0 / 4.0 * q5 * q5 * q2 * rt2 + 3.0 / 2.0 * q3 * rt2 * q5 * q4
                                       + q4 * q4 * q1 * rt6 / 4.0 - 3.0 / 4.0 * q4 * q4 * q2 * rt2 + q1 * q1 * q1 * rt6 / 6.0) / 3.0
                                  + C * (R * R) / 4.0;
                Fth += mul * (Fth_elem - f0);
            }//end for igp, loop through gauss points
        }//end if domain1 element
    }//end for loop through each element
    //
    //  END BULK ENERGY
    //
    //if (L1){
    //    fprintf(fid,"%e\t%e\t%e\t%e;\n", simu->getCurrentTime(), Fd, Fth, Fe);
    //}
    //else{
    fprintf(fid, "%e\t%e\t%e\t%e\t%e\t%e;\n", currentTime, F11, F22, F33 , Fth , Fe);
    //}
    printf("OK\n");
}

void closeEnergyFile(FILE *fid, Simu &simu) {
    if (simu.getOutputEnergy() == 1) {
        fprintf(fid, "];");
        fclose(fid);
    }
}
