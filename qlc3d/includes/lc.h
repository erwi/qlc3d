#ifndef LC_H
#define LC_H
#include <stdio.h>

class LC {
    /*!Liquid crystal material parameters class*/
public:
    // Decalre default material parameters
    const static double DEFAULT_K11;
    const static double DEFAULT_K22;
    const static double DEFAULT_K33;
    const static double DEFAULT_P0;
    const static double DEFAULT_A;
    const static double DEFAULT_B;
    const static double DEFAULT_C;
    const static double DEFAULT_EPS_PAR;
    const static double DEFAULT_EPS_PER;
    const static double DEFAULT_E1;
    const static double DEFAULT_E3;
    const static double DEFAULT_GAMMA1;
    // elastic coefficients
    double K11, K22, K33;
    double L1, L2, L3, L4, L5, L6;
    double p0;
    // thermotropic coefficients
    double A, B, C;
    // equilibrium order parameter
    double S0;
    // permittivity
    double eps_par;
    double eps_per;
    //flexoelectricity
    double e11;
    double e33;
    //viscosities
    double u1;

    LC();
    double getS0();
    double gamma1;
    void convert_params_n2Q();
};

#endif

