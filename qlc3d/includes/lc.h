#ifndef LC_H
#define LC_H
#include <stdio.h>
#include <string>
#include <algorithm>
class LC {
    /*!Liquid crystal material parameters*/
public:
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
    double u2;
    double gamma1;
    double gamma2;
    double alpha1;
    double alpha4;
    double alpha5;
    double alpha6;
    //
    LC();
    void printLC();
    void WriteLC(FILE *fid);
    void convert_params_n2Q();
    double getS0();
};

#endif

