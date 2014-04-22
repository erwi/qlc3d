#ifndef LC_H
#define LC_H
#include <stdio.h>
class Reader; // forward declaration
class LC {
    /*!Liquid crystal material parameters*/
    double gamma1;
    void convert_params_n2Q();
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

    LC();
    void printLC();
    void WriteLC(FILE *fid);

    double getS0();
    void readSettingsFile(Reader &reader);
};

#endif

