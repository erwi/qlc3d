#include <lc.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <reader.h>
LC::LC() {
    // Load default parameter values by default
    A = -1.2e5;
    B = -2.1333e6;
    C =  1.7333e6;
    // elastic coefficients
    K11 = 10e-12;
    K22 = 10e-12;
    K33 = 10e-12;
    p0 = 0.0;
    // permittivity
    eps_par = 18.5;
    eps_per = 7.0;
    //flexoelectricity
    e11 = 0;
    e33 = 0;
    // rotational viscosity
    gamma1 = 0.0777;
    convert_params_n2Q(); // converts from "vector" to "tensor" parameter values
}

void LC::convert_params_n2Q() {
    /*!converts paraemters defined in Oseen-Frank vecotr model values
       to Q-tensor values that take into the account order parameter*/
    S0 = (-B+sqrt(B*B-24*A*C))/(6*C);
    if ( (S0<0) || (S0 > 1) ) {
        std::cerr << "ERROR!, unusual equilibrium order parameter value: "<< S0 << std::endl;
        std::cerr	<< "Check Thermotropic coefficients, A: " << A
                    << " B: " << B
                    << " C: " << C << std::endl;
        std::exit(1);
    }
    // WHERE DOES EXTRA FACTOR OF 2 COME IN L-TERMS?
    L1=2.0*(K33-K11+3.0*K22)/(S0*S0*27.0);
    L2=4.0*(K11-K22)/(9.0*S0*S0);
    L3 = 0;
    L4 = 0.0;
    if (p0!=0) {
        double q0 = 2*3.14159265/p0;
        L4 = ( 8.0 * q0 * K22 ) / ( S0 * S0 * 9.0 );
    }
    L5 = 0;
    L6=4.0*(K33-K11)/(S0*S0*S0*27.0);
    // VISCOSITIES
    u1 = 4*gamma1 / (9 * S0);
}

double LC::getS0() {
    return S0;
}

