#include <lc.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <reader.h>

//<editor-fold desc=LC>
double LC::calculateL1(double K11, double K22, double K33, double A, double B, double C) {
    double S0 = LC::calculateS0(A, B, C);
    return 2.0 * (K33 - K11 + 3.0 * K22) / (S0 * S0 * 27.0);
}

double LC::calculateL2(double K11, double K22, double A, double B, double C) {
    double S0 = LC::calculateS0(A, B, C);
    return 4.0 * (K11 - K22) / (9.0 * S0 * S0);
}

double LC::calculateL4(double p0, double K22, double A, double B, double C) {
    if (p0 == 0.) {
        return 0;
    }
    double S0 = LC::calculateS0(A, B, C);
    double q0 = 2 * M_PI / p0;
    return (8.0 * q0 * K22) / (S0 * S0 * 9.0);
}

double LC::calculateL6(double K11, double K33, double A, double B, double C) {
    double S0 = LC::calculateS0(A, B, C);
    return 4.0 * (K33 - K11) / (S0 * S0 * S0 * 27.0);
}

double LC::calculateU1(double gamma1, double A, double B, double C) {
    double S0 = LC::calculateS0(A, B, C);
    return 2 * gamma1 / (9 * S0 * S0);
}

double LC::calculateS0(double A, double B, double C) {
    double S0 = (-B + sqrt(B * B - 24 * A * C)) / (6 * C);
    if (S0 < 0 || S0 > 1 || std::isnan(S0)) {
        std::string msg = "Unwonted value of S0=" + std::to_string(S0) + ". Expected value in range 0 to 1.";
        msg += "Check thermotropic coefficients A, B, C=" + std::to_string(A)
               + ", " + std::to_string(B) + ", " + std::to_string(C);
        throw std::invalid_argument(msg);
    }
    return S0;
}
//</editor-fold>

//<editor-fold desc=LCBuilder>
LCBuilder &LCBuilder::K11(double K11) {
    K11_ = K11;
    return *this;
}

LCBuilder &LCBuilder::K22(double K22) {
    K22_ = K22;
    return *this;
}

LCBuilder &LCBuilder::K33(double K33) {
    K33_ = K33;
    return *this;
}

LCBuilder &LCBuilder::p0(double p0) {
    p0_ = p0;
    return *this;
}

LCBuilder &LCBuilder::A(double A) {
    A_ = A;
    return *this;
}

LCBuilder &LCBuilder::B(double B) {
    B_ = B;
    return *this;
}

LCBuilder &LCBuilder::C(double C) {
    C_ = C;
    return *this;
}

LCBuilder &LCBuilder::eps_par(double eps_par) {
    eps_par_ = eps_par;
    return *this;
}

LCBuilder &LCBuilder::eps_per(double eps_per) {
    eps_per_ = eps_per;
    return *this;
}

LCBuilder &LCBuilder::e1(double e1) {
    e1_ = e1;
    return *this;
}

LCBuilder &LCBuilder::e3(double e3) {
    e3_ = e3;
    return *this;
}

LCBuilder &LCBuilder::gamma1(double gamma1) {
    gamma1_ = gamma1;
    return *this;
}
//</editor-fold>
