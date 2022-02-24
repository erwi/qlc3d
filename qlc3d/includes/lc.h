#ifndef LC_H
#define LC_H
#include <stdio.h>

/*!
 * Liquid crystal material parameters class
 */
class LC {
    const double K11_, K22_, K33_;
    double p0_;
    double A_, B_, C_;
    double eps_par_;
    double eps_per_;
    double e1_;
    double e3_;
    double gamma1_;
    // implicit values, calculated from the explicitly set ones
    double L1_, L2_, L3_, L4_, L5_, L6_;
    double S0_;
    double u1_;

    double static calculateS0(double A, double B, double C);
    double static calculateL1(double K11, double K22, double K33, double A, double B, double C);
    double static calculateL2(double K11, double K22, double A, double B, double C);
    double static calculateL3(double K24, double A, double B, double C);
    double static calculateL4(double p0, double K22, double A, double B, double C);
    double static calculateL6(double K11, double K33, double A, double B, double C);
    double static calculateU1(double gamma1, double A, double B, double C);

public:
    LC(double K11, double K22, double K33, double p0, double A, double B, double C, double eps_par, double eps_per,
       double e1, double e3, double gamma1):
    K11_{ K11 }, K22_{ K22 }, K33_{ K33 }, p0_{ p0 }, A_{ A }, B_{ B }, C_{ C }, eps_par_{ eps_par }, eps_per_{ eps_per },
    e1_{ e1 }, e3_{ e3 }, gamma1_{ gamma1 },
    S0_{ LC::calculateS0(A, B, C) },
    L1_{ LC::calculateL1(K11, K22, K33, A, B, C) },
    L2_{ LC::calculateL2(K11, K22, A, B, C) },
    L3_{ 0 },
    L4_{ calculateL4(p0, K22, A, B, C) },
    L5_{ 0 },
    L6_{ calculateL6(K11, K33, A, B, C) },
    u1_{ calculateU1(gamma1, A, B, C) }
    {}

    [[nodiscard]] const double & K11() const { return K11_; }
    [[nodiscard]] const double & K22() const { return K22_; }
    [[nodiscard]] const double & K33() const { return K33_; }
    [[nodiscard]] const double & p0() const { return p0_; }
    [[nodiscard]] const double & A() const { return A_; }
    [[nodiscard]] const double & B() const { return B_; }
    [[nodiscard]] const double & C() const { return C_; }

    [[nodiscard]] const double & eps_par() const { return eps_par_; }
    [[nodiscard]] const double & eps_per() const { return eps_per_; }

    [[nodiscard]] const double & e1() const { return e1_; }
    [[nodiscard]] const double & e3() const { return e3_; }
    [[nodiscard]] const double & gamma1() const { return gamma1_; }

    [[nodiscard]] const double & u1() const { return u1_; }

    [[nodiscard]] const double & L1() const { return L1_; }
    [[nodiscard]] const double & L2() const { return L2_; }
    [[nodiscard]] const double & L3() const { return L3_; }
    [[nodiscard]] const double & L4() const { return L4_; }
    [[nodiscard]] const double & L5() const { return L5_; }
    [[nodiscard]] const double & L6() const { return L6_; }

    [[nodiscard]] const double & S0() const { return S0_; };
};

class LCBuilder {
    double K11_ = 10e-12;
    double K22_ = 10e-12;
    double K33_ = 10e-12;
    double K24_ = 0;
    double p0_ = 0;
    double A_ = -1.2e5;
    double B_ = -2.1333e6;
    double C_ = 1.7333e6;
    double eps_par_ = 18.5;
    double eps_per_ = 7.0;
    double e1_ = 0;
    double e3_ = 0;
    double gamma1_ = 0.0777;

public:
    LCBuilder() { }

    LCBuilder &K11(double K11);
    LCBuilder &K22(double K22);
    LCBuilder &K33(double K33);
    LCBuilder &K24(double K24); // TODO: Actually read this from settings file and set.
    LCBuilder &p0(double p0);
    LCBuilder &A(double A);
    LCBuilder &B(double B);
    LCBuilder &C(double C);
    LCBuilder &eps_par(double eps_par);
    LCBuilder &eps_per(double eps_per);
    LCBuilder &e1(double e1);
    LCBuilder &e3(double e3);
    LCBuilder &gamma1(double gamma1);

    LC* build() const {
        return new LC { K11_, K22_, K33_, p0_, A_, B_, C_, eps_par_, eps_per_, e1_, e3_, gamma1_ };
    }
};



#endif

