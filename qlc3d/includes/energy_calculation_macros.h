#ifndef ENERGY_CALCULATION_MACROS_H
#define ENERGY_CALCULATION_MACROS_H



#define ENERGY_ELASTIC_2K(E) {\
    double t2 = 1.0/(rt2*rt2);\
    double t3 = 1.0/rt2;\
    double t4 = 1.0/rt6;\
    double t5 = 1.0/(rt6*rt6);\
    double t6 = q2x*q2x;\
    double t7 = t2*t6*(1.0/2.0);\
    double t8 = q3x*q3x;\
    double t9 = t2*t8*(1.0/2.0);\
    double t10 = q1x*q1x;\
    double t11 = q5x*q5x;\
    double t12 = t2*t11*(1.0/2.0);\
    double t13 = q2y*q2y;\
    double t14 = t2*t13*(1.0/2.0);\
    double t15 = q3y*q3y;\
    double t16 = t2*t15*(1.0/2.0);\
    double t17 = q4y*q4y;\
    double t18 = t2*t17*(1.0/2.0);\
    double t19 = q1y*q1y;\
    double t20 = q4z*q4z;\
    double t21 = t2*t20*(1.0/2.0);\
    double t22 = q1z*q1z;\
    double t23 = q5z*q5z;\
    double t24 = t2*t23*(1.0/2.0);\
    double t25 = q1x*q2x*t3*t4;\
    E = K1*(t7+t9+t12+t14+t16+t18+t21+t24-t25+t5*t10*(1.0/2.0)+t5*t19*(1.0/2.0)+t5*t22*2.0+q1y*q2y*t3*t4)+K2*(t7+t9+t12+t14+t16+t18+t21+t24+t25+t5*t10*(5.0/2.0)+t5*t19*(5.0/2.0)+t5*t22+(q4x*q4x)*t2+(q5y*q5y)*t2+(q2z*q2z)*t2+(q3z*q3z)*t2+q2x*q3y*t2-q3x*q2y*t2-q4x*q5y*t2-q4x*q3z*t2-q5x*q2z*t2+q4y*q2z*t2-q5y*q3z*t2+q1x*q3y*t3*t4+q3x*q1y*t3*t4-q1y*q2y*t3*t4-q1x*q5z*t3*t4*2.0+q5x*q1z*t3*t4-q1y*q4z*t3*t4*2.0+q4y*q1z*t3*t4)+K2*(q0*q0)*((q2*q2)*t2*4.0+(q3*q3)*t2*4.0+(q1*q1)*t5*1.2E1+(q4*q4)*t2*4.0+(q5*q5)*t2*4.0)+K2*q0*(q2*q4x*t2*2.0-q4*q2x*t2*2.0-q3*q5x*t2*2.0+q5*q3x*t2*2.0+q2*q5y*t2*2.0+q3*q4y*t2*2.0-q4*q3y*t2*2.0-q5*q2y*t2*2.0-q2*q3z*t2*4.0+q3*q2z*t2*4.0+q4*q5z*t2*2.0-q5*q4z*t2*2.0+q1*q4x*t3*t4*6.0-q4*q1x*t3*t4*6.0-q1*q5y*t3*t4*6.0+q5*q1y*t3*t4*6.0);\
}\



#endif // ENERGY_CALCULATION_MACROS_H
