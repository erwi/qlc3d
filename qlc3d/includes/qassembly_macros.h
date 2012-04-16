#ifndef QASSEMBLY_MACROS_H
#define QASSEMBLY_MACROS_H

// MACROS FOR ADDING THERMOTROPIC TERMS TO LOCAL MATRIX OUTSIDE OF
// NESTED ROW COLUMN i/j ASSEMBLY LOOPS

// ADDS DIAGONAL ENTRIES TO LOCAL MATRIX
// a = OFFSET (0,4,8,12,16)
#define THERMOdiag(a,Taa)  { Taa*=mul;\
    for (int i = 0 ; i < 4 ; ++i){\
    double t1 = Taa*Sh[i];\
        for (int j = 0 ; j < 4 ; ++j){\
            lK[(a)+i][(a)+j]+= Sh[j]*t1;\
        }\
    }\
}

// ADDS OFF-DIAGONAL ENTRIES TO LOCAL MATRIX
// a AND b ARE OFFSETS (0,4,8,12,16)
#define THERMO(a,b,Tab) {Tab*=mul;\
    for (int i = 0 ; i < 4 ; ++i){\
        double t1 = Tab*Sh[i];\
        for (int j = 0 ; j < 4 ; ++j){\
            double t2 = t1*Sh[j];\
            lK[(a)+i][(b)+j] += t2;\
            lK[(b)+i][(a)+j] += t2;\
        }\
    }\
}

// ADDS R.H.S BULK TERMS TO LOCAL ELEMENT R.H.S VECTOR
#define ADD_RHS_BULK_TERMS(offset, T) {\
    T*=mul;\
    for (int i = 0 ; i < 4 ; ++i){\
        lL[(offset)+i] += T*Sh[i];\
    }\
}


//--------------------------------------------------
//        THERMOTROPIC ENERGY COMPONENTS
//--------------------------------------------------

//  5 R.H.S. VECTOR COMPONENTS RHS_THERMO1 -> RHS_THERMO5
#define RHS_THERMO1 (A*q1    +   D3*0.5*B*(q5*q5*rt6*0.5    -   q3*q3*rt6   -   rt6*q2*q2   +   q1*q1*rt6   +   q4*q4*rt6*0.5)  +   C*R*q1)
#define RHS_THERMO2 (A*q2    +   D3*B*(0.75*q5*q5*rt2    -   q1*rt6*q2   -   0.75*q4*q4*rt2)   +   C*R*q2)
#define RHS_THERMO3 (A*q3    +   D3*B*(-q3*q1*rt6   +   1.5*rt2*q5*q4 )     +   C*R*q3)
#define RHS_THERMO4 (A*q4    +   D3*0.5*B*(3.0*q3*rt2*q5    +   q4*q1*rt6   -   3.0*q4*q2*rt2)   +   C*R*q4)
#define RHS_THERMO5 (A*q5    +   D3*0.5*B*(q5*q1*rt6    +   3.0*q5*q2*rt2   +   3.0*q3*rt2*q4)   +   C*R*q5)

// MATRIX THERMOTROPIC TERMS
#define MATRIX_THERMO11 (A  +   D3*B*q1*rt6 +   2.0*C*q1*q1 + C*R)
#define MATRIX_THERMO12 (-D3*B*rt6*q2   +       2.0*C*q2*q1)
#define MATRIX_THERMO13 (-B*q3*rt6*D3   +       2.0*C*q3*q1)
#define MATRIX_THERMO14 (B*q4*rt6*D6    +       2.0*C*q4*q1)
#define MATRIX_THERMO15 (B*q5*rt6*D6    +       2.0*C*q5*q1)

#define MATRIX_THERMO22 (A  -   D3*B*q1*rt6 +   2.0*C*q2*q2+C*R)
#define MATRIX_THERMO23 (2.0*C*q3*q2)
#define MATRIX_THERMO24 (-B*q4*rt2*D2   +   2.0*C*q4*q2)
#define MATRIX_THERMO25 (B*q5*rt2*D2    +   2.0*C*q5*q2)

#define MATRIX_THERMO33 (A  -   B*q1*rt6*D3 +   2.0*C*q3*q3+C*R)
#define MATRIX_THERMO34 (B*q5*rt2*D2    +   2.0*C*q4*q3)
#define MATRIX_THERMO35 (B*q4*rt2*D2    +   2.0*C*q5*q3)

#define MATRIX_THERMO44 (A  +   B*(q1*rt6-3.0*q2*rt2)*D6    +   2.0*C*q4*q4+C*R)
#define MATRIX_THERMO45 (B*q3*rt2*D2    +   2.0*C*q5*q4)

#define MATRIX_THERMO55 (A  +   B*(q1*rt6+3.0*q2*rt2)*D6    +   2.0*C*q5*q5+C*R)


#endif // QASSEMBLY_MACROS_H
