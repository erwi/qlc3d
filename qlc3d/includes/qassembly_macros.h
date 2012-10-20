#ifndef QASSEMBLY_MACROS_H
#define QASSEMBLY_MACROS_H

#define EPS0 8.8541878176e-12
// MACROS FOR ADDING THERMOTROPIC TERMS TO LOCAL MATRIX OUTSIDE OF
// NESTED ROW COLUMN i/j ASSEMBLY LOOPS

// ADDS DIAGONAL ENTRIES TO LOCAL MATRIX
// a = OFFSET (0,4,8,12,16)
/*
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
*/
// ADDS R.H.S BULK TERMS TO LOCAL ELEMENT R.H.S VECTOR
#define ADD_RHS_BULK_TERMS(lL, L){\
for( int i = 0 ; i < 4 ; ++i){\
lL[i + 0] += L[0]*mul*Sh[i];\
lL[i + 4] += L[1]*mul*Sh[i];\
lL[i + 8] += L[2]*mul*Sh[i];\
lL[i +12] += L[3]*mul*Sh[i];\
lL[i +16] += L[4]*mul*Sh[i];}\
}


//--------------------------------------------------
//        THERMOTROPIC ENERGY COMPONENTS
//--------------------------------------------------

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

//==============================================
// THERMOTROPIC ENERGY MATRIX ASSEMBLY
//==============================================
#define MATRIX_THERMOTROPIC(K){\
    const double D2 = 0.5;\
    const double D3 = 0.33333333333333333333333333333;\
    const double D6 = 0.16666666666666666666666666667;\
\
const double t11 = MATRIX_THERMO11*mul;\
const double t22 = MATRIX_THERMO22*mul;\
const double t33 = MATRIX_THERMO33*mul;\
const double t44 = MATRIX_THERMO44*mul;\
const double t55 = MATRIX_THERMO55*mul;\
const double t12 = MATRIX_THERMO12*mul;\
const double t13 = MATRIX_THERMO13*mul;\
const double t14 = MATRIX_THERMO14*mul;\
const double t15 = MATRIX_THERMO15*mul;\
const double t23 = MATRIX_THERMO23*mul;\
const double t24 = MATRIX_THERMO24*mul;\
const double t25 = MATRIX_THERMO25*mul;\
const double t34 = MATRIX_THERMO34*mul;\
const double t35 = MATRIX_THERMO35*mul;\
const double t45 = MATRIX_THERMO45*mul;\
\
    for (int i = 0 ; i < 4 ; ++i){\
    K[i+ 0][i+ 0] += t11*Sh[i];\
    K[i+ 4][i+ 4] += t22*Sh[i];\
    K[i+ 8][i+ 8] += t33*Sh[i];\
    K[i+12][i+12] += t44*Sh[i];\
    K[i+16][i+16] += t55*Sh[i];\
    for (int j = 0 ; j < 4 ; ++ j){\
    const double ShRC = Sh[i]*Sh[j];\
    K[i+0][j+4] += t12*ShRC;\
    K[i+4][j+0] += t12*ShRC;\
    K[i+0][j+8] += t13*ShRC;\
    K[i+8][j+0] += t13*ShRC;\
    K[i+0 ][j+12] += t14*ShRC;\
    K[i+12][j+0 ] += t14*ShRC;\
    K[i+0 ][j+16]   += t15*ShRC;\
    K[i+16][j+ 0]   += t15*ShRC;\
    \
    K[i +4][j + 8] += t23*ShRC;\
    K[i +8][j + 4] += t23*ShRC;\
    K[i +4][j +12] += t24*ShRC;\
    K[i+12][j + 4] += t24*ShRC;\
    K[i+4 ][j +16] += t25*ShRC;\
    K[i+16][j + 4] += t25*ShRC;\
    \
    K[i+ 8][j +12] += t34*ShRC;\
    K[i+12][j + 8] += t34*ShRC;\
    K[i+ 8][j +16] += t35*ShRC;\
    K[i+16][j + 8] += t35*ShRC;\
    \
    K[i+12][j+16] += t45*ShRC;\
    K[i+16][j+12] += t45*ShRC;\
    }\
    }\
    }


//==============================================
// THERMOTROPIC ENERGY RHS TERMS
//==============================================
#define RHS_THERMOTROPIC(L){\
    const double D3 = 0.33333333333333333333333333333;\
    L[0]+= (A*q1    +   D3*0.5*B*(q5*q5*rt6*0.5    -   q3*q3*rt6   -   rt6*q2*q2   +   q1*q1*rt6   +   q4*q4*rt6*0.5)  +   C*R*q1);\
    L[1]+= (A*q2    +   D3*B*(0.75*q5*q5*rt2    -   q1*rt6*q2   -   0.75*q4*q4*rt2)   +   C*R*q2);\
    L[2]+= (A*q3    +   D3*B*(-q3*q1*rt6   +   1.5*rt2*q5*q4 )     +   C*R*q3);\
    L[3]+= (A*q4    +   D3*0.5*B*(3.0*q3*rt2*q5    +   q4*q1*rt6   -   3.0*q4*q2*rt2)   +   C*R*q4);\
    L[4]+= (A*q5    +   D3*0.5*B*(q5*q1*rt6    +   3.0*q5*q2*rt2   +   3.0*q3*rt2*q4)   +   C*R*q5);\
}

//==============================================
//  DIELECTRIC ENERGY RHS TERMS
//==============================================
#define RHS_DIELECTRIC(L){\
    const double D3 = 0.33333333333333333333333333333;\
    const double D6 = 0.16666666666666666666666666667;\
    L[0] += rt6*(Vx*Vx + Vy*Vy-2.0*Vz*Vz)*deleps*D3*D6*EPS0;\
    L[1] += -rt2*(Vx*Vx - Vy*Vy)*deleps*D6*EPS0;\
    L[2] += -rt2*Vx*Vy*deleps*D3*EPS0;\
    L[3] += -rt2*Vz*Vy*deleps*D3*EPS0;\
    L[4] += -rt2*Vx*Vz*deleps*D3*EPS0;\
}

//==============================================
//  ELASTIC ENERGY RHS VECTOR FOR 2K FORMULATION
//==============================================

#define RHS_ELASTIC_2K_FORMULATION(L){\
  const double t2 = 1.0/rt2;\
  const double t3 = 1.0/rt6;\
  const double t4 = 1.0/(6.0);\
  const double t5 = 0.5;\
  const double t6 = q0*q0;\
  const double t7 = K1*q0*q5;\
  const double t8 = K1*q5y*t5;\
  const double t9 = K1*q2z*t5;\
  const double t10 = K1*q0*q3;\
  const double t11 = K1*q1z*t2*t3;\
  const double t12 = K1*q4x*t5;\
  const double t13 = K1*q0*q2;\
  const double t14 = K1*q0*q1*t2*t3*6.0;\
  L[i + 0] += ShRx*(K1*q1x*t4*6.0+K1*q3y*t2*t3-K1*q5z*t2*t3*2.0-K1*q0*q4*t2*t3*6.0)+ShRy*(K1*q1y*t4*6.0+K1*q3x*t2*t3-K1*q4z*t2*t3*2.0+K1*q0*q5*t2*t3*6.0)+ShRz*(K1*q1z*t4*6.0+K1*q5x*t2*t3+K1*q4y*t2*t3)+ShR*(K1*q1*t4*t6*2.4E1+K1*q0*q4x*t2*t3*6.0-K1*q0*q5y*t2*t3*6.0);\
  L[i + 4] += ShR*(K1*q0*q4x*t5*2.0+K1*q0*q5y*t5*2.0-K1*q0*q3z*t5*4.0+K1*q2*t5*t6*8.0)+ShRx*(K1*q2x*t5*2.0+K1*q3y*t5-K1*q0*q4*t5*2.0)-ShRy*(t7+K1*q3x*t5-K1*q2y*t5*2.0)+ShRz*(-K1*q5x*t5+K1*q4y*t5+K1*q2z*t5*2.0+K1*q0*q3*t5*4.0);\
  L[i + 8] += ShR*(K1*q0*q5x*t5*-2.0+K1*q0*q4y*t5*2.0+K1*q0*q2z*t5*4.0+K1*q3*t5*t6*8.0)+ShRy*(K1*q2x*t5+K1*q3y*t5*2.0-K1*q0*q4*t5*2.0+K1*q1x*t2*t3)-ShRz*(t8+t12-K1*q3z*t5*2.0+K1*q0*q2*t5*4.0)+ShRx*(t7+K1*q3x*t5*2.0-K1*q2y*t5+K1*q1y*t2*t3);\
  L[i +12] += -ShRz*(t7-K1*q4z*t5*2.0+K1*q1y*t2*t3*2.0)-ShR*(K1*q0*q2x*t5*2.0+K1*q0*q3y*t5*2.0-K1*q0*q5z*t5*2.0-K1*q4*t5*t6*8.0+K1*q0*q1x*t2*t3*6.0)+ShRx*(-t8+t13+t14+K1*q4x*t5*2.0-K1*q3z*t5)+ShRy*(t9+t10+t11+K1*q4y*t5*2.0);\
  L[i +16] += -ShRx*(t9+t10-t11-K1*q5x*t5*2.0)+ShR*(K1*q0*q3x*t5*2.0-K1*q0*q2y*t5*2.0-K1*q0*q4z*t5*2.0+K1*q5*t5*t6*8.0+K1*q0*q1y*t2*t3*6.0)+ShRz*(K1*q5z*t5*2.0+K1*q0*q4*t5*2.0-K1*q1x*t2*t3*2.0)-ShRy*(t12-t13+t14-K1*q5y*t5*2.0+K1*q3z*t5);\
}

//=============================================
// EALSTIC ENERGY MATRIX FOR 2K FORMULATION
//=============================================
#define MATRIX_ELASTIC_2K_FORMULATION(M){\
const double t16 = 1.0/(6.0);\
const double t17 = 1.0/rt2;\
const double t18 = 1.0/rt6;\
const double t19 = 0.5;\
const double t20 = q0*q0;\
const double t21 = K1*ShCx*ShRy*t17*t18;\
const double t22 = K1*ShCy*ShRx*t17*t18;\
const double t23 = t21+t22;\
const double t24 = K1*ShCy*ShRx*t19;\
const double t25 = K1*ShC*ShRz*q0*t19*4.0;\
const double t26 = K1*ShCx*ShRx*t19*2.0;\
const double t27 = K1*ShCy*ShRy*t19*2.0;\
const double t28 = K1*ShCz*ShRz*t19*2.0;\
const double t29 = K1*ShC*ShR*t19*t20*8.0;\
const double t30 = t26+t27+t28+t29;\
const double t31 = K1*ShCy*ShR*q0*t19*2.0;\
const double t37 = K1*ShC*ShRy*q0*t19*2.0;\
const double t32 = t31-t37-K1*ShCx*ShRz*t19;\
const double t33 = K1*ShCy*ShRz*t19;\
const double t34 = K1*ShCx*ShR*q0*t19*2.0;\
const double t35 = K1*ShCx*ShR*q0*t17*t18*6.0;\
const double t36 = K1*ShC*ShRx*q0*t19*2.0;\
const double t38 = K1*ShC*ShRy*q0*t17*t18*6.0;\
const double t39 = -t31+t37-K1*ShCz*ShRx*t19;\
const double t40 = K1*ShCz*ShRy*t19;\
const double t41 = K1*ShCx*ShRy*t19;\
const double t42 = K1*ShCz*ShR*q0*t19*2.0;\
  M[i + 0][j + 0] += K1*ShCx*ShRx*t16*6.0+K1*ShCy*ShRy*t16*6.0+K1*ShCz*ShRz*t16*6.0+K1*ShC*ShR*t16*t20*2.4E1;\
  M[i + 0][j + 8] += t23;\
  M[i + 0][j +12] += t35+K1*ShCy*ShRz*t17*t18-K1*ShCz*ShRy*t17*t18*2.0-K1*ShC*ShRx*q0*t17*t18*6.0;\
  M[i + 0][j +16] += t38+K1*ShCx*ShRz*t17*t18-K1*ShCz*ShRx*t17*t18*2.0-K1*ShCy*ShR*q0*t17*t18*6.0;\
  M[i + 4][j + 4] += t30;\
  M[i + 4][j + 8] += t24+t25-K1*ShCx*ShRy*t19-K1*ShCz*ShR*q0*t19*4.0;\
  M[i + 4][j +12] += t33+t34-K1*ShC*ShRx*q0*t19*2.0;\
  M[i + 4][j +16] += t32;\
  M[i + 8][j + 0] += t23;\
  M[i + 8][j + 4] += -t24-t25+t41+K1*ShCz*ShR*q0*t19*4.0;\
  M[i + 8][j + 8] += t30;\
  M[i + 8][j +12] += t32;\
  M[i + 8][j +16] += -t33-t34+t36;\
  M[i +12][j + 0] += -t35-K1*ShCy*ShRz*t17*t18*2.0+K1*ShCz*ShRy*t17*t18+K1*ShC*ShRx*q0*t17*t18*6.0;\
  M[i +12][j + 4] += -t34+t36+t40;\
  M[i +12][j + 8] += t39;\
  M[i +12][j +12] += t30;\
  M[i +12][j +16] += -t24+t42-K1*ShC*ShRz*q0*t19*2.0;\
  M[i +16][j + 0] += -t38-K1*ShCx*ShRz*t17*t18*2.0+K1*ShCz*ShRx*t17*t18+K1*ShCy*ShR*q0*t17*t18*6.0;\
  M[i +16][j + 4] += t39;\
  M[i +16][j + 8] += t34-t36-t40;\
  M[i +16][j + 12] += -t41-t42+K1*ShC*ShRz*q0*t19*2.0;\
  M[i +16][j + 16] += t30;\
}


//==============================================
//  ELASTIC ENERGY RHS FOR SINGLE K
//==============================================
#define RHS_ELASTIC_SINGLE_K(lL){\
    lL[i+0]  += (ShRx*q1x+ShRy*q1y+ShRz*q1z)*L1; \
    lL[i+4]  += (ShRx*q2x+ShRy*q2y+ShRz*q2z)*L1; \
    lL[i+8]  += (ShRx*q3x+ShRy*q3y+ShRz*q3z)*L1; \
    lL[i+12] += (ShRx*q4x+ShRy*q4y+ShRz*q4z)*L1;\
    lL[i+16] += (ShRx*q5x+ShRy*q5y+ShRz*q5z)*L1;\
}
//================================================
//  ELASTIC ENERGY MATRIX FOR SINGLE K
//================================================
#define MATRIX_ELASTIC_SINGLE_K(lK){\
double dot = L1*mul*(dSh[i][0]*dSh[j][0]+dSh[i][1]*dSh[j][1]+dSh[i][2]*dSh[j][2]);\
    lK[i   ][j   ] += dot;\
    lK[i+4 ][j+4 ] += dot;\
    lK[i+8 ][j+8 ] += dot;\
    lK[i+12][j+12] += dot;\
    lK[i+16][j+16] += dot;\
}

#endif // QASSEMBLY_MACROS_H
