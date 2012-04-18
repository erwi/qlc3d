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


//==============================================
//  ELASTIC ENERGY RHS VECTOR FOR 2K FORMULATION
//==============================================

#define RHS_ELASTIC_2K_FORMULATION(lL){\
  double t2 = 1.0/(rt6*rt6); \
  double t3 = 1.0/rt2; \
  double t4 = 1.0/rt6; \
  double t5 = 1.0/(rt2*rt2);\
  double t6 = q0*q0;\
  double t7 = K2*q2x*t5;\
  double t8 = K2*q3y*t5;\
  double t9 = K2*q1x*t3*t4;\
  double t10 = K2*q3x*t5;\
  double t11 = K2*q0*q5*t5*2.0;\
  double t12 = K2*q1y*t3*t4;\
  double t13 = K2*q5y*t5;\
  double t14 = K2*q4y*t5;\
  double t15 = K2*q4x*t5;\
  double t16 = K2*q0*q2*t5*2.0;\
  double t17 = K2*q0*q1*t3*t4*6.0;\
  double t18 = K2*q2z*t5;\
  double t19 = K2*q0*q3*t5*2.0;\
  double t20 = K2*q1z*t3*t4;\
 lL[i + 0] += ShRz*(K1*q1z*t2*4.0+K2*q1z*t2*2.0+K2*q5x*t3*t4+K2*q4y*t3*t4)+ShRx*(K1*q1x*t2+K2*q1x*t2*5.0-K1*q2x*t3*t4+K2*q2x*t3*t4+K2*q3y*t3*t4-K2*q5z*t3*t4*2.0-K2*q0*q4*t3*t4*6.0)+ShRy*(K1*q1y*t2+K2*q1y*t2*5.0+K2*q3x*t3*t4+K1*q2y*t3*t4-K2*q2y*t3*t4-K2*q4z*t3*t4*2.0+K2*q0*q5*t3*t4*6.0)+ShR*(K2*q1*t2*t6*2.4E1+K2*q0*q4x*t3*t4*6.0-K2*q0*q5y*t3*t4*6.0);\
 lL[i + 4] += ShRx*(t7+t8+t9+K1*q2x*t5-K2*q0*q4*t5*2.0-K1*q1x*t3*t4)+ShR*(K2*q0*q4x*t5*2.0+K2*q0*q5y*t5*2.0-K2*q0*q3z*t5*4.0+K2*q2*t5*t6*8.0)-ShRy*(t10+t11+t12-K1*q2y*t5-K2*q2y*t5-K1*q1y*t3*t4)+ShRz*(t14-K2*q5x*t5+K2*q2z*t5*2.0+K2*q0*q3*t5*4.0);\
 lL[i + 8] += ShR*(K2*q0*q5x*t5*-2.0+K2*q0*q4y*t5*2.0+K2*q0*q2z*t5*4.0+K2*q3*t5*t6*8.0)-ShRz*(t13+t15-K2*q3z*t5*2.0+K2*q0*q2*t5*4.0)+ShRx*(t10+t11+t12+K1*q3x*t5-K2*q2y*t5)+ShRy*(t7+t8+t9+K1*q3y*t5-K2*q0*q4*t5*2.0);\
 lL[i +12] += ShRy*(t14+t18+t19+t20+K1*q4y*t5)-ShR*(K2*q0*q2x*t5*2.0+K2*q0*q3y*t5*2.0-K2*q0*q5z*t5*2.0-K2*q4*t5*t6*8.0+K2*q0*q1x*t3*t4*6.0)-ShRz*(t11-K1*q4z*t5-K2*q4z*t5+K2*q1y*t3*t4*2.0)+ShRx*(-t13+t16+t17+K2*q4x*t5*2.0-K2*q3z*t5);\
 lL[i +16] += ShRz*(K1*q5z*t5+K2*q5z*t5+K2*q0*q4*t5*2.0-K2*q1x*t3*t4*2.0)+ShRx*(-t18-t19+t20+K1*q5x*t5+K2*q5x*t5)+ShR*(K2*q0*q3x*t5*2.0-K2*q0*q2y*t5*2.0-K2*q0*q4z*t5*2.0+K2*q5*t5*t6*8.0+K2*q0*q1y*t3*t4*6.0)-ShRy*(t15-t16+t17-K2*q5y*t5*2.0+K2*q3z*t5);\
}\

//=============================================
// EALSTIC ENERGY MATRIX FOR 2K FORMULATION
//=============================================
#define MATRIX_ELASTIC_2K_FORMULATION(M){\
    double t22 = 1.0/(rt6*rt6);\
     double  t23 = K1*t22;\
      double t24 = K2*t22*5.0;\
      double t25 = t23+t24;\
      double t26 = 1.0/rt2;\
      double t27 = 1.0/rt6;\
      double t28 = K1*t26*t27;\
      double t30 = K2*t26*t27;\
      double t29 = t28-t30;\
      double t31 = ShCy*ShRy*t29;\
      double t32 = t31-ShCx*ShRx*t29;\
      double t33 = 1.0/(rt2*rt2);\
      double t34 = K1*t33;\
      double t35 = K2*t33;\
      double t36 = t34+t35;\
      double t37 = q0*q0;\
      double t38 = K2*ShCx*ShRy*t26*t27;\
      double t39 = K2*ShCy*ShRx*t26*t27;\
      double t40 = t38+t39;\
      double t41 = K2*ShCy*ShRx*t33;\
      double t42 = K2*ShC*ShRz*q0*t33*4.0;\
      double t43 = ShCx*ShRx*t36;\
      double t44 = ShCy*ShRy*t36;\
      double t45 = K2*ShCz*ShRz*t33*2.0;\
      double t46 = K2*ShC*ShR*t33*t37*8.0;\
      double t47 = t43+t44+t45+t46;\
      double t48 = K2*ShCy*ShR*q0*t33*2.0;\
      double t54 = K2*ShC*ShRy*q0*t33*2.0;\
      double t49 = t48-t54-K2*ShCx*ShRz*t33;\
      double t50 = K2*ShCy*ShRz*t33;\
      double t51 = K2*ShCx*ShR*q0*t33*2.0;\
      double t52 = K2*ShCx*ShR*q0*t26*t27*6.0;\
      double t53 = K2*ShC*ShRx*q0*t33*2.0;\
      double t55 = K2*ShC*ShRy*q0*t26*t27*6.0;\
      double t56 = -t48+t54-K2*ShCz*ShRx*t33;\
      double t57 = K2*ShCz*ShRy*t33;\
      double t58 = K2*ShCx*ShRy*t33;\
      double t59 = K2*ShCz*ShR*q0*t33*2.0;\
      double t60 = ShCz*ShRz*t36;\
      M[i+0][j+0]  += ShCx*ShRx*t25+ShCy*ShRy*t25+ShCz*ShRz*(K1*t22*4.0+K2*t22*2.0)+K2*ShC*ShR*t22*t37*2.4E1;\
      M[i+0][j+4]  += t32;\
      M[i+0][j+8]  += t40;\
      M[i+0][j+12] += t52+K2*ShCy*ShRz*t26*t27-K2*ShCz*ShRy*t26*t27*2.0-K2*ShC*ShRx*q0*t26*t27*6.0;\
      M[i+0][j+16] += t55+K2*ShCx*ShRz*t26*t27-K2*ShCz*ShRx*t26*t27*2.0-K2*ShCy*ShR*q0*t26*t27*6.0;\
      M[i+4][j+0]  += t32;\
      M[i+4][j+4]  += t47;\
      M[i+4][j+8]  += t41+t42-K2*ShCx*ShRy*t33-K2*ShCz*ShR*q0*t33*4.0;\
      M[i+4][j+12] += t50+t51-K2*ShC*ShRx*q0*t33*2.0;\
      M[i+4][j+16] += t49;\
      M[i+8][j+0]  += t40;\
      M[i+8][j+4]  += -t41-t42+t58+K2*ShCz*ShR*q0*t33*4.0;\
      M[i+8][j+8]  += t47;\
      M[i+8][j+12] += t49;\
      M[i+8][j+16] += -t50-t51+t53;\
      M[i+12][j+0] += -t52-K2*ShCy*ShRz*t26*t27*2.0+K2*ShCz*ShRy*t26*t27+K2*ShC*ShRx*q0*t26*t27*6.0;\
      M[i+12][j+4] += -t51+t53+t57;\
      M[i+12][j+8] += t56;\
      M[i+12][j+12] += t44+t46+t60+K2*ShCx*ShRx*t33*2.0;\
      M[i+12][j+16] += -t41+t59-K2*ShC*ShRz*q0*t33*2.0;\
      M[i+16][j+0] += -t55-K2*ShCx*ShRz*t26*t27*2.0+K2*ShCz*ShRx*t26*t27+K2*ShCy*ShR*q0*t26*t27*6.0;\
      M[i+16][j+4] += t56;\
      M[i+16][j+8] += t51-t53-t57;\
      M[i+16][j+12] += -t58-t59+K2*ShC*ShRz*q0*t33*2.0;\
      M[i+16][j+16] += t43+t46+t60+K2*ShCy*ShRy*t33*2.0;\
}\


//==============================================
//  ELASTIC ENERGY RHS FOR SINGLE K
//==============================================
#define RHS_ELASTIC_SINGLE_K(lL){\
    lL[i+0]  += (ShRx*q1x+ShRy*q1y+ShRz*q1z)*L1; \
    lL[i+4]  += (ShRx*q2x+ShRy*q2y+ShRz*q2z)*L1; \
    lL[i+8]  += (ShRx*q3x+ShRy*q3y+ShRz*q3z)*L1; \
    lL[i+12] += (ShRx*q4x+ShRy*q4y+ShRz*q4z)*L1;\
    lL[i+16] += (ShRx*q5x+ShRy*q5y+ShRz*q5z)*L1;\
}\
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
}\

#endif // QASSEMBLY_MACROS_H
