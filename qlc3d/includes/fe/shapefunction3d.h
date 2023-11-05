#include <vector>
#include "geom/vec3.h"
#ifndef SHAPEFUNCTION3D_H
#define SHAPEFUNCTION3D_H

// forward declarations
namespace qlc3d {
  class TTensor;
  class DielectricPermittivity;
}


// Table 10.4 Quadrature for unit tetrahedra in Zienkiewicz
#define NGP4 11
#define W4_11 -74.0/5625.0
#define W4_12 343.0/45000.0
#define W4_13 56.0/2250.0

#define A4  (1+sqrt(5.0/14.0))/4.0
#define B4  (1-sqrt(5.0/14.0))/4.0
    const double W4[NGP4] = {W4_11, W4_12, W4_12, W4_12, W4_12, W4_13, W4_13, W4_13, W4_13, W4_13, W4_13};
    const double gp[NGP4][4]={	{0.25	  , 0.25	,	0.25	,0.25}
                            ,{11.0/14.0     ,	1.0/14.0	,	1.0/14.0	,1.0/14.0}
                            ,{1.0/14.0      ,	11.0/14.0	,	1.0/14.0	,1.0/14.0}
                            ,{1.0/14.0	  ,	1.0/14.0	,	11.0/14.0   ,1.0/14.0},
                            {1.0/14.0  , 1.0/14.0	,	1.0/14.0	,11.0/14.0},
                            {A4	  , A4	,  B4 , B4},
                            {A4	  , B4	,  A4 , B4},
                            {A4   , B4  ,  B4 , A4},
                            {B4	  , A4  ,  A4 , B4},
                            {B4	  , A4  ,  B4 , A4},
                            {B4	  , B4  ,  A4 , A4}};

class Shape4thOrder {
    public:
    const unsigned int ngp;
    double w[NGP4];
    double sh1[NGP4][4];
    double sh1r[NGP4][4];
    double sh1s[NGP4][4];
    double sh1t[NGP4][4];



    Shape4thOrder(): ngp(NGP4) {
        for (unsigned int i = 0 ; i < ngp ; i++) {
            w[i] = W4[i];
            // P1 Shape functions
            sh1[i][0]=1-gp[i][0]-gp[i][1]-gp[i][2];
            sh1[i][1]=gp[i][0];
            sh1[i][2]=gp[i][1];
            sh1[i][3]=gp[i][2];
            // P1 Shape functions r-derivatives
            sh1r[i][0]=-1.0;
            sh1r[i][1]=1.0;
            sh1r[i][2]=0.0;
            sh1r[i][3]=0.0;
            // P1 Shape functions s-derivatives
            sh1s[i][0]=-1.0;
            sh1s[i][1]=0.0;
            sh1s[i][2]=1.0;
            sh1s[i][3]=0.0;
            // P1 Shape functions t-derivatives
            sh1t[i][0]=-1.0;
            sh1t[i][1]=0.0;
            sh1t[i][2]=0.0;
            sh1t[i][3]=1.0;
        }
    }
};

#endif // SHAPEFUNCTION3D_H
#undef NGP4
#undef W4_11
#undef W4_12
#undef W4_13
#undef A4
#undef B4
