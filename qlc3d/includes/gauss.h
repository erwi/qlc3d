#ifndef GAUSS_H
#define GAUSS_H

#include <math.h>

// Gauss integration
// ---------------------------------------------------------
//     3D Gauss-Legendre weights for N = 11, D = 4
// ---------------------------------------------------------

const int ngp=11;
const double a=(1+sqrt(5.0/14.0))/4.0;
const double b=(1-sqrt(5.0/14.0))/4.0;

static double gp[ngp][4]={
        {0.25	  , 0.25	,	0.25	,0.25},
        {11.0/14.0     ,	1.0/14.0	,	1.0/14.0	,1.0/14.0},
        {1.0/14.0      ,	11.0/14.0	,	1.0/14.0	,1.0/14.0},
        {1.0/14.0	  ,	1.0/14.0	,	11.0/14.0   ,1.0/14.0},
        {1.0/14.0	  , 1.0/14.0	,	1.0/14.0	,11.0/14.0},
        {a, a, b, b},
        {a, b, a, b},
        {a, b, b, a},
        {b, a, a, b},
        {b, a, b, a},
        {b, b, a, a}};

const double w11 = -74.0/5625.0;
const double w12 = 343.0/45000.0;
const double w13 = 56.0/2250.0;
static double w[ngp]={w11,w12,w12,w12,w12,w13,w13,w13,w13,w13,w13};

static double sh1[ngp][4]; // P1 Shape functions
static double sh1r[ngp][4]; // P1 Shape functions r-derivatives
static double sh1s[ngp][4]; // P1 Shape functions s-derivatives
static double sh1t[ngp][4]; //P1 shape functions t-derivative

//GLOBALLY USED VALUES
//static double rt2 = sqrt(2.0);
//static double rt3 = sqrt(3.0);
//static double rt6 = sqrt(6.0);


#endif
