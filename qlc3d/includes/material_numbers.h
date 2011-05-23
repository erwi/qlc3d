
// Define Material nubers
#ifndef MATERIAL_NUMBERS_H
#define MATERIAL_NUMBERS_H

#define MAT_DOMAIN1 4
#define MAT_DOMAIN2 8
#define MAT_DOMAIN3 12
#define MAT_DOMAIN4 16
#define MAT_DOMAIN5 20
#define MAT_DOMAIN6 24
#define MAT_DOMAIN7 28

#define MAT_DIELECTRIC1 32
#define MAT_DIELECTRIC2 36
#define MAT_DIELECTRIC3 40
#define MAT_DIELECTRIC4 44
#define MAT_DIELECTRIC5 48
#define MAT_DIELECTRIC6 52
#define MAT_DIELECTRIC7 56

#define MAT_ELECTRODE1  64
#define MAT_ELECTRODE2  128
#define MAT_ELECTRODE3  192
#define MAT_ELECTRODE4  256
#define MAT_ELECTRODE5  320
#define MAT_ELECTRODE6  384
#define MAT_ELECTRODE7  448
#define MAT_ELECTRODE8  518
#define MAT_ELECTRODE9  576

#define MAT_FIXLC1 2048
#define MAT_FIXLC2 4096
#define MAT_FIXLC3 6144
#define MAT_FIXLC4 8192
#define MAT_FIXLC5 10240
#define MAT_FIXLC6 12288
#define MAT_FIXLC7 14336
#define MAT_FIXLC8 16384
#define MAT_FIXLC9 18432

#define MAT_PERIODIC 3
#define MAT_NEUMANN 2

inline
int FIXLCN_TO_MATNUM(const int& n){
    /*! converts fixlc number to material number. e.g. 1 -> 2048 etc.*/
    return n*MAT_FIXLC1;
}





#endif




