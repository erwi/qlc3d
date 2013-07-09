#include <math.h>

#include <omp.h>
#include <time.h>
#include <qlc3d.h>
#include <shapefunction3d.h>
#include <qassembly_macros.h>
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>

#define	BIGNUM 2e16
const int	npt = 4; //Number of Points per Tetrahedra


// ---------------------------------------------------------
//     2D Gauss-Legendre weights for N = 27, D = 11
// ---------------------------------------------------------
const int ngps=27;
const double wsurf[ngps]={ 0.006829866, 0.006829866, 0.006829866, 0.01809227, 0.01809227, 0.01809227, 0.0004635032, 0.0004635032,
                           0.0004635032,0.02966149 , 0.02966149 , 0.02966149, 0.03857477, 0.03857477, 0.03857477, 0.02616856,
                           0.02616856  ,0.02616856 , 0.02616856 , 0.02616856, 0.02616856, 0.01035383, 0.01035383, 0.01035383,
                           0.01035383  ,0.01035383 , 0.01035383 };
const double sgp[ngps][2] =
{
    {0.9352701	,0.03236495},
    {0.03236495	,0.9352701 },
    {0.03236495	,0.03236495},
    {0.7612982	,0.1193509},
    {0.1193509	,0.7612982},
    {0.1193509	,0.1193509},
    {0.0692221	,0.534611} ,
    {0.534611	,-0.0692221},
    {0.534611	,0.534611},
    {0.5933802	,0.2033099},
    {0.2033099	,0.5933802},
    {0.2033099	,0.2033099},
    {0.2020614	,0.3989693},
    {0.3989693	,0.2020614},
    {0.3989693	,0.3989693},
    {0.05017814	,0.5932012},
    {0.05017814	,0.3566206},
    {0.5932012	,0.05017814},
    {0.5932012	,0.3566206},
    {0.3566206	,0.05017814},
    {0.3566206	,0.5932012},
    {0.02102202	,0.807489},
    {0.02102202	,0.171489},
    {0.807489	,0.02102202},
    {0.807489	,0.171489},
    {0.171489	,0.02102202},
    {0.171489	,0.807489}
};

// using ngps for volume shapes too, since ngps>gps
double sh1[ngps][4]; // P1 Shape functions
double sh1r[ngps][4]; // P1 Shape functions r-derivatives
double sh1s[ngps][4]; // P1 Shape functions s-derivatives
double sh1t[ngps][4]; //P1 shape functions t-derivative
double ssh1[ngps][3];	//SURFACE term P1 shape function

const double rt2 = sqrt(2.0);
const double rt6 = sqrt(6.0);
const double rt3 = sqrt(3.0);

static double S0;
static double L1, L2, L4, L6;
static double A, B, C;
static double deleps, efe, efe2;
static unsigned int npLC;
void init_globals(LC& mat_par, SolutionVector& q){
    S0	= mat_par.S0;
    L1 = mat_par.L1 ;
    L2 = mat_par.L2 ;
    L4 = mat_par.L4 ;
    L6 = mat_par.L6 ;
    A = mat_par.A;
    B = mat_par.B;
    C = mat_par.C;
    npLC = (unsigned int) q.getnDoF();
    deleps = ( mat_par.eps_par - mat_par.eps_per ) / S0;
    efe  = (2.0/S0/3.0)*(mat_par.e11+2*mat_par.e33);
    efe2 = (4.0/S0/9.0)*(mat_par.e11 - mat_par.e33);
}


void init_shape_N() // initialise Neumann and surface derivative shape functions
{
    //for (int i=0; i<ngp; i++)
    for (int i=0; i<ngps; i++){
        // P1 Shape functions
        sh1[i][0]=1-sgp[i][0]-sgp[i][1];
        sh1[i][1]=sgp[i][0];
        sh1[i][2]=sgp[i][1];
        sh1[i][3]=0;//gps[i][2];	// this is always zero
        // P1 Shape functions r-derivatives
        sh1r[i][0]=-1;
        sh1r[i][1]=1;
        sh1r[i][2]=0;
        sh1r[i][3]=0;
        // P1 Shape functions s-derivatives
        sh1s[i][0]=-1;
        sh1s[i][1]=0;
        sh1s[i][2]=1;
        sh1s[i][3]=0;
        //P1 Shape functions t-derivatives
        sh1t[i][0]=-1;
        sh1t[i][1]=0;
        sh1t[i][2]=0;
        sh1t[i][3]=1;
    }
}				

void init_shape_surf(){
    for (int i=0; i<ngps;i++){
        ssh1[i][0]=1-sgp[i][0]-sgp[i][1];
        ssh1[i][1]=sgp[i][0];
        ssh1[i][2]=sgp[i][1];
    }
}//end void init_shape_surf


inline void localKL(double* p,Mesh* t,
                    idx element_num, SolutionVector* q ,
                    SolutionVector* v, double lK[20][20],
double lL[20],LC* mat_par, Simu* simu,
const Shape4thOrder& shapes)
{

    memset(lK,0,20*20*sizeof(double));
    memset(lL,0,20*sizeof(double));
    double lI[20][20]; // LOCAL IDENTITY MATRIX
    memset(lI,0,20*20*sizeof(double));


    // LOCAL COPY OF ELEMENT
    idx tt[4] = {t->getNode(element_num,0),
                 t->getNode(element_num,1),
                 t->getNode(element_num,2),
                 t->getNode(element_num,3)};
    double Jdet = t->getDeterminant(element_num)*1e18 ; // SCALE BACK TO METRES FOR NOW...

    if (Jdet < 0) {
        printf("Warning, Jdet < 0\n");
        Jdet = -Jdet;
    }

    bool three_elastic_constants = false;
    if ((L2!=0)&&(L6!=0)) three_elastic_constants = true;

    //1. Calculate Inverse Jacobian - for 1st order elements can be done outside integration loop -> igp = 0
    double xr(0),xs(0),xt(0),yr(0),ys(0),yt(0),zr(0),zs(0),zt(0);
    const double pp[4][3] ={ {p[tt[0]*3] , p[tt[0]*3+1] , p[tt[0]*3+2]} ,
                             {p[tt[1]*3] , p[tt[1]*3+1] , p[tt[1]*3+2]} ,
                             {p[tt[2]*3] , p[tt[2]*3+1] , p[tt[2]*3+2]} ,
                             {p[tt[3]*3] , p[tt[3]*3+1] , p[tt[3]*3+2]} };

    for (int i=0;i<4;++i){
        xr+=shapes.sh1r[0][i]*pp[i][0];
        xs+=shapes.sh1s[0][i]*pp[i][0];
        xt+=shapes.sh1t[0][i]*pp[i][0];

        yr+=shapes.sh1r[0][i]*pp[i][1];
        ys+=shapes.sh1s[0][i]*pp[i][1];
        yt+=shapes.sh1t[0][i]*pp[i][1];

        zr+=shapes.sh1r[0][i]*pp[i][2];
        zs+=shapes.sh1s[0][i]*pp[i][2];
        zt+=shapes.sh1t[0][i]*pp[i][2];
    }//end for i
    //Inverse Jacobian

    // SCALING TO MICRONS HERE: (1e-6*1e-6)/1e-18 = x1e6
    const double invJdet = 1e6/Jdet;
    const double Jinv[3][3]={{ (zt*ys-yt*zs)*invJdet ,(xt*zs-zt*xs)*invJdet, (xs*yt-ys*xt)*invJdet}
                             ,{(yt*zr-zt*yr)*invJdet ,(zt*xr-xt*zr)*invJdet, (xt*yr-yt*xr)*invJdet}
                             ,{(yr*zs-ys*zr)*invJdet ,(xs*zr-xr*zs)*invJdet, (ys*xr-xs*yr)*invJdet}};

    Jdet *= 1e-18; // NEEDS TO BE microns^3 LATER IN mul = ...
    // shape function derivatives
    double dSh[4][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    for(int i=0;i<4;++i){
        dSh[i][0]=shapes.sh1r[0][i]*Jinv[0][0]
                + shapes.sh1s[0][i]*Jinv[1][0]
                + shapes.sh1t[0][i]*Jinv[2][0];

        dSh[i][1]=shapes.sh1r[0][i]*Jinv[0][1]
                + shapes.sh1s[0][i]*Jinv[1][1]
                + shapes.sh1t[0][i]*Jinv[2][1];

        dSh[i][2]=shapes.sh1r[0][i]*Jinv[0][2]
                + shapes.sh1s[0][i]*Jinv[1][2]
                + shapes.sh1t[0][i]*Jinv[2][2];
    }//end for i

    // LOCAL COPY OF Q-TENSOR AND POTENTIAL VARIABLES AT THE FOUR NODES OF THE CURRENT ELEMENT
    // THIS SEEMS TO SPEED UP THE EXECUTION OF THIS FUNCTION
    const double qbuff[6*4] = {q->getValue(tt[0],0), q->getValue(tt[0],1), q->getValue(tt[0],2), q->getValue(tt[0],3), q->getValue(tt[0],4),
                               q->getValue(tt[1],0), q->getValue(tt[1],1), q->getValue(tt[1],2), q->getValue(tt[1],3), q->getValue(tt[1],4),
                               q->getValue(tt[2],0), q->getValue(tt[2],1), q->getValue(tt[2],2), q->getValue(tt[2],3), q->getValue(tt[2],4),
                               q->getValue(tt[3],0), q->getValue(tt[3],1), q->getValue(tt[3],2), q->getValue(tt[3],3), q->getValue(tt[3],4),
                               v->getValue(tt[0])  , v->getValue(tt[1]),   v->getValue(tt[2]),   v->getValue(tt[3])
                              };
    for (unsigned int igp = 0 ; igp < shapes.ngp ; ++igp){
        //shape function
        const double Sh[4]= {   shapes.sh1[igp][0], shapes.sh1[igp][1],
                                shapes.sh1[igp][2], shapes.sh1[igp][3]};

        // Function variables and derivatives
        double q1=0, q2=0, q3=0, q4=0, q5=0;
        double q1x=0,q2x=0,q3x=0,q4x=0,q5x=0;
        double q1y=0,q2y=0,q3y=0,q4y=0,q5y=0;
        double q1z=0,q2z=0,q3z=0,q4z=0,q5z=0;
        double Vx=0,Vy=0,Vz=0;

        // Solution and derivatives
#define IND(i,j) 5*(i) + (j)
        for(int i=0;i<4;++i){
            q1+=Sh[i]*qbuff[IND(i,0)];// OPTIMIZE BY PREFETCHING Q AND V TO LOCAL BUFFER AT START OF FUNCTION
            q2+=Sh[i]*qbuff[IND(i,1)];
            q3+=Sh[i]*qbuff[IND(i,2)];
            q4+=Sh[i]*qbuff[IND(i,3)];
            q5+=Sh[i]*qbuff[IND(i,4)];

            q1x+=dSh[i][0]*qbuff[IND(i,0)];
            q2x+=dSh[i][0]*qbuff[IND(i,1)];
            q3x+=dSh[i][0]*qbuff[IND(i,2)];
            q4x+=dSh[i][0]*qbuff[IND(i,3)];
            q5x+=dSh[i][0]*qbuff[IND(i,4)];

            q1y+=dSh[i][1]*qbuff[IND(i,0)];
            q2y+=dSh[i][1]*qbuff[IND(i,1)];
            q3y+=dSh[i][1]*qbuff[IND(i,2)];
            q4y+=dSh[i][1]*qbuff[IND(i,3)];
            q5y+=dSh[i][1]*qbuff[IND(i,4)];

            q1z+=dSh[i][2]*qbuff[IND(i,0)];
            q2z+=dSh[i][2]*qbuff[IND(i,1)];
            q3z+=dSh[i][2]*qbuff[IND(i,2)];
            q4z+=dSh[i][2]*qbuff[IND(i,3)];
            q5z+=dSh[i][2]*qbuff[IND(i,4)];
            // voltages
            Vx+=dSh[i][0]*qbuff[20+i];
            Vy+=dSh[i][1]*qbuff[20+i];
            Vz+=dSh[i][2]*qbuff[20+i];
        }//end for i

        const double mul=shapes.w[igp]*Jdet;
        const double R=q1*q1+q2*q2+q3*q3+q5*q5+q4*q4; // frequently reoccurring term
        //const double D2 = 0.5;
        //const double D3 = 0.33333333333333333333333333333;
        //const double D6 = 0.16666666666666666666666666667;


        // BULK RHS TERMS
        double L[5] = {0,0,0,0,0};
        RHS_THERMOTROPIC(L);
        RHS_DIELECTRIC(L);
        ADD_RHS_BULK_TERMS(lL,L);

        MATRIX_THERMOTROPIC(lK);
        /*
        // BULK THERMOTROPIC MATRIX TERMS
        {
            double Th;

            Th = MATRIX_THERMO11;//(A  +   D3*B*q1*rt6 +   2.0*C*q1*q1 + C*R); // T11
            THERMOdiag(0,Th);

            Th = MATRIX_THERMO12;//(-D3*B*rt6*q2   +       2.0*C*q2*q1);    // T12,T21
            THERMO(0,4,Th);

            Th = MATRIX_THERMO13;//(-B*q3*rt6*D3   +       2.0*C*q3*q1);    // T13,T31
            THERMO(0,8,Th);

            Th = MATRIX_THERMO14;//(B*q4*rt6*D6    +       2.0*C*q4*q1);    // T14, T41
            THERMO(0,12,Th);

            Th = MATRIX_THERMO15;//(B*q5*rt6*D6    +       2.0*C*q5*q1);    // T15, T51
            THERMO(0,16,Th);

            Th = MATRIX_THERMO22;//(A  -   D3*B*q1*rt6 +   2.0*C*q2*q2+C*R);// T22
            THERMOdiag(4,Th);

            Th = MATRIX_THERMO23;//(2.0*C*q3*q2);   // T23, T32
            THERMO(4,8,Th);

            Th = MATRIX_THERMO24;//(-B*q4*rt2*D2   +   2.0*C*q4*q2);    // T24, T42
            THERMO(4,12,Th);

            Th = MATRIX_THERMO25;//(B*q5*rt2*D2    +   2.0*C*q5*q2);    // T25, T52
            THERMO(4,16,Th);

            Th = MATRIX_THERMO33;//(A  -   B*q1*rt6*D3 +   2.0*C*q3*q3+C*R);    // T33
            THERMOdiag(8,Th);

            Th = MATRIX_THERMO34;//(B*q5*rt2*D2    +   2.0*C*q4*q3);    // T34, T43
            THERMO(8,12,Th);

            Th = MATRIX_THERMO35;//(B*q4*rt2*D2    +   2.0*C*q5*q3);    // T35, T53
            THERMO(8,16,Th);

            Th = MATRIX_THERMO44;//(A  +   B*(q1*rt6-3.0*q2*rt2)*D6    +   2.0*C*q4*q4+C*R); // T44
            THERMOdiag(12,Th);

            Th = MATRIX_THERMO45;//(B*q3*rt2*D2    +   2.0*C*q5*q4);    // T45, T54
            THERMO(12,16,Th);

            Th = MATRIX_THERMO55;//(A  +   B*(q1*rt6+3.0*q2*rt2)*D6    +   2.0*C*q5*q5+C*R); // T55
            THERMOdiag(16,Th);
        }// END BULK THERMO SCOPE
*/
        for (int i=0;i<4;++i) // matrix rows
        {
            const double ShRx=mul*dSh[i][0];//including weight and jacobian in trial function
            const double ShRy=mul*dSh[i][1];
            const double ShRz=mul*dSh[i][2];
            const double ShR =mul*Sh[i];


            // ADD THERMOTROPIC, ELASTIC (SINGLE K) AND POTENTIAL TERMS TO RHS VECTOR
            // THESE ARE ALWAYS PRESENT
            double L1_term;
            L1_term=(ShRx*q1x+ShRy*q1y+ShRz*q1z)*L1;
            lL[i+0] +=L1_term;
            L1_term=(ShRx*q2x+ShRy*q2y+ShRz*q2z)*L1;
            lL[i+4] += L1_term;
            L1_term=(ShRx*q3x+ShRy*q3y+ShRz*q3z)*L1;
            lL[i+8] += L1_term;
            L1_term=(ShRx*q4x+ShRy*q4y+ShRz*q4z)*L1;
            lL[i+12] += L1_term;
            L1_term=(ShRx*q5x+ShRy*q5y+ShRz*q5z)*L1;
            lL[i+16] += L1_term;

            //Chiral term
            if (L4!=0)
            {
                double Lc;//[5];
                double D4 = 0.25;
                Lc = (   (q5y-q4x)*ShR   +   ShRx*q4  -   ShRy*q5)*L4*D4*rt3;
                lL[i+0] += Lc;
                Lc = (   (2.0*q3z -q5y-q4x)*ShR   +   ShRx*q4 + ShRy*q5 - ShRz*q3*2.0 ) *L4*D4;
                lL[i+4] += Lc;
                Lc = (   (-q2z*2.0*D4-q4y*D4+q5x*D4)*ShR  +   ShRy*q4*D4 +   ShRz*q2*2.0*D4 -   ShRx*q5*D4)*L4;
                lL[i+8] += Lc;
                Lc = ((-q5z*D4+q1x*rt3*D4+q3y*D4+q2x*D4)*ShR-ShRx*q1*rt3*D4-ShRx*q2*D4-ShRy*q3*D4+ShRz*q5*D4)*L4;
                lL[i+12]+= Lc;
                Lc = ((q4z*D4-q1y*rt3*D4+q2y*D4-q3x*D4)*ShR+ShRx*q3*D4+ShRy*q1*rt3*D4-ShRy*q2*D4-ShRz*q4*D4)*L4;
                lL[i+16]+= Lc;
            }
            if (three_elastic_constants)
            {
                const double D2 = 0.5;\
                const double D3 = 0.33333333333333333333333333333;\
                const double D6 = 0.16666666666666666666666666667;\
                double temp(0);
                const double D4 = 0.25;
                const double rt23 = rt2*rt3;
                //L2 and L6 q1 terms
                temp = (ShRx*q1x*D6-ShRx*rt3*q2x*D6-ShRx*rt3*q3y*D6-ShRx*rt3*q5z*D6-ShRy*q3x*rt3*D6+ShRy*q1y*D6+ShRy*rt3*q2y*D6-ShRy*rt3*q4z*D6+ShRz*q5x*rt3*D3+ShRz*q4y*rt3*D3+2.0*D3*ShRz*q1z)	*L2;
                // L6 q1 term
                temp += (-ShRx*q1x*q1*rt23*D6-ShRy*q1y*q1*rt23*D6+ShRz*q1*rt23*q1z*D3-ShR*rt23*q5y*q5y*D2*D6+ShRx*q5*rt2*q1z*D2+ShRy*q3*rt2*q1x*D2+ShRy*q4*rt2*q1z*D2+ShRz*q5*rt2*q1x*D2-ShR*rt23*q4x*q4x*D6*D2-ShR*rt23*q4y*q4y*D6*D2-ShR*rt23*q3x*q3x*D2*D6-ShR*rt23*q2y*q2y*D2*D6+ShR*rt23*q2z*q2z*D6+ShR*rt23*q3z*q3z*D6+ShR*rt23*q5z*q5z*D6+ShR*rt23*q4z*q4z*D6+ShR*rt23*q1z*q1z*D6+ShRz*q4*rt2*q1y*D2-ShR*rt23*q1x*q1x*D2*D6-ShR*rt23*q3y*q3y*D2*D6-ShR*rt23*q2x*q2x*D2*D6-ShR*rt23*q1y*q1y*D2*D6+ShRx*q1x*q2*rt2*D2+ShRx*q3*rt2*q1y*D2-ShRy*q1y*q2*rt2*D2-ShR*rt23*q5x*q5x*D2*D6)*L6;
                lL[i+0] += temp;

                // L2 and L6 q2 term
                temp = (-ShRx*q1x*rt3*D6+q2x*ShRx*D2+q3y*ShRx*D2+q5z*ShRx*D2-q3x*ShRy*D2+ShRy*q1y*rt3*D6+q2y*ShRy*D2-q4z*ShRy*D2)*L2;
                temp += (-ShR*rt2*q5y*q5y*D4-ShR*rt2*q1y*q1y*D4+ShR*rt2*q4x*q4x*D4+ShR*rt2*q1x*q1x*D4-ShR*rt2*q3y*q3y*D4+ShR*rt2*q2x*q2x*D4+ShR*rt2*q5x*q5x*D4+ShR*rt2*q3x*q3x*D4-ShR*rt2*q4y*q4y*D4-ShR*rt2*q2y*q2y*D4-ShRx*q1*rt23*q2x*D6+ShRx*rt2*q2*q2x*D2+ShRx*q3*q2y*rt2*D2+ShRx*q5*q2z*rt2*D2-ShRy*q1*rt23*q2y*D6-ShRy*rt2*q2*q2y*D2+ShRy*q3*q2x*rt2*D2+ShRy*q4*q2z*rt2*D2+ShRz*q5*q2x*rt2*D2+ShRz*q4*q2y*rt2*D2+ShRz*q1*rt23*q2z*D3)*L6;
                lL[i+4] += temp;

                // L2 and L6 q3 terms
                temp = (ShRx*q3x*D2-ShRx*q1y*rt3*D6-ShRx*q2y*D2+ShRx*q4z*D2-ShRy*q1x*rt3*D6+ShRy*q2x*D2+ShRy*q3y*D2+ShRy*q5z*D2)*L2;
                temp+= (ShR*rt2*q1x*q1y*D2+ShR*rt2*q2x*q2y*D2+ShR*rt2*q3x*q3y*D2+ShR*rt2*q5x*q5y*D2+ShR*rt2*q4x*q4y*D2-ShRx*q3x*q1*rt23*D6+ShRx*q3x*q2*rt2*D2+ShRx*q3*rt2*q3y*D2+ShRx*q5*rt2*q3z*D2-ShRy*q3y*q1*rt23*D6-ShRy*q3y*q2*rt2*D2+ShRy*q3*rt2*q3x*D2+ShRy*q4*rt2*q3z*D2+ShRz*q5*rt2*q3x*D2+ShRz*q4*rt2*q3y*D2+ShRz*q1*rt23*q3z*D3)*L6;
                lL[i+8] += temp;

                // L2 and L6 q4 terms
                temp = (ShRy*q5x*D2+ShRy*q4y*D2+ShRy*q1z*rt3*D3+ShRz*q3x*D2-ShRz*q1y*rt3*D6-ShRz*q2y*D2+ShRz*q4z*D2)*L2;
                temp += (ShR*rt2*q1y*q1z*D2+ShR*rt2*q2y*q2z*D2+ShR*rt2*q3y*q3z*D2+ShR*rt2*q5y*q5z*D2+ShR*rt2*q4y*q4z*D2-ShRx*q4x*q1*rt23*D6+ShRx*q4x*q2*rt2*D2+ShRx*q3*rt2*q4y*D2+ShRx*q5*rt2*q4z*D2-ShRy*q4y*q1*rt23*D6-ShRy*q4y*q2*rt2*D2+ShRy*q3*rt2*q4x*D2+ShRy*q4*rt2*q4z*D2+ShRz*q5*rt2*q4x*D2+ShRz*q4*rt2*q4y*D2+ShRz*q1*rt23*q4z*D3)*L6;
                lL[i+12] += temp;

                // L2 and L6 q5 terms
                temp = (ShRx*q5x*D2+ShRx*q4y*D2+ShRx*q1z*rt3*D3-ShRz*q1x*rt3*D6+ShRz*q2x*D2+ShRz*q3y*D2+ShRz*q5z*D2)*L2;
                temp += (ShR*rt2*q1x*q1z*D2+ShR*rt2*q2x*q2z*D2+ShR*rt2*q3x*q3z*D2+ShR*rt2*q5x*q5z*D2+ShR*rt2*q4x*q4z*D2-ShRx*q5x*q1*rt23*D6+ShRx*q5x*q2*rt2*D2+ShRx*q3*rt2*q5y*D2+ShRx*q5*rt2*q5z*D2-ShRy*q5y*q1*rt23*D6-ShRy*q5y*q2*rt2*D2+ShRy*q3*rt2*q5x*D2+ShRy*q4*rt2*q5z*D2+ShRz*q5*rt2*q5x*D2+ShRz*q4*rt2*q5y*D2+ShRz*q1*rt23*q5z*D3)*L6;
                lL[i+16] += temp;

            } // end if 3 elastic contants

            if ((efe!=0.0)||(efe2!=0.0)) // IF FLEXOELECTRIC COEFFICIENTS ARN'T 0
            {
                double flexo = rt6*(Vx*ShRx+Vy*ShRy-2.0*Vz*ShRz)*efe*0.16666666666666666666667;
                lL[i+0] += flexo;
                flexo = -0.5*rt2*(Vx*ShRx-Vy*ShRy)*efe;
                lL[i+4] += flexo;
                flexo = -0.5*rt2*(Vy*ShRx+Vx*ShRy)*efe;
                lL[i+8] += flexo;
                flexo = -0.5*rt2*(Vz*ShRy+Vy*ShRz)*efe;
                lL[i+12] += flexo;
                flexo = -0.5*rt2*(Vz*ShRx+Vx*ShRz)*efe;
                lL[i+16] += flexo;
            }



            // ADD AND ASSEMBLE MATRIX TERMS
            for (int j=0 ; j<4 ; ++j)
            {
                const double ShCx=dSh[j][0];
                const double ShCy=dSh[j][1];
                const double ShCz=dSh[j][2];
                const double ShC =Sh[j];
                const double ShRC=ShR*Sh[j];

                //L1- matrix term 'dot' only appears on diagonal
                const double dot=L1*mul*(dSh[i][0]*dSh[j][0]+dSh[i][1]*dSh[j][1]+dSh[i][2]*dSh[j][2]);

                lK[i   ][j   ] += dot;
                lK[i+4 ][j+4 ] += dot;
                lK[i+8 ][j+8 ] += dot;
                lK[i+12][j+12] += dot;
                lK[i+16][j+16] += dot;

                //Local identity matrix

                lI[i   ][j   ]+=ShRC;
                lI[i+4 ][j+4 ]+=ShRC;
                lI[i+8 ][j+8 ]+=ShRC;
                lI[i+12][j+12]+=ShRC;
                lI[i+16][j+16]+=ShRC;

                if (three_elastic_constants)//if three elastic constants used
                {
                    double temp(0);
                    const double rt2L = rt2;
                    const double rt23 = rt2*rt3;
                    const double D2 = 0.5;\
                    const double D3 = 0.33333333333333333333333333333;\
                    const double D6 = 0.16666666666666666666666666667;\
                    // dlL[0]/dq1 -> dlL[0]/dq5 L2 and L6 terms
                    temp =  (ShRx*ShCx*D6+ShRy*ShCy*D6+2.0*D3*ShRz*ShCz)*L2;
                    temp += (-ShC*ShRx*q1x*rt23*D6-ShC*ShRy*q1y*rt23*D6+ShC*ShRz*rt23*q1z*D3-ShCx*ShRx*q1*rt23*D6+ShCx*ShRy*q3*rt2L/2.0+ShCx*ShRz*q5*rt2L*D2-ShCx*ShR*rt23*q1x*D6+ShCx*ShRx*q2*rt2L*D2-ShCy*ShRy*q1*rt23*D6+ShCy*ShRz*q4*rt2L*D2-ShCy*ShR*rt23*q1y*D6+ShCy*ShRx*q3*rt2L*D2-ShCy*ShRy*q2*rt2L*D2+ShCz*ShRz*q1*rt23*D3+ShCz*ShRx*q5*rt2L*D2+ShCz*ShRy*q4*rt2L*D2+ShCz*ShR*rt23*q1z*D3)*L6;
                    lK[i][j   ] += temp;

                    temp = (-ShRx*rt3*ShCx*D6+ShRy*rt3*ShCy*D6)*L2;
                    temp += (ShC*ShRx*q1x*rt2L*D2-ShC*ShRy*q1y*rt2L*D2-ShR*rt23*q2x*ShCx*D6-ShR*rt23*q2y*ShCy*D6+ShR*rt23*q2z*ShCz*D3)*L6;;
                    lK[i  ][j+4] += temp;
                    lK[i+4][j  ] += temp;

                    temp = (-ShRy*rt3*ShCx*D6-ShRx*rt3*ShCy*D6)*L2;
                    temp += (ShC*ShRy*rt2L*q1x*D2+ShC*ShRx*rt2L*q1y*D2-ShR*rt23*q3x*ShCx*D6-ShR*rt23*q3y*ShCy*D6+ShR*rt23*q3z*ShCz*D3)*L6;
                    lK[i  ][j+8] += temp;
                    lK[i+8][j  ] += temp;

                    temp = (ShRz*rt3*ShCy*D3-ShRy*rt3*ShCz*D6)*L2;
                    temp += (ShC*ShRy*rt2L*q1z*D2+ShC*ShRz*rt2L*q1y*D2-ShR*rt23*q4x*ShCx*D6-ShR*rt23*q4y*ShCy*D6+ShR*rt23*q4z*ShCz*D3)*L6;
                    lK[i   ][j+12] += temp;
                    lK[i+12][j   ] += temp;

                    temp = (ShRz*rt3*ShCx*D3-ShRx*rt3*ShCz*D6)*L2;
                    temp += (ShC*ShRx*rt2L*q1z*D2+ShC*ShRz*rt2L*q1x*D2-ShR*rt23*q5x*ShCx*D6-ShR*rt23*q5y*ShCy*D6+ShR*rt23*q5z*ShCz*D3)*L6;
                    lK[i   ][j+16] += temp;
                    lK[i+16][j   ] += temp;

                    // dlL[4]/dq2 -> dlL[4]/dq5 L2 and L6 terms
                    temp = 0.5*(ShRx*ShCx+ShRy*ShCy)*L2;
                    temp += (ShC*ShRx*rt2L*q2x*D2-ShC*ShRy*rt2L*q2y*D2+ShCx*ShR*rt2L*q2x*D2-ShCx*ShRx*q1*rt23*D6+ShCx*ShRx*q2*rt2L*D2+ShCx*ShRy*q3*rt2L*D2+ShCx*ShRz*q5*rt2L*D2-ShCy*ShR*rt2L*q2y*D2+ShCy*ShRx*q3*rt2L*D2-ShCy*ShRy*q1*rt23*D6-ShCy*ShRy*q2*rt2L*D2+ShCy*ShRz*q4*rt2L*D2+ShCz*ShRx*q5*rt2L*D2+ShCz*ShRy*q4*rt2L*D2+ShCz*ShRz*q1*rt23*D3)*L6;
                    lK[i+4][j+4] += temp;
                    temp = 0.5*(-ShRy*ShCx+ShRx*ShCy)*L2;
                    temp += (ShC*ShRx*q2y*rt2L*D2+ShC*ShRy*q2x*rt2L*D2+ShR*rt2L*q3x*ShCx*D2-ShR*rt2L*q3y*ShCy*D2)*L6;
                    lK[i+4][j+8] += temp;
                    lK[i+8][j+4] += temp;

                    temp = -0.5*ShCz*L2*ShRy;
                    temp += (ShC*ShRy*q2z*rt2L*D2+ShC*ShRz*q2y*rt2L*D2+ShR*rt2L*q4x*ShCx*D2-ShR*rt2L*q4y*ShCy*D2)*L6;
                    lK[i+4 ][j+12] += temp;
                    lK[i+12][j+4 ] += temp;

                    temp = 0.5*ShCz*L2*ShRx;
                    temp += (ShC*ShRx*q2z*rt2L*D2+ShC*ShRz*q2x*rt2L*D2+ShR*rt2L*q5x*ShCx*D2-ShR*rt2L*q5y*ShCy*D2)*L6;
                    lK[i+4 ][j+16] += temp;
                    lK[i+16][j+4 ] += temp;

                    // dlL[8] / dq3 -> dlL[8]/dq5 L2 and L6 terms
                    temp = 0.5*(ShRx*ShCx+ShRy*ShCy)*L2;
                    temp += (ShC*ShRx*rt2L*q3y*D2+ShC*ShRy*rt2L*q3x*D2+ShCx*ShR*rt2L*q3y*D2-ShCx*ShRx*q1*rt23*D6+ShCx*ShRx*q2*rt2L*D2+ShCx*ShRy*q3*rt2L*D2+ShCx*ShRz*q5*rt2L*D2+ShCy*ShR*rt2L*q3x*D2+ShCy*ShRx*q3*rt2L*D2-ShCy*ShRy*q1*rt23*D6-ShCy*ShRy*q2*rt2L*D2+ShCy*ShRz*q4*rt2L*D2+ShCz*ShRx*q5*rt2L*D2+ShCz*ShRy*q4*rt2L*D2+ShCz*ShRz*q1*rt23*D3)*L6;
                    lK[i+8][j+8]+= temp;
                    temp = ShCz*L2*ShRx*0.5;
                    temp += (ShC*ShRy*rt2L*q3z*D2+ShC*ShRz*rt2L*q3y*D2+ShR*rt2L*q4y*ShCx*D2+ShR*rt2L*q4x*ShCy*D2)*L6;
                    lK[i+8 ][j+12] += temp;
                    lK[i+12][j+8 ] += temp;
                    temp = 0.5*ShCz*L2*ShRy;
                    temp += (ShC*ShRx*rt2L*q3z*D2+ShC*ShRz*rt2L*q3x*D2+ShR*rt2L*q5y*ShCx*D2+ShR*rt2L*q5x*ShCy*D2)*L6;
                    lK[i+8 ][j+16] += temp;
                    lK[i+16][j+8 ] += temp;

                    // dlL[12] / dq4 -> dlL[8]/dq5 L2 and L6 terms
                    temp = 0.5*(ShRy*ShCy+ShRz*ShCz)*L2;
                    temp += (ShC*ShRy*rt2L*q4z*D2+ShC*ShRz*rt2L*q4y*D2-ShCx*ShRx*q1*rt23*D6+ShCx*ShRx*q2*rt2L*D2+ShCx*ShRy*q3*rt2L*D2+ShCx*ShRz*q5*rt2L*D2+ShCy*ShR*rt2L*q4z*D2+ShCy*ShRx*q3*rt2L*D2-ShCy*ShRy*q1*rt23*D6-ShCy*ShRy*q2*rt2L*D2+ShCy*ShRz*q4*rt2L*D2+ShCz*ShR*rt2L*q4y*D2+ShCz*ShRx*q5*rt2L*D2+ShCz*ShRy*q4*rt2L*D2+ShCz*ShRz*q1*rt23*D3)*L6;
                    lK[i+12][j+12] += temp;
                    temp = 0.5*ShCx*L2*ShRy;
                    temp += (ShC*ShRx*rt2L*q4z*D2+ShC*ShRz*rt2L*q4x*D2+ShR*rt2L*q5z*ShCy*D2+ShR*rt2L*q5y*ShCz*D2)*L6;
                    lK[i+12][j+16] += temp;
                    lK[i+16][j+12] += temp;
                    // dlL[12] / dq5 L2 and L6 terms
                    temp = 0.5*(ShRx*ShCx+ShRz*ShCz)*L2;
                    temp += (ShC*ShRx*rt2L*q5z*D2+ShC*ShRz*rt2L*q5x*D2+ShCx*ShR*rt2L*q5z*D2-ShCx*ShRx*q1*rt23*D6+ShCx*ShRx*q2*rt2L*D2+ShCx*ShRy*q3*rt2L*D2+ShCx*ShRz*q5*rt2L*D2+ShCy*ShRx*q3*rt2L*D2-ShCy*ShRy*q1*rt23*D6-ShCy*ShRy*q2*rt2L*D2+ShCy*ShRz*q4*rt2L*D2+ShCz*ShR*rt2L*q5x*D2+ShCz*ShRx*q5*rt2L*D2+ShCz*ShRy*q4*rt2L*D2+ShCz*ShRz*q1*rt23*D3)*L6;
                    lK[i+16][j+16] += temp;

                } // end if 3 elastic contants

                // chirality terms
                if (L4!=0)
                {
                    double temp(0);
                    const double D4 = 0.25;
                    //Kc[0][3] = (ShRx*rt3*ShC/4.0-ShR*rt3*ShCx/4.0)*L4;
                    temp = (ShRx*ShC - ShR*ShCx)*L4*D4*rt3;
                    lK[i   ][j+12] += temp;
                    lK[i+12][j+0 ] -= temp;
                    //Kc[0][4] = (-ShRy*rt3*ShC/4.0+ShR*rt3*ShCy/4.0)*L4;
                    temp = (-ShRy*ShC    +   ShR*ShCy)*L4*D4*rt3;
                    lK[i   ][j+16] += temp;
                    lK[i+16][j   ] -= temp;
                    //Kc[1][2] = (-ShRz*ShC/2.0+ShR*ShCz/2.0)*L4;
                    temp = (-ShRz*ShC   +   ShR*ShCz)*L4*0.5;
                    lK[i+4][j+8] += temp;
                    lK[i+8][j+4] -= temp;
                    //Kc[1][3] = (ShC*ShRx/4.0-ShCx*ShR/4.0)*L4;
                    temp = (ShC*ShRx -   ShCx*ShR)*L4*D4;
                    lK[i+4 ][j+12] += temp;
                    lK[i+12][j+4 ] -= temp;
                    //Kc[1][4] = (ShC*ShRy/4.0-ShCy*ShR/4.0)*L4;    // <----SAME A
                    temp = (ShC*ShRy -   ShCy*ShR)*L4*D4;
                    lK[i+4 ][j+16] += temp;
                    lK[i+16][j+4 ] -= temp;
                    //Kc[2][3] = (ShC*ShRy/4.0-ShCy*ShR/4.0)*L4;    // <----SAME A
                    temp = (ShC*ShRy    -   ShCy*ShR)*L4*D4;
                    lK[i+8 ][j+12] += temp;
                    lK[i+12][j+8 ] -= temp;
                    //Kc[2][4] = (-ShC*ShRx/4.0+ShCx*ShR/4.0)*L4;
                    temp = (-ShC*ShRx   +   ShCx*ShR)*L4*D4;
                    lK[i+8 ][j+16] += temp;
                    lK[i+16][j+8 ] -= temp;
                    //Kc[3][4] = (ShRz*ShC/4.0-ShR*ShCz/4.0)*L4;
                    temp = (ShRz*ShC    -   ShR*ShCz)*L4*D4;
                    lK[i+12][j+16] += temp;
                    lK[i+16][j+12] -= temp;

                }


            }//end for j

        } // end for i  - rows
    }//end for igp

    // IF CRANK-NICHOLSON
    if(simu->dt!=0)// Crank-Nicolson time stepping
    {
        // %0 calculates the steady state
        // the product of the mass matrix and the various q vectors
        double Mq[20];	// M * q
        memset(Mq,0,20*sizeof(double));

        for (int i=0;i<4;i++) 		//each node row
        {
            for (int j=0;j<4;j++) 	//each node column
            {
                Mq[4*0+i]+=lI[i][j]*qbuff[IND(j,0)];
                Mq[4*1+i]+=lI[i][j]*qbuff[IND(j,1)];
                Mq[4*2+i]+=lI[i][j]*qbuff[IND(j,2)];
                Mq[4*3+i]+=lI[i][j]*qbuff[IND(j,3)];
                Mq[4*4+i]+=lI[i][j]*qbuff[IND(j,4)];
            }
        }

        const double temp = mat_par->u1 / simu->dt;
        for (int i=0;i<20;i++)
        {
            lL[i] =  0.5*lL[i] + Mq[i]*temp;    // current RHS
            for (int j=0 ; j<20 ; j++)
            {
                lK[i][j] = 0.5*lK[i][j] + lI[i][j]*temp;
            }
        }
    }//if(dt!=0)
}// end void localKL

void localKL_NQ(
        double* p,
        idx* tt,
        double lL[20],
idx it,
idx index_to_Neumann,
Mesh*  mesh,
Mesh* surf_mesh,
SolutionVector* v){

    int i;


    memset(lL,0,20*sizeof(double));
    double n[3];
    surf_mesh->CopySurfaceNormal(it,&n[0]);

    double eDet = surf_mesh->getDeterminant(it);
    double Jdet = mesh->getDeterminant(index_to_Neumann);

    // Jacobian
    double xr,xs,xt,yr,ys,yt,zr,zs,zt;
    xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;
    for (i=0; i<4; i++) {
        xr+=sh1r[0][i]*p[(tt[i])*3+0]*1e-6; //  <- tt is reordered volume element
        xs+=sh1s[0][i]*p[(tt[i])*3+0]*1e-6;
        xt+=sh1t[0][i]*p[(tt[i])*3+0]*1e-6;

        yr+=sh1r[0][i]*p[(tt[i])*3+1]*1e-6;
        ys+=sh1s[0][i]*p[(tt[i])*3+1]*1e-6;
        yt+=sh1t[0][i]*p[(tt[i])*3+1]*1e-6;

        zr+=sh1r[0][i]*p[(tt[i])*3+2]*1e-6;
        zs+=sh1s[0][i]*p[(tt[i])*3+2]*1e-6;
        zt+=sh1t[0][i]*p[(tt[i])*3+2]*1e-6;
    }//end for i

    double Jinv[3][3]={{ (zt*ys-yt*zs)/Jdet , (xt*zs-zt*xs)/Jdet , (xs*yt-ys*xt)/Jdet}
                       ,{ (yt*zr-zt*yr)/Jdet , (zt*xr-xt*zr)/Jdet , (xt*yr-yt*xr)/Jdet}
                       ,{ (yr*zs-ys*zr)/Jdet , (xs*zr-xr*zs)/Jdet , (ys*xr-xs*yr)/Jdet}};

    double dSh[4][3];
    for(i=0;i<4;i++){
        dSh[i][0]=sh1r[0][i]*Jinv[0][0]+sh1s[0][i]*Jinv[1][0]+sh1t[0][i]*Jinv[2][0];
        dSh[i][1]=sh1r[0][i]*Jinv[0][1]+sh1s[0][i]*Jinv[1][1]+sh1t[0][i]*Jinv[2][1];
        dSh[i][2]=sh1r[0][i]*Jinv[0][2]+sh1s[0][i]*Jinv[1][2]+sh1t[0][i]*Jinv[2][2];
    }//end for i

    for (int igp=0; igp<ngps; igp++) {
        //double Sh[4];
        //Sh[0] = sh1[igp][0];
        //Sh[1] = sh1[igp][1];
        //Sh[2] = sh1[igp][2];
        //Sh[3] = sh1[igp][3];

        double Vx=0,Vy=0,Vz=0;
        for(i=0;i<4;i++){
            // voltages
            Vx+=dSh[i][0]*v->getValue(tt[i]);
            Vy+=dSh[i][1]*v->getValue(tt[i]);
            Vz+=dSh[i][2]*v->getValue(tt[i]);
        }//end for i

        double mul=wsurf[igp]*eDet;
        double ShR;
        for (i=0; i<4; i++) {
            ShR = sh1[igp][i]*mul;
            lL[i+0] += (rt6* ShR*(Vx*n[0]+Vy*n[1]-2.0*Vz*n[2])*efe/6.0);
            lL[i+4] += (-rt2*ShR*(Vx*n[0]-Vy*n[1])*efe/2.0);
            lL[i+8] += (-rt2*ShR*(Vy*n[0]+Vx*n[1])*efe/2.0);
            lL[i+12]+= (-rt2*ShR*(Vz*n[1]+Vy*n[2])*efe/2.0);
            lL[i+16]+= (-rt2*ShR*(Vz*n[0]+Vx*n[2])*efe/2.0);
        }//end for i
    }//end for igp
}
// end void localKL_NQ

// ASSEMBLE LOCAL WEAK ANCHORING SURFACES
void wk_localKL(
        Mesh* e,
        idx element_num ,
        SolutionVector* q,
        double lL[15],
double lK[15][15],
idx FixLCNumber,
Alignment* alignment,
LC* lc ,
double* NodeNormals)
{		
    double Ss=lc->S0;
    double 	W = alignment->getStrength(FixLCNumber);
    W = W / (3*Ss);// SCALE TO RP ANCHORING

    double K1 	= alignment->getK1(FixLCNumber);
    double K2   = alignment->getK2(FixLCNumber);
    if ( alignment->getUsesSurfaceNormal(FixLCNumber) ) // if degenerate -> ignore K2 value
    {
        K1 = 1.0;
        K2 = 0;
        if (W<0) // if weak homeotropic
            Ss = -0.5*Ss;
    }

    double A = (K1+K2)/(Ss*6.0);
    double v1[3]= {0,0,0};
    double v2[3]= {0,0,0};

    if (! alignment->getUsesSurfaceNormal(FixLCNumber) ) // if a non-degenerate surface
    {
        v1[0] = * ( alignment->getPtrTov1(FixLCNumber) + 0);
        v1[1] = * ( alignment->getPtrTov1(FixLCNumber) + 1);
        v1[2] = * ( alignment->getPtrTov1(FixLCNumber) + 2);
        v2[0] = * ( alignment->getPtrTov2(FixLCNumber) + 0);
        v2[1] = * ( alignment->getPtrTov2(FixLCNumber) + 1);
        v2[2] = * ( alignment->getPtrTov2(FixLCNumber) + 2);

        //	printf("v1= [%f,%f,%f], v2 =[%f,%f,%f]\n", v1[0], v1[1], v1[2], v2[0], v2[1],v2[2] );
        //	printf("Easy = [%f,%f,%f]\n", alignment->Easy[0], alignment->Easy[1], alignment->Easy[2]);
    }

    memset(lK,0,15*15*sizeof(double));	//SET LOCAL MATRICES TO ZERO
    memset(lL,0,15*sizeof(double));


    for (int igp=0;igp<ngps;igp++){
        double q1,q2,q3,q4,q5,v1x,v1y,v1z,v2x,v2y,v2z;
        q1=q2=q3=q4=q5=v1x=v1y=v1z=v2x=v2y=v2z=0;
        for (int i = 0 ;i < 3 ; i++){
            // Q-tensor components with shape functions
            q1+=ssh1[igp][i]*q->getValue(e->getNode( element_num , i ) , 0);
            q2+=ssh1[igp][i]*q->getValue(e->getNode( element_num , i ) , 1);
            q3+=ssh1[igp][i]*q->getValue(e->getNode( element_num , i ) , 2);
            q4+=ssh1[igp][i]*q->getValue(e->getNode( element_num , i ) , 3);
            q5+=ssh1[igp][i]*q->getValue(e->getNode( element_num , i ) , 4);

            // vector components with shape functions
            if (! alignment->getUsesSurfaceNormal(FixLCNumber) ){ // if a non-degenerate surface
                v1x+=ssh1[igp][i]*v1[0];// v1---- x,y and z-components
                v1y+=ssh1[igp][i]*v1[1];
                v1z+=ssh1[igp][i]*v1[2];

                v2x+=ssh1[igp][i]*v2[0];// v2 ---- x,y and z-components
                v2y+=ssh1[igp][i]*v2[1];
                v2z+=ssh1[igp][i]*v2[2];
            }
            else // use node normals
            {
                v1x+=ssh1[igp][i]*NodeNormals[e->getNode(element_num,i)*3 + 0];
                v1y+=ssh1[igp][i]*NodeNormals[e->getNode(element_num,i)*3 + 1];
                v1z+=ssh1[igp][i]*NodeNormals[e->getNode(element_num,i)*3 + 2];
            }
        }//end for i

        //Tii = thermotropic stiffness term
        double	Tii=2*A*W;
        //Ti = thermotropic RHS vector terms
        double	T1=Tii*q1;
        double  T2=Tii*q2;
        double 	T3=Tii*q3;
        double 	T4=Tii*q4;
        double 	T5=Tii*q5;

        double S1v1, S2v1, S3v1, S4v1, S5v1;
        //vector v1 - terms
        S1v1 = (-v1x*v1x*rt6  - v1y*v1y*rt6 + 2*v1z*v1z*rt6)*W*K1/6.0;
        S2v1 = (v1x*v1x*rt2   - v1y*v1y*rt2)*W*K1/2.0;
        S3v1 = v1x * rt2 * v1y*W*K1;
        S4v1 = v1y * rt2 * v1z*W*K1;
        S5v1 = v1x * rt2 * v1z*W*K1;
        //vector v2 - terms
        double S1v2=0, S2v2=0, S3v2=0, S4v2=0, S5v2=0;
        if ( !alignment->getUsesSurfaceNormal(FixLCNumber) ){ // only needed for non-degenerate surfaces
            S1v2 = (-v2x*v2x*rt6 - v2y*v2y*rt6 +  2*v2z*v2z*rt6)*W*K2/6.0;
            S2v2 = (v2x*v2x*rt2  - v2y*v2y*rt2)*W*K2/2.0;
            S3v2 = v2x * rt2 * v2y *W*K2;
            S4v2 = v2y * rt2 * v2z *W*K2;
            S5v2 = v2x * rt2 * v2z *W*K2;
        }

        double		mul=wsurf[igp]*e->getDeterminant(element_num);

        for (int i=0;i<3;i++){
            //LOCAL RHS
            lL[i+0] +=ssh1[igp][i]*mul*(S1v1  + S1v2 + T1);//
            lL[i+3] +=ssh1[igp][i]*mul*(S2v1  + S2v2 + T2);//
            lL[i+6] +=ssh1[igp][i]*mul*(S3v1  + S3v2 + T3);//
            lL[i+9] +=ssh1[igp][i]*mul*(S4v1  + S4v2 + T4);//
            lL[i+12]+=ssh1[igp][i]*mul*(S5v1  + S5v2 + T5);//

            for ( int j=0;j<3;j++){
                double Sh2=ssh1[igp][i]*ssh1[igp][j]*mul;
                //K MATRIX
                lK[i][j   ]   += Sh2*Tii;
                lK[i+3][j+3 ] += Sh2*Tii;
                lK[i+6][j+6]  += Sh2*Tii;
                lK[i+9][j+9]  += Sh2*Tii;
                lK[i+12][j+12]+= Sh2*Tii;
            }//end for j
        }//end for i
    }//end for igp
}//end void wk_localKL


void assemble_volumes(
        SpaMtrix::IRCMatrix &K,
        SpaMtrix::Vector &L,
        SolutionVector* q,
        SolutionVector* v,
        Mesh* t, double* p,
        LC* mat_par,
        Simu* simu){

    idx npLC = q->getnDoF();
    //init_shape();
    Shape4thOrder shapes;
#ifndef DEBUG
#pragma omp parallel for
#endif
    // LOOP OVER EACH ELEMENT  it
    for (idx it= 0 ; it < t->getnElements () ;it++){
        // IF THIS ELEMENT IS LC ELEMENT, ASSEMBLE LOCAL MATRIX
        if( t->getMaterialNumber(it) != MAT_DOMAIN1 ){ // if LC element
            continue;
        }
        double lK[20][20]; 	// local element matrix
        double lL[20];		// local RHS vector
        localKL(p,t,it,q,v,lK,lL,mat_par, simu, shapes);

        // ADD LOCAL MATRIX TO GLOBAL MATRIX
        for (int i=0;i<20;i++){  // LOOP OVER ROWS
            int ri = t->getNode(it,i%4) + npLC*(i/4);   // LOCAL TO GLOBAL
            idx eqr = q->getEquNode(ri);    // eqr IS MAPPED INDEX TO GLOBAL MATRIX ROW
            if ( eqr != NOT_AN_INDEX ){ // ONLY FOR NON-FIXED NODES
#ifndef DEBUG
#pragma omp atomic
#endif
                L[eqr]+= lL[i]*BIGNUM;
                for (int j = 0 ; j < 20 ; j++){ // LOOP OVER COLUMNS
                    idx rj = t->getNode(it,j%4) + npLC*(j/4);
                    idx eqc = q->getEquNode( rj );
                    if ( eqc != NOT_AN_INDEX ){ // IF NOT FIXED
                        K.sparse_add(eqr,eqc,lK[i][j]*BIGNUM);
                    }
                }
            }// end NON-FIXED NODE
        }//end for i
    }//end fr it
}
// end void assemble_volumes

void assemble_Neumann_surfaces(
        SpaMtrix::Vector &L,
        SolutionVector* q,
        SolutionVector* v,
        Mesh* mesh,
        Mesh* surf_mesh,
        double* p)
{
    int npLC = q->getnDoF();
    init_shape_N();
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (idx it=0; it< surf_mesh->getnElements(); it++) // LOOP OVER EVERY SURFACE ELEMENT
    {
        // ONLY TRIS CONNECTED TO LC TETS ARE ASSEMBLED
        int index_to_Neumann = surf_mesh->getConnectedVolume(it);
        if ( (index_to_Neumann > -1) &&                             // IF CONNECTED TO LC TET
             (surf_mesh->getMaterialNumber(it) != MAT_PERIODIC))    // IF THIS IS NOT A PERIODIC TRIANGLE
        {
            double lL[20];

            // ELEMENT NODE NUMBERS ARE RE-ORDERED SO THAT t[4] IS NOT PART OF TRI ELEMENT
            idx ee[3] = {   surf_mesh->getNode(it,0) ,
                            surf_mesh->getNode(it,1) ,
                            surf_mesh->getNode(it,2) } ;

            idx tt[4] = {   mesh->getNode(index_to_Neumann,0),
                            mesh->getNode(index_to_Neumann,1),
                            mesh->getNode(index_to_Neumann,2),
                            mesh->getNode(index_to_Neumann,3)};

            idx intr=-1;//find  index to internal node
            for (idx i=0;i<4;i++)
            {
                if ( (tt[i]!= ee[0]) && (tt[i]!= ee[1]) && (tt[i]!= ee[2]) )
                {
                    intr = i;
                    break;
                }
            }
            idx ti[4] = { ee[0], ee[1], ee[2], tt[intr] }; // REORDER LOCAL TET ELEMENT
            // NODE-NUMBERING SO THAT
            // INTERNAL NODE IS ALWAYS LAST

            localKL_NQ(p, tt, lL , it , index_to_Neumann,mesh, surf_mesh, v);

            for (unsigned int i=0; i<20; i++) // LOOP OVER ROWS
            {
                idx ri = ti[i%4] + npLC*(i/4);
                idx eqr = q->getEquNode(ri);


                if (eqr != NOT_AN_INDEX ) // IF NOT FIXED
                {
#ifndef DEBUG
#pragma omp atomic
#endif
                    L[eqr]+=lL[i]*BIGNUM;
                }
            }//end for i
        }//end if LC
    }//end for it
}
//end void assemble_Neumann


// ASSEMBLE WEAK ANCHORING SURFACES
void assemble_surfaces(
        SpaMtrix::IRCMatrix &K ,
        SpaMtrix::Vector &L ,
        SolutionVector* q ,
        Mesh* e ,
        LC* lc ,
        Alignment* alignment,
        double* NodeNormals){

    init_shape_surf();
    int npLC = q->getnDoF();
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (idx ie = 0 ; ie < e->getnElements() ; ie ++){
        int FixLCNum = e->getFixLCNumber(ie); // gets FixLC number for surface element ie
        if ((FixLCNum > 0) && ( !alignment->IsStrong(FixLCNum-1) ) ){ // if alignment surface
            double lK[15][15];
            double lL[15];




            wk_localKL( e , ie , q , lL , lK , FixLCNum , alignment, lc , NodeNormals);
            for (unsigned int i=0;i<15;i++){// LOOP OVER ROWS
                idx ri 	= e->getNode(ie,i%3) + npLC*(i/3);
                idx eqr = q->getEquNode(ri);
                if (eqr != NOT_AN_INDEX ) // IF NOT FIXED
                {
#ifndef DEBUG
#pragma omp atomic
#endif

                    L[eqr]+=lL[i]*2e16;

                    for (unsigned int j=0;j<15;j++) // LOOP OVER COLUMNS
                    {

                        idx rj  = e->getNode(ie,j%3) + npLC*(j/3);
                        idx eqc = q->getEquNode(rj);

                        if ( eqc != NOT_AN_INDEX )// IF NOT FIXED
                        {
                            int ii = i; // SURFACE CONTRIBUTION MATRIX IS
                            int jj = j; // SYMMETRIC -> ONLY UPPER DIAGONAL
                            if (j<i)    // IS ASSEMBLED
                            {
                                ii = j;
                                jj = i;
                            }

                            K.sparse_add(eqr, eqc, lK[ii][jj]*BIGNUM);
                        } // end if j not fixed
                    }// end for j
                }// end if i not fixed
            }//end for i
        }// end if alignment surfce
    }// end for ie, loop over surface elements
}
//*/
void assembleQ(
        SpaMtrix::IRCMatrix &K,
        SpaMtrix::Vector &L,  // current RHS
        SolutionVector *q,  // current Q-Tensor
        SolutionVector* v,
        Mesh* t,
        Mesh* e,
        double* p,
        LC* mat_par,
        Simu* simu,
        Alignment* alignment,
        double* NodeNormals)
{
    S0	= mat_par->S0;
    L1 = mat_par->L1 ;
    L2 = mat_par->L2 ;
    L4 = mat_par->L4 ;
    L6 = mat_par->L6 ;
    A = mat_par->A;
    B = mat_par->B;
    C = mat_par->C;
    npLC = (unsigned int) q->getnDoF();
    deleps = ( mat_par->eps_par - mat_par->eps_per ) / S0;
    efe  = (2.0/S0/3.0)*(mat_par->e11+2*mat_par->e33);
    efe2 = (4.0/S0/9.0)*(mat_par->e11 - mat_par->e33);






    assemble_volumes(K, L, q,  v, t, p, mat_par, simu);
    //SHOULD ADD CHECK TO WHETHER NEUMANN SURFACES ACTUALLY EXIST
    assemble_Neumann_surfaces( L, q, v, t, e, p);
    if ( alignment->WeakSurfacesExist() ) // if weak anchoring surfaces exist
        assemble_surfaces(K , L , q ,  e , mat_par ,  alignment, NodeNormals);


}
// end void assembleQ

// Assembles previous time step part of RHS when doing non-linear Crank-Nicholson
/*================================================================*/
/*  ASSEMBLES LOCAL ELEMENT CONTRIBUTIONS FROM PREVIOUS TIME-STEP */
/*================================================================*/

inline void assemble_local_prev_volumes(   double lL[20],
SolutionVector& q,  SolutionVector& v,
Mesh& t, double* p,  idx element_num,
LC& mat_par, Simu& simu,
const Shape4thOrder& shapes)
{


    memset(lL,0,20*sizeof(double));
    double lI[20][20];
    memset(lI,0,20*20*sizeof(double));

    idx tt[4] = {   t.getNode(element_num,0), t.getNode(element_num,1),
                    t.getNode(element_num,2), t.getNode(element_num,3)};

    double Jdet = t.getDeterminant(element_num);

    bool three_elastic_constants = false;
    if ((L2!=0)&&(L6!=0)) three_elastic_constants = true;



    double Kc[5][5];
    memset(Kc,0,5*5*sizeof(double));



    //1. Calculate Inverse Jacobian - for 1st order elements can be done outside integration loop -> igp = 0
    double xr,xs,xt,yr,ys,yt,zr,zs,zt;
    xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;
    for (int i=0;i<4;i++) {
        xr+=shapes.sh1r[0][i]*p[(tt[i])*3+0]*1e-6;
        xs+=shapes.sh1s[0][i]*p[(tt[i])*3+0]*1e-6;
        xt+=shapes.sh1t[0][i]*p[(tt[i])*3+0]*1e-6;

        yr+=shapes.sh1r[0][i]*p[(tt[i])*3+1]*1e-6;
        ys+=shapes.sh1s[0][i]*p[(tt[i])*3+1]*1e-6;
        yt+=shapes.sh1t[0][i]*p[(tt[i])*3+1]*1e-6;

        zr+=shapes.sh1r[0][i]*p[(tt[i])*3+2]*1e-6;
        zs+=shapes.sh1s[0][i]*p[(tt[i])*3+2]*1e-6;
        zt+=shapes.sh1t[0][i]*p[(tt[i])*3+2]*1e-6;
    }//end for i
    //Inverse Jacobian
    Jdet = fabs(xr*ys*zt-xr*zs*yt+xs*yt*zr-xs*yr*zt+xt*yr*zs-xt*ys*zr);
    const double Jinv[3][3]={{(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet}
                             ,{(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet}
                             ,{(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet}};


    // shape function derivatives
    double dSh[4][3] = {{0}};
    for(int i=0;i<4;i++){
        dSh[i][0]=sh1r[0][i]*Jinv[0][0]+sh1s[0][i]*Jinv[1][0]+sh1t[0][i]*Jinv[2][0];
        dSh[i][1]=sh1r[0][i]*Jinv[0][1]+sh1s[0][i]*Jinv[1][1]+sh1t[0][i]*Jinv[2][1];
        dSh[i][2]=sh1r[0][i]*Jinv[0][2]+sh1s[0][i]*Jinv[1][2]+sh1t[0][i]*Jinv[2][2];
    }//end for i
    const double qbuff[6*4] = { q.getValue(tt[0],0), q.getValue(tt[0],1), q.getValue(tt[0],2), q.getValue(tt[0],3), q.getValue(tt[0],4),
                                q.getValue(tt[1],0), q.getValue(tt[1],1), q.getValue(tt[1],2), q.getValue(tt[1],3), q.getValue(tt[1],4),
                                q.getValue(tt[2],0), q.getValue(tt[2],1), q.getValue(tt[2],2), q.getValue(tt[2],3), q.getValue(tt[2],4),
                                q.getValue(tt[3],0), q.getValue(tt[3],1), q.getValue(tt[3],2), q.getValue(tt[3],3), q.getValue(tt[3],4),
                                v.getValue(tt[0])  , v.getValue(tt[1])  , v.getValue(tt[2])  , v.getValue(tt[3])};
    for (unsigned int igp = 0 ; igp < shapes.ngp ; igp ++){
        //shape function
        const double Sh[4]={  shapes.sh1[igp][0], shapes.sh1[igp][1],
                              shapes.sh1[igp][2], shapes.sh1[igp][3] };

        // Function variables and derivatives
        double q1=0, q2=0, q3=0, q4=0, q5=0;
        double q1x=0,q2x=0,q3x=0,q4x=0,q5x=0;
        double q1y=0,q2y=0,q3y=0,q4y=0,q5y=0;
        double q1z=0,q2z=0,q3z=0,q4z=0,q5z=0;
        double Vx=0,Vy=0,Vz=0;

        // Solution and derivatives
#define IND(i,j) 5*(i) + (j)
        for(int i=0;i<4;i++) {
            q1+=Sh[i]*qbuff[IND(i,0)];// PREFETCHED Q AND V TO LOCAL BUFFER OUTSIDE LOOP
            q2+=Sh[i]*qbuff[IND(i,1)];
            q3+=Sh[i]*qbuff[IND(i,2)];
            q4+=Sh[i]*qbuff[IND(i,3)];
            q5+=Sh[i]*qbuff[IND(i,4)];

            q1x+=dSh[i][0]*qbuff[IND(i,0)];
            q2x+=dSh[i][0]*qbuff[IND(i,1)];
            q3x+=dSh[i][0]*qbuff[IND(i,2)];
            q4x+=dSh[i][0]*qbuff[IND(i,3)];
            q5x+=dSh[i][0]*qbuff[IND(i,4)];

            q1y+=dSh[i][1]*qbuff[IND(i,0)];
            q2y+=dSh[i][1]*qbuff[IND(i,1)];
            q3y+=dSh[i][1]*qbuff[IND(i,2)];
            q4y+=dSh[i][1]*qbuff[IND(i,3)];
            q5y+=dSh[i][1]*qbuff[IND(i,4)];

            q1z+=dSh[i][2]*qbuff[IND(i,0)];
            q2z+=dSh[i][2]*qbuff[IND(i,1)];
            q3z+=dSh[i][2]*qbuff[IND(i,2)];
            q4z+=dSh[i][2]*qbuff[IND(i,3)];
            q5z+=dSh[i][2]*qbuff[IND(i,4)];
            // voltages
            Vx+=dSh[i][0]*qbuff[20+i];
            Vy+=dSh[i][1]*qbuff[20+i];
            Vz+=dSh[i][2]*qbuff[20+i];
        }//end for i
#undef IND
        const double R=q1*q1+q2*q2+q3*q3+q5*q5+q4*q4; // frequently reoccurring term
        const double mul=shapes.w[igp]*Jdet;

        // BULK RHS TERMS
        double L[5] = {0,0,0,0,0};
        RHS_THERMOTROPIC(L);
        RHS_DIELECTRIC(L);
        ADD_RHS_BULK_TERMS(lL,L);

        for (int i=0;i<4;i++){ // matrix rows
            const double ShRx=mul*dSh[i][0];//including weight and jacobian in trial function
            const double ShRy=mul*dSh[i][1];
            const double ShRz=mul*dSh[i][2];
            const double ShR =mul*Sh[i];
            // L1 ELASTIC TERM
            lL[i   ] += (ShRx*q1x+ShRy*q1y+ShRz*q1z)*L1;
            lL[i+4 ] += (ShRx*q2x+ShRy*q2y+ShRz*q2z)*L1;
            lL[i+8 ] += (ShRx*q3x+ShRy*q3y+ShRz*q3z)*L1;
            lL[i+12] += (ShRx*q4x+ShRy*q4y+ShRz*q4z)*L1;
            lL[i+16] += (ShRx*q5x+ShRy*q5y+ShRz*q5z)*L1;

            //Chiral term
            if (L4!=0) {
                lL[i   ] += ((rt3*q5y/4.0-rt3*q4x/4.0)*ShR+ShRx*q4*rt3/4.0-ShRy*q5*rt3/4.0)*L4;
                lL[i+4 ] += ((q3z/2.0-q5y/4.0-q4x/4.0)*ShR+ShRx*q4/4.0+ShRy*q5/4.0-ShRz*q3/2.0)*L4;
                lL[i+8 ] += ((-q2z/2.0-q4y/4.0+q5x/4.0)*ShR+ShRy*q4/4.0+ShRz*q2/2.0-ShRx*q5/4.0)*L4;
                lL[i+12] += ((-q5z/4.0+q1x*rt3/4.0+q3y/4.0+q2x/4.0)*ShR-ShRx*q1*rt3/4.0-ShRx*q2/4.0-ShRy*q3/4.0+ShRz*q5/4.0)*L4;
                lL[i+16] += ((q4z/4.0-q1y*rt3/4.0+q2y/4.0-q3x/4.0)*ShR+ShRx*q3/4.0+ShRy*q1*rt3/4.0-ShRy*q2/4.0-ShRz*q4/4.0)*L4;
            }
            if (three_elastic_constants){
                double L2_1,L2_2,L2_3,L2_4,L2_5;//L2 elastic RHS  terms;
                double L6_1,L6_2,L6_3,L6_4,L6_5;//L6 elastic RHS  terms;
                L2_1 = (ShRx*q1x/6.0-ShRx*rt3*q2x/6.0-ShRx*rt3*q3y/6.0-ShRx*rt3*q5z/6.0-ShRy*q3x*rt3/6.0+ShRy*q1y/6.0+ShRy*rt3*q2y/6.0-	ShRy*rt3*q4z/6.0+ShRz*q5x*rt3/3.0	+ShRz*q4y*rt3/3.0+2.0/3.0*ShRz*q1z)	*L2;
                L2_2 = (-ShRx*q1x*rt3/6.0+q2x*ShRx/2.0+q3y*ShRx/2.0+q5z*ShRx/2.0-q3x*ShRy/2.0+ShRy*q1y*rt3/6.0+q2y*ShRy/2.0-q4z*ShRy/2.0)*L2;
                L2_3 = (ShRx*q3x/2.0-ShRx*q1y*rt3/6.0-ShRx*q2y/2.0+ShRx*q4z/2.0-ShRy*q1x*rt3/6.0+ShRy*q2x/2.0+ShRy*q3y/2.0+ShRy*q5z/2.0)*L2;
                L2_4 = (ShRy*q5x/2.0+ShRy*q4y/2.0+ShRy*q1z*rt3/3.0+ShRz*q3x/2.0-	ShRz*q1y*rt3/6.0-ShRz*q2y/2.0+ShRz*q4z/2.0)*L2;
                L2_5 = (ShRx*q5x/2.0+ShRx*q4y/2.0+ShRx*q1z*rt3/3.0-ShRz*q1x*rt3/6.0	+ShRz*q2x/2.0+ShRz*q3y/2.0+ShRz*q5z/2.0)*L2;
                L6_1 = (-ShRx*q1x*q1*rt2*rt3/6.0-ShRy*q1y*q1*rt2*rt3/6.0+ShRz*q1*rt2*rt3*q1z/3.0-ShR*rt2*rt3*q5y*q5y/12.0+ShRx*q5*rt2*q1z/2.0+ShRy*q3*rt2*q1x/2.0+ShRy*q4*rt2*q1z/2.0+ShRz*q5*rt2*q1x/2.0-ShR*rt2*rt3*q4x*q4x/12.0-ShR*rt2*rt3*q4y*q4y/12.0-ShR*rt2*rt3*q3x*q3x/12.0-ShR*rt2*rt3*q2y*q2y/12.0+ShR*rt2*rt3*q2z*q2z/6.0+ShR*rt2*rt3*q3z*q3z/6.0+ShR*rt2*rt3*q5z*q5z/6.0+ShR*rt2*rt3*q4z*q4z/6.0+ShR*rt2*rt3*q1z*q1z/6.0+ShRz*q4*rt2*q1y/2.0-ShR*rt2*rt3*q1x*q1x/12.0-ShR*rt2*rt3*q3y*q3y/12.0-ShR*rt2*rt3*q2x*q2x/12.0-ShR*rt2*rt3*q1y*q1y/12.0+ShRx*q1x*q2*rt2/2.0+ShRx*q3*rt2*q1y/2.0-ShRy*q1y*q2*rt2/2.0-ShR*rt2*rt3*q5x*q5x/12.0)*L6;
                L6_2 = (-ShR*rt2*q5y*q5y/4.0-ShR*rt2*q1y*q1y/4.0+ShR*rt2*q4x*q4x/4.0+ShR*rt2*q1x*q1x/4.0-ShR*rt2*q3y*q3y/4.0+ShR*rt2*q2x*q2x/4.0+ShR*rt2*q5x*q5x/4.0+ShR*rt2*q3x*q3x/4.0-ShR*rt2*q4y*q4y/4.0-ShR*rt2*q2y*q2y/4.0-ShRx*q1*rt2*rt3*q2x/6.0+ShRx*rt2*q2*q2x/2.0+ShRx*q3*q2y*rt2/2.0+ShRx*q5*q2z*rt2/2.0-ShRy*q1*rt2*rt3*q2y/6.0-ShRy*rt2*q2*q2y/2.0+ShRy*q3*q2x*rt2/2.0+ShRy*q4*q2z*rt2/2.0+ShRz*q5*q2x*rt2/2.0+ShRz*q4*q2y*rt2/2.0+ShRz*q1*rt2*rt3*q2z/3.0)*L6;
                L6_3 = (ShR*rt2*q1x*q1y/2.0+ShR*rt2*q2x*q2y/2.0+ShR*rt2*q3x*q3y/2.0+ShR*rt2*q5x*q5y/2.0+ShR*rt2*q4x*q4y/2.0-ShRx*q3x*q1*rt2*rt3/6.0+ShRx*q3x*q2*rt2/2.0+ShRx*q3*rt2*q3y/2.0+ShRx*q5*rt2*q3z/2.0-ShRy*q3y*q1*rt2*rt3/6.0-ShRy*q3y*q2*rt2/2.0+ShRy*q3*rt2*q3x/2.0+ShRy*q4*rt2*q3z/2.0+ShRz*q5*rt2*q3x/2.0+ShRz*q4*rt2*q3y/2.0+ShRz*q1*rt2*rt3*q3z/3.0)*L6;
                L6_4 = (ShR*rt2*q1y*q1z/2.0+ShR*rt2*q2y*q2z/2.0+ShR*rt2*q3y*q3z/2.0+ShR*rt2*q5y*q5z/2.0+ShR*rt2*q4y*q4z/2.0-ShRx*q4x*q1*rt2*rt3/6.0+ShRx*q4x*q2*rt2/2.0+ShRx*q3*rt2*q4y/2.0+ShRx*q5*rt2*q4z/2.0-ShRy*q4y*q1*rt2*rt3/6.0-ShRy*q4y*q2*rt2/2.0+ShRy*q3*rt2*q4x/2.0+ShRy*q4*rt2*q4z/2.0+ShRz*q5*rt2*q4x/2.0+ShRz*q4*rt2*q4y/2.0+ShRz*q1*rt2*rt3*q4z/3.0)*L6;
                L6_5 = (ShR*rt2*q1x*q1z/2.0+ShR*rt2*q2x*q2z/2.0+ShR*rt2*q3x*q3z/2.0+ShR*rt2*q5x*q5z/2.0+ShR*rt2*q4x*q4z/2.0-ShRx*q5x*q1*rt2*rt3/6.0+ShRx*q5x*q2*rt2/2.0+ShRx*q3*rt2*q5y/2.0+ShRx*q5*rt2*q5z/2.0-ShRy*q5y*q1*rt2*rt3/6.0-ShRy*q5y*q2*rt2/2.0+ShRy*q3*rt2*q5x/2.0+ShRy*q4*rt2*q5z/2.0+ShRz*q5*rt2*q5x/2.0+ShRz*q4*rt2*q5y/2.0+ShRz*q1*rt2*rt3*q5z/3.0)*L6;
                lL[i+0]  +=  L2_1 + L6_1;
                lL[i+4]  +=  L2_2 + L6_2;
                lL[i+8]  +=  L2_3 + L6_3;
                lL[i+12] +=  L2_4 + L6_4;
                lL[i+16] +=  L2_5 + L6_5;

            }
            if ((efe!=0.0)||(efe2!=0.0)){ // IF FLEXOELECTRIC COEFFICIENTS ARN'T 0
                lL[i+0]  += (rt6*(Vx*ShRx+Vy*ShRy-2.0*Vz*ShRz)*efe/6.0);
                lL[i+4]  += (-rt2*(Vx*ShRx-Vy*ShRy)*efe/2.0);
                lL[i+8]  += (-rt2*(Vy*ShRx+Vx*ShRy)*efe/2.0);
                lL[i+12] += (-rt2*(Vz*ShRy+Vy*ShRz)*efe/2.0);
                lL[i+16] += (-rt2*(Vz*ShRx+Vx*ShRz)*efe/2.0);
            }

            // ASSEMBLE LOCAL MASS MATRIX
            for (int j = 0 ; j < 4 ; j++){
                double ShRC = ShR*Sh[j];
                lI[i   ][j   ]+=ShRC;
                lI[i+4 ][j+4 ]+=ShRC;
                lI[i+8 ][j+8 ]+=ShRC;
                lI[i+12][j+12]+=ShRC;
                lI[i+16][j+16]+=ShRC;
            }
        } // end for i  - rows
    }//end for igp

    // MAKE CRANK-NICHOLSON RHS TERM
    // the product of the mass matrix and the various q vectors
    double Mq[20]={0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0};	// M * q

    for (int i=0;i<4;i++) {		//each node row
        for (int j=0;j<4;j++) {	//each node column
            for (int k=0;k<5;k++){//each component
                Mq[4*k+i]+=lI[i][j]*q.Values[npLC*k+tt[j]];
            }
        }
    }
    double dt = simu.dt;
    double u1 = mat_par.u1;
    for (int i=0;i<20;i++){
        /* SEGFAULT WITH OPENMP HERE*/
        lL[i] =   ( lL[i] / 2.0 ) -  ( Mq[i]*(u1 / dt) ) ;   // ORIGNAL
        // lL[i] =   ( lL[i] / 2.0 ) +  ( Mq[i]*(u1 / dt) ) ;   // M*current Q
    }
}// end local rhs prev


/*=====================================================*/
/*  ASSEMBLES RHS CONTRIBUTION FROM PREVIOUS TIME-STEP */
/*=====================================================*/

void assemble_prev_rhs(SpaMtrix::Vector &Ln,
                       SolutionVector& qn,
                       SolutionVector& v,
                       LC& mat_par,
                       Simu& simu,
                       Geometry& geom
                       )
{
    init_globals(mat_par, qn);
    Shape4thOrder shapes;
    unsigned int elem_cnt = geom.t->getnElements();//unsigned int) t.getnElements();

    // OPENMP LOOP COMPILED WITH -march=native an -O3 RESULTS IN SEGFAULT ON
    // WINXP32, COMPILED WITHMinGW. THIS IS NOT A PROBLEM WITH UBUNTU,
    // NO PROBLEMS FOUND WITH gdb / valgrind. SGFAULTING LINE MARKED IN FUNTION
    // assemble_local_prev_volumes. ENABLING OPENMP ONLY FOR LINUX (02/12/2010)
    // 07/04/2012 compiling with TDM gcc4.6.1 on win7-64 -> no problems
    //#ifndef __WIN32__
    //#pragma omp parallel for // PARALLEL LOOP IN LINUX
    //#endif
    //int th = 0; // debug thread number
    Mesh& t = *geom.t;
    double* p = geom.getPtrTop();
#ifndef DEBUG
#pragma omp parallel for
#endif
    for ( idx it = 0 ; it < elem_cnt ; it++){
        // IF THIS ELEMENT IS LC ELEMENT, ASSEMBLE LOCAL MATRIX
        if( t.getMaterialNumber(it) == MAT_DOMAIN1 ){// if LC element
            idx eqr;
            double lL[20];		// local RHS vector
            assemble_local_prev_volumes(lL,
                                        qn, v ,
                                        t , p , it,
                                        mat_par , simu,
                                        shapes);
            // ADD LOCAL MATRIX TO GLOBAL MATRIX
            for (unsigned int i=0;i<20;i++){
                int ri = t.getNode(it,i%4) + npLC*(i/4);
                eqr = qn.getEquNode(ri);
                if (eqr != NOT_AN_INDEX ){ // IF NOT FIXED
#ifndef DEBUG
#pragma omp atomic
#endif
                    Ln[eqr] += lL[i]*BIGNUM;
                }
            }// end for i
        }//end if LC material
    }//end for it
}// end function


