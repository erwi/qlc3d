#include <math.h>

#include <omp.h>
#include <time.h>
#include <qlc3d.h>
#include <sparsematrix.h>
//#include "gauss.h"

#define	BIGNUM 2e16
const int	npt = 4; //Number of Points per Tetrahedra
// Gauss integration
// ---------------------------------------------------------
//     3D Gauss-Legendre weights for N = 11, D = 4
// ---------------------------------------------------------

const int ngp = 11;
const double a=(1+sqrt(5.0/14.0))/4.0;
const double b=(1-sqrt(5.0/14.0))/4.0;

static double gp[ngp][4]={	{0.25	  , 0.25	,	0.25	,0.25}
	,{11.0/14.0     ,	1.0/14.0	,	1.0/14.0	,1.0/14.0}
	,{1.0/14.0      ,	11.0/14.0	,	1.0/14.0	,1.0/14.0}
	,{1.0/14.0	  ,	1.0/14.0	,	11.0/14.0   ,1.0/14.0},
	{1.0/14.0  , 1.0/14.0	,	1.0/14.0	,11.0/14.0},
	{a	  , a	,  b , b},
	{a	  , b	,  a , b},
	{a        , b   ,  b , a},
	{b	  , a   ,  a , b},
	{b	  , a   ,  b , a},
	{b	  , b   ,  a , a}};
const double w11 = -74.0/5625.0;
const double w12 = 343.0/45000.0;
const double w13 = 56.0/2250.0;
static double w[ngp]={w11,w12,w12,w12,w12,w13,w13,w13,w13,w13,w13};

// ---------------------------------------------------------
//     2D Gauss-Legendre weights for N = 27, D = 11
// ---------------------------------------------------------
		const int ngps=27;
		static double wsurf[ngps]={ 0.006829866, 0.006829866, 0.006829866, 0.01809227, 0.01809227, 0.01809227, 0.0004635032, 0.0004635032,
									0.0004635032,0.02966149 , 0.02966149 , 0.02966149, 0.03857477, 0.03857477, 0.03857477, 0.02616856,
									0.02616856  ,0.02616856 , 0.02616856 , 0.02616856, 0.02616856, 0.01035383, 0.01035383, 0.01035383,
									0.01035383  ,0.01035383 , 0.01035383 };
static double sgp[ngps][2] = 
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
static double sh1[ngps][4]; // P1 Shape functions
static double sh1r[ngps][4]; // P1 Shape functions r-derivatives
static double sh1s[ngps][4]; // P1 Shape functions s-derivatives
static double sh1t[ngps][4]; //P1 shape functions t-derivative
static double ssh1[ngps][3];	//SURFACE term P1 shape function
	
static double rt2 = sqrt(2.0);
static double rt6 = sqrt(6.0);
static double rt3 = sqrt(3.0);

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

void init_shape()
{
    memset(sh1,0,ngps*4*sizeof(double));
    memset(sh1r,0,ngps*4*sizeof(double));
    memset(sh1s,0,ngps*4*sizeof(double));
    memset(sh1t,0,ngps*4*sizeof(double));
    for (int i=0; i<ngp; i++) {
              //  cout<<"i=" << i << endl;
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
	//printf("sh1[%i][1] = %f\n",i,sh1[i][1]);
	}
}					
void init_shape_N() // initialise Neumann and surface derivative shape functions
{
  //  memset(sh1,0,ngps*4*sizeof(double));
  //  memset(sh1r,0,ngps*4*sizeof(double));
  //  memset(sh1s,0,ngps*4*sizeof(double));
  //  memset(sh1t,0,ngps*4*sizeof(double));


    for (int i=0; i<ngp; i++) { // use surface shape functions ( <-!! ngps = 27, but array defined for gps = 11. WHY NO ERROR?!) FIX!!!
	// P1 Shape functions
        //cout << "i ="<< i << endl;
        cout<<"";
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


void localKL(double* p,Mesh* t, int element_num, SolutionVector* q ,
			SolutionVector* v, double lK[20][20], double lL[20],LC* mat_par, Simu* simu)
{
        //cout << "dt(localKL) = " << simu->getdt() << endl;
        memset(lK,0,20*20*sizeof(double));
	memset(lL,0,20*sizeof(double));
	double lI[20][20];
	memset(lI,0,20*20*sizeof(double));
	
	//printf("tid %i, elem# = %i\n", tid, element_num);
	
	int tt[4] = {0,0,0,0};
	tt[0] = t->getNode(element_num,0);
	tt[1] = t->getNode(element_num,1);
	tt[2] = t->getNode(element_num,2);
	tt[3] = t->getNode(element_num,3);
	
	double Jdet = t->getDeterminant(element_num);
		
	if (Jdet < 0) Jdet = -Jdet;
	
	
	bool three_elastic_constants = false;
	if ((L2!=0)&&(L6!=0)) three_elastic_constants = true;

	
	double L1_1,L1_2,L1_3,L1_4,L1_5;//L1 elastic RHS  terms;
	double L2_1,L2_2,L2_3,L2_4,L2_5;//L2 elastic RHS  terms;
	double L6_1,L6_2,L6_3,L6_4,L6_5;//L6 elastic RHS  terms;
	double Lc[5] = {0,0,0,0,0};
	double V1,V2,V3,V4,V5;	
	
	double T1,T2,T3,T4,T5,T11,T12,T13,T14,T15,T22,T23,T24,T25,T33,T34,T35,T44,T45,T55;
	double L2_K11,L2_K12,L2_K13,L2_K14,L2_K15,L2_K22,L2_K23,L2_K24,L2_K25,
		L2_K33,L2_K34,L2_K35,L2_K44,L2_K45,L2_K55;
	double L6_K11,L6_K12,L6_K13,L6_K14,L6_K15,L6_K22,L6_K23,L6_K24,L6_K25,
		L6_K33,L6_K34,L6_K35,L6_K44,L6_K45,L6_K55;	
	double Kc[5][5];
	memset(Kc,0,5*5*sizeof(double));
	
	
	
	//1. Calculate Inverse Jacobian - for 1st order elements can be done outside integration loop -> igp = 0
		double xr,xs,xt,yr,ys,yt,zr,zs,zt;
		xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;
		for (int i=0;i<4;i++) {
			xr+=sh1r[0][i]*p[(tt[i])*3+0]*1e-6;
			xs+=sh1s[0][i]*p[(tt[i])*3+0]*1e-6;
			xt+=sh1t[0][i]*p[(tt[i])*3+0]*1e-6;

			yr+=sh1r[0][i]*p[(tt[i])*3+1]*1e-6;
			ys+=sh1s[0][i]*p[(tt[i])*3+1]*1e-6;
			yt+=sh1t[0][i]*p[(tt[i])*3+1]*1e-6;

			zr+=sh1r[0][i]*p[(tt[i])*3+2]*1e-6;
			zs+=sh1s[0][i]*p[(tt[i])*3+2]*1e-6;
			zt+=sh1t[0][i]*p[(tt[i])*3+2]*1e-6;
		}//end for i
		//Inverse Jacobian
		Jdet = fabs(xr*ys*zt-xr*zs*yt+xs*yt*zr-xs*yr*zt+xt*yr*zs-xt*ys*zr);
		double Jinv[3][3]={{(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet}
						  ,{(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet}
						  ,{(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet}};
	
	
	// shape function derivatives
	double dSh[4][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};	
	for(int i=0;i<4;i++){
	    dSh[i][0]=sh1r[0][i]*Jinv[0][0]+sh1s[0][i]*Jinv[1][0]+sh1t[0][i]*Jinv[2][0];
	    dSh[i][1]=sh1r[0][i]*Jinv[0][1]+sh1s[0][i]*Jinv[1][1]+sh1t[0][i]*Jinv[2][1];
	    dSh[i][2]=sh1r[0][i]*Jinv[0][2]+sh1s[0][i]*Jinv[1][2]+sh1t[0][i]*Jinv[2][2];
	}//end for i

    for (int igp = 0 ; igp < ngp ; igp ++){
	//shape function
	double Sh[4];
	Sh[0] = sh1[igp][0];
	Sh[1] = sh1[igp][1];
	Sh[2] = sh1[igp][2];
	Sh[3] = sh1[igp][3];
		
	// Function variables and derivatives
	double q1=0, q2=0, q3=0, q4=0, q5=0;
	double q1x=0,q2x=0,q3x=0,q4x=0,q5x=0;
	double q1y=0,q2y=0,q3y=0,q4y=0,q5y=0;
	double q1z=0,q2z=0,q3z=0,q4z=0,q5z=0;
	double Vx=0,Vy=0,Vz=0;

	// Solution and derivatives
	for(int i=0;i<4;i++){
	    q1+=Sh[i]*q->getValue(tt[i],0);				// q1i * Ni  = A1
	    q2+=Sh[i]*q->getValue(tt[i],1);
	    q3+=Sh[i]*q->getValue(tt[i],2);
	    q4+=Sh[i]*q->getValue(tt[i],3);
	    q5+=Sh[i]*q->getValue(tt[i],4);
		  
	    // voltages
	    Vx+=dSh[i][0]*v->getValue(tt[i]);
	    Vy+=dSh[i][1]*v->getValue(tt[i]);
	    Vz+=dSh[i][2]*v->getValue(tt[i]);
	        
	    q1x+=dSh[i][0]*q->getValue(tt[i],0);
	    q2x+=dSh[i][0]*q->getValue(tt[i],1);
	    q3x+=dSh[i][0]*q->getValue(tt[i],2);
	    q4x+=dSh[i][0]*q->getValue(tt[i],3);
	    q5x+=dSh[i][0]*q->getValue(tt[i],4);
			
	    q1y+=dSh[i][1]*q->getValue(tt[i],0);
	    q2y+=dSh[i][1]*q->getValue(tt[i],1);
	    q3y+=dSh[i][1]*q->getValue(tt[i],2);
	    q4y+=dSh[i][1]*q->getValue(tt[i],3);
	    q5y+=dSh[i][1]*q->getValue(tt[i],4);

	    q1z+=dSh[i][2]*q->getValue(tt[i],0);
	    q2z+=dSh[i][2]*q->getValue(tt[i],1);
	    q3z+=dSh[i][2]*q->getValue(tt[i],2);
	    q4z+=dSh[i][2]*q->getValue(tt[i],3);
	    q5z+=dSh[i][2]*q->getValue(tt[i],4);
	}//end for i

	double R=q1*q1+q2*q2+q3*q3+q5*q5+q4*q4; // frequently reoccurring term
    	double mul=w[igp]*Jdet;
	
	T1=((A*q1+B*(q5*q5*rt6/4.0-q3*q3*rt6/2.0-rt6*q2*q2/2.0+q1*q1*rt6/2.0+q4*q4*rt6/4.0)/3.0)+C*R*q1);
	T2=(A*q2+B*(3.0/4.0*q5*q5*rt2-q1*rt6*q2-3.0/4.0*q4*q4*rt2)/3.0+C*R*q2);
	T3=(A*q3+B*(-q3*q1*rt6+3.0/2.0*rt2*q5*q4)/3.0+C*R*q3);
	T4=(A*q4+B*(3.0/2.0*q3*rt2*q5+q4*q1*rt6/2.0-3.0/2.0*q4*q2*rt2)/3.0+C*R*q4);
	T5=(A*q5+B*(q5*q1*rt6/2.0+3.0/2.0*q5*q2*rt2+3.0/2.0*q3*rt2*q4)/3.0+C*R*q5);
		
	V1 =  rt6*(Vx*Vx + Vy*Vy-2.0*Vz*Vz)*deleps/18.0*eps0;
	V2 = -rt2*(Vx*Vx - Vy*Vy)*deleps/6.0*eps0;
	V3 = -rt2*Vx*Vy*deleps/3.0*eps0;
	V4 = -rt2*Vz*Vy*deleps/3.0*eps0;
	V5 = -rt2*Vx*Vz*deleps/3.0*eps0;
		
	T11=(A+B*q1*rt6/3.0+2.0*C*q1*q1+C*R);
	T12=(-B*rt6*q2/3.0+2.0*C*q2*q1);
	T13=(-B*q3*rt6/3.0+2.0*C*q3*q1);
	T14=(B*q4*rt6/6.0+2.0*C*q4*q1);
	T15=(B*q5*rt6/6.0+2.0*C*q5*q1);

	T22=(A-B*q1*rt6/3.0+2.0*C*q2*q2+C*R);
	T23=(2.0*C*q3*q2);
	T24=(-B*q4*rt2/2.0+2.0*C*q4*q2);
	T25=(B*q5*rt2/2.0+2.0*C*q5*q2);

	T33=(A-B*q1*rt6/3.0+2.0*C*q3*q3+C*R);
	T34=(B*q5*rt2/2.0+2.0*C*q4*q3);
	T35=(B*q4*rt2/2.0+2.0*C*q5*q3);

	T44=(A+B*(q1*rt6-3.0*q2*rt2)/6.0+2.0*C*q4*q4+C*R);
	T45=(B*q3*rt2/2.0+2.0*C*q5*q4);
	T55=(A+B*(q1*rt6+3.0*q2*rt2)/6.0+2.0*C*q5*q5+C*R);
		
	double Lflexo1, Lflexo2, Lflexo3, Lflexo4, Lflexo5;
	Lflexo1 = 0; Lflexo2 =0; Lflexo3=0; Lflexo4=0; Lflexo5=0;
		
	for (int i=0;i<4;i++){ // matrix rows
	    double ShRx=mul*dSh[i][0];//including weight and jacobian in trial function
	    double ShRy=mul*dSh[i][1];
	    double ShRz=mul*dSh[i][2];
	    double ShR =mul*Sh[i];
			
	    L1_1=(ShRx*q1x+ShRy*q1y+ShRz*q1z)*L1;
	    L1_2=(ShRx*q2x+ShRy*q2y+ShRz*q2z)*L1;
	    L1_3=(ShRx*q3x+ShRy*q3y+ShRz*q3z)*L1;
	    L1_4=(ShRx*q4x+ShRy*q4y+ShRz*q4z)*L1;
	    L1_5=(ShRx*q5x+ShRy*q5y+ShRz*q5z)*L1;

	    //Chiral term
	    if (L4!=0)
            {
			Lc[0] = ((rt3*q5y/4.0-rt3*q4x/4.0)*ShR+ShRx*q4*rt3/4.0-ShRy*q5*rt3/4.0)*L4;
			Lc[1] = ((q3z/2.0-q5y/4.0-q4x/4.0)*ShR+ShRx*q4/4.0+ShRy*q5/4.0-ShRz*q3/2.0)*L4;
			Lc[2] = ((-q2z/2.0-q4y/4.0+q5x/4.0)*ShR+ShRy*q4/4.0+ShRz*q2/2.0-ShRx*q5/4.0)*L4;
			Lc[3] = ((-q5z/4.0+q1x*rt3/4.0+q3y/4.0+q2x/4.0)*ShR-ShRx*q1*rt3/4.0-ShRx*q2/4.0-ShRy*q3/4.0+ShRz*q5/4.0)*L4;
			Lc[4] = ((q4z/4.0-q1y*rt3/4.0+q2y/4.0-q3x/4.0)*ShR+ShRx*q3/4.0+ShRy*q1*rt3/4.0-ShRy*q2/4.0-ShRz*q4/4.0)*L4;
	    }
	    if (three_elastic_constants){
		L2_1= (ShRx*q1x/6.0-ShRx*rt3*q2x/6.0-ShRx*rt3*q3y/6.0-ShRx*rt3*q5z/6.0-ShRy*q3x*rt3/6.0+ShRy*q1y/6.0+ShRy*rt3*q2y/6.0-	ShRy*rt3*q4z/6.0+ShRz*q5x*rt3/3.0	+ShRz*q4y*rt3/3.0+2.0/3.0*ShRz*q1z)	*L2;
		L2_2 = (-ShRx*q1x*rt3/6.0+q2x*ShRx/2.0+q3y*ShRx/2.0+q5z*ShRx/2.0-q3x*ShRy/2.0+ShRy*q1y*rt3/6.0+q2y*ShRy/2.0-q4z*ShRy/2.0)*L2;
		L2_3 = (ShRx*q3x/2.0-ShRx*q1y*rt3/6.0-ShRx*q2y/2.0+ShRx*q4z/2.0-ShRy*q1x*rt3/6.0+ShRy*q2x/2.0+ShRy*q3y/2.0+ShRy*q5z/2.0)*L2;
		L2_4= (ShRy*q5x/2.0+ShRy*q4y/2.0+ShRy*q1z*rt3/3.0+ShRz*q3x/2.0-	ShRz*q1y*rt3/6.0-ShRz*q2y/2.0+ShRz*q4z/2.0)*L2;
		L2_5 = (ShRx*q5x/2.0+ShRx*q4y/2.0+ShRx*q1z*rt3/3.0-ShRz*q1x*rt3/6.0	+ShRz*q2x/2.0+ShRz*q3y/2.0+ShRz*q5z/2.0)*L2;
		L6_1 = (-ShRx*q1x*q1*rt2*rt3/6.0-ShRy*q1y*q1*rt2*rt3/6.0+ShRz*q1*rt2*rt3*q1z/3.0-ShR*rt2*rt3*q5y*q5y/12.0+ShRx*q5*rt2*q1z/2.0+ShRy*q3*rt2*q1x/2.0+ShRy*q4*rt2*q1z/2.0+ShRz*q5*rt2*q1x/2.0-ShR*rt2*rt3*q4x*q4x/12.0-ShR*rt2*rt3*q4y*q4y/12.0-ShR*rt2*rt3*q3x*q3x/12.0-ShR*rt2*rt3*q2y*q2y/12.0+ShR*rt2*rt3*q2z*q2z/6.0+ShR*rt2*rt3*q3z*q3z/6.0+ShR*rt2*rt3*q5z*q5z/6.0+ShR*rt2*rt3*q4z*q4z/6.0+ShR*rt2*rt3*q1z*q1z/6.0+ShRz*q4*rt2*q1y/2.0-ShR*rt2*rt3*q1x*q1x/12.0-ShR*rt2*rt3*q3y*q3y/12.0-ShR*rt2*rt3*q2x*q2x/12.0-ShR*rt2*rt3*q1y*q1y/12.0+ShRx*q1x*q2*rt2/2.0+ShRx*q3*rt2*q1y/2.0-ShRy*q1y*q2*rt2/2.0-ShR*rt2*rt3*q5x*q5x/12.0)*L6;
		L6_2 = (-ShR*rt2*q5y*q5y/4.0-ShR*rt2*q1y*q1y/4.0+ShR*rt2*q4x*q4x/4.0+ShR*rt2*q1x*q1x/4.0-ShR*rt2*q3y*q3y/4.0+ShR*rt2*q2x*q2x/4.0+ShR*rt2*q5x*q5x/4.0+ShR*rt2*q3x*q3x/4.0-ShR*rt2*q4y*q4y/4.0-ShR*rt2*q2y*q2y/4.0-ShRx*q1*rt2*rt3*q2x/6.0+ShRx*rt2*q2*q2x/2.0+ShRx*q3*q2y*rt2/2.0+ShRx*q5*q2z*rt2/2.0-ShRy*q1*rt2*rt3*q2y/6.0-ShRy*rt2*q2*q2y/2.0+ShRy*q3*q2x*rt2/2.0+ShRy*q4*q2z*rt2/2.0+ShRz*q5*q2x*rt2/2.0+ShRz*q4*q2y*rt2/2.0+ShRz*q1*rt2*rt3*q2z/3.0)*L6;
		L6_3 = (ShR*rt2*q1x*q1y/2.0+ShR*rt2*q2x*q2y/2.0+ShR*rt2*q3x*q3y/2.0+ShR*rt2*q5x*q5y/2.0+ShR*rt2*q4x*q4y/2.0-ShRx*q3x*q1*rt2*rt3/6.0+ShRx*q3x*q2*rt2/2.0+ShRx*q3*rt2*q3y/2.0+ShRx*q5*rt2*q3z/2.0-ShRy*q3y*q1*rt2*rt3/6.0-ShRy*q3y*q2*rt2/2.0+ShRy*q3*rt2*q3x/2.0+ShRy*q4*rt2*q3z/2.0+ShRz*q5*rt2*q3x/2.0+ShRz*q4*rt2*q3y/2.0+ShRz*q1*rt2*rt3*q3z/3.0)*L6;
		L6_4 = (ShR*rt2*q1y*q1z/2.0+ShR*rt2*q2y*q2z/2.0+ShR*rt2*q3y*q3z/2.0+ShR*rt2*q5y*q5z/2.0+ShR*rt2*q4y*q4z/2.0-ShRx*q4x*q1*rt2*rt3/6.0+ShRx*q4x*q2*rt2/2.0+ShRx*q3*rt2*q4y/2.0+ShRx*q5*rt2*q4z/2.0-ShRy*q4y*q1*rt2*rt3/6.0-ShRy*q4y*q2*rt2/2.0+ShRy*q3*rt2*q4x/2.0+ShRy*q4*rt2*q4z/2.0+ShRz*q5*rt2*q4x/2.0+ShRz*q4*rt2*q4y/2.0+ShRz*q1*rt2*rt3*q4z/3.0)*L6;
		L6_5 = (ShR*rt2*q1x*q1z/2.0+ShR*rt2*q2x*q2z/2.0+ShR*rt2*q3x*q3z/2.0+ShR*rt2*q5x*q5z/2.0+ShR*rt2*q4x*q4z/2.0-ShRx*q5x*q1*rt2*rt3/6.0+ShRx*q5x*q2*rt2/2.0+ShRx*q3*rt2*q5y/2.0+ShRx*q5*rt2*q5z/2.0-ShRy*q5y*q1*rt2*rt3/6.0-ShRy*q5y*q2*rt2/2.0+ShRy*q3*rt2*q5x/2.0+ShRy*q4*rt2*q5z/2.0+ShRz*q5*rt2*q5x/2.0+ShRz*q4*rt2*q5y/2.0+ShRz*q1*rt2*rt3*q5z/3.0)*L6;
	    }else{L2_1=L2_2=L2_3=L2_4=L2_5=L6_1=L6_2=L6_3=L6_4=L6_5=0.0;}//end if 3 elestic constants
	    if ((efe!=0.0)||(efe2!=0.0)){ // IF FLEXOELECTRIC COEFFICIENTS ARN'T 0
		Lflexo1 = (rt6*(Vx*ShRx+Vy*ShRy-2.0*Vz*ShRz)*efe/6.0);
		Lflexo2 = (-rt2*(Vx*ShRx-Vy*ShRy)*efe/2.0);
		Lflexo3 = (-rt2*(Vy*ShRx+Vx*ShRy)*efe/2.0);
		Lflexo4 = (-rt2*(Vz*ShRy+Vy*ShRz)*efe/2.0);
		Lflexo5 = (-rt2*(Vz*ShRx+Vx*ShRz)*efe/2.0);
	    }
		
	    //ADD TERMS TO RHS VECTOR
	    lL[i+0]  += T1*ShR + L1_1 + L2_1 + L6_1 + Lc[0] + V1*ShR + Lflexo1;
	    lL[i+4]  += T2*ShR + L1_2 + L2_2 + L6_2 + Lc[1] + V2*ShR + Lflexo2;
	    lL[i+8]  += T3*ShR + L1_3 + L2_3 + L6_3 + Lc[2] + V3*ShR + Lflexo3;
	    lL[i+12] += T4*ShR + L1_4 + L2_4 + L6_4 + Lc[3] + V4*ShR + Lflexo4;
	    lL[i+16] += T5*ShR + L1_5 + L2_5 + L6_5 + Lc[4] + V5*ShR + Lflexo5;
			
	    if ( simu->IsAssembleMatrix() ){
		for (int j=0 ; j<4 ; j++){
		    double ShCx=dSh[j][0];
		    double ShCy=dSh[j][1];
		    double ShCz=dSh[j][2];
		    double ShC =Sh[j];
		    double ShRC=ShR*Sh[j];

		    //L1- matrix terms
		    double dot=L1*mul*(dSh[i][0]*dSh[j][0]+dSh[i][1]*dSh[j][1]+dSh[i][2]*dSh[j][2]);
		    if (three_elastic_constants){//if three elastic constants used
			L2_K11 = (ShRx*ShCx/6.0+ShRy*ShCy/6.0+2.0/3.0*ShRz*ShCz)*L2;
			L2_K12 = (-ShRx*rt3*ShCx/6.0+ShRy*rt3*ShCy/6.0)*L2;
			L2_K13 = (-ShRy*rt3*ShCx/6.0-ShRx*rt3*ShCy/6.0)*L2;
			L2_K14 = (ShRz*rt3*ShCy/3.0-ShRy*rt3*ShCz/6.0)*L2;
			L2_K15 = (ShRz*rt3*ShCx/3.0-ShRx*rt3*ShCz/6.0)*L2;

			L2_K22 = (ShRx*ShCx/2.0+ShRy*ShCy/2.0)*L2;
			L2_K23 = (-ShRy*ShCx/2.0+ShRx*ShCy/2.0)*L2;
			L2_K24 = -ShCz*L2*ShRy/2.0;
			L2_K25 = ShCz*L2*ShRx/2.0;

			L2_K33 = (ShRx*ShCx/2.0+ShRy*ShCy/2.0)*L2;
			L2_K34 = ShCz*L2*ShRx/2.0;
			L2_K35 = ShCz*L2*ShRy/2.0;

			L2_K44 = (ShRy*ShCy/2.0+ShRz*ShCz/2.0)*L2;
			L2_K45 = ShCx*L2*ShRy/2.0;
			L2_K55 = (ShRx*ShCx/2.0+ShRz*ShCz/2.0)*L2;

			//L6 terms
			L6_K11 = (-ShC*ShRx*q1x*rt2*rt3/6.0-ShC*ShRy*q1y*rt2*rt3/6.0+ShC*ShRz*rt2*rt3*q1z/3.0-ShCx*ShRx*q1*rt2*rt3/6.0+ShCx*ShRy*q3*rt2/2.0+ShCx*ShRz*q5*rt2/2.0-ShCx*ShR*rt2*rt3*q1x/6.0+ShCx*ShRx*q2*rt2/2.0-ShCy*ShRy*q1*rt2*rt3/6.0+ShCy*ShRz*q4*rt2/2.0-ShCy*ShR*rt2*rt3*q1y/6.0+ShCy*ShRx*q3*rt2/2.0-ShCy*ShRy*q2*rt2/2.0+ShCz*ShRz*q1*rt2*rt3/3.0+ShCz*ShRx*q5*rt2/2.0+ShCz*ShRy*q4*rt2/2.0+ShCz*ShR*rt2*rt3*q1z/3.0)*L6;
			L6_K12 = (ShC*ShRx*q1x*rt2/2.0-ShC*ShRy*q1y*rt2/2.0-ShR*rt2*rt3*q2x*ShCx/6.0-ShR*rt2*rt3*q2y*ShCy/6.0+ShR*rt2*rt3*q2z*ShCz/3.0)*L6;
			L6_K13 = (ShC*ShRy*rt2*q1x/2.0+ShC*ShRx*rt2*q1y/2.0-ShR*rt2*rt3*q3x*ShCx/6.0-ShR*rt2*rt3*q3y*ShCy/6.0+ShR*rt2*rt3*q3z*ShCz/3.0)*L6;
			L6_K14 = (ShC*ShRy*rt2*q1z/2.0+ShC*ShRz*rt2*q1y/2.0-ShR*rt2*rt3*q4x*ShCx/6.0-ShR*rt2*rt3*q4y*ShCy/6.0+ShR*rt2*rt3*q4z*ShCz/3.0)*L6;
			L6_K15 = (ShC*ShRx*rt2*q1z/2.0+ShC*ShRz*rt2*q1x/2.0-ShR*rt2*rt3*q5x*ShCx/6.0-ShR*rt2*rt3*q5y*ShCy/6.0+ShR*rt2*rt3*q5z*ShCz/3.0)*L6;
      
			L6_K22 = (ShC*ShRx*rt2*q2x/2.0-ShC*ShRy*rt2*q2y/2.0+ShCx*ShR*rt2*q2x/2.0-ShCx*ShRx*q1*rt2*rt3/6.0+ShCx*ShRx*q2*rt2/2.0+ShCx*ShRy*q3*rt2/2.0+ShCx*ShRz*q5*rt2/2.0-ShCy*ShR*rt2*q2y/2.0+ShCy*ShRx*q3*rt2/2.0-ShCy*ShRy*q1*rt2*rt3/6.0-ShCy*ShRy*q2*rt2/2.0+ShCy*ShRz*q4*rt2/2.0+ShCz*ShRx*q5*rt2/2.0+ShCz*ShRy*q4*rt2/2.0+ShCz*ShRz*q1*rt2*rt3/3.0)*L6;
			L6_K23 = (ShC*ShRx*q2y*rt2/2.0+ShC*ShRy*q2x*rt2/2.0+ShR*rt2*q3x*ShCx/2.0-ShR*rt2*q3y*ShCy/2.0)*L6;
			L6_K24 = (ShC*ShRy*q2z*rt2/2.0+ShC*ShRz*q2y*rt2/2.0+ShR*rt2*q4x*ShCx/2.0-ShR*rt2*q4y*ShCy/2.0)*L6;
			L6_K25 = (ShC*ShRx*q2z*rt2/2.0+ShC*ShRz*q2x*rt2/2.0+ShR*rt2*q5x*ShCx/2.0-ShR*rt2*q5y*ShCy/2.0)*L6;
      
			L6_K33 = (ShC*ShRx*rt2*q3y/2.0+ShC*ShRy*rt2*q3x/2.0+ShCx*ShR*rt2*q3y/2.0-ShCx*ShRx*q1*rt2*rt3/6.0+ShCx*ShRx*q2*rt2/2.0+ShCx*ShRy*q3*rt2/2.0+ShCx*ShRz*q5*rt2/2.0+ShCy*ShR*rt2*q3x/2.0+ShCy*ShRx*q3*rt2/2.0-ShCy*ShRy*q1*rt2*rt3/6.0-ShCy*ShRy*q2*rt2/2.0+ShCy*ShRz*q4*rt2/2.0+ShCz*ShRx*q5*rt2/2.0+ShCz*ShRy*q4*rt2/2.0+ShCz*ShRz*q1*rt2*rt3/3.0)*L6;
			L6_K34 = (ShC*ShRy*rt2*q3z/2.0+ShC*ShRz*rt2*q3y/2.0+ShR*rt2*q4y*ShCx/2.0+ShR*rt2*q4x*ShCy/2.0)*L6;
			L6_K35 = (ShC*ShRx*rt2*q3z/2.0+ShC*ShRz*rt2*q3x/2.0+ShR*rt2*q5y*ShCx/2.0+ShR*rt2*q5x*ShCy/2.0)*L6;

			L6_K44 = (ShC*ShRy*rt2*q4z/2.0+ShC*ShRz*rt2*q4y/2.0-ShCx*ShRx*q1*rt2*rt3/6.0+ShCx*ShRx*q2*rt2/2.0+ShCx*ShRy*q3*rt2/2.0+ShCx*ShRz*q5*rt2/2.0+ShCy*ShR*rt2*q4z/2.0+ShCy*ShRx*q3*rt2/2.0-ShCy*ShRy*q1*rt2*rt3/6.0-ShCy*ShRy*q2*rt2/2.0+ShCy*ShRz*q4*rt2/2.0+ShCz*ShR*rt2*q4y/2.0+ShCz*ShRx*q5*rt2/2.0+ShCz*ShRy*q4*rt2/2.0+ShCz*ShRz*q1*rt2*rt3/3.0)*L6;
			L6_K45 = (ShC*ShRx*rt2*q4z/2.0+ShC*ShRz*rt2*q4x/2.0+ShR*rt2*q5z*ShCy/2.0+ShR*rt2*q5y*ShCz/2.0)*L6;

			L6_K55 = (ShC*ShRx*rt2*q5z/2.0+ShC*ShRz*rt2*q5x/2.0+ShCx*ShR*rt2*q5z/2.0-ShCx*ShRx*q1*rt2*rt3/6.0+ShCx*ShRx*q2*rt2/2.0+ShCx*ShRy*q3*rt2/2.0+ShCx*ShRz*q5*rt2/2.0+ShCy*ShRx*q3*rt2/2.0-ShCy*ShRy*q1*rt2*rt3/6.0-ShCy*ShRy*q2*rt2/2.0+ShCy*ShRz*q4*rt2/2.0+ShCz*ShR*rt2*q5x/2.0+ShCz*ShRx*q5*rt2/2.0+ShCz*ShRy*q4*rt2/2.0+ShCz*ShRz*q1*rt2*rt3/3.0)*L6;
		    }
		    else{
			L2_K11=L2_K12=L2_K13=L2_K14=L2_K15=L2_K22=L2_K23=L2_K24=L2_K25=L2_K33=L2_K34=L2_K35=L2_K44=L2_K45=L2_K55=  L6_K11=L6_K12=L6_K13=L6_K14=L6_K15=L6_K22=L6_K23=L6_K24=L6_K25=L6_K33=L6_K34=L6_K35=L6_K44=L6_K45=L6_K55=0.0;
		    }//end if three elastic constants
				
		    // chirality terms
		    if (L4!=0){
			Kc[0][3] = (ShRx*rt3*ShC/4.0-ShR*rt3*ShCx/4.0)*L4;
			Kc[0][4] = (-ShRy*rt3*ShC/4.0+ShR*rt3*ShCy/4.0)*L4;
				
			Kc[1][2] = (-ShRz*ShC/2.0+ShR*ShCz/2.0)*L4;
			Kc[1][3] = (ShC*ShRx/4.0-ShCx*ShR/4.0)*L4;
			Kc[1][4] = (ShC*ShRy/4.0-ShCy*ShR/4.0)*L4;
				
			Kc[2][1] = (ShRz*ShC/2.0-ShR*ShCz/2.0)*L4;
			Kc[2][3] = (ShC*ShRy/4.0-ShCy*ShR/4.0)*L4;
			Kc[2][4] = (-ShC*ShRx/4.0+ShCx*ShR/4.0)*L4;
				
			Kc[3][0] = (-ShRx*rt3*ShC/4.0+ShR*rt3*ShCx/4.0)*L4;
			Kc[3][1] = (-ShC*ShRx/4.0+ShCx*ShR/4.0)*L4;
			Kc[3][2] = (-ShC*ShRy/4.0+ShCy*ShR/4.0)*L4;
			Kc[3][4] = (ShRz*ShC/4.0-ShR*ShCz/4.0)*L4;
				
			Kc[4][0] = (ShRy*rt3*ShC/4.0-ShR*rt3*ShCy/4.0)*L4;
			Kc[4][1] = (-ShC*ShRy/4.0+ShCy*ShR/4.0)*L4;
			Kc[4][2] = (ShC*ShRx/4.0-ShCx*ShR/4.0)*L4;
			Kc[4][3] = (-ShRz*ShC/4.0+ShR*ShCz/4.0)*L4;
		    }

		    lK[i][j   ] +=ShRC*T11	+ dot	+L2_K11	+L6_K11;
		    lK[i][j+4 ] +=ShRC*T12	        +L2_K12	+L6_K12;
		    lK[i][j+8 ] +=ShRC*T13          +L2_K13	+L6_K13;
		    lK[i][j+12] +=ShRC*T14          +L2_K14	+L6_K14   +Kc[0][3];
		    lK[i][j+16] +=ShRC*T15          +L2_K15	+L6_K15   +Kc[0][4];

		    lK[i+4][j+0 ] +=ShRC*T12        + L2_K12 +L6_K12;
		    lK[i+4][j+4 ] +=ShRC*T22 +	dot	+ L2_K22 +L6_K22;
		    lK[i+4][j+8 ] +=ShRC*T23        + L2_K23 +L6_K23 +Kc[1][2];
		    lK[i+4][j+12] +=ShRC*T24        + L2_K24 +L6_K24 +Kc[1][3];
		    lK[i+4][j+16] +=ShRC*T25        + L2_K25 +L6_K25 +Kc[1][4];

		    lK[i+8][j+0] +=ShRC*T13         + L2_K13 +L6_K13;
		    lK[i+8][j+4] +=ShRC*T23         + L2_K23 +L6_K23 +Kc[2][1];
		    lK[i+8][j+8]  +=ShRC*T33 + dot	+ L2_K33 +L6_K33 ;
		    lK[i+8][j+12] +=ShRC*T34        + L2_K34 +L6_K34 +Kc[2][3];
		    lK[i+8][j+16] +=ShRC*T35        + L2_K35 +L6_K35 +Kc[2][4];
			
		    lK[i+12][j+0 ]+=ShRC*T14 +		  L2_K14 +L6_K14 + Kc[3][0];
		    lK[i+12][j+4 ]+=ShRC*T24 +        L2_K24 +L6_K24 + Kc[3][1];
		    lK[i+12][j+8 ]+=ShRC*T34 +        L2_K34 +L6_K34 + Kc[3][2];
		    lK[i+12][j+12]+=ShRC*T44 + dot	+ L2_K44 +L6_K44;
		    lK[i+12][j+16]+=ShRC*T45 +        L2_K45 +L6_K45 + Kc[3][4];
				
		    lK[i+16][j+0 ]+=ShRC*T15        + L2_K15 +L6_K15 +Kc[4][0];
		    lK[i+16][j+4] +=ShRC*T25        + L2_K25 +L6_K25 +Kc[4][1];
		    lK[i+16][j+8] +=ShRC*T35        + L2_K35 +L6_K35 +Kc[4][2];
		    lK[i+16][j+12]+=ShRC*T45	    + L2_K45 +L6_K45 +Kc[4][3];
		    lK[i+16][j+16]+=ShRC*T55 + dot	+ L2_K55 +L6_K55;
				
		    //Local identity matrix
		    lI[i   ][j   ]+=ShRC;
		    lI[i+4 ][j+4 ]+=ShRC;
		    lI[i+8 ][j+8 ]+=ShRC;
		    lI[i+12][j+12]+=ShRC;
		    lI[i+16][j+16]+=ShRC;

		}//end for j
	    } // end if assemble matrix
	} // end for i  - rows
    }//end for igp

    if(simu->dt!=0){// Crank-Nicolson time stepping
	// %0 calculates the steady state
	// the product of the mass matrix and the various q vectors
	double Mq[20];	// M * q

	memset(Mq,0,20*sizeof(double));

	for (int i=0;i<4;i++) {		//each node row
	    for (int j=0;j<4;j++) {	//each node column
            for (int k=0;k<5;k++){//each component
                Mq[4*k+i]+=lI[i][j]*q->Values[npLC*k+tt[j]];
            }
	    }
	}
	double dt = simu->dt;
	double u1 = mat_par->u1;
	for (int i=0;i<20;i++){
	    lL[i] =  ( lL[i] / 2.0 ) + + Mq[i]*(u1 / dt);    // current RHS


		for (int j=0 ; j<20 ; j++){
             lK[i][j] = lK[i][j]/2.0 +lI[i][j]*(u1/dt);
		}

	}
    }//if(dt!=0)
}// end void localKL

void localKL_NQ(
	double* p,
	int* tt,
	double lL[20],
	int it,
	int index_to_Neumann,
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
	double Sh[4];
	Sh[0] = sh1[igp][0];
	Sh[1] = sh1[igp][1];
	Sh[2] = sh1[igp][2];
	Sh[3] = sh1[igp][3];
		
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
	int element_num ,
	SolutionVector* q,
	double lL[15],
	double lK[15][15],
	int FixLCNumber,
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
    SparseMatrix* K,
    double* L,
    SolutionVector* q,
    SolutionVector* v,
    Mesh* t, double* p,
    LC* mat_par,
    Simu* simu,
    Settings* settings)
{

    int npLC = q->getnDoF();
    init_shape();
    omp_set_num_threads( settings->getnThreads() ); // number of threads used
    #pragma omp parallel for
    // LOOP OVER EACH ELEMENT  it
    for (int it= 0 ; it < t->getnElements () ;it++)
    {
	double lK[20][20]; 	// local element matrix
	double lL[20];		// local RHS vector
	int eqr,eqc;
        // IF THIS ELEMENT IS LC ELEMENT, ASSEMBLE LOCAL MATRIX

        if( t->getMaterialNumber(it) == MAT_DOMAIN1 ) // if LC element
        {
            localKL(p,t,it,q,v,lK,lL,mat_par, simu);

            // ADD LOCAL MATRIX TO GLOBAL MATRIX
            for (int i=0;i<20;i++)  // LOOP OVER ROWS
            {
                //int ri = tt[i%4]+npLC*(i/4); //
                int ri = t->getNode(it,i%4) + npLC*(i/4);   // LOCAL TO GLOBAL

                eqr = q->getEquNode(ri);    // eqr IS MAPPED INDEX TO GLOBAL MATRIX ROW

                if ( eqr != SolutionVector::FIXED_NODE ) // ONLY FOR NON-FIXED NODES
                {
                    #pragma omp atomic
                    L[eqr]+= lL[i]*BIGNUM;

                    for (int j = 0 ; j < 20 ; j++) // LOOP OVER COLUMNS
                    {
                        int rj = t->getNode(it,j%4) + npLC*(j/4);
                        eqc = q->getEquNode( rj );

                        if ( eqc != SolutionVector::FIXED_NODE )
                        {
                            K->sparse_add(eqr,eqc,lK[i][j]*BIGNUM);
                        }
                    }
                }// end NON-FIXED NODE


	    }//end for i
        }//end for if tmat
    }//end fr it
}
// end void assemble_volumes

void assemble_Neumann_surfaces(
	double* L,
	SolutionVector* q,
	SolutionVector* v,
	Mesh* mesh,
	Mesh* surf_mesh,
	double* p){

    int npLC = q->getnDoF();
    init_shape_N();
#pragma omp parallel for

    for (int it=0; it< surf_mesh->getnElements(); it++) // LOOP OVER EVERY SURFACE ELEMENT
    {
        // ONLY TRIS CONNECTED TO LC TETS ARE ASSEMBLED
        int index_to_Neumann = surf_mesh->getConnectedVolume(it);
        if ( (index_to_Neumann > -1) &&                             // IF CONNECTED TO LC TET
             (surf_mesh->getMaterialNumber(it) != MAT_PERIODIC))    // IF THIS IS NOT A PERIODIC TRIANGLE
        {
            double lL[20];

            // ELEMENT NODE NUMBERS ARE RE-ORDERED SO THAT t[4] IS NOT PART OF TRI ELEMENT
	    int ee[3] = {   surf_mesh->getNode(it,0) ,
			    surf_mesh->getNode(it,1) ,
			    surf_mesh->getNode(it,2) } ;

	    int tt[4] = {   mesh->getNode(index_to_Neumann,0),
			    mesh->getNode(index_to_Neumann,1),
			    mesh->getNode(index_to_Neumann,2),
			    mesh->getNode(index_to_Neumann,3)};

	    int intr=-1;//find  index to internal node
            for (int i=0;i<4;i++)
            {
                if ( (tt[i]!= ee[0]) && (tt[i]!= ee[1]) && (tt[i]!= ee[2]) )
                {
                    intr = i;
		    break;
                }
	    }
            int ti[4] = { ee[0], ee[1], ee[2], tt[intr] }; // REORDER LOCAL TET ELEMENT
                                                           // NODE-NUMBERING SO THAT
                                                           // INTERNAL NODE IS ALWAYS LAST

	    localKL_NQ(p, tt, lL , it , index_to_Neumann,mesh, surf_mesh, v);

            for (int i=0; i<20; i++) // LOOP OVER ROWS
            {
                int ri = ti[i%4] + npLC*(i/4);
                int eqr = q->getEquNode(ri);

                if (eqr != SolutionVector::FIXED_NODE )
                {
                    #pragma omp atomic
                    L[eqr]+=lL[i]*BIGNUM;
                }
            }//end for i
	}//end if LC
    }//end for it
}
//end void assemble_Neumann


// ASSEMBLE WEAK ANCHORING SURFACES
void assemble_surfaces(
	SparseMatrix* K ,
	double* L ,
	SolutionVector* q ,
	Mesh* e ,
	LC* lc ,
	Alignment* alignment,
	double* NodeNormals){

    init_shape_surf();
    int npLC = q->getnDoF();

    #pragma omp parallel for
    for (int ie = 0 ; ie < e->getnElements() ; ie ++)
    {
        int FixLCNum = e->getFixLCNumber(ie); // gets FixLC number for surface element ie
        if ((FixLCNum > 0) && ( !alignment->IsStrong(FixLCNum-1) ) ) // if alignment surface
        {
            double lK[15][15];
            double lL[15];
            wk_localKL( e , ie , q , lL , lK , FixLCNum , alignment, lc , NodeNormals);

            for (int i=0;i<15;i++)// LOOP OVER ROWS
            {
                int ri 	= e->getNode(ie,i%3) + npLC*(i/3);
                int eqr = q->getEquNode(ri);
                if (eqr != SolutionVector::FIXED_NODE )
                {
                    #pragma omp atomic
                    L[eqr]+=lL[i]*2e16;

                    for (int j=0;j<15;j++) // LOOP OVER COLUMNS
                    {
                        int rj  = e->getNode(ie,j%3) + npLC*(j/3);
                        int eqc = q->getEquNode(rj);

                        if ( eqc != SolutionVector::FIXED_NODE )
                        {
                            int ii = i; // SURFACE CONTRIBUTION MATRIX IS
                            int jj = j; // SYMMETRIC -> ONLY UPPER DIAGONAL
                            if (j<i)    // IS ASSEMBLED
                            {
                                ii = j;
                                jj = i;
                            }

                            K->sparse_add(eqr, eqc, lK[ii][jj]*BIGNUM);
                        } // end if j not fixed
                    }// end for j
                }// end if i not fixed
            }//end for i
        }// end if alignment surfce
    }// end for ie, loop over surface elements
}
//*/
void assembleQ(
            SparseMatrix* K,
            double* L,  // current RHS
            SolutionVector *q,  // current Q-Tensor
            SolutionVector* v,
            Mesh* t,
            Mesh* e,
            double* p,
            LC* mat_par,
            Simu* simu,
            Settings* settings,
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


    assemble_volumes(K, L, q,  v, t, p, mat_par, simu, settings);

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

void assemble_local_prev_volumes(double lL[20],
				 SolutionVector& q,  SolutionVector& v,
				 Mesh& t, double* p,  unsigned int& element_num,
				 LC& mat_par, Simu& simu, int th  ){

    // th = current thread number for debugging
    if (th){} // no warnings of unused variables

    memset(lL,0,20*sizeof(double));
    double lI[20][20];
    memset(lI,0,20*20*sizeof(double));



    int tt[4] = {0,0,0,0};
    tt[0] = t.getNode(element_num,0);
    tt[1] = t.getNode(element_num,1);
    tt[2] = t.getNode(element_num,2);
    tt[3] = t.getNode(element_num,3);

    double Jdet = t.getDeterminant(element_num);

    if (Jdet < 0) Jdet = -Jdet;


    bool three_elastic_constants = false;
    if ((L2!=0)&&(L6!=0)) three_elastic_constants = true;


    double L1_1,L1_2,L1_3,L1_4,L1_5;//L1 elastic RHS  terms;
    double L2_1,L2_2,L2_3,L2_4,L2_5;//L2 elastic RHS  terms;
    double L6_1,L6_2,L6_3,L6_4,L6_5;//L6 elastic RHS  terms;
    double Lc[5] = {0,0,0,0,0};
    double V1,V2,V3,V4,V5;
    double T1,T2,T3,T4,T5;
    double Kc[5][5];
    memset(Kc,0,5*5*sizeof(double));



    //1. Calculate Inverse Jacobian - for 1st order elements can be done outside integration loop -> igp = 0
	    double xr,xs,xt,yr,ys,yt,zr,zs,zt;
	    xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;
	    for (int i=0;i<4;i++) {
		    xr+=sh1r[0][i]*p[(tt[i])*3+0]*1e-6;
		    xs+=sh1s[0][i]*p[(tt[i])*3+0]*1e-6;
		    xt+=sh1t[0][i]*p[(tt[i])*3+0]*1e-6;

		    yr+=sh1r[0][i]*p[(tt[i])*3+1]*1e-6;
		    ys+=sh1s[0][i]*p[(tt[i])*3+1]*1e-6;
		    yt+=sh1t[0][i]*p[(tt[i])*3+1]*1e-6;

		    zr+=sh1r[0][i]*p[(tt[i])*3+2]*1e-6;
		    zs+=sh1s[0][i]*p[(tt[i])*3+2]*1e-6;
		    zt+=sh1t[0][i]*p[(tt[i])*3+2]*1e-6;
	    }//end for i
	    //Inverse Jacobian
	    Jdet = fabs(xr*ys*zt-xr*zs*yt+xs*yt*zr-xs*yr*zt+xt*yr*zs-xt*ys*zr);
	    double Jinv[3][3]={{(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet}
					      ,{(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet}
					      ,{(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet}};


    // shape function derivatives
    double dSh[4][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    for(int i=0;i<4;i++){
		dSh[i][0]=sh1r[0][i]*Jinv[0][0]+sh1s[0][i]*Jinv[1][0]+sh1t[0][i]*Jinv[2][0];
		dSh[i][1]=sh1r[0][i]*Jinv[0][1]+sh1s[0][i]*Jinv[1][1]+sh1t[0][i]*Jinv[2][1];
		dSh[i][2]=sh1r[0][i]*Jinv[0][2]+sh1s[0][i]*Jinv[1][2]+sh1t[0][i]*Jinv[2][2];
    }//end for i

    for (int igp = 0 ; igp < ngp ; igp ++){
    //shape function
    double Sh[4];
    Sh[0] = sh1[igp][0];
    Sh[1] = sh1[igp][1];
    Sh[2] = sh1[igp][2];
    Sh[3] = sh1[igp][3];

    // Function variables and derivatives
    double q1=0, q2=0, q3=0, q4=0, q5=0;
    double q1x=0,q2x=0,q3x=0,q4x=0,q5x=0;
    double q1y=0,q2y=0,q3y=0,q4y=0,q5y=0;
    double q1z=0,q2z=0,q3z=0,q4z=0,q5z=0;
    double Vx=0,Vy=0,Vz=0;

    // Solution and derivatives
    for(int i=0;i<4;i++){
	q1+=Sh[i]*q.getValue(tt[i],0);				// q1i * Ni  = A1
	q2+=Sh[i]*q.getValue(tt[i],1);
	q3+=Sh[i]*q.getValue(tt[i],2);
	q4+=Sh[i]*q.getValue(tt[i],3);
	q5+=Sh[i]*q.getValue(tt[i],4);

	// voltages
	Vx+=dSh[i][0]*v.getValue(tt[i]);
	Vy+=dSh[i][1]*v.getValue(tt[i]);
	Vz+=dSh[i][2]*v.getValue(tt[i]);

	q1x+=dSh[i][0]*q.getValue(tt[i],0);
	q2x+=dSh[i][0]*q.getValue(tt[i],1);
	q3x+=dSh[i][0]*q.getValue(tt[i],2);
	q4x+=dSh[i][0]*q.getValue(tt[i],3);
	q5x+=dSh[i][0]*q.getValue(tt[i],4);

	q1y+=dSh[i][1]*q.getValue(tt[i],0);
	q2y+=dSh[i][1]*q.getValue(tt[i],1);
	q3y+=dSh[i][1]*q.getValue(tt[i],2);
	q4y+=dSh[i][1]*q.getValue(tt[i],3);
	q5y+=dSh[i][1]*q.getValue(tt[i],4);

	q1z+=dSh[i][2]*q.getValue(tt[i],0);
	q2z+=dSh[i][2]*q.getValue(tt[i],1);
	q3z+=dSh[i][2]*q.getValue(tt[i],2);
	q4z+=dSh[i][2]*q.getValue(tt[i],3);
	q5z+=dSh[i][2]*q.getValue(tt[i],4);
    }//end for i

    double R=q1*q1+q2*q2+q3*q3+q5*q5+q4*q4; // frequently reoccurring term
    double mul=w[igp]*Jdet;

    T1=((A*q1+B*(q5*q5*rt6/4.0-q3*q3*rt6/2.0-rt6*q2*q2/2.0+q1*q1*rt6/2.0+q4*q4*rt6/4.0)/3.0)+C*R*q1);
    T2=(A*q2+B*(3.0/4.0*q5*q5*rt2-q1*rt6*q2-3.0/4.0*q4*q4*rt2)/3.0+C*R*q2);
    T3=(A*q3+B*(-q3*q1*rt6+3.0/2.0*rt2*q5*q4)/3.0+C*R*q3);
    T4=(A*q4+B*(3.0/2.0*q3*rt2*q5+q4*q1*rt6/2.0-3.0/2.0*q4*q2*rt2)/3.0+C*R*q4);
    T5=(A*q5+B*(q5*q1*rt6/2.0+3.0/2.0*q5*q2*rt2+3.0/2.0*q3*rt2*q4)/3.0+C*R*q5);

    V1 =  rt6*(Vx*Vx + Vy*Vy-2.0*Vz*Vz)*deleps/18.0*eps0;
    V2 = -rt2*(Vx*Vx - Vy*Vy)*deleps/6.0*eps0;
    V3 = -rt2*Vx*Vy*deleps/3.0*eps0;
    V4 = -rt2*Vz*Vy*deleps/3.0*eps0;
    V5 = -rt2*Vx*Vz*deleps/3.0*eps0;

    double Lflexo1, Lflexo2, Lflexo3, Lflexo4, Lflexo5;
    Lflexo1 = 0; Lflexo2 =0; Lflexo3=0; Lflexo4=0; Lflexo5=0;

    for (int i=0;i<4;i++){ // matrix rows
		double ShRx=mul*dSh[i][0];//including weight and jacobian in trial function
		double ShRy=mul*dSh[i][1];
		double ShRz=mul*dSh[i][2];
		double ShR =mul*Sh[i];

		L1_1=(ShRx*q1x+ShRy*q1y+ShRz*q1z)*L1;
		L1_2=(ShRx*q2x+ShRy*q2y+ShRz*q2z)*L1;
		L1_3=(ShRx*q3x+ShRy*q3y+ShRz*q3z)*L1;
		L1_4=(ShRx*q4x+ShRy*q4y+ShRz*q4z)*L1;
		L1_5=(ShRx*q5x+ShRy*q5y+ShRz*q5z)*L1;

		//Chiral term
		if (L4!=0){
			Lc[0] = ((rt3*q5y/4.0-rt3*q4x/4.0)*ShR+ShRx*q4*rt3/4.0-ShRy*q5*rt3/4.0)*L4;
			Lc[1] = ((q3z/2.0-q5y/4.0-q4x/4.0)*ShR+ShRx*q4/4.0+ShRy*q5/4.0-ShRz*q3/2.0)*L4;
			Lc[2] = ((-q2z/2.0-q4y/4.0+q5x/4.0)*ShR+ShRy*q4/4.0+ShRz*q2/2.0-ShRx*q5/4.0)*L4;
			Lc[3] = ((-q5z/4.0+q1x*rt3/4.0+q3y/4.0+q2x/4.0)*ShR-ShRx*q1*rt3/4.0-ShRx*q2/4.0-ShRy*q3/4.0+ShRz*q5/4.0)*L4;
			Lc[4] = ((q4z/4.0-q1y*rt3/4.0+q2y/4.0-q3x/4.0)*ShR+ShRx*q3/4.0+ShRy*q1*rt3/4.0-ShRy*q2/4.0-ShRz*q4/4.0)*L4;
			}
		if (three_elastic_constants){
			L2_1= (ShRx*q1x/6.0-ShRx*rt3*q2x/6.0-ShRx*rt3*q3y/6.0-ShRx*rt3*q5z/6.0-ShRy*q3x*rt3/6.0+ShRy*q1y/6.0+ShRy*rt3*q2y/6.0-	ShRy*rt3*q4z/6.0+ShRz*q5x*rt3/3.0	+ShRz*q4y*rt3/3.0+2.0/3.0*ShRz*q1z)	*L2;
			L2_2 = (-ShRx*q1x*rt3/6.0+q2x*ShRx/2.0+q3y*ShRx/2.0+q5z*ShRx/2.0-q3x*ShRy/2.0+ShRy*q1y*rt3/6.0+q2y*ShRy/2.0-q4z*ShRy/2.0)*L2;
			L2_3 = (ShRx*q3x/2.0-ShRx*q1y*rt3/6.0-ShRx*q2y/2.0+ShRx*q4z/2.0-ShRy*q1x*rt3/6.0+ShRy*q2x/2.0+ShRy*q3y/2.0+ShRy*q5z/2.0)*L2;
			L2_4= (ShRy*q5x/2.0+ShRy*q4y/2.0+ShRy*q1z*rt3/3.0+ShRz*q3x/2.0-	ShRz*q1y*rt3/6.0-ShRz*q2y/2.0+ShRz*q4z/2.0)*L2;
			L2_5 = (ShRx*q5x/2.0+ShRx*q4y/2.0+ShRx*q1z*rt3/3.0-ShRz*q1x*rt3/6.0	+ShRz*q2x/2.0+ShRz*q3y/2.0+ShRz*q5z/2.0)*L2;
			L6_1 = (-ShRx*q1x*q1*rt2*rt3/6.0-ShRy*q1y*q1*rt2*rt3/6.0+ShRz*q1*rt2*rt3*q1z/3.0-ShR*rt2*rt3*q5y*q5y/12.0+ShRx*q5*rt2*q1z/2.0+ShRy*q3*rt2*q1x/2.0+ShRy*q4*rt2*q1z/2.0+ShRz*q5*rt2*q1x/2.0-ShR*rt2*rt3*q4x*q4x/12.0-ShR*rt2*rt3*q4y*q4y/12.0-ShR*rt2*rt3*q3x*q3x/12.0-ShR*rt2*rt3*q2y*q2y/12.0+ShR*rt2*rt3*q2z*q2z/6.0+ShR*rt2*rt3*q3z*q3z/6.0+ShR*rt2*rt3*q5z*q5z/6.0+ShR*rt2*rt3*q4z*q4z/6.0+ShR*rt2*rt3*q1z*q1z/6.0+ShRz*q4*rt2*q1y/2.0-ShR*rt2*rt3*q1x*q1x/12.0-ShR*rt2*rt3*q3y*q3y/12.0-ShR*rt2*rt3*q2x*q2x/12.0-ShR*rt2*rt3*q1y*q1y/12.0+ShRx*q1x*q2*rt2/2.0+ShRx*q3*rt2*q1y/2.0-ShRy*q1y*q2*rt2/2.0-ShR*rt2*rt3*q5x*q5x/12.0)*L6;
			L6_2 = (-ShR*rt2*q5y*q5y/4.0-ShR*rt2*q1y*q1y/4.0+ShR*rt2*q4x*q4x/4.0+ShR*rt2*q1x*q1x/4.0-ShR*rt2*q3y*q3y/4.0+ShR*rt2*q2x*q2x/4.0+ShR*rt2*q5x*q5x/4.0+ShR*rt2*q3x*q3x/4.0-ShR*rt2*q4y*q4y/4.0-ShR*rt2*q2y*q2y/4.0-ShRx*q1*rt2*rt3*q2x/6.0+ShRx*rt2*q2*q2x/2.0+ShRx*q3*q2y*rt2/2.0+ShRx*q5*q2z*rt2/2.0-ShRy*q1*rt2*rt3*q2y/6.0-ShRy*rt2*q2*q2y/2.0+ShRy*q3*q2x*rt2/2.0+ShRy*q4*q2z*rt2/2.0+ShRz*q5*q2x*rt2/2.0+ShRz*q4*q2y*rt2/2.0+ShRz*q1*rt2*rt3*q2z/3.0)*L6;
			L6_3 = (ShR*rt2*q1x*q1y/2.0+ShR*rt2*q2x*q2y/2.0+ShR*rt2*q3x*q3y/2.0+ShR*rt2*q5x*q5y/2.0+ShR*rt2*q4x*q4y/2.0-ShRx*q3x*q1*rt2*rt3/6.0+ShRx*q3x*q2*rt2/2.0+ShRx*q3*rt2*q3y/2.0+ShRx*q5*rt2*q3z/2.0-ShRy*q3y*q1*rt2*rt3/6.0-ShRy*q3y*q2*rt2/2.0+ShRy*q3*rt2*q3x/2.0+ShRy*q4*rt2*q3z/2.0+ShRz*q5*rt2*q3x/2.0+ShRz*q4*rt2*q3y/2.0+ShRz*q1*rt2*rt3*q3z/3.0)*L6;
			L6_4 = (ShR*rt2*q1y*q1z/2.0+ShR*rt2*q2y*q2z/2.0+ShR*rt2*q3y*q3z/2.0+ShR*rt2*q5y*q5z/2.0+ShR*rt2*q4y*q4z/2.0-ShRx*q4x*q1*rt2*rt3/6.0+ShRx*q4x*q2*rt2/2.0+ShRx*q3*rt2*q4y/2.0+ShRx*q5*rt2*q4z/2.0-ShRy*q4y*q1*rt2*rt3/6.0-ShRy*q4y*q2*rt2/2.0+ShRy*q3*rt2*q4x/2.0+ShRy*q4*rt2*q4z/2.0+ShRz*q5*rt2*q4x/2.0+ShRz*q4*rt2*q4y/2.0+ShRz*q1*rt2*rt3*q4z/3.0)*L6;
			L6_5 = (ShR*rt2*q1x*q1z/2.0+ShR*rt2*q2x*q2z/2.0+ShR*rt2*q3x*q3z/2.0+ShR*rt2*q5x*q5z/2.0+ShR*rt2*q4x*q4z/2.0-ShRx*q5x*q1*rt2*rt3/6.0+ShRx*q5x*q2*rt2/2.0+ShRx*q3*rt2*q5y/2.0+ShRx*q5*rt2*q5z/2.0-ShRy*q5y*q1*rt2*rt3/6.0-ShRy*q5y*q2*rt2/2.0+ShRy*q3*rt2*q5x/2.0+ShRy*q4*rt2*q5z/2.0+ShRz*q5*rt2*q5x/2.0+ShRz*q4*rt2*q5y/2.0+ShRz*q1*rt2*rt3*q5z/3.0)*L6;
		}else{L2_1=L2_2=L2_3=L2_4=L2_5=L6_1=L6_2=L6_3=L6_4=L6_5=0.0;}//end if 3 elestic constants
		if ((efe!=0.0)||(efe2!=0.0)){ // IF FLEXOELECTRIC COEFFICIENTS ARN'T 0
	    Lflexo1 = (rt6*(Vx*ShRx+Vy*ShRy-2.0*Vz*ShRz)*efe/6.0);
	    Lflexo2 = (-rt2*(Vx*ShRx-Vy*ShRy)*efe/2.0);
	    Lflexo3 = (-rt2*(Vy*ShRx+Vx*ShRy)*efe/2.0);
	    Lflexo4 = (-rt2*(Vz*ShRy+Vy*ShRz)*efe/2.0);
	    Lflexo5 = (-rt2*(Vz*ShRx+Vx*ShRz)*efe/2.0);
	}

	//ADD TERMS TO RHS VECTOR
	lL[i+0]  += T1*ShR + L1_1 + L2_1 + L6_1 + Lc[0] + V1*ShR + Lflexo1;
	lL[i+4]  += T2*ShR + L1_2 + L2_2 + L6_2 + Lc[1] + V2*ShR + Lflexo2;
	lL[i+8]  += T3*ShR + L1_3 + L2_3 + L6_3 + Lc[2] + V3*ShR + Lflexo3;
	lL[i+12] += T4*ShR + L1_4 + L2_4 + L6_4 + Lc[3] + V4*ShR + Lflexo4;
	lL[i+16] += T5*ShR + L1_5 + L2_5 + L6_5 + Lc[4] + V5*ShR + Lflexo5;

	for (int j = 0 ; j < 4 ; j++){
	    double ShRC = ShR*Sh[j];
	    //Local identity matrix
	     lI[i   ][j   ]+=ShRC;
	     lI[i+4 ][j+4 ]+=ShRC;
	     lI[i+8 ][j+8 ]+=ShRC;
	     lI[i+12][j+12]+=ShRC;
	     lI[i+16][j+16]+=ShRC;
	}//

    } // end for i  - rows
}//end for igp

    // MAKE CRANK-NICHOLSON RHS TERM
    // the product of the mass matrix and the various q vectors
    double Mq[20]={0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0};	// M * q
    //memset(Mq,0,20*sizeof(double));
    for (int i=0;i<4;i++) {		//each node row
		for (int j=0;j<4;j++) {	//each node column
			for (int k=0;k<5;k++){//each component
                                Mq[4*k+i]+=lI[i][j]*q.Values[npLC*k+tt[j]];
				}
			}
    }

    double dt = simu.dt;
    double u1 = mat_par.u1;


    for (size_t i=0;i<20;i++){
    /* SEGFAULT WITH OPENMP HERE*/
        lL[i] =   ( lL[i] / 2.0 ) -  ( Mq[i]*(u1 / dt) ) ;   // M*current Q

    }


}// end local rhs prev


/*=====================================================*/
/*  ASSEMBLES RHS CONTRIBUTION FROM PREVIOUS TIME-STEP */
/*=====================================================*/

void assemble_prev_rhs(double* Ln,
		       SolutionVector& qn,
		       SolutionVector& v,
                       //Mesh& t,
                       //Mesh& e,
                       //double* p ,
		       LC& mat_par,
                       Simu& simu,
                       Geometry& geom
                       )
{
    //if (e.getnElements()){} // no warning of unused variables. WRITE fUNCTION FOR WEAK ANCHORING CONTRIBUTIONS

    init_globals(mat_par, qn);
    init_shape();
    unsigned int elem_cnt = geom.t->getnElements();//unsigned int) t.getnElements();

    // OPENMP LOOP COMPILED WITH -march=native an -O3 RESULTS IN SEGFAULT ON
    // WINXP32, COMPILED WITHMinGW. THIS IS NOT A PROBLEM WITH UBUNTU,
    // NO PROBLEMS FOUND WITH gdb / valgrind. SGFAULTING LINE MARKED IN FUNTION
    // assemble_local_prev_volumes. ENABLING OPENMP ONLY FOR LINUX (02/12/2010)
#ifndef __WIN32__
#pragma omp parallel for // PARALLEL LOOP IN LINUX
#endif
    int th = 0; // debug thread number
    Mesh& t = *geom.t;
    double* p = geom.getPtrTop();
    for ( unsigned int it = 0 ; it < elem_cnt ; it++)
    {

        int eqr;
        double lL[20];		// local RHS vector

        // IF THIS ELEMENT IS LC ELEMENT, ASSEMBLE LOCAL MATRIX
        if( t.getMaterialNumber(it) == MAT_DOMAIN1 )// if LC element
        {
            assemble_local_prev_volumes(lL,
                                        qn, v ,
                                        t , p , it,
                                        mat_par , simu, th );

	    // ADD LOCAL MATRIX TO GLOBAL MATRIX

            for (int i=0;i<20;i++)
            {
                int ri = t.getNode(it,i%4) + npLC*(i/4);
                eqr = qn.getEquNode(ri);
                if (eqr != SolutionVector::FIXED_NODE )
                {
                    #pragma omp atomic
                    Ln[eqr] += lL[i]*BIGNUM;
                }
            }// end for i
        }//end if LC material

    }//end for it

}// end function


