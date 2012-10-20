#include <assembleq2k.h>
#include <globals.h>
#include <shapefunction3d.h>
#include <material_numbers.h>
#include <qassembly_macros.h>
#include <ircmatrix.h>
#define	BIGNUM 2e16

// LOCAL ELEMENT MATRIX ASSEMBLY
inline void localKL_2K( double* p,                      // COORINDATES
                        const Mesh& t,                  // TETS MESH
                        const idx it,                   // CURRENT ELEMENT IDX
                        const SolutionVector& q,        // Q-TENSOR
                        const SolutionVector& v,        // POTENTIAL
                        const LC& mat_par,              // MATERIAL PARAMS
                        const Simu& simu,               //
                        const Shape4thOrder& shapes,    // SHAPE FUNCTIONS
                        double lK[20][20],              // LOCAL MATRIX
                        double* lL )                 // LOCAL R.H.S.
{
    memset(lK, 0 , 20*20*sizeof(double) );
    memset(lL, 0 , 20*sizeof(double) );

    double lI[20][20]; // LOCAL IDENTITY MATRIX
    memset(lI, 0 , 20*20*sizeof(double) );

    // LOCAL CACHED COPY OF CURRENT ELEMENT AND ITS COORDINATES
    idx tt[4] = {t.getNode(it,0), t.getNode(it,1) ,
                 t.getNode(it,2), t.getNode(it,3) };

    double pp[4][3] = { {p[tt[0]*3] , p[tt[0]*3+1] , p[tt[0]*3+2]} ,
                        {p[tt[1]*3] , p[tt[1]*3+1] , p[tt[1]*3+2]} ,
                        {p[tt[2]*3] , p[tt[2]*3+1] , p[tt[2]*3+2]} ,
                        {p[tt[3]*3] , p[tt[3]*3+1] , p[tt[3]*3+2]} };

    // CALCULATE INVERSE JACOBIAN MATRIX FOR CURRENT ELEMENT
    double xr(0), xs(0), xt(0), yr(0), ys(0), yt(0), zr(0), zs(0), zt(0);
    for (int i=0;i<4;++i)
    {
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
    double Jdet = t.getDeterminant(it)*1e18;
    double invJdet = 1e6 / Jdet;
    Jdet*=1e-18;
    // ACTUAL INVERSE JACOBIAN
    double Jinv[3][3]={{ (zt*ys-yt*zs)*invJdet ,(xt*zs-zt*xs)*invJdet, (xs*yt-ys*xt)*invJdet}
                       ,{(yt*zr-zt*yr)*invJdet ,(zt*xr-xt*zr)*invJdet, (xt*yr-yt*xr)*invJdet}
                       ,{(yr*zs-ys*zr)*invJdet ,(xs*zr-xr*zs)*invJdet, (ys*xr-xs*yr)*invJdet}};

    // CREATE LOCAL SHAPE FUNCTION DERIVATIVES
    double dSh[4][3];
    for (int i = 0 ; i < 4 ; i++)
    {
        dSh[i][0]=shapes.sh1r[0][i]*Jinv[0][0]      // X-DERIVATIVES
                + shapes.sh1s[0][i]*Jinv[1][0]
                + shapes.sh1t[0][i]*Jinv[2][0];

        dSh[i][1]=shapes.sh1r[0][i]*Jinv[0][1]      // Y
                + shapes.sh1s[0][i]*Jinv[1][1]
                + shapes.sh1t[0][i]*Jinv[2][1];

        dSh[i][2]=shapes.sh1r[0][i]*Jinv[0][2]      // Z
                + shapes.sh1s[0][i]*Jinv[1][2]
                + shapes.sh1t[0][i]*Jinv[2][2];
    }

    // CREATE CACHED COPIES OF VARIABLES Q AND V
    const double vars[6*4] = {q.getValue(tt[0],0)  , q.getValue(tt[0],1), q.getValue(tt[0],2), q.getValue(tt[0],3), q.getValue(tt[0],4),
                              q.getValue(tt[1],0), q.getValue(tt[1],1), q.getValue(tt[1],2), q.getValue(tt[1],3), q.getValue(tt[1],4),
                              q.getValue(tt[2],0), q.getValue(tt[2],1), q.getValue(tt[2],2), q.getValue(tt[2],3), q.getValue(tt[2],4),
                              q.getValue(tt[3],0), q.getValue(tt[3],1), q.getValue(tt[3],2), q.getValue(tt[3],3), q.getValue(tt[3],4),
                              v.getValue(tt[0])  , v.getValue(tt[1])  , v.getValue(tt[2])  , v.getValue(tt[3]   , 0.0                )
                             };


    const double A = mat_par.A;
    const double B = mat_par.B;
    const double C = mat_par.C;
    const double rt6 = sqrt(6.0);
    const double rt2 = sqrt(2.0);
    const double deleps = mat_par.eps_par - mat_par.eps_per;

    double q0(0);
    if (mat_par.p0 != 0.0 )
    {q0 = -2*PI/mat_par.p0;} // '-' SIGN RESULTS IN LEFT HADED CHEIRALITY FOR p0 > 0

    const double K1 = mat_par.K11;
    //const double K2 = mat_par.K22;
    //const double L1 = mat_par.L1;
    // FOR EACH GAUSS POINT
    for (unsigned int igp = 0 ; igp < shapes.ngp; ++igp)
    {
        // LOCAL SHAPE FUNCTION
        const double Sh[4] = {shapes.sh1[igp][0], shapes.sh1[igp][1],
                              shapes.sh1[igp][2], shapes.sh1[igp][3]};

        // FUNCTION VARIABLES FOR THIS GAUSS POINT
        double q1(0), q2(0), q3(0), q4(0), q5(0);
        double q1x(0), q2x(0), q3x(0), q4x(0), q5x(0);
        double q1y(0), q2y(0), q3y(0), q4y(0), q5y(0);
        double q1z(0), q2z(0), q3z(0), q4z(0), q5z(0);
        double Vx(0),  Vy(0) , Vz(0);
        for(int i=0;i<4;++i)
        {
#define IND(i,j) 5*(i) + (j)
            q1+=Sh[i]*vars[IND(i,0)];// OPTIMIZE BY PREFETCHING Q AND V TO LOCAL BUFFER AT START OF FUNCTION
            q2+=Sh[i]*vars[IND(i,1)];
            q3+=Sh[i]*vars[IND(i,2)];
            q4+=Sh[i]*vars[IND(i,3)];
            q5+=Sh[i]*vars[IND(i,4)];

            q1x+=dSh[i][0]*vars[IND(i,0)];
            q2x+=dSh[i][0]*vars[IND(i,1)];
            q3x+=dSh[i][0]*vars[IND(i,2)];
            q4x+=dSh[i][0]*vars[IND(i,3)];
            q5x+=dSh[i][0]*vars[IND(i,4)];

            q1y+=dSh[i][1]*vars[IND(i,0)];
            q2y+=dSh[i][1]*vars[IND(i,1)];
            q3y+=dSh[i][1]*vars[IND(i,2)];
            q4y+=dSh[i][1]*vars[IND(i,3)];
            q5y+=dSh[i][1]*vars[IND(i,4)];

            q1z+=dSh[i][2]*vars[IND(i,0)];
            q2z+=dSh[i][2]*vars[IND(i,1)];
            q3z+=dSh[i][2]*vars[IND(i,2)];
            q4z+=dSh[i][2]*vars[IND(i,3)];
            q5z+=dSh[i][2]*vars[IND(i,4)];
            // voltages
            Vx+=dSh[i][0]*vars[20+i];
            Vy+=dSh[i][1]*vars[20+i];
            Vz+=dSh[i][2]*vars[20+i];
        }//end for i
        // printf("haa %u\n",igp);
        const double mul=shapes.w[igp]*Jdet;
        const double R=q1*q1+q2*q2+q3*q3+q5*q5+q4*q4; // frequently reoccurring term
        //const double D2 = 0.5;
        //const double D3 = 0.33333333333333333333333333333;
        //const double D6 = 0.16666666666666666666666666667;

        double L[5] = {0,0,0,0,0};
        RHS_THERMOTROPIC(L);
        RHS_DIELECTRIC(L);
        ADD_RHS_BULK_TERMS(lL,L);
        MATRIX_THERMOTROPIC(lK);

        /*
        // BULK RHS TERMS
        {

            double T;
            // Q1
            T = RHS_THERMO1;
            //ADD_RHS_BULK_TERMS(0,T);

            // Q2
            T = RHS_THERMO2;
            //ADD_RHS_BULK_TERMS(4,T);

            // Q3
            T = RHS_THERMO3;
            //ADD_RHS_BULK_TERMS(8,T);

            // Q4
            T = RHS_THERMO4;
            //ADD_RHS_BULK_TERMS(12,T);

            // Q5
            T = RHS_THERMO5;
            //ADD_RHS_BULK_TERMS(16,T);

        } // END BULK RHS TERMS
        // BULK THERMOTROPIC MATRIX TERMS
*/
        /*
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
//*/
        // FOR ROWS i
        for (unsigned int i = 0 ; i < 4 ; ++i)
        {
            const double ShRx=mul*dSh[i][0];//including weight and jacobian in trial function
            const double ShRy=mul*dSh[i][1];
            const double ShRz=mul*dSh[i][2];
            const double ShR =mul*Sh[i];

            RHS_ELASTIC_2K_FORMULATION(lL);   // CALCULATES AND ADDS TO RHS
            //RHS_ELASTIC_SINGLE_K(lL);
            // FOR COLUMNS j
            for (unsigned int j = 0 ; j < 4 ; ++j)
            {
                const double ShCx=dSh[j][0];
                const double ShCy=dSh[j][1];
                const double ShCz=dSh[j][2];
                const double ShC =Sh[j];
                const double ShRC=ShR*Sh[j];

                MATRIX_ELASTIC_2K_FORMULATION(lK);
                //MATRIX_ELASTIC_SINGLE_K(lK);

                // LOCAL IDENTITY MATRIX, NEEDED FOR CRANK-NICHOLSON
                lI[i   ][j   ]+=ShRC;
                lI[i+4 ][j+4 ]+=ShRC;
                lI[i+8 ][j+8 ]+=ShRC;
                lI[i+12][j+12]+=ShRC;
                lI[i+16][j+16]+=ShRC;


            }// END FOR COLUMNS j

        }// END FOR ROWS i
    }// end for igp


    // IF CRANK-NICHOLSON
    //*
    if(simu.dt>0)// Crank-Nicolson time stepping
    {
        // the product of the mass matrix and the various q vectors
        double Mq[20];	// M * q
        memset(Mq,0,20*sizeof(double));

        for (int i=0;i<4;i++) 		//each node row
        {
            for (int j=0;j<4;j++) 	//each node column
            {
                Mq[4*0+i]+=lI[i][j]*vars[IND(j,0)];
                Mq[4*1+i]+=lI[i][j]*vars[IND(j,1)];
                Mq[4*2+i]+=lI[i][j]*vars[IND(j,2)];
                Mq[4*3+i]+=lI[i][j]*vars[IND(j,3)];
                Mq[4*4+i]+=lI[i][j]*vars[IND(j,4)];
            }
        }

        const double temp = mat_par.u1 / simu.dt;
        for (int i=0;i<20;i++)
        {
            lL[i] =  0.5*lL[i] + Mq[i]*temp;    // current RHS
            for (int j=0 ; j<20 ; j++)
            {
                lK[i][j] = 0.5*lK[i][j] + lI[i][j]*temp;
            }
        }
    }//if(dt!=0)

    //*/
    // printf("elem %u done\n", it);
}

void assemble_volumes2K(IRCMatrix &K,
                        double* L,
                        const SolutionVector& q,
                        const SolutionVector& v,
                        const LC& mat_par,
                        const Simu& simu,
                        const Mesh& t,
                        double* p
                        )
{
    idx npLC = q.getnDoF();
    const Shape4thOrder shapes;

#pragma omp parallel for
    for (idx it = 0 ; it < t.getnElements() ; it++)
    {
        if (t.getMaterialNumber(it) != MAT_DOMAIN1 )    // ONLY LC ELEMENTS ARE ASSEMBLED HERE
            continue;

        double lK[20][20];  // LOCAL MATRIX
        double lL[20];      // LOCAL R.H.S. VECTOR


        localKL_2K(p, t, it, q, v, mat_par, simu, shapes, lK, lL );
        //printf("ok\n");

        // ADD LOCAL MATRIX AND R.H.S. VECTOR CONTRIBUTIONS TO
        // GLOBAL MATRIX AND R.H.S.
        // TAKE INTO ACCOUNT PERIODIC NODES AND FIXED NODES

        // FOR ROWS i
        for ( unsigned int i = 0 ; i < 20 ; i++)
        {
            // LOCAL TO GLOBAL ROW INDEX NUMBER CONVERSION
            idx ri = t.getNode(it, i%4 ) + npLC*(i/4);  // GLOBAL
            ri = q.getEquNode(ri); // GLOBAL, TAKING INTO ACCOUNT PERIODICITY

            if (ri == NOT_AN_INDEX ) // IF FIXED NODE, SKIP GLOBAL ASSEMBLY
                continue;
#pragma omp atomic
            L[ri] += lL[i]*BIGNUM;

            // FOR COLUMNS j
            for ( unsigned int j = 0 ; j < 20 ; j++)
            {
                idx rj = t.getNode(it, j%4) + npLC*(j/4);
                rj = q.getEquNode(rj);

                // ADD GLOBAL MATRIX, IF NOT FIXED NODE j
                if (rj != NOT_AN_INDEX)
                    K.sparse_add(ri,rj,lK[i][j]*BIGNUM);
            }
        }

    }// end for elements
}




void assemble_volumes2K_previous(double lL[20],
                                 const SolutionVector& q,
                                 const SolutionVector& v,
                                 const Mesh& t,
                                 double* p,
                                 const unsigned int it,
                                 const LC& mat_par,
                                 const Simu& simu,
                                 const Shape4thOrder& shapes
                                 )
{

    memset(lL, 0 , 20*sizeof(double) );

    double lI[20][20]; // LOCAL IDENTITY MATRIX
    memset(lI, 0 , 20*20*sizeof(double) );

    // LOCAL CACHED COPY OF CURRENT ELEMENT AND ITS COORDINATES
    idx tt[4] = {t.getNode(it,0), t.getNode(it,1) ,
                 t.getNode(it,2), t.getNode(it,3) };

    double pp[4][3] = { {p[tt[0]*3] , p[tt[0]*3+1] , p[tt[0]*3+2]} ,
                        {p[tt[1]*3] , p[tt[1]*3+1] , p[tt[1]*3+2]} ,
                        {p[tt[2]*3] , p[tt[2]*3+1] , p[tt[2]*3+2]} ,
                        {p[tt[3]*3] , p[tt[3]*3+1] , p[tt[3]*3+2]} };

    // CALCULATE INVERSE JACOBIAN MATRIX FOR CURRENT ELEMENT
    double xr(0), xs(0), xt(0), yr(0), ys(0), yt(0), zr(0), zs(0), zt(0);
    for (int i=0;i<4;++i)
    {
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
    double Jdet = t.getDeterminant(it)*1e18;
    double invJdet = 1e6 / Jdet;
    Jdet*=1e-18;
    // ACTUAL INVERSE JACOBIAN
    double Jinv[3][3]={{ (zt*ys-yt*zs)*invJdet ,(xt*zs-zt*xs)*invJdet, (xs*yt-ys*xt)*invJdet}
                       ,{(yt*zr-zt*yr)*invJdet ,(zt*xr-xt*zr)*invJdet, (xt*yr-yt*xr)*invJdet}
                       ,{(yr*zs-ys*zr)*invJdet ,(xs*zr-xr*zs)*invJdet, (ys*xr-xs*yr)*invJdet}};

    // CREATE LOCAL SHAPE FUNCTION DERIVATIVES
    double dSh[4][3];
    for (int i = 0 ; i < 4 ; i++)
    {
        dSh[i][0]=shapes.sh1r[0][i]*Jinv[0][0]      // X-DERIVATIVES
                + shapes.sh1s[0][i]*Jinv[1][0]
                + shapes.sh1t[0][i]*Jinv[2][0];

        dSh[i][1]=shapes.sh1r[0][i]*Jinv[0][1]      // Y
                + shapes.sh1s[0][i]*Jinv[1][1]
                + shapes.sh1t[0][i]*Jinv[2][1];

        dSh[i][2]=shapes.sh1r[0][i]*Jinv[0][2]      // Z
                + shapes.sh1s[0][i]*Jinv[1][2]
                + shapes.sh1t[0][i]*Jinv[2][2];
    }

    // CREATE CACHED COPIES OF VARIABLES Q AND V
    const double vars[6*4] = {q.getValue(tt[0],0)  , q.getValue(tt[0],1), q.getValue(tt[0],2), q.getValue(tt[0],3), q.getValue(tt[0],4),
                              q.getValue(tt[1],0), q.getValue(tt[1],1), q.getValue(tt[1],2), q.getValue(tt[1],3), q.getValue(tt[1],4),
                              q.getValue(tt[2],0), q.getValue(tt[2],1), q.getValue(tt[2],2), q.getValue(tt[2],3), q.getValue(tt[2],4),
                              q.getValue(tt[3],0), q.getValue(tt[3],1), q.getValue(tt[3],2), q.getValue(tt[3],3), q.getValue(tt[3],4),
                              v.getValue(tt[0])  , v.getValue(tt[1])  , v.getValue(tt[2])  , v.getValue(tt[3]   , 0.0                )
                             };


    const double A = mat_par.A;
    const double B = mat_par.B;
    const double C = mat_par.C;
    const double rt6 = sqrt(6.0);
    const double rt2 = sqrt(2.0);
    const double deleps = mat_par.eps_par - mat_par.eps_per;

    double q0(0);
    if (mat_par.p0 != 0.0 )
    {q0 = -2*PI/mat_par.p0;} // '-' SIGN RESULTS IN LEFT HANDED CHIRALITY FOR p0 > 0

    const double K1 = mat_par.K11;
    //const double K2 = mat_par.K22;
    //const double L1 = mat_par.L1;

    // FOR EACH GAUSS POINT
    for (unsigned int igp = 0 ; igp < shapes.ngp; ++igp)
    {
        // LOCAL SHAPE FUNCTION
        const double Sh[4] = {shapes.sh1[igp][0], shapes.sh1[igp][1],
                              shapes.sh1[igp][2], shapes.sh1[igp][3]};

        // FUNCTION VARIABLES FOR THIS GAUSS POINT
        double q1(0), q2(0), q3(0), q4(0), q5(0);
        double q1x(0), q2x(0), q3x(0), q4x(0), q5x(0);
        double q1y(0), q2y(0), q3y(0), q4y(0), q5y(0);
        double q1z(0), q2z(0), q3z(0), q4z(0), q5z(0);
        double Vx(0),  Vy(0) , Vz(0);
        for(int i=0;i<4;++i)
        {
#define IND(i,j) 5*(i) + (j)
            q1+=Sh[i]*vars[IND(i,0)];// OPTIMIZE BY PREFETCHING Q AND V TO LOCAL BUFFER AT START OF FUNCTION
            q2+=Sh[i]*vars[IND(i,1)];
            q3+=Sh[i]*vars[IND(i,2)];
            q4+=Sh[i]*vars[IND(i,3)];
            q5+=Sh[i]*vars[IND(i,4)];

            q1x+=dSh[i][0]*vars[IND(i,0)];
            q2x+=dSh[i][0]*vars[IND(i,1)];
            q3x+=dSh[i][0]*vars[IND(i,2)];
            q4x+=dSh[i][0]*vars[IND(i,3)];
            q5x+=dSh[i][0]*vars[IND(i,4)];

            q1y+=dSh[i][1]*vars[IND(i,0)];
            q2y+=dSh[i][1]*vars[IND(i,1)];
            q3y+=dSh[i][1]*vars[IND(i,2)];
            q4y+=dSh[i][1]*vars[IND(i,3)];
            q5y+=dSh[i][1]*vars[IND(i,4)];

            q1z+=dSh[i][2]*vars[IND(i,0)];
            q2z+=dSh[i][2]*vars[IND(i,1)];
            q3z+=dSh[i][2]*vars[IND(i,2)];
            q4z+=dSh[i][2]*vars[IND(i,3)];
            q5z+=dSh[i][2]*vars[IND(i,4)];
            // voltages
            Vx+=dSh[i][0]*vars[20+i];
            Vy+=dSh[i][1]*vars[20+i];
            Vz+=dSh[i][2]*vars[20+i];
        }//end for i
        // printf("haa %u\n",igp);
        const double mul=shapes.w[igp]*Jdet;
        const double R=q1*q1+q2*q2+q3*q3+q5*q5+q4*q4; // frequently reoccurring term

        // BULK RHS TERMS
        double L[5] = {0,0,0,0,0};
        RHS_THERMOTROPIC(L);
        RHS_DIELECTRIC(L);
        ADD_RHS_BULK_TERMS(lL,L);
        // FOR ROWS i
        for (unsigned int i = 0 ; i < 4 ; ++i)
        {
            const double ShRx=mul*dSh[i][0];//including weight and jacobian in trial function
            const double ShRy=mul*dSh[i][1];
            const double ShRz=mul*dSh[i][2];
            const double ShR =mul*Sh[i];

            RHS_ELASTIC_2K_FORMULATION(lL);   // CALCULATES AND ADDS TO RHS
            //RHS_ELASTIC_SINGLE_K(lL);



            // FOR COLUMNS j
            for (unsigned int j = 0 ; j < 4 ; ++j)
            {
                const double ShRC=ShR*Sh[j];
                // LOCAL IDENTITY MATRIX, NEEDED FOR CRANK-NICHOLSON
                lI[i   ][j   ]+=ShRC;
                lI[i+4 ][j+4 ]+=ShRC;
                lI[i+8 ][j+8 ]+=ShRC;
                lI[i+12][j+12]+=ShRC;
                lI[i+16][j+16]+=ShRC;
            }// END FOR COLUMNS j

        }// END FOR ROWS i
    }// end for igp


    // IF CRANK-NICHOLSON
    //*
    if(simu.dt!=0)// Crank-Nicolson time stepping
    {

        // the product of the mass matrix and the various q vectors
        double Mq[20];	// M * q
        memset(Mq,0,20*sizeof(double));

        for (int i=0;i<4;i++) 		//each node row
        {
            for (int j=0;j<4;j++) 	//each node column
            {
                Mq[4*0+i]+=lI[i][j]*vars[IND(j,0)];
                Mq[4*1+i]+=lI[i][j]*vars[IND(j,1)];
                Mq[4*2+i]+=lI[i][j]*vars[IND(j,2)];
                Mq[4*3+i]+=lI[i][j]*vars[IND(j,3)];
                Mq[4*4+i]+=lI[i][j]*vars[IND(j,4)];
            }
        }

        const double temp = mat_par.u1 / simu.dt;
        for (int i=0;i<20;i++)
        {
            lL[i] =  0.5*lL[i] - Mq[i]*temp;    // current RHS
        }
    }//if(dt!=0)
}



void assemble_prev_rhs_K2(double* Ln,
                          SolutionVector& qn,
                          SolutionVector& v,
                          LC& mat_par,
                          Simu& simu,
                          Geometry& geom)
{
    Shape4thOrder shapes;
    unsigned int elem_cnt = geom.t->getnElements();
    Mesh& t = *geom.t;
    double* p = geom.getPtrTop();
    unsigned int npLC = qn.getnDoF();
#ifndef DEBUG
#pragma opm parallel for
#endif
    for (idx it = 0 ; it < elem_cnt ; it++)
    {
        // ONLY ASSEMBLE LC ELEMENTS
        if (t.getMaterialNumber(it) != MAT_DOMAIN1 )
            continue;

        double lL[20];
        assemble_volumes2K_previous(lL, qn, v, t, p, it, mat_par, simu, shapes );

        // ADD TO GLOBAL RHS VECTOR
        for (unsigned int i = 0 ; i < 20 ; i++)
        {
            // ROW INDEX (GLOBAL)
            idx ri = t.getNode( it,i%4 ) + npLC*( i/4 );
            ri = qn.getEquNode( ri );
            if (ri != NOT_AN_INDEX) // IF NOT FIXED
            {
#pragma omp atomic
                Ln[ri] += lL[i]*BIGNUM;
            }
        }//end for i


    }// end for it


}


