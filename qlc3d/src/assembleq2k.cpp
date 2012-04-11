#include <assembleq2k.h>
#include <globals.h>
#include <shapefunction3d.h>
#include <material_numbers.h>
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
                        double lL[20] )                 // LOCAL R.H.S.
{
    memset(lK, 0 , 20*20*sizeof(double) );
    memset(lL, 0 , 20*sizeof(double) );

    double lI[20][20]; // LOCAL IDENTITY MATRIX
    memset(lL, 0 , 20*sizeof(double) );

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
    double invJdet = 1.0 / t.getDeterminant(it)*1e12;
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


    }// end for igp

}

void assemble_volumes2K(SparseMatrix& K,
                        double* L,
                        const SolutionVector& q,
                        SolutionVector& v,
                        LC& mat_par,
                        Simu& simu,
                        const Mesh& t,
                        double* p
                        )
{
    idx npLC = q.getnDoF();
    Shape4thOrder shapes;

    for (idx it = 0 ; it < t.getnElements() ; it++)
    {
        if (t.getMaterialNumber(it) != MAT_DOMAIN1 )    // ONLY LC ELEMENTS ARE ASSEMBLED HERE
            continue;

        double lK[20][20];  // LOCAL MATRIX
        double lL[20];      // LOCAL R.H.S. VECTOR

        localKL_2K(p, t, it, q, v, mat_par, simu, shapes, lK, lL );


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


