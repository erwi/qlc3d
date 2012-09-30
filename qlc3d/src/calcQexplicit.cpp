#include <qlc3d.h>
#include <globals.h>
#include <shapefunction3d.h>
#include <qassembly_macros.h>
#define	BIGNUM 2e16
//#include <sparsematrix.h>
//#include <solutionvector.h>


void assembleLocalMassMatrix(double lK[4][4],
                             Geometry &geom,
                             const idx it,
                             const Shape4thOrder &shapes)
{
    // ELEMENT MASS MATRIX ASSEMBLY. MATRIX IS FILLED WITH
    // INTEGRAL OF ROW*COL ELEMENTS

    memset(lK,0,4*4*sizeof(double) );
    double Jdet = geom.t->getDeterminant(it);
    for (idx igp = 0 ; igp < shapes.ngp ; ++igp)
    {
        const double mul = shapes.w[igp]*Jdet;
        for (idx i = 0 ; i < 4 ; ++i)
        {
            for (idx j = 0 ; j < 4 ; ++j)
                lK[i][j] += mul*shapes.sh1[igp][i]*shapes.sh1[igp][j];
        }// end for i
    }// end for igp
}

void fillInMassMatrix(SparseMatrix &K,
                      Geometry &geom,
                      SolutionVector &q // NEEDED ONLY FOR EquNodes
                      )
{
    idx N = q.getnFreeNodes();
    Shape4thOrder shapes;
    // FOR EACH TETRAHEDRAL ELEMENT
    for (idx it = 0; it < geom.t->getnElements(); ++it)
    {
        // IF NOT LC ELEMENT, IGNORE IT
        if (geom.t->getMaterialNumber(it) != MAT_DOMAIN1)
            continue;

        double lK[4][4]; // LOCAL ELEMENT MATRIX
        assembleLocalMassMatrix(lK, geom, it, shapes);

        for (int i = 0 ; i < 4 ; ++i)
        {
            const idx ri = geom.t->getNode(it, i);
            const idx eqr = q.getEquNode(ri);

            // IF FIXED NODE, DON'T INSERT TO MATRIX
            if ( eqr == NOT_AN_INDEX )
                continue;


            for (int j = 0 ; j < 4 ; ++j)
            {
                const idx rj = geom.t->getNode(it, j);
                const idx eqc = q.getEquNode(rj);
                if (eqc == NOT_AN_INDEX)
                    continue;
                K.sparse_add(eqr,eqc, lK[i][j]*BIGNUM );            // q1
                K.sparse_add(eqr+N, eqc+N, lK[i][j]*BIGNUM);        // q2
                K.sparse_add(eqr+2*N, eqc+2*N, lK[i][j]*BIGNUM);    // q3
                K.sparse_add(eqr+3*N, eqc+3*N, lK[i][j]*BIGNUM);    // q4
                K.sparse_add(eqr+4*N, eqc+4*N, lK[i][j]*BIGNUM);    // q5

            }// end for j
        }// end for i
    }// end for all elements
}

void assembleLocalRHS(SolutionVector &q,
                      SolutionVector &v,
                      LC &lc,
                      Geometry &geom,
                      const idx &it,
                      const Shape4thOrder shapes,
                      double lL[20]
                      )
{
    memset( lL,0, 20*sizeof(double) );


    double A = lc.A;
    double B = lc.B;
    double C = lc.C;
    double deleps = lc.eps_par-lc.eps_per;
    double rt6 = sqrt(6.0);
    double rt2 = sqrt(2.0);

    // LOCAL COPY OF ELEMENT
    idx tt[4] = {geom.t->getNode(it,0),
                 geom.t->getNode(it,1),
                 geom.t->getNode(it,2),
                 geom.t->getNode(it,3)};
    double Jdet = geom.t->getDeterminant(it)*1e18 ; // SCALE BACK TO METRES FOR NOW...
    //1. Calculate Inverse Jacobian - for 1st order elements can be done outside integration loop -> igp = 0
    double xr(0),xs(0),xt(0),yr(0),ys(0),yt(0),zr(0),zs(0),zt(0);
    double *p = geom.getPtrTop();
    const double pp[4][3] ={ {p[tt[0]*3] , p[tt[0]*3+1] , p[tt[0]*3+2]} ,
                             {p[tt[1]*3] , p[tt[1]*3+1] , p[tt[1]*3+2]} ,
                             {p[tt[2]*3] , p[tt[2]*3+1] , p[tt[2]*3+2]} ,
                             {p[tt[3]*3] , p[tt[3]*3+1] , p[tt[3]*3+2]} };

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

    //Inverse Jacobian
    //Jdet = fabs(xr*ys*zt-xr*zs*yt+xs*yt*zr-xs*yr*zt+xt*yr*zs-xt*ys*zr);
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
    const double qbuff[6*4] = {q.getValue(tt[0],0), q.getValue(tt[0],1), q.getValue(tt[0],2), q.getValue(tt[0],3), q.getValue(tt[0],4),
                               q.getValue(tt[1],0), q.getValue(tt[1],1), q.getValue(tt[1],2), q.getValue(tt[1],3), q.getValue(tt[1],4),
                               q.getValue(tt[2],0), q.getValue(tt[2],1), q.getValue(tt[2],2), q.getValue(tt[2],3), q.getValue(tt[2],4),
                               q.getValue(tt[3],0), q.getValue(tt[3],1), q.getValue(tt[3],2), q.getValue(tt[3],3), q.getValue(tt[3],4),
                               v.getValue(tt[0])  , v.getValue(tt[1]),   v.getValue(tt[2]),   v.getValue(tt[3])
                              };
    for (idx igp = 0; igp < shapes.ngp ; ++igp)
    {
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
        for(int i=0;i<4;++i)
        {
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
        // BULK RHS TERMS
        double L[5] = {0,0,0,0,0};
        RHS_THERMOTROPIC(L);
        RHS_DIELECTRIC(L);
        ADD_RHS_BULK_TERMS(lL,L);
    }//end for igp

}

void assembleRHS(SolutionVector &q,
                 SolutionVector &v,
                 Geometry &geom,
                 LC &lc,
                 Alignment &Alignment,
                 std::vector<double> &L
                 )
{
    idx npLC = q.getnDoF();
    const Shape4thOrder shapes;
    for (idx it = 0 ; it < geom.t->getnElements() ; it++)
    {
        double lL[20];
        // ONLY LC ELEMENTS
        if (geom.t->getMaterialNumber(it) != MAT_DOMAIN1 )
            continue;

        assembleLocalRHS(q,v,lc,geom,it, shapes, lL);

        //ASSEBLE LOCAL RHS ELEMENT CONRIBUTION
        for (idx i = 0 ; i < 20 ; ++i)
        {
            const idx ri = geom.t->getNode(it, i%4) + npLC*(i/4); // LOCAL TO GLOBAL
            const idx eqr = q.getEquNode(ri);
            if (eqr == NOT_AN_INDEX)
                continue;
            L[eqr] += lL[i]*BIGNUM;
        }// end for i
    }// end fot it


}


double calcQExplicit(SolutionVector &q,
                     SolutionVector &v,
                     SparseMatrix &K,
                     Geometry &geom,
                     LC &lc,
                     Alignment &alignment,
                     Simu &simu,
                     Settings &settings)
{

    // ON FIRST ITERATION OR IF MESH HAS BEEN MODFIED
    // VALUES MUST BE FILLED INTO K
    if ( ( simu.IsMeshModified() ) ||
         ( simu.getCurrentIteration() == 1 ) )
        fillInMassMatrix(K, geom, q);

    idx nd = K.rows;
    std::vector<double> M(nd,0);

    double *dq = new double [nd];
    double *qtemp = new double[nd];

    std::vector<double> L(nd,0);


    K.lump(M);

    idx iter, maxiter;
    double maxdq = 0;
    maxiter = 1000;
    for (iter = 0 ; iter < maxiter ; ++iter)
    {
        maxdq = 0;
        assembleRHS(q,v,geom,lc,alignment,L);


        for (idx i = 0 ; i < nd ; ++i)
        {
            dq[i] =-1.0*simu.getdt() *  ( L[i] / M[i] );
        }
        idx npLC = q.getnDoF();


        for (idx d = 0 ; d < 5 ; ++d)
        {
            for (idx i = 0 ; i < npLC ; ++i)
            {
                const idx n = i +d*npLC;
                const idx effDoF = q.getEquNode(n);

                if (effDoF == NOT_AN_INDEX)
                    continue;

                const double dqj = dq[effDoF];
                q.Values[n] += dqj;

                if (fabs(dqj) > maxdq ) maxdq = fabs(dqj);

            }
        }

        //printf("iter : %i , maxdq : %f\n", iter, maxdq);
    }
    printf("maxdq = %f", maxdq);
    // SOLVE: K*dq=L
    //UPDATE: Q += dq

    //K.PrintDiagonal();
    //K.PrintMatrix();
    //exit(1);


    delete [] dq;

    return maxdq;
}
