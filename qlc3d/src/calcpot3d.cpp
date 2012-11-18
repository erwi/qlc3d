#include <calcpot3d.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#include <solutionvector.h>
#include <lc.h>
#include <sparsematrix.h>
#include <settings.h>
#include <geometry.h>
#include <shapefunction3d.h>
#include <shapefunction2d.h>
#ifndef COMPLEX
#define COMPLEX std::complex<double>
#endif
//#include <compcol_double.h>	// compressed column matrix
//#include <cg.h>
//#include <gmres.h>
//#include <icpre_double.h>
//#include <diagpre_double.h>
//#include <ilupre_double.h>

//#include MATRIX_H

// SPAMTRIX INCLUDES
//#include <setup.h>
#include <ircmatrix.h>
#include <vector.h>
#include <iterativesolvers.h>
#include <luincpreconditioner.h>






const idx npt = 4; //Number of Points per Tetrahedra

double rt2 = sqrt(2.0);
double rt6 = sqrt(6.0);

//void init_shapes_surf();
void potasm(SolutionVector *v, SolutionVector* q, LC* lc,Mesh *mesh, double *p, SparseMatrix *K, double* L);


void assemble_volume(double *p,
                     SolutionVector *v,
                     SolutionVector* q,
                     LC* lc,
                     Mesh *mesh,
                     SpaMtrix::IRCMatrix &K,
                     SpaMtrix::Vector &L,
                     Electrodes* electrodes);
void assemble_Neumann(double *p,
                      SolutionVector *v,
                      SolutionVector* q,
                      LC* lc,
                      Mesh *mesh,
                      Mesh* surf_mesh,
                      SpaMtrix::IRCMatrix &K,
                      SpaMtrix::Vector &L);


void Pot_GMRES(SpaMtrix::IRCMatrix &K,
               SpaMtrix::Vector &B,
               SpaMtrix::Vector &X,
               Settings* settings );


// DEBUG FUNCTION FOR PRINTING VALUES OF LOCAL MATRIX
void printlK(double* lK , int size){
    for (int r = 0 ; r < size ; r++){
        for (int c = 0 ; c < size ; c++){
            cout << lK[r + size*c] << " ";
        }
	cout << endl;
    } //end for r
}


// DECLARATION ONLY
void setUniformEField( Electrodes& electrodes, SolutionVector& v, double* p);

void calcpot3d(
        SpaMtrix::IRCMatrix &K,
        SolutionVector *v,
        SolutionVector* q,
        LC* lc,
        Geometry& geom,
        Settings* settings,
        Electrodes* electrodes)
{
    // First check whether potential calculation is actually needed...

    // NO NEED TO CALCULATE POTENTIAL IF...
    if ( (v->getnFixed() == 0 ) ||      // no fixed potential nodes OR
         (!electrodes->getCalcPot() ) ) // no need to calculate potential
    {
        v->setValuesTo(0.0); // if no potential calculation, set all values to zero
        if ( electrodes->isEField() )
        {
            setUniformEField( *electrodes, *v, geom.getPtrTop());
        }
        return;
    }
    K = 0.0; // clears values but keeps sparsity structure

    SpaMtrix::Vector L(v->getnFreeNodes());
    SpaMtrix::Vector V(v->getnFreeNodes());
    // Assemble system
    assemble_volume(geom.getPtrTop(),v,q,lc,geom.t, K , L, electrodes);
    assemble_Neumann(geom.getPtrTop() , v , q , lc , geom.t , geom.e , K , L);

//#ifdef DEBUG
//    K->DetectZeroDiagonals();
//#endif
    // GMRES SOLVER
    Pot_GMRES(K,L,V, settings);

    // COPY NON-FIXED VALUES BACK TO SOLUTIONVECTOR
    for (idx i = 0 ; i < v->getnDoF() ; i++)
    {
        idx ind = v->getEquNode(i);
        // EQU NODES OF FIXED DOFS ARE ALL "NOT_AN_INDEX"
        if (ind != NOT_AN_INDEX)
            v->setValue(i,0, V[ind] );
    }// end for i
}
//end calcpot3d

inline void localKL(
    double *p,
    int *tt,
    double lK[npt][npt],
    double lL[npt],
    int it,
    Mesh* mesh,
    SolutionVector* q,
    LC* lc,
    Electrodes* electrodes,
    const Shape4thOrder& shapes)

{
    idx i,j;
    double eper, deleps;
    double S0 = lc->getS0();

    eper 	= 0;
    deleps	= 0;
    if (mesh->getMaterialNumber(it) == MAT_DOMAIN1){ // if LC element
        eper = lc->eps_per / S0;
        deleps = (lc->eps_par - lc->eps_per) /S0;
    }
    else{ // otherwise dielectric
        idx ind_de = mesh->getDielectricNumber(it) - 1; // -1 for 0 indexing
        eper = electrodes->getDielectricPermittivity(ind_de);
    }

    //printf("deleps = %e\n", deleps);
    memset(lK,0,4*4*sizeof(double));
    memset(lL,0,4*sizeof(double));

    //cout << "element:"<< it << endl;
    //printlK( &lK[0][0] , 4);
    //printf("gettin")
    double Jdet= mesh->getDeterminant(it);
    // Jacobian
    double xr,xs,xt,yr,ys,yt,zr,zs,zt;
    xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;

    for (i=0; i<npt; i++)
    {
        xr+=shapes.sh1r[0][i]*p[ tt[i]*3 + 0 ]*1e-6;
        xs+=shapes.sh1s[0][i]*p[ tt[i]*3 + 0 ]*1e-6;
        xt+=shapes.sh1t[0][i]*p[ tt[i]*3 + 0 ]*1e-6;

        yr+=shapes.sh1r[0][i]*p[ tt[i]*3 + 1 ]*1e-6;
        ys+=shapes.sh1s[0][i]*p[ tt[i]*3 + 1 ]*1e-6;
        yt+=shapes.sh1t[0][i]*p[ tt[i]*3 + 1 ]*1e-6;

        zr+=shapes.sh1r[0][i]*p[ tt[i]*3 + 2 ]*1e-6;
        zs+=shapes.sh1s[0][i]*p[ tt[i]*3 + 2 ]*1e-6;
        zt+=shapes.sh1t[0][i]*p[ tt[i]*3 + 2 ]*1e-6;
    }//end for i

    if (Jdet<0) Jdet = -Jdet;

    double Jinv[3][3]={{(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet}
                       ,{(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet}
                       ,{(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet}};

    double Sh[4],dSh[4][3];
    // x,y,z derivatives of shape functions
    for(i=0;i<4;i++){
        dSh[i][0]=shapes.sh1r[0][i]*Jinv[0][0]+
                shapes.sh1s[0][i]*Jinv[1][0]+
                shapes.sh1t[0][i]*Jinv[2][0];

        dSh[i][1]=shapes.sh1r[0][i]*Jinv[0][1]+
                shapes.sh1s[0][i]*Jinv[1][1]+
                shapes.sh1t[0][i]*Jinv[2][1];

        dSh[i][2]=shapes.sh1r[0][i]*Jinv[0][2]+
                shapes.sh1s[0][i]*Jinv[1][2]+
                shapes.sh1t[0][i]*Jinv[2][2];
    }//end for i

    for (unsigned int igp=0; igp<shapes.ngp; igp++)
    {
        Sh[0] = shapes.sh1[igp][0];
        Sh[1] = shapes.sh1[igp][1];
        Sh[2] = shapes.sh1[igp][2];
        Sh[3] = shapes.sh1[igp][3];

        double e11,e22,e33,e12,e13,e23;
        e11=0;e22=0;e33=0;e12=0;e13=0;e23=0;
        //if (1){//
        if ( mesh->getMaterialNumber(it) != MAT_DOMAIN1){ // if this element is not LC
            // only set diagonal permittivities to non-zero
            for ( i = 0 ; i<npt ;++i){
                e11+= Sh[i]*eper;
                e22+= Sh[i]*eper;
                e33+= Sh[i]*eper;
            }

        }
        else{// otherwise LC ->
            for(i=0;i<npt;i++){
                e11+= Sh[i]*(((2.0/3.0/S0)*(-q->getValue(tt[i],0) /rt6 + q->getValue(tt[i],1)/rt2)+(1.0/3.0))*deleps + eper);	//~nx*nx
                e22+= Sh[i]*(((2.0/3.0/S0)*(-q->getValue(tt[i],0)/rt6 - q->getValue(tt[i],1)/rt2)+(1.0/3.0))*deleps + eper);	//~ny*ny
                e33+= Sh[i]*(((2.0/3.0/S0)*(2.0*q->getValue(tt[i],0)/rt6)	     +(1.0/3.0))*deleps + eper);			//~nz*nz
                e12+= Sh[i]*(2.0/3.0/S0)*(q->getValue(tt[i],2)/rt2)			*deleps;						//~nx*ny
                e13+= Sh[i]*(2.0/3.0/S0)*(q->getValue(tt[i],4)/rt2)			*deleps;						//~nx*nz
                e23+= Sh[i]*(2.0/3.0/S0)*(q->getValue(tt[i],3)/rt2)			*deleps;
            }
            //cout << "e =" << e11<<","<< e22<<","<< e33 <<","<< e12 <<","<< e13 <<","<< e23 << endl;
        }
        // Local K and L
        double mul=shapes.w[igp]*Jdet;

        for (i=0; i<4; i++)
        {
            for (j=0; j<4; j++)
            {

                lK[i][j]+=mul*(
                            dSh[i][0]*dSh[j][0]*e11+
                            dSh[i][1]*dSh[j][1]*e22+
                            dSh[i][2]*dSh[j][2]*e33+

                            dSh[i][0]*dSh[j][1]*(e12)+
                            dSh[i][1]*dSh[j][0]*(e12)+
                            dSh[i][1]*dSh[j][2]*(e23)+
                            dSh[i][2]*dSh[j][1]*(e23)+
                            dSh[i][0]*dSh[j][2]*(e13)+
                            dSh[i][2]*dSh[j][0]*(e13)
                            );
            }//end for j
        }//end for i
    }//end for igp


}// end void localKL
void localKL_N(
    double* p,
    int* tt,
    double lK[npt][npt],
    double lL[npt],
    int it,
    int index_to_Neumann,
    Mesh*  mesh,
    Mesh* surf_mesh,
    SolutionVector* q,
    LC* lc,
    const ShapeSurf4thOrder& shapes)
{
    int i,j;
    double S0=lc->getS0();

    double eper   = 4.5;
    double deleps = 0.0;	// default for dielectric material

    if (mesh->getMaterialNumber(index_to_Neumann) == MAT_DOMAIN1)
    {
        eper   = lc->eps_per / S0;
        deleps = (lc->eps_par - lc->eps_per) / S0;
    }
    else
    {
        printf("NON LC NEUMANN\n"); // THIS SHOULD NEVER HAPPEN
    }
    memset(lK,0,npt*npt*sizeof(double));
    memset(lL,0,4*sizeof(double));
    double n[3];

    surf_mesh->CopySurfaceNormal(it,n);

#ifdef DEBUG
    if ( abs(n[0]*n[0] + n[1]*n[1] + n[2]*n[2] -1.0 ) > 0.01 )
    {
        printf("bad triangle normal in neumann tri %i ", it);
        exit(1);
    }
#endif
    double eDet = surf_mesh->getDeterminant(it);
    double Jdet = mesh->getDeterminant(index_to_Neumann);

#ifdef DEBUG
    if (eDet <= 0 )
        printf("ZERO EDET\n");
    if (Jdet <= 0 )
        printf("ZERO JDET\n");
#endif

    double xr,xs,xt,yr,ys,yt,zr,zs,zt;
    xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;

    for (i=0; i<4; i++){
        xr+=shapes.sh1r[0][i]*p[ tt[i] *3+0]*1e-6; //  <- tt is reordered volume element
        xs+=shapes.sh1s[0][i]*p[ tt[i] *3+0]*1e-6;
        xt+=shapes.sh1t[0][i]*p[ tt[i] *3+0]*1e-6;

        yr+=shapes.sh1r[0][i]*p[ tt[i] *3+1]*1e-6;
        ys+=shapes.sh1s[0][i]*p[ tt[i] *3+1]*1e-6;
        yt+=shapes.sh1t[0][i]*p[ tt[i] *3+1]*1e-6;

        zr+=shapes.sh1r[0][i]*p[ tt[i] *3+2]*1e-6;
        zs+=shapes.sh1s[0][i]*p[ tt[i] *3+2]*1e-6;
        zt+=shapes.sh1t[0][i]*p[ tt[i] *3+2]*1e-6;
    }//end for i

    double Jinv[3][3]={{ (zt*ys-yt*zs)/Jdet , (xt*zs-zt*xs)/Jdet , (xs*yt-ys*xt)/Jdet}
                       ,{ (yt*zr-zt*yr)/Jdet , (zt*xr-xt*zr)/Jdet , (xt*yr-yt*xr)/Jdet}
                       ,{ (yr*zs-ys*zr)/Jdet , (xs*zr-xr*zs)/Jdet , (ys*xr-xs*yr)/Jdet}};

    double Sh[4],dSh[4][3];
    for(i=0;i<4;i++){
    dSh[i][0]=shapes.sh1r[0][i]*Jinv[0][0]+
            shapes.sh1s[0][i]*Jinv[1][0]+
            shapes.sh1t[0][i]*Jinv[2][0];
    dSh[i][1]=shapes.sh1r[0][i]*Jinv[0][1]+
            shapes.sh1s[0][i]*Jinv[1][1]+
            shapes.sh1t[0][i]*Jinv[2][1];
    dSh[i][2]=shapes.sh1r[0][i]*Jinv[0][2]+
            shapes.sh1s[0][i]*Jinv[1][2]+
            shapes.sh1t[0][i]*Jinv[2][2];
    }//end for i



    // Jacobian
    for (unsigned int igp=0; igp<shapes.ngps; igp++)
    {
        Sh[0] = shapes.sh1[igp][0];
        Sh[1] = shapes.sh1[igp][1];
        Sh[2] = shapes.sh1[igp][2];
        Sh[3] = shapes.sh1[igp][3];

        double e11=0,e22=0,e33=0,e12=0,e23=0,e13=0;
        for(i=0;i<4;i++){
            e11+= shapes.sh1[igp][i]*(((2.0/3.0/S0)*(-q->getValue(tt[i],0) /rt6 + q->getValue(tt[i],1)/rt2)+(1.0/3.0))*deleps + eper);	//~nx*nx
            e22+= shapes.sh1[igp][i]*(((2.0/3.0/S0)*(-q->getValue(tt[i],0)/rt6 - q->getValue(tt[i],1)/rt2)+(1.0/3.0))*deleps + eper);	//~ny*ny
            e33+= shapes.sh1[igp][i]*(((2.0/3.0/S0)*(2.0*q->getValue(tt[i],0)/rt6)	     +(1.0/3.0))*deleps + eper);			//~nz*nz
            e12+= shapes.sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],2)/rt2)			*deleps;						//~nx*ny
            e13+= shapes.sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],4)/rt2)			*deleps;						//~nx*nz
            e23+= shapes.sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],3)/rt2)			*deleps;					//~ny*nz
        }//end for i
        double mul=shapes.w[igp]*eDet;
        for (i=0; i<4; i++){
            for (j=0; j<4; j++){
                lK[i][j]+=mul*Sh[i]*(((e11-1)*dSh[j][0] + e12*dSh[j][1] + e13*dSh[j][2])*n[0]
                                     +(e12*dSh[j][0] + (e22-1)*dSh[j][1] + e23*dSh[j][2])*n[1]
                                     +(e13*dSh[j][0] + e23*dSh[j][1] + (e33-1)*dSh[j][2])*n[2] );
            }//end for j
        }//end for i
    }//end for igp
}
// end void localKL




void assemble_volume(
    double *p,
    SolutionVector *v,
    SolutionVector* q,
    LC* lc,
    Mesh *mesh,
    SpaMtrix::IRCMatrix &K,
    SpaMtrix::Vector &L,
    Electrodes* electrodes)
{
    idx it;
    Shape4thOrder shapes;
    #pragma omp parallel for
    for (it=0; it< mesh->getnElements(); it++)
    {
        double lK[npt][npt];
        double lL[npt];
        int tmat;
        int t[4] = {0,0,0,0};
        t[0] =  mesh->getNode(it,0);
        t[1] =  mesh->getNode(it,1);
        t[2] =  mesh->getNode(it,2);
        t[3] =  mesh->getNode(it,3);
        tmat = mesh->getMaterialNumber(it);

        localKL(p,t,lK,lL,it, mesh,q,lc,electrodes, shapes);
        //printlK( &lK[0][0] , npt);
        for (idx i=0; i<npt; i++) // FOR ROWS
        {
            idx ri=v->getEquNode(t[i]);

            // RHS FIXED NODE HANDLING
            if ( ri == NOT_AN_INDEX ) // IF THIS NODE IS FIXED
            {
                for (int j = 0; j<4 ; j++)// SET CONTRIBUTION TO CONNECTED *FREE* NODES
                {
                    idx nc = v->getEquNode( t[j] ); // INDEX TO CONNECTED NODE DEGREE OF FREEDOM POSITION
                    if (nc != NOT_AN_INDEX )
                    {
#ifndef DEBUG
#pragma omp atomic
#endif
                        L[ nc ] -= lK[i][j]*v->getValue(t[i]); // L = -K*v
                    }
                }

            }// END IF ROW NODE IS FIXED

            if ( ri != NOT_AN_INDEX )
                for (idx j=0; j<npt  ; j++) // FOR COLUMNS
                {
                    idx rj=v->getEquNode(t[j]);
                    if (rj != NOT_AN_INDEX )
                    {
                        K.sparse_add(rj,ri,lK[j][i]);
                    }
                }//end for j
        }//end for i
    }//end for it
}// end void assemble_volume

////////////////////////////////////////////////////////////////////////////////
// assemble Neumann boundaries
void assemble_Neumann(
    double *p,
    SolutionVector *v,
    SolutionVector* q,
    LC* lc,
    Mesh *mesh,
    Mesh* surf_mesh,
    SpaMtrix::IRCMatrix &K,
    SpaMtrix::Vector &L)
{
    idx it = 0;
    ShapeSurf4thOrder shapes;
#ifndef DEBUG
#pragma omp parallel for
#endif
    for (it=0; it < surf_mesh->getnElements(); it++)
    {
        int index_to_Neumann = surf_mesh->getConnectedVolume(it);

        if ( (index_to_Neumann > -1) && (surf_mesh->getMaterialNumber(it) == MAT_NEUMANN))// if connected to LC tet
        {

            double lK[4][4];
            double lL[4];
            int ee[3] = {   surf_mesh->getNode(it,0) ,
                            surf_mesh->getNode(it,1) ,
                            surf_mesh->getNode(it,2) } ;

            int tt[4] = { tt[0] =   mesh->getNode(index_to_Neumann,0),
                          mesh->getNode(index_to_Neumann,1),
                          mesh->getNode(index_to_Neumann,2),
                          mesh->getNode(index_to_Neumann,3)};

            int intr = 0;//find  index to internal node

            for (int i=0;i<4;i++){
                if ( (tt[i]!= ee[0]) && (tt[i]!= ee[1]) && (tt[i]!= ee[2]) ){
                    intr = i;
                    break;
                }
            }

            int ti[4] = { ee[0], ee[1], ee[2], tt[intr] }; // reordered local element, internal node is always last

            localKL_N(p,&ti[0], lK , lL, it,index_to_Neumann, mesh, surf_mesh, q, lc, shapes);

            for (int i=0; i<4; i++)
            {
                idx ri=v->getEquNode(ti[i]);
                if (ri == NOT_AN_INDEX ) // HANDLE FIXED NODE
                {
                    for (int j = 0; j<4 ; j++)
                    {
                        idx cr = v->getEquNode( ti[j] ); // CONNECTED NODE DOF ORDER
                        if (cr != NOT_AN_INDEX )
                        {
#ifndef DEBUG
#pragma omp atomic
#endif
                            L[cr ] -= lK[j][i]*v->getValue(ti[i]);
                        }
                    }
                }// END HANDLE FIXED NODE
                else // HANDLE FREE NODE
                {
                    for (int j=0; j<4; j++) // FOR COLUMNS
                    {
                        idx rj=v->getEquNode(ti[j]);
                        if (rj != NOT_AN_INDEX )// NON-FIXED NODE
                        {
                            K.sparse_add( rj,ri,lK[j][i] );
                        }
                    }//end for j
                }// END HANDLE FREE NODES
            }//end for i
        }//end if LC
    }//end for it
}//end void assemble_Neumann


void setUniformEField(Electrodes &electrodes, SolutionVector &v, double *p)
{
    // SETS POTENTIAL VALUES IN V SUCH THA A UNIFORM E-FIELD IS CREATED
    // CENTRE OF STRUCTURE IS ASSUMED TO BE AT 0V, VOLTAGE VALUES FOR NODES
    // IN THE DIRECTION OF THE EFIELD ARE SET ACCORDING TO THEIR DISTANCE
    // TO THE CENTRE


    // GET E-FIELD DIRECTION VECTOR AND MAGNITUDE
    double*E = &(electrodes.EField[0]);
    double Emag = sqrt( E[0]*E[0] + E[1]*E[1] + E[2]*E[2] );
    double Ehat[3] = { E[0] / Emag,
                       E[1] / Emag,
                       E[2] / Emag};

    //1. CALCULATE CENTRE OF STRUCTURE
    int np = v.getnDoF();
    double xmax = 0.0;
    double ymax = 0.0;
    double zmax = 0.0;
    for (int i = 0 ; i < np ; i++)
    {
        xmax = p[3*i + 0 ] > xmax? p[3*i + 0] : xmax;
        ymax = p[3*i + 1 ] > ymax? p[3*i + 1] : ymax;
        zmax = p[3*i + 2 ] > zmax? p[3*i + 2] : zmax;
    }

    double centre[3] = {xmax/2.0 , ymax/2.0 , zmax / 2.0}; // centre of structure

    //2. LOOP OVER EACH NODE AND CACLCULATE ITS DISTANCE TO CENTRE

    for (int i = 0 ; i < np ; i++)
    {
        double* pos = &p[i*3]; // shortcut to this node
        double vec[3] = { centre[0] - pos[0],
                          centre[1] - pos[1],
                          centre[2] - pos[2] }; // vector from centre to node i

        // want distance along EField, i.e. dot product
        double dist = vec[0]*Ehat[0] + vec[1]*Ehat[1] + vec[2]*Ehat[2];


        // set potential value as distance*magnitude
        v.setValue(i,0, dist*Emag + v.getValue(i) );
    }



}

void Pot_GMRES(SpaMtrix::IRCMatrix &K, SpaMtrix::Vector &B, SpaMtrix::Vector &X, Settings* settings )
{
    /*!
        Solves the Linear simulatenous equation Ax=b using the GMRES method
    */


    idx size = K.getNumRows();
    idx maxiter 	= settings->getV_GMRES_Maxiter();
    idx restart 	= settings->getV_GMRES_Restart();

    maxiter = maxiter < size ? maxiter:size;
    restart = restart < maxiter ? restart :maxiter;
    double toler 	= settings->getV_GMRES_Toler();

    SpaMtrix::LUIncPreconditioner LU(K); // DOES THIS HAVE TO BE RECOMPUTED EACH TIME??
    SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);
    if (!solver.gmres(K, X, B, LU) )
        printf("GMRES did not converge in %i iterations \nTolerance achieved is %f\n",solver.maxIter,solver.toler);

}

