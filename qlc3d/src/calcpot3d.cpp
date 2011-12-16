# include <math.h>
# include <math.h>
# include <omp.h>
#include <time.h>

#include <solutionvector.h>
#include <lc.h>
#include <sparsematrix.h>
#include <settings.h>
#include <geometry.h>
#ifndef COMPLEX
    #define COMPLEX std::complex<double>
#endif
#include <compcol_double.h>	// compressed column matrix
#include <cg.h>
#include <gmres.h>
#include <icpre_double.h>
#include <diagpre_double.h>
#include <ilupre_double.h>

#include MATRIX_H




const int	npt = 4; //Number of Points per Tetrahedra

double rt2 = sqrt(2.0);
double rt6 = sqrt(6.0);

void init_shapes_surf();
void potasm(SolutionVector *v, SolutionVector* q, LC* lc,Mesh *mesh, double *p, SparseMatrix *K, double* L);


void assemble_volume(double *p,SolutionVector *v,SolutionVector* q, LC* lc, Mesh *mesh, SparseMatrix *K, double* L, Electrodes* electrodes);
void assemble_Neumann(double *p,SolutionVector *v, SolutionVector* q, LC* lc,Mesh *mesh, Mesh* surf_mesh, SparseMatrix *K, double* L);

//void Pot_PCG(SparseMatrix *K, double *b, SolutionVector* sv, Settings* settings );
void Pot_PCG(SparseMatrix *K, double *b, double* V, Settings* settings );
void Pot_GMRES(SparseMatrix *K, double *b, double* V, Settings* settings );
void Pot_SuperLU(SparseMatrix *K, double *b, SolutionVector* sv, Settings* settings );



// DEBUG FUNCTION FOR PRINTING VALUES OF LOCAL MATRIX
void printlK(double* lK , int size){
   for (int r = 0 ; r < size ; r++){
       for (int c = 0 ; c < size ; c++){
	   cout << lK[r + size*c] << " ";
       }
	cout << endl;
   } //end for r
}

// Gauss integration
// ---------------------------------------------------------
//     3D Gauss-Legendre weights for N = 11, D = 4
// ---------------------------------------------------------

const int ngp=11;
const double a=(1+sqrt(5.0/14.0))/4.0;
const double b=(1-sqrt(5.0/14.0))/4.0;

const double gp[ngp][4]={
	{0.25	  , 0.25	,	0.25	,0.25},
	{11.0/14.0     ,	1.0/14.0	,	1.0/14.0	,1.0/14.0},
	{1.0/14.0      ,	11.0/14.0	,	1.0/14.0	,1.0/14.0},
	{1.0/14.0	  ,	1.0/14.0	,	11.0/14.0   ,1.0/14.0},
	{1.0/14.0	  , 1.0/14.0	,	1.0/14.0	,11.0/14.0},
	{a		  ,	a		,	b       ,b},
	{a		  , b		,   a       ,b},
	{a        , b       ,   b	    ,a},
	{b		  , a       ,   a       ,b},
	{b		  , a      	,   b		,a},
	{b		  , b		,   a		,a}};
const double w11 = -74.0/5625.0;
const double w12 = 343.0/45000.0;
const double w13 = 56.0/2250.0;
static double w[ngp]={w11,w12,w12,w12,w12,w13,w13,w13,w13,w13,w13};

// shape functions used to be 'static', but this caused run-time error when compiled with MinGW g++ 4.4
// and -march=native flag.
double sh1[ngp][4]; // P1 Shape functions
double sh1r[ngp][4]; // P1 Shape functions r-derivatives
double sh1s[ngp][4]; // P1 Shape functions s-derivatives
double sh1t[ngp][4]; //P1 shape functions t-derivative

void init_shapes()
{
    // TERAHEDRA SHAPE FUCNTIONS AND ITS DERIVATIVES
    for (int i=0; i<ngp; i++) {
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
	 //   cout << i << endl;
    }
}


// ---------------------------------------------------------
//     2D Gauss-Legendre weights for Neumann Boundary integrals N = 6, D = 4
// ---------------------------------------------------------
const int ngps=6;
const double gps[ngps][2] ={
		{0.8168476, 0.09157621},
		{0.09157621,0.8168476},
		{0.09157621,0.09157621},
		{0.1081030, 0.4459485},
		{0.4459485, 0.1081030},
		{0.4459485, 0.4459485}};
static double wsurf[ngps]={ 0.05497587, 0.05497587, 0.05497587, 0.1116908, 0.1116908, 0.1116908};

void init_shapes_surf() // surface integral shape functions
{
	memset(sh1,0,ngps*4*sizeof(double));
	memset(sh1r,0,ngps*4*sizeof(double));
	memset(sh1s,0,ngps*4*sizeof(double));
	memset(sh1t,0,ngps*4*sizeof(double));
	for (int i=0; i<ngps; i++) {
                //cout << "shapes_surf" << i << endl;
            // P1 Shape functions
		sh1[i][0]=1-gps[i][0]-gps[i][1];
		sh1[i][1]=gps[i][0];
		sh1[i][2]=gps[i][1];
		sh1[i][3]=0;//gps[i][2];
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
void setUniformEField( Electrodes& electrodes, SolutionVector& v, double* p);

void calcpot3d(
        SparseMatrix* K,
        SolutionVector *v,
        SolutionVector* q,
        LC* lc,
        //Mesh *mesh,
        //Mesh* surf_mesh,
        //double *p,
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

    K->setAllValuesTo(0); // clears values but keeps sparsity structure

    double* L = (double*) malloc(v->getnFreeNodes()*sizeof(double) );
    double* V = (double*) malloc(v->getnFreeNodes()*sizeof(double) );
    if ( (L == NULL) || (V == NULL) )
    {
        printf("error in %s, malloc returned NULL - bye!\n", __func__);
        exit(1);
    }

    memset(L,0,v->getnFreeNodes()*sizeof(double) );
    memset(V,0,v->getnFreeNodes()*sizeof(double) );
    // Assemble system

    init_shapes();
    assemble_volume(geom.getPtrTop(),v,q,lc,geom.t, K , L, electrodes);


    init_shapes_surf();
    assemble_Neumann(geom.getPtrTop() , v , q , lc , geom.t , geom.e , K , L);

#ifdef DEBUG
    K->DetectZeroDiagonals();
#endif

    // Solve System
    if (settings->getV_Solver() == V_SOLVER_PCG)
        Pot_PCG(K,L,V, settings);
    else if (settings->getV_Solver() == V_SOLVER_GMRES)
        Pot_GMRES(K,L,V, settings);
    else
    {
        printf("error - potential solver is set to %i\n", settings->getV_Solver());
        exit(1);
    }
    free(L);

    for (int i = 0 ; i < v->getnDoF() ; i++)
    {
        int ind = v->getEquNode(i);
        if (ind != SolutionVector::FIXED_NODE)
        {
           v->setValue(i,0, V[ind] );
        }
    }

    free(V);

    // UNIFORM FIELD
    if ( electrodes->isEField() )
        setUniformEField( *electrodes, *v, geom.getPtrTop() );

}
//end calcpot3d

void localKL(
	double *p,
	int *tt,
	double lK[npt][npt],
	double lL[npt],
	int it,
	Mesh* mesh,
	SolutionVector* q,
	LC* lc,
	Electrodes* electrodes){
    int i,j;
    double eper, deleps;
    double S0 = lc->getS0();

    eper 	= 0;
    deleps	= 0;
    if (mesh->getMaterialNumber(it) == MAT_DOMAIN1){ // if LC element
        eper = lc->eps_per / S0;
        deleps = (lc->eps_par - lc->eps_per) /S0;
    }
    else{ // otherwise dielectric
       int ind_de = mesh->getDielectricNumber(it) - 1; // -1 for 0 indexing
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

		for (i=0; i<npt; i++) {
			xr+=sh1r[0][i]*p[ tt[i]*3 + 0 ]*1e-6;
			xs+=sh1s[0][i]*p[ tt[i]*3 + 0 ]*1e-6;
			xt+=sh1t[0][i]*p[ tt[i]*3 + 0 ]*1e-6;

			yr+=sh1r[0][i]*p[ tt[i]*3 + 1 ]*1e-6;
			ys+=sh1s[0][i]*p[ tt[i]*3 + 1 ]*1e-6;
			yt+=sh1t[0][i]*p[ tt[i]*3 + 1 ]*1e-6;

			zr+=sh1r[0][i]*p[ tt[i]*3 + 2 ]*1e-6;
			zs+=sh1s[0][i]*p[ tt[i]*3 + 2 ]*1e-6;
			zt+=sh1t[0][i]*p[ tt[i]*3 + 2 ]*1e-6;
		}//end for i
		//(xr*ys*zt-xr*zs*yt+xs*yt*zr-xs*yr*zt+xt*yr*zs-xt*ys*zr);
		if (Jdet<0) Jdet = -Jdet;// printf("negative jacobian!\n");
		//cout << "Jdet = " << Jdet << endl;
		double Jinv[3][3]={{(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet}
						  ,{(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet}
						  ,{(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet}};

       double Sh[4],dSh[4][3];
	   // x,y,z derivatives of shape functions
	for(i=0;i<4;i++){
	    dSh[i][0]=sh1r[0][i]*Jinv[0][0]+sh1s[0][i]*Jinv[1][0]+sh1t[0][i]*Jinv[2][0];
            dSh[i][1]=sh1r[0][i]*Jinv[0][1]+sh1s[0][i]*Jinv[1][1]+sh1t[0][i]*Jinv[2][1];
	    dSh[i][2]=sh1r[0][i]*Jinv[0][2]+sh1s[0][i]*Jinv[1][2]+sh1t[0][i]*Jinv[2][2];
	}//end for i

	for (int igp=0; igp<ngp; igp++){
	   Sh[0] = sh1[igp][0];
	   Sh[1] = sh1[igp][1];
	   Sh[2] = sh1[igp][2];
	   Sh[3] = sh1[igp][3];

	   double e11,e22,e33,e12,e13,e23;
	   e11=0;e22=0;e33=0;e12=0;e13=0;e23=0;
          //if (1){//
           if ( mesh->getMaterialNumber(it) != MAT_DOMAIN1){ // if this element is not LC
               // only set diagonal permittivities to non-zero
               for ( i = 0 ; i<npt ;++i){
                   e11+= sh1[igp][i]*eper;
                   e22+= sh1[igp][i]*eper;
                   e33+= sh1[igp][i]*eper;
               }

           }
           else{// otherwise LC ->
               for(i=0;i<npt;i++){
                   e11+= sh1[igp][i]*(((2.0/3.0/S0)*(-q->getValue(tt[i],0) /rt6 + q->getValue(tt[i],1)/rt2)+(1.0/3.0))*deleps + eper);	//~nx*nx
                   e22+= sh1[igp][i]*(((2.0/3.0/S0)*(-q->getValue(tt[i],0)/rt6 - q->getValue(tt[i],1)/rt2)+(1.0/3.0))*deleps + eper);	//~ny*ny
                   e33+= sh1[igp][i]*(((2.0/3.0/S0)*(2.0*q->getValue(tt[i],0)/rt6)	     +(1.0/3.0))*deleps + eper);			//~nz*nz
                   e12+= sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],2)/rt2)			*deleps;						//~nx*ny
                   e13+= sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],4)/rt2)			*deleps;						//~nx*nz
                   e23+= sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],3)/rt2)			*deleps;
	       }
	       //cout << "e =" << e11<<","<< e22<<","<< e33 <<","<< e12 <<","<< e13 <<","<< e23 << endl;
	   }
       // Local K and L
       double mul=w[igp]*Jdet;

       for (i=0; i<4; i++) {
	   for (j=0; j<4; j++) {

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
  // printlK( &lK[0][0] , 4);

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
	LC* lc){
    int i,j;
    double S0=lc->getS0();

    double eper   = 4.5;
    double deleps = 0.0;	// default for dielectric material

    if (mesh->getMaterialNumber(index_to_Neumann) == MAT_DOMAIN1){
	eper   = lc->eps_per / S0;
	deleps = (lc->eps_par - lc->eps_per) / S0;
    }
    memset(lK,0,npt*npt*sizeof(double));
    memset(lL,0,4*sizeof(double));
    double n[3];

    surf_mesh->CopySurfaceNormal(it,n);
    if ( n[0]*n[0] + n[1]*n[1] + n[2]*n[2]  < 0.1 )
    {
	printf("zero normal in neumann tri %i ", it);
	exit(1);
    }


   // if ( ( p[ tt[0]*3 + 2 ] <= 0.01)&&
	// ( p[ tt[1]*3 + 2 ] <= 0.01)&&
	 //( p[ tt[2]*3 + 2 ] <= 0.01) ){

//	printf("normal[%i] = %f,%f,%f\n",it, n[0], n[1] , n[2]);

  //  }
//*/


    double eDet = surf_mesh->getDeterminant(it);
    double Jdet = mesh->getDeterminant(index_to_Neumann);

    double xr,xs,xt,yr,ys,yt,zr,zs,zt;
    xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;




    for (i=0; i<4; i++){
	xr+=sh1r[0][i]*p[ tt[i] *3+0]*1e-6; //  <- tt is reordered volume element
	xs+=sh1s[0][i]*p[ tt[i] *3+0]*1e-6;
	xt+=sh1t[0][i]*p[ tt[i] *3+0]*1e-6;

	yr+=sh1r[0][i]*p[ tt[i] *3+1]*1e-6;
	ys+=sh1s[0][i]*p[ tt[i] *3+1]*1e-6;
	yt+=sh1t[0][i]*p[ tt[i] *3+1]*1e-6;

	zr+=sh1r[0][i]*p[ tt[i] *3+2]*1e-6;
	zs+=sh1s[0][i]*p[ tt[i] *3+2]*1e-6;
	zt+=sh1t[0][i]*p[ tt[i] *3+2]*1e-6;
    }//end for i

    double Jinv[3][3]={{ (zt*ys-yt*zs)/Jdet , (xt*zs-zt*xs)/Jdet , (xs*yt-ys*xt)/Jdet}
			 ,{ (yt*zr-zt*yr)/Jdet , (zt*xr-xt*zr)/Jdet , (xt*yr-yt*xr)/Jdet}
			 ,{ (yr*zs-ys*zr)/Jdet , (xs*zr-xr*zs)/Jdet , (ys*xr-xs*yr)/Jdet}};

    double Sh[4],dSh[4][3];
    for(i=0;i<4;i++){
	dSh[i][0]=sh1r[0][i]*Jinv[0][0]+sh1s[0][i]*Jinv[1][0]+sh1t[0][i]*Jinv[2][0];
	dSh[i][1]=sh1r[0][i]*Jinv[0][1]+sh1s[0][i]*Jinv[1][1]+sh1t[0][i]*Jinv[2][1];
	dSh[i][2]=sh1r[0][i]*Jinv[0][2]+sh1s[0][i]*Jinv[1][2]+sh1t[0][i]*Jinv[2][2];
    }//end for i



    // Jacobian
    for (int igp=0; igp<ngps; igp++) {
	Sh[0] = sh1[igp][0];
	Sh[1] = sh1[igp][1];
	Sh[2] = sh1[igp][2];
	Sh[3] = sh1[igp][3];

	double e11=0,e22=0,e33=0,e12=0,e23=0,e13=0;
	for(i=0;i<4;i++){
	    e11+= sh1[igp][i]*(((2.0/3.0/S0)*(-q->getValue(tt[i],0) /rt6 + q->getValue(tt[i],1)/rt2)+(1.0/3.0))*deleps + eper);	//~nx*nx
	    e22+= sh1[igp][i]*(((2.0/3.0/S0)*(-q->getValue(tt[i],0)/rt6 - q->getValue(tt[i],1)/rt2)+(1.0/3.0))*deleps + eper);	//~ny*ny
	    e33+= sh1[igp][i]*(((2.0/3.0/S0)*(2.0*q->getValue(tt[i],0)/rt6)	     +(1.0/3.0))*deleps + eper);			//~nz*nz
	    e12+= sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],2)/rt2)			*deleps;						//~nx*ny
	    e13+= sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],4)/rt2)			*deleps;						//~nx*nz
	    e23+= sh1[igp][i]*(2.0/3.0/S0)*(q->getValue(tt[i],3)/rt2)			*deleps;					//~ny*nz
	 }//end for i
	double mul=wsurf[igp]*eDet;
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
	SparseMatrix *K,
	double* L,
	Electrodes* electrodes){
    int it;
    //#pragma omp parallel for



    for (it=0; it< mesh->getnElements(); it++){
	double lK[npt][npt];
	double lL[npt];
	int tmat;
	int t[4] = {0,0,0,0};
	t[0] =  mesh->getNode(it,0);
	t[1] =  mesh->getNode(it,1);
	t[2] =  mesh->getNode(it,2);
	t[3] =  mesh->getNode(it,3);
	tmat = mesh->getMaterialNumber(it);

	localKL(p,t,lK,lL,it, mesh,q,lc,electrodes);
	//printlK( &lK[0][0] , npt);
        for (int i=0; i<npt; i++) // FOR ROWS
        {
	    int ri=v->getEquNode(t[i]);

            // RHS FIXED NODE HANDLING
            if ( ri == SolutionVector::FIXED_NODE ) // IF THIS NODE IS FIXED
            {
                for (int j = 0; j<4 ; j++)// SET CONTRIBUTION TO CONNECTED *FREE* NODES
                {
                    int nc = v->getEquNode( t[j] ); // INDEX TO CONNECTED NODE DEGREE OF FREEDOM POSITION
                    if (nc != SolutionVector::FIXED_NODE )
                    {
                        #pragma omp atomic
                        L[ nc ] -= lK[i][j]*v->getValue(t[i]); // L = -K*v
                    }
		}
                //L[ri] = 0; // if not dirichlet, RHS = 0: L = 0
            }// END IF ROW NODE IS FIXED

            if ( ri!=SolutionVector::FIXED_NODE )
                for (int j=0; j<npt  ; j++) // FOR COLUMNS
                {
                    int rj=v->getEquNode(t[j]);
                    if (rj != SolutionVector::FIXED_NODE )
                    {
                        K->sparse_add(rj,ri,lK[j][i]);
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
       SparseMatrix *K,
       double* L){
   int it = 0;
#pragma omp parallel for
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

	   localKL_N(p,&ti[0], lK , lL, it,index_to_Neumann, mesh, surf_mesh, q, lc);

	   for (int i=0; i<4; i++){
	       int ri=v->getEquNode(ti[i]);

               if (ri == SolutionVector::FIXED_NODE ) // HANDLE FIXED NODE
               {
		   for (int j = 0; j<4 ; j++)
                   {
                       int cr = v->getEquNode( ti[j] ); // CONNECTED NODE DOF ORDER
                       if (cr != SolutionVector::FIXED_NODE )
                       {
                            #pragma omp atomic
                            L[cr ] -= lK[j][i]*v->getValue(ti[i]);
                       }
                   }
               }// END HANDLE FIXED NODE
               else // HANDLE FREE NODE
               {
                   for (int j=0; j<4; j++) // FOR COLUMNS
                   {
                       int rj=v->getEquNode(ti[j]);

                       if (rj != SolutionVector::FIXED_NODE )
                       {
                           K->sparse_add( rj,ri,lK[j][i] );
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

void Pot_PCG(SparseMatrix *K, double *b, double *V, Settings* settings )
{
    int nnz = K->nnz;
    int size = K->rows;

       //for (int ee = 0 ; ee < size; ee++)
       //	printf("b[%i] =%f\n", ee, b[ee]*1e18);

       // Create SparseLib++ data structures
    CompCol_Mat_double A;
    A.point_to(size, nnz, K->P, K->I, K->J);


    //convert solution vector and RHS vector to SparseLib++
    VECTOR_double X;// = VECTOR_double(V,A.dim(0));// cannot use "point_to" with potential due to reoredring of values.
    X.point_to(V, size);
    VECTOR_double B;// = VECTOR_double(b,A.dim(0));
    B.point_to(b , size );
       // PCG settings...
    int return_flag =10;
    int maxiter =settings->getV_PCG_Maxiter();

    double toler = settings->getV_PCG_Toler();

    // Solves with different preconditioners...
    if (settings->getV_PCG_Preconditioner() == DIAG_PRECONDITIONER )
    {
        DiagPreconditioner_double D(A); // diagonal preconditioning, ~+3 times faster than cholesky
       return_flag = CG(A,X,B,D,maxiter,toler);
    }
    else if (settings->getV_PCG_Preconditioner() == IC_PRECONDITIONER )
    {
        ICPreconditioner_double D(A);
        return_flag = CG(A,X,B,D,maxiter,toler);
    }
    else if (settings->getV_PCG_Preconditioner() == ILU_PRECONDITIONER )
    {
       CompCol_ILUPreconditioner_double D(A); // compressed column format ILU
       return_flag = CG(A,X,B,D,maxiter,toler);
    }

    if (return_flag == 1) // if no convergence, print warni
        printf("PCG did not converge in %i iterations \nTolerance achieved is %f\n",maxiter,toler);


       //copy solution back to solution vector
       //for (int i = 0; i < sv->getnDoF() ; i++)
//	       sv->setValue(i , 0 , -X(sv->getEquNode(i)));

}
void Pot_GMRES(SparseMatrix *K, double *b, double* V, Settings* settings ){
/*! Solves the Linear simulatenous equation Ax=b using the GMRES method*/

    int nnz = K->nnz;
    int size = K->rows;
    // Create SparseLib++ data structures
    CompCol_Mat_double A;
    A.point_to(size , nnz, K->P, K->I, K->J);

   //Convert solution vector and RHS vector to SparseLib++
    VECTOR_double X; //=  VECTOR_double(sv->Values,A.dim(0)); // cannot use "point_to" with potential values because reoredring is done after calculation
    VECTOR_double B;
    X.point_to( V , size );
    B.point_to( b , size );

    // GMRES settings...
    int return_flag =10;
    int maxiter 	= settings->getV_GMRES_Maxiter();
    int restart 	= settings->getV_GMRES_Restart();
    double toler 	= settings->getV_GMRES_Toler();

    MATRIX_double H(restart+1, restart, 0.0);	// storage for upper Hessenberg H

    // Solves with different preconditioners...
    if (settings->getV_GMRES_Preconditioner() == DIAG_PRECONDITIONER )
    {
        DiagPreconditioner_double D(A);
        return_flag = GMRES(A,X,B,D,H,restart,maxiter,toler);
    }
    else if (settings->getV_GMRES_Preconditioner() == IC_PRECONDITIONER )
    {
        ICPreconditioner_double D(A);
        return_flag = GMRES(A,X,B,D,H,restart,maxiter,toler);
    }
    else if (settings->getV_GMRES_Preconditioner() == ILU_PRECONDITIONER )
    {
        CompCol_ILUPreconditioner_double D(A); // compressed column format ILU
        return_flag = GMRES(A,X,B,D,H,restart,maxiter,toler);
    }

    if (return_flag == 1)
        printf("GMRES did not converge in %i iterations \nTolerance achieved is %f\n",maxiter,toler);


    //copy solution back to solution vector
//#pragma omp parallel for
//    for (int i = 0; i < sv->getnDoF() ; i++)
//        sv->setValue(i , 0 , -X(sv->getEquNode(i)) );

}

