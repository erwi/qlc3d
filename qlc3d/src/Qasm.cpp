//#include "mex.h"
#include <time.h>
#include <cmath>
#include <cstring>

#include "sparsematrix.h"

struct spm{				// define sparse matrix link
		int row;
		spm * next;
		spm * prev;
	};
//#include <algo.h>
////////////////////////////////////////////////////////////////////////////////
// Global Functions
//mxArray* sparse_create(int *t, int *tmat,int nt,int npLC, int np,int ndof, int &nnz);

mxArray* sparsecreatell(int *t, int *tmat, int rt, int ct,int dof_per_node,int *equ_n, int *nnz,int npLC);
void assemble(	double *p,int *t, double *volume, int *tmat, double *q0,double *u0,double *qn, double *qd,double *v, int nt,int npLC, int np,
				double *Pr,int *Ir,int *Jc,double *L, double *params,int *per);

void pl_assemble(double *e,int e_length ,double *p,int npLC, double *q, double *qn, double *qd, int *pl_surface_index,
			  int pl_n_surfaces, double *pl_W, double *params, double *L, double *Pr, int *Ir, int *Jc,int *per, double *snorm);

void wk_assemble(double *e,int e_length,double *p,int npLC,double *q0, double* qn,double *qd,
				 int *wk_surface_index,int wk_n_surfaces,double *wk_vector1,double *wk_vector2,
				 double *params, double *L, double *Pr,int *Ir, int *Jc,int *per);

void assemble_surfaces(double *p, int *t,int *tmat, double *e, double *volume,double *v, int nt,
						int np, double *Pr, int *Ir, int *Jc, double *L, double *params, int npLC,
						int *e_to_t,int n_e, int *per);
				 
void sort(int *ps,int*pe);
void FlowReOrientation(double K[20][20],double L[20],double *q0,double *u0,int *tt, double sh[4],double dsh[4][3],double mul,int npLC, double *params);
void sp_elem(spm *col, int r, int c, int *nnz);
////////////////////////////////////////////////////////////////////////////////
// MAIN mexFunction
void Qasm()
{

// ___________Get data from Matlab__________________________
	double *p=mxGetPr(prhs[0]);
	int np=mxGetN(prhs[0]);
	//double *perm=(double*)mxGetPr(prhs[1]);
	double *q0	=mxGetPr(prhs[1]);
	double *qn	=mxGetPr(prhs[2]);
	double *qd	=mxGetPr(prhs[3]);
	double *v	=mxGetPr(prhs[4]);
	
	int *t=(int*)mxGetPr(prhs[5]);
	int nt=mxGetN(prhs[5]);
	int *tmat=(int*)mxGetPr(prhs[6]);
	
	double *volume = mxGetPr(prhs[7]);

	int ndof=(int)mxGetScalar(prhs[8]);
	//int ndof=5*npLC;//+npLC;	//q1->q5
	
	double *params=mxGetPr(prhs[9]);
	double *e=mxGetPr(prhs[10]);
	int e_length=mxGetN(prhs[10]);
//____________ANCHORING DATA_____________________
//planar degenerate anchoring data
	//index to planar surfaces
	int *pl_surface_index=(int*)mxGetPr(prhs[11]); //index to pl surfcs
	int pl_n_surfaces=mxGetM(prhs[11]);				//number of pl surfcs
	double *pl_W=mxGetPr(prhs[12]);					//pl ancoring strength
	


//weak surface anchoring data	
	int *wk_surface_index=(int*)mxGetPr(prhs[13]);	//index to weak surfaces
	int wk_n_surfaces=mxGetM(prhs[13]);				//number of weak surfaces
	double *wk_vector1=mxGetPr(prhs[14]);			//[strength,relativestrength,v1x,v1y,v1z;....]
	double *wk_vector2=mxGetPr(prhs[15]);			//[strength,relativestrength,v2x,v2y,v2z;....]
//___________________________________________


	int *per=(int*)mxGetPr(prhs[16]);
	int npLC=(int)mxGetScalar(prhs[17]);
    double *snorm = mxGetPr(prhs[18]);
	double *u0 = mxGetPr(prhs[19]);
	
	int *e_to_t =(int*) mxGetPr(prhs[20]);  // links between surfaces and tets
	
	
//______________________________________________________________

	// Create sparse matrix K
	int nnz=0;
	plhs[0]= sparsecreatell(t,tmat,4,nt,5,per,&nnz,npLC);
	
	double *Pr	= mxGetPr(plhs[0]);
	int *Ir		= mxGetIr(plhs[0]);
	int *Jc		= mxGetJc(plhs[0]);
	//	sparsecreatell(t,)
	int k,nnznew;
	// Create residual vector L
	plhs[1]=mxCreateDoubleMatrix(ndof,1,mxREAL);
	double *L=mxGetPr(plhs[1]);
	
	assemble(p,t,volume,tmat,q0,u0,qn,qd,v,nt,npLC,np,Pr,Ir,Jc,L,params,per);
	//assemble_surfaces(p,t,tmat,e,volume,v,nt,np,Pr,Ir,Jc,L,params,npLC,e_to_t,e_length,per);
	
	//if (pl_n_surfaces>0)
	//	pl_assemble(e,e_length,p,npLC,q0,qn,qd,pl_surface_index,pl_n_surfaces,pl_W,params,L,Pr,Ir,Jc,per,snorm);
	//if (wk_n_surfaces>0)
	//	wk_assemble(e,e_length,p,npLC,q0,qn,qd,wk_surface_index,wk_n_surfaces,wk_vector1,wk_vector2,params,L,Pr,Ir,Jc,per);
	
	//free the zero entries in the sparse array
	nnznew=0;
	for(k=0; k<nnz; k++){
			if(Pr[k]!=0.0) 	nnznew=nnznew+1;
	}
	


}
////////////////////////////////////////////////////////////////////////////////
// Gauss integration

// ---------------------------------------------------------
//     3D Gauss-Legendre weights for N = 11, D = 4
// ---------------------------------------------------------

const int ngp=11;
const double a=(1+sqrt(5.0/14.0))/4.0;
const double b=(1-sqrt(5.0/14.0))/4.0;

static double gp[ngp][4]={
	0.25	  , 0.25	,	0.25	,0.25,
	11.0/14.0     ,	1.0/14.0	,	1.0/14.0	,1.0/14.0,
	1.0/14.0      ,	11.0/14.0	,	1.0/14.0	,1.0/14.0,
	1.0/14.0	  ,	1.0/14.0	,	11.0/14.0   ,1.0/14.0,	
	1.0/14.0	  , 1.0/14.0	,	1.0/14.0	,11.0/14.0,
	a		  ,	a		,	b       ,b,
	a		  , b		,   a       ,b,
	a         , b       ,   b       ,a,
	b		  , a       ,   a       ,b,
	b		  , a       ,   b		,a,
	b		  , b		,   a		,a};

const double w11 = -74.0/5625.0;
const double w12 = 343.0/45000.0;
const double w13 = 56.0/2250.0;
static double w[ngp]={w11,w12,w12,w12,w12,w13,w13,w13,w13,w13,w13};

//-------------------------------------------------------
//	3D Gauss-Legendre weights for N = 1, D = 1
//-------------------------------------------------------
const int ngp2 = 1;
static double gp2[ngp2][3] = {0.25 , 0.25 , 0.25};
static double w2[ngp2] = {1.0 / 6.0};
// ---------------------------------------------------------
//     2D Gauss-Legendre weights for N = 6, D = 4
// ---------------------------------------------------------
const int sngp=6;
static double sgp[sngp][2] ={
		0.8168476, 0.09157621,
		0.09157621,0.8168476,
        0.09157621,0.09157621,
        0.1081030, 0.4459485,
        0.4459485, 0.1081030,
        0.4459485, 0.4459485};
static double sw[sngp]={ 0.05497587, 0.05497587, 0.05497587, 0.1116908, 0.1116908, 0.1116908};

////////////////////////////////////////////////////////////////////////////////
// Global data
static double ssh1[sngp][3];	//SURFACE term P1 shape function

static double sh1[ngp][4]; // P1 Shape functions
static double sh1r[ngp][4]; // P1 Shape functions r-derivatives
static double sh1s[ngp][4]; // P1 Shape functions s-derivatives
static double sh1t[ngp][4]; //P1 shape functions t-derivative

//static double sh2[ngp2][4]; // low order shape function
//static double sh2r[ngp2][4];
//static double sh2s[ngp2][4];
//static double sh2t[ngp2][4];

//GLOBALLY USED VALUES
static double rt2 = sqrt(2.0);
static double rt3 = sqrt(3.0);
static double rt6 = sqrt(6.0);


////////////////////////////////////////////////////////////////////////////////
// init_shape
void init_shape()
{
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
	}


}
 void init_shape_surf(){
	for (int i=0; i<sngp;i++){
	ssh1[i][0]=1-sgp[i][0]-sgp[i][1];
	ssh1[i][1]=sgp[i][0];
	ssh1[i][2]=sgp[i][1];
	}
}//end void init_shape
void init_shape_N()
{
memset(sh1,0,ngp*4*sizeof(double));
memset(sh1r,0,ngp*4*sizeof(double));
memset(sh1s,0,ngp*4*sizeof(double));
memset(sh1t,0,ngp*4*sizeof(double));
for (int i=0; i<sngp; i++) { // use surface shape functions
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
//
//------------------------------------------------------
//
mxArray* sparsecreatell(int *t, int *tmat, int rt, int ct, int dof_per_node,int *equ_n,int *nnz, int npLC)
{
	int maxp	= 0;// number of columns 
	int i=0,j=0,n1=0,n2=0,x=0; // counters & temporary indicies
	int *tt;		// temporary element

// find matrix size, taking into account periodic nodes
// determine matrix size -find largest node number
	
	for (i=0;i<ct;i++) //loop columns
	{
		if (tmat[i]==4)
		{
			for(j=0;j<rt;j++) //loop rows
			{
			if ((equ_n[t[i*rt+j]-1]+1)>maxp) maxp=equ_n[t[i*rt+j]-1]+1;
			}
		}
	}
//mexPrintf("maxp = %i\n",maxp);
// allocate memory for columns pointer array
	spm *col =(spm*) mxMalloc(dof_per_node*maxp*sizeof(spm)); // columns array
	for (i=0;i<maxp*dof_per_node;i++)  
	{
		col[i].next = NULL;  // initialize with no elements 
		col[i].prev = NULL;
		col[i].row  = 0;
	}


// fill in first row of blocks... column blocks are added later
	int r,c;
	tt=t;
	for (i=0;i<ct;i++,tt+=rt)
	{
		if(tmat[i]==4) // if LC tet
		{
		//mexPrintf("element %i\n",i);
		for (x=0;x<dof_per_node;x++){
		for (n1=0;n1<rt;n1++)
		{
			r = equ_n[tt[n1]-1] + x*maxp + 1;
		for (n2=n1+1;n2<rt;n2++)
		{
			c = equ_n[tt[n2]-1] + x*maxp + 1;
			sp_elem(col,1+equ_n[tt[n1]-1],r,nnz); // add entries to linked lists
			sp_elem(col,1+equ_n[tt[n1]-1],c,nnz); // row-column pairs	
			sp_elem(col,1+equ_n[tt[n2]-1],r,nnz);  
			sp_elem(col,1+equ_n[tt[n2]-1],c,nnz);
		}//end for n2
		}//end for n1
		}//end  for x
		}
	}//end for i

// extend linked lists, column by column to make matrix square
	if (dof_per_node>1) // only needed if multiple variables/node
	{
		spm *head, *tail, *temp;
		int rows=0;
		nnz[0] = dof_per_node*nnz[0]; // must increase number of nonzeros
		for (i=0;i<maxp*dof_per_node;i++)
		{
		head = &col[i];
		if (head->next) // if column exists
		{
			head = head->next; // first valid node
			tail = head;
			rows=1;
			while (tail->next!=NULL) // find end of list 
			{
				tail=tail->next; 
				rows++;
			}
			
			for (j=1;j<dof_per_node;j++) // create of_per_node copies of existing LL
			{
				temp = head;
				for (x=0;x<rows;x++)	// make a copy of LL to the end 
				{
					tail->next = (spm*)mxMalloc(sizeof(spm));	
					tail->next->prev = tail;				
					tail=tail->next;						
					tail->next = NULL;
					tail->row = (temp->row) + j*maxp; 
					temp=temp->next;
				}//end for x
			}//end for j
	}//end if column exists
	}//end for i
	}//end if 

// ALLOCATE MATLAB SPARSE MATRIX
	mxArray *K=mxCreateSparse(dof_per_node*maxp,dof_per_node*maxp,nnz[0],mxREAL);
	int *Ir=mxGetIr(K);
	int *Jc=mxGetJc(K);
	double *Pr = mxGetPr(K);
	spm *node;
	spm *temp;

// FILL IN CONNECTIVITY INFORMATION TO Ir AND Jc ARRAYS FROM LINKED LISTS
// Linked list is freed as it is copied into arrays
	int ir = 0;
	int jc = 0;
	nnz[0]=0; // counts number of entries

	for (i=0;i<maxp*dof_per_node;i++)// loop through col[i];
	{
		node= &col[i];
		node = node->next; // first node in list
		Jc[jc]=nnz[0];		   // column index	
		while (node!=NULL) // follow list untill NULL
		{
			Ir[ir] = (node->row) -1 ;	     // set row number in Ir 	
			nnz[0]++;
			ir++; //next position in Ir array

			//Cleanup, free linked list on the go...
			temp=node;
			node = node->next;
			if (!temp) mxFree(temp); 
					
		}//end while
		
		if (!node) mxFree(node);		  // free linked list, so this doesnt need to be done later	
		jc++; // next position in Jc array
	}// end for i
	Jc[jc]=nnz[0]; // last value in Jc is nnz
	mxFree(col);

	for (i=0;i<nnz[0];i++) Pr[i]=0.0;
	return K;

}//end void sprsecreatell


// void sp_elem inserts nodes into linked lists acording to row and column indicies
void sp_elem(spm *col, int r, int c, int *nnz)
{
//	mexPrintf("[%i,%i]",r,c);
	spm *node = &col[c-1];
	
	if (node->next == NULL)//first link in list
	{
		node->next = (spm*)mxMalloc(sizeof(spm));
		node->next->prev = node;
		node->next->row = r;
		node->next->next = NULL;
		node=node->next;
		nnz[0] ++;
	}
	else{
	node = node->next; // new db
	while ((node->row <r) && (node->next!=NULL))
	{
		node = node->next;
	}

	if (node->row==r)
	{
		
	//	mexPrintf("[%i,%i] -> do nothing\n",r,c);
	}
	 //do nothing 
	
	else if (node->row>r) // must insert new node in middle of list
	{
	//	mexPrintf("[%i,%i] -> insert in between\n",r,c);
		//spm *temp = new spm;		// new temp node
		spm *temp = (spm*)mxMalloc(sizeof(spm));

		node->prev->next = temp;	// insert temp node in between list << insertion fails
		temp->next=node;			//
		temp->prev=node->prev;
		node->prev=temp;			//
		temp->row = r;				// set row value
		nnz[0]++;						// add one to nonzero entries
	}




	else if(node->next==NULL) // found end of column, add new at end
	{
	//	mexPrintf("[%i,%i] -> add to end\n",r,c);
		//node->next = new spm;
		node->next=(spm*)mxMalloc(sizeof(spm));
		node->next->prev = node;	// link backwards
		node->next->next = NULL;	// link forwards = NULL, because end of list
		node->next->row = r;		// set row value
		nnz[0]++;					// add one to nonzero entries
	}
	
	else
	{mexPrintf("mysparsecreate.cpp -- !! unknown !!\n");}
	}//else not first link
}


 ////////////////////////////////////////////////////////////////////////////////
// localKL
void localKL(double *p,int *tt,double *volume,double *q0,double *u0, double *qn, double *qd,double *v, int npLC, int np,
			 double lK[20][20],double lL[20],double *params,int element_num)
{
	double lI[20][20];
	double e0=8.85e-12;
	double pi=3.14159265358979;
	double L1,L2,L6,dt,A,B,C,k11,k22,k33,deleps,u1,S0,P,qc,efe;
	bool bIsFlow = (params[11]==1);
	dt=params[0];
	A=params[1];
	B=params[2];
	C=params[3];
	k11=params[4];
	k22=params[5];
	k33=params[6];
	deleps=params[7];
	u1=params[9];
	S0=params[12];
	P=params[13];
	double e11=params[14]; //flexoelectric coeff.
	double e33=params[15];
	qc=0;
	int i,j,k;
	efe  = (2.0/S0/3.0)*(e11+2*e33);
	double efe2 = (4.0/S0/9.0)*(e11 - e33);

	deleps = deleps/S0;
	
	if(P!=0) qc=8.0*pi*k22/(9.0*S0*S0*P); //chirality
	//mexPrintf("qc = %f\n",qc);

	L1=2.0*(k33-k11+3.0*k22)/(S0*S0*27.0);
	L2=4.0*(k11-k22)/(9.0*S0*S0);
	L6=4.0*(k33-k11)/(S0*S0*S0*27.0);
	bool three_elastic_constant =false;
//mexPrintf("L1 %f, deleps %f\n",L1*1e12,deleps);
	if ((k11!=k22)||(k11!=k33)||(k22!=k33)) three_elastic_constant = true;
	//if(P!=0) qc=-2*pi*L1*4/P;

	memset(lK,0,20*20*sizeof(double)); 
	memset(lI,0,20*20*sizeof(double));
	memset(lL,0,20*sizeof(double));
double T1,T2,T3,T4,T5;//thermotropic 
double L1_1,L1_2,L1_3,L1_4,L1_5;//L1 elastic terms
double L2_1,L2_2,L2_3,L2_4,L2_5;//L2 elastic terms
double L6_1,L6_2,L6_3,L6_4,L6_5;//L6 elastic terms
double V1,V2,V3,V4,V5;			//Electric energy terms

double T11,T12,T13,T14,T15,T22,T23,T24,T25,T33,T34,T35,T44,T45,T55;
double L2_K11,L2_K12,L2_K13,L2_K14,L2_K15,L2_K22,L2_K23,L2_K24,L2_K25,
		L2_K33,L2_K34,L2_K35,L2_K44,L2_K45,L2_K55;
double L6_K11,L6_K12,L6_K13,L6_K14,L6_K15,L6_K22,L6_K23,L6_K24,L6_K25,
		L6_K33,L6_K34,L6_K35,L6_K44,L6_K45,L6_K55;	
double Lefe_1, Lefe_2, Lefe_3, Lefe_4, Lefe_5;	//Flexoelectric energy residual terms
double Kefe11=0, Kefe12=0, Kefe13=0, Kefe14=0, Kefe15=0, Kefe22=0, Kefe23=0, Kefe24=0, Kefe25=0,
	   Kefe33=0, Kefe34=0, Kefe35=0, Kefe44=0, Kefe45=0, Kefe55=0;
double Lc[5] = {0,0,0,0,0};
double Kc[5][5] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};

for (int igp=0;igp<ngp;igp++) {

//1. Calculate Inverse Jacobian
	double xr,xs,xt,yr,ys,yt,zr,zs,zt;
	xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;
		
		for (i=0;i<4;i++) {
			xr+=sh1r[igp][i]*p[(tt[i]-1)*3+0]*1e-6;
			xs+=sh1s[igp][i]*p[(tt[i]-1)*3+0]*1e-6;
			xt+=sh1t[igp][i]*p[(tt[i]-1)*3+0]*1e-6;

			yr+=sh1r[igp][i]*p[(tt[i]-1)*3+1]*1e-6;
			ys+=sh1s[igp][i]*p[(tt[i]-1)*3+1]*1e-6;
			yt+=sh1t[igp][i]*p[(tt[i]-1)*3+1]*1e-6;

			zr+=sh1r[igp][i]*p[(tt[i]-1)*3+2]*1e-6;
			zs+=sh1s[igp][i]*p[(tt[i]-1)*3+2]*1e-6;
			zt+=sh1t[igp][i]*p[(tt[i]-1)*3+2]*1e-6;
		}//end for i
		//Inverse Jacobian
		double Jdet = volume[0];
		double Jinv[3][3]={(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet
						  ,(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet
						  ,(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet};
//2. Shape function derivatives
		
		double dSh[4][3]; 
		double Sh[4];
		for(i=0;i<4;i++){ 
			Sh[i] = sh1[igp][i];
			dSh[i][0]=sh1r[igp][i]*Jinv[0][0]+sh1s[igp][i]*Jinv[1][0]+sh1t[igp][i]*Jinv[2][0];  
			dSh[i][1]=sh1r[igp][i]*Jinv[0][1]+sh1s[igp][i]*Jinv[1][1]+sh1t[igp][i]*Jinv[2][1];
			dSh[i][2]=sh1r[igp][i]*Jinv[0][2]+sh1s[igp][i]*Jinv[1][2]+sh1t[igp][i]*Jinv[2][2];
		}//end for i
//3. Function variables and derivatives
		//Construct variable/shape function integrals
		double q1=0, q2=0, q3=0, q4=0, q5=0;
		double q1x=0,q2x=0,q3x=0,q4x=0,q5x=0;
		double q1y=0,q2y=0,q3y=0,q4y=0,q5y=0;
		double q1z=0,q2z=0,q3z=0,q4z=0,q5z=0;
		double Vx=0,Vy=0,Vz=0;
			   
			
		for(i=0;i<4;i++){
			// Solution and derivatives
			q1+=Sh[i]*q0[tt[i]-1];				// q1i * Ni  = A1
			q2+=Sh[i]*q0[npLC+tt[i]-1];
			q3+=Sh[i]*q0[2*npLC+tt[i]-1];
			q4+=Sh[i]*q0[3*npLC+tt[i]-1];
			q5+=Sh[i]*q0[4*npLC+tt[i]-1];
	        //if (element_num==2)mexPrintf("q1-q5 = [%f,%f,%f,%f,%f]\n",q0[tt[i]-1],q0[npLC+tt[i]-1],q0[2*npLC+tt[i]-1],q0[3*npLC+tt[i]-1],q0[4*npLC+tt[i]-1]);
			// voltages
			Vx+=dSh[i][0]*v[tt[i]-1];
			Vy+=dSh[i][1]*v[tt[i]-1];
			Vz+=dSh[i][2]*v[tt[i]-1];
	        //mexPrintf(" node potential : %f \n",v[tt[i]-1]); 
			q1x+=dSh[i][0]*q0[tt[i]-1];
			q2x+=dSh[i][0]*q0[npLC+tt[i]-1];
			q3x+=dSh[i][0]*q0[npLC*2+tt[i]-1];
			q4x+=dSh[i][0]*q0[npLC*3+tt[i]-1];
			q5x+=dSh[i][0]*q0[npLC*4+tt[i]-1];
			
			q1y+=dSh[i][1]*q0[tt[i]-1];
			q2y+=dSh[i][1]*q0[npLC+tt[i]-1];
			q3y+=dSh[i][1]*q0[npLC*2+tt[i]-1];
			q4y+=dSh[i][1]*q0[npLC*3+tt[i]-1];
			q5y+=dSh[i][1]*q0[npLC*4+tt[i]-1];

			q1z+=dSh[i][2]*q0[tt[i]-1];
			q2z+=dSh[i][2]*q0[npLC+tt[i]-1];
			q3z+=dSh[i][2]*q0[npLC*2+tt[i]-1];
			q4z+=dSh[i][2]*q0[npLC*3+tt[i]-1];
			q5z+=dSh[i][2]*q0[npLC*4+tt[i]-1];
		}//end for i
		//mexPrintf("Vx %f, Vy %f, Vz %f\n",Vx*10,Vy*10,Vz*10);		
		double R=q1*q1+q2*q2+q3*q3+q5*q5+q4*q4;


//4. Calculate Residual and Stiffness matrix
	
	double mul=w[igp]*volume[0];
	
	

// ADD EFFECT OF FLOW	
	if (bIsFlow) FlowReOrientation(lK,lL,q0,u0,tt,Sh,dSh,mul,npLC,params);
	// i = local row index
	for (i=0;i<4;i++) {
		
		double ShRx=mul*dSh[i][0];//including weight and jacobian in trial function
		double ShRy=mul*dSh[i][1];
		double ShRz=mul*dSh[i][2];
		double ShR =mul*Sh[i];

		//Vi = Electric enery terms
		//V1= (deleps*e0*rt6/12.0*(Vx*Vx+Vy*Vy-Vz*Vz*2.0))*ShR;
		//V2= (-deleps*e0*rt2/4.0*(Vx*Vx-Vy*Vy))*ShR;
		//V3= (-deleps*e0*Vx*Vy*rt2/2.0)*ShR;
		//V4= (-deleps*e0*Vy*Vz*rt2/2.0)*ShR;
		//V5= (-deleps*e0*Vx*Vz*rt2/2.0)*ShR;
		
			V1 =  rt6*(Vx*Vx + Vy*Vy-2.0*Vz*Vz)*deleps/18.0*ShR*e0;
			V2 = -rt2*(Vx*Vx - Vy*Vy)*deleps/6.0*ShR*e0;
			V3 = -rt2*Vx*Vy*deleps/3.0*ShR*e0;
			V4 = -rt2*Vz*Vy*deleps/3.0*ShR*e0;
			V5 = -rt2*Vx*Vz*deleps/3.0*ShR*e0;
		//mexPrintf("deleps %f \n",deleps);
		//IF flexoelectricity considered		
	if ((efe!=0.0)||(efe2!=0.0)){
		Lefe_1 = (rt6*(Vx*ShRx+Vy*ShRy-2.0*Vz*ShRz)*efe/6.0);
		Lefe_2 = (-rt2*(Vx*ShRx-Vy*ShRy)*efe/2.0);
		Lefe_3 = (-rt2*(Vy*ShRx+Vx*ShRy)*efe/2.0);
		Lefe_4 = (-rt2*(Vz*ShRy+Vy*ShRz)*efe/2.0);
		Lefe_5 = (-rt2*(Vz*ShRx+Vx*ShRz)*efe/2.0);
		
		//double 	MapleGenVar1 = -ShRy*Vx*efe2*q3*rt2*rt6/12.0-ShR*Vz*efe2*rt6*q4y*rt2/6.0+ShR*Vx*efe2*rt6*q2x*rt2/12.0+ShR*Vx*efe2*q3y*rt2*rt6/12.0+ShR*Vx*efe2*q5z*rt2*rt6/12.0+ShR*Vy*efe2*q3x*rt2*rt6/12.0-ShR*Vy*efe2*rt6*q2y*rt2/12.0+ShR*Vy*efe2*q4z*rt2*rt6/12.0-ShR*Vz*efe2*rt6*q5x*rt2/6.0-ShRx*Vx*efe2*q2*rt2*rt6/12.0-ShRx*Vy*efe2*q3*rt2*rt6/12.0-ShRx*Vz*efe2*q5*rt2*rt6/12.0;
		//Lefe_1	= MapleGenVar1+ShRy*Vy*efe2*q2*rt2*rt6/12.0-ShRy*Vz*efe2*q4*rt2*rt6/12.0+ShRz*Vx*efe2*q5*rt2*rt6/6.0+ShRz*Vy*efe2*q4*rt2*rt6/6.0-efe*Vy*rt6*ShRy/6.0+efe*Vz*rt6*ShRz/3.0-ShR*Vx*efe2*q1x/6.0-ShR*Vy*efe2*q1y/6.0-2.0/3.0*ShR*Vz*efe2*q1z+ShRx*Vx*efe2*q1/6.0+ShRy*Vy*efe2*q1/6.0+2.0/3.0*ShRz*Vz*efe2*q1-efe*Vx*rt6*ShRx/6.0;
		//Lefe_2	= efe*Vx*rt2*ShRx/2.0-efe*Vy*rt2*ShRy/2.0-ShR*Vx*efe2*q2x/2.0-ShR*Vx*efe2*q3y/2.0-ShR*Vx*efe2*q5z/2.0+ShR*Vy*efe2*q3x/2.0-ShR*Vy*efe2*q2y/2.0+ShR*Vy*efe2*q4z/2.0+ShR*Vx*efe2*rt2*q1x*rt6/12.0-ShR*Vy*efe2*rt2*q1y*rt6/12.0+ShRx*Vx*efe2*q2/2.0+ShRx*Vy*efe2*q3/2.0+ShRx*Vz*efe2*q5/2.0-ShRx*Vx*efe2*q1*rt6*rt2/12.0-ShRy*Vx*efe2*q3/2.0+ShRy*Vy*efe2*q2/2.0-ShRy*Vz*efe2*q4/2.0+ShRy*Vy*efe2*q1*rt6*rt2/12.0;
		//Lefe_3  = efe*Vy*rt2*ShRx/2.0+efe*Vx*rt2*ShRy/2.0+ShR*Vx*efe2*q2y/2.0-ShR*Vy*efe2*q2x/2.0+ShR*Vx*efe2*rt2*q1y*rt6/12.0+ShR*Vy*efe2*rt2*q1x*rt6/12.0-ShR*Vx*efe2*q3x/2.0-ShR*Vx*efe2*q4z/2.0-ShR*Vy*efe2*q3y/2.0-ShR*Vy*efe2*q5z/2.0-ShRx*Vy*efe2*q2/2.0-ShRx*Vy*efe2*q1*rt6*rt2/12.0+ShRx*Vx*efe2*q3/2.0+ShRx*Vz*efe2*q4/2.0+ShRy*Vx*efe2*q2/2.0+ShRy*Vy*efe2*q3/2.0+ShRy*Vz*efe2*q5/2.0-ShRy*Vx*efe2*q1*rt6*rt2/12.0;
		//Lefe_4	= efe*Vz*rt2*ShRy/2.0+efe*Vy*rt2*ShRz/2.0+ShR*Vz*efe2*q2y/2.0-ShR*Vy*efe2*rt2*q1z*rt6/6.0+ShR*Vz*efe2*rt2*q1y*rt6/12.0-ShR*Vy*efe2*q5x/2.0-ShR*Vy*efe2*q4y/2.0-ShR*Vz*efe2*q3x/2.0-ShR*Vz*efe2*q4z/2.0+ShRy*Vz*efe2*q1*rt6*rt2/6.0+ShRy*Vx*efe2*q5/2.0+ShRy*Vy*efe2*q4/2.0-ShRz*Vy*efe2*q2/2.0-ShRz*Vy*efe2*q1*rt6*rt2/12.0+ShRz*Vx*efe2*q3/2.0+ShRz*Vz*efe2*q4/2.0;
		//Lefe_5	= efe*Vz*rt2*ShRx/2.0+efe*Vx*rt2*ShRz/2.0-	ShR*Vz*efe2*q2x/2.0-ShR*Vx*efe2*rt2*q1z*rt6/6.0+ShR*Vz*efe2*rt2*q1x*rt6/12.0-ShR*Vx*efe2*q5x/2.0-ShR*Vx*efe2*q4y/2.0-ShR*Vz*efe2*q3y/2.0-ShR*Vz*efe2*q5z/2.0+ShRx*Vz*efe2*q1*rt6*rt2/6.0+ShRx*Vx*efe2*q5/2.0+ShRx*Vy*efe2*q4/2.0+ShRz*Vx*efe2*q2/2.0+ShRz*Vy*efe2*q3/2.0+ShRz*Vz*efe2*q5/2.0-ShRz*Vx*efe2*q1*rt6*rt2/12.0;
	}
	else
		{
			Lefe_1=Lefe_2=Lefe_3=Lefe_4=Lefe_5=0.0;	
		}//end if flexoelectricity
	if (qc!=0)
	{
		memset(&Lc,0,5*sizeof(double));
		Lc[0] = -ShR*qc*rt3*q4x/2.0+ShR*qc*rt3*q5y/2.0-ShRx*qc*q4*rt3/2.0+ShRy*qc*q5*rt3/2.0;
		Lc[1] = ShR*qc*q3z-q5y*qc*ShR/2.0-q4x*qc*ShR/2.0-q4*qc*ShRx/2.0-q5*qc*ShRy/2.0+ShRz*qc*q3;
		Lc[2] = -ShR*qc*q4y/2.0+ShR*qc*q5x/2.0-ShR*qc*q2z+ShRx*qc*q5/2.0-ShRy*qc*q4/2.0-ShRz*qc*q2;
		Lc[3] = ShR*qc*q3y/2.0-ShR*qc*q5z/2.0+ShR*qc*q2x/2.0+ShR*qc*q1x*rt3/2.0+ShRx*qc*q2/2.0+ShRx*qc*q1*rt3/2.0+ShRy*qc*q3/2.0-ShRz*qc*q5/2.0;
		Lc[4] = ShR*qc*q4z/2.0-ShR*qc*q3x/2.0+ShR*qc*q2y/2.0-ShR*qc*q1y*rt3/2.0-ShRx*qc*q3/2.0+ShRy*qc*q2/2.0-ShRy*qc*q1*rt3/2.0+ShRz*qc*q4/2.0;
	}
	
	//L1_i = L1-Elastic terms
		L1_1=(ShRx*q1x+ShRy*q1y+ShRz*q1z)*L1;
		L1_2=(ShRx*q2x+ShRy*q2y+ShRz*q2z)*L1;
		L1_3=(ShRx*q3x+ShRy*q3y+ShRz*q3z)*L1;
		L1_4=(ShRx*q4x+ShRy*q4y+ShRz*q4z)*L1;
		L1_5=(ShRx*q5x+ShRy*q5y+ShRz*q5z)*L1;
	//L2_i = L2-Elastic terms
	if (three_elastic_constant){
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
				
	//Ti = thermotropic residual terms
		T1=ShR*(A*q1+B*(q5*q5*rt6/4.0-q3*q3*rt6/2.0-rt6*q2*q2/2.0+q1*q1*rt6/2.0+q4*q4*rt6/4.0)/3.0+C*R*q1);
		T2=ShR*(A*q2+B*(3.0/4.0*q5*q5*rt2-q1*rt6*q2-3.0/4.0*q4*q4*rt2)/3.0+C*R*q2);
		T3=ShR*(A*q3+B*(-q3*q1*rt6+3.0/2.0*rt2*q5*q4)/3.0+C*R*q3);
		T4=ShR*(A*q4+B*(3.0/2.0*q3*rt2*q5+q4*q1*rt6/2.0-3.0/2.0*q4*q2*rt2)/3.0+C*R*q4);
		T5=ShR*(A*q5+B*(q5*q1*rt6/2.0+3.0/2.0*q5*q2*rt2+3.0/2.0*q3*rt2*q4)/3.0+C*R*q5);

	
	
	
	
	//Add energy terms to residual vector
		//mexPrintf("L2_1, %f, L2_3 %f, L2_4 %f, L2_5 %f\n",L6_1*1e10,L6_2*1e10,L6_3*1e10,L6_4*1e10,L6_5*1e10);
		lL[i+0]  += L1_1 + L2_1 + L6_1 + V1 + Lefe_1 + T1 +Lc[0];
	    lL[i+4]  += L1_2 + L2_2 + L6_2 + V2 + Lefe_2 + T2 +Lc[1];
		lL[i+8]  += L1_3 + L2_3 + L6_3 + V3 + Lefe_3 + T3 +Lc[2];
		lL[i+12] += L1_4 + L2_4 + L6_4 + V4 + Lefe_4 + T4 +Lc[3];
		lL[i+16] += L1_5 + L2_5 + L6_5 + V5 + Lefe_5 + T5 +Lc[4];
			
				
	// j = local column index			
	for (int j=0;j<4;j++){
		double ShCx=dSh[j][0];
		double ShCy=dSh[j][1];
		double ShCz=dSh[j][2];	
		double ShC =Sh[j];
		double ShRC=ShR*Sh[j];

	// Calculate stiffnes matrix terms
	//Tij = thermotropic stiffness matrix terms
		T11=ShRC*(A+B*q1*rt6/3.0+2.0*C*q1*q1+C*R);
		T12=ShRC*(-B*rt6*q2/3.0+2.0*C*q2*q1);
		T13=ShRC*(-B*q3*rt6/3.0+2.0*C*q3*q1);
		T14=ShRC*(B*q4*rt6/6.0+2.0*C*q4*q1);
		T15=ShRC*(B*q5*rt6/6.0+2.0*C*q5*q1);

		T22=ShRC*(A-B*q1*rt6/3.0+2.0*C*q2*q2+C*R);
		T23=ShRC*(2.0*C*q3*q2);
		T24=ShRC*(-B*q4*rt2/2.0+2.0*C*q4*q2);
		T25=ShRC*(B*q5*rt2/2.0+2.0*C*q5*q2);

		T33=ShRC*(A-B*q1*rt6/3.0+2.0*C*q3*q3+C*R);
		T34=ShRC*(B*q5*rt2/2.0+2.0*C*q4*q3);
		T35=ShRC*(B*q4*rt2/2.0+2.0*C*q5*q3);

		T44=ShRC*(A+B*(q1*rt6-3.0*q2*rt2)/6.0+2.0*C*q4*q4+C*R);
		T45=ShRC*(B*q3*rt2/2.0+2.0*C*q5*q4);
		T55=ShRC*(A+B*(q1*rt6+3.0*q2*rt2)/6.0+2.0*C*q5*q5+C*R);
	//flexoelectric stiffnes matrix terms
	/*
		if((efe!=0.0)||(efe2!=0.0))
		{
			
			Kefe11 =ShC*ShRx*Vx*efe2/6.0+ShC*ShRy*Vy*efe2/6.0+2.0/3.0*ShC*ShRz*Vz*efe-ShR*Vx*efe2*ShCx/6.0-2.0/3.0*ShR*Vz*efe2*ShCz;
			Kefe12 =-ShC*ShRx*Vx*efe2*rt2*rt6/12.0+ShC*ShRy*Vy*efe2*rt2*rt6/12.0+ShR*Vx*efe2*rt6*rt2*ShCx/12.0;
			Kefe13 =-ShC*ShRy*Vx*efe2*rt2*rt6/12.0-	ShC*ShRx*Vy*efe2*rt2*rt6/12.0+ShR*Vy*efe2*rt2*rt6*ShCx/12.0;
			Kefe14 =-ShC*ShRy*Vz*efe2*rt2*rt6/12.0+ShC*ShRz*Vy*efe2*rt2*rt6/6.0+ShR*Vy*efe2*rt2*rt6*ShCz/12.0;
			Kefe15  =-ShC*ShRx*Vz*efe2*rt2*rt6/12.0+ShC*ShRz*Vx*efe2*rt2*rt6/6.0-ShR*Vz*efe2*rt6*rt2*ShCx/6.0+ShR*Vx*efe2*rt6*rt2*ShCz/12.0;

			Kefe22 = ShC*ShRx*Vx*efe2/2.0+ShC*ShRy*Vy*efe2/2.0-ShR*Vx*efe2*ShCx/2.0;
			Kefe23 = ShRx*Vy*efe2*ShC/2.0-ShRy*Vx*efe2*ShC/2.0+ShR*Vy*efe2*ShCx/2.0;
			Kefe24 = -ShRy*Vz*efe2*ShC/2.0+ShR*Vy*efe2*ShCz/2.0;
			Kefe25 = ShRx*Vz*efe2*ShC/2.0-ShR*Vx*efe2*ShCz/2.0;

			Kefe33 = ShC*ShRx*Vx*efe2/2.0+ShC*ShRy*Vy*efe2/2.0-	ShR*Vx*efe2*ShCx/2.0;
			Kefe34 = ShRx*Vz*efe2*ShC/2.0-ShR*Vx*efe2*ShCz/2.0;
			Kefe35 = ShRy*Vz*efe2*ShC/2.0-ShR*Vy*efe2*ShCz/2.0;

			Kefe44 = ShC*ShRy*Vy*efe2/2.0+ShC*ShRz*Vz*efe2/2.0-	ShR*Vz*efe2*ShCz/2.0;
			Kefe45 = ShRy*Vx*efe2*ShC/2.0-ShR*Vy*efe2*ShCx/2.0;

			Kefe55 = ShC*ShRx*Vx*efe2/2.0+ShC*ShRz*Vz*efe2/2.0-	ShR*Vx*efe2*ShCx/2.0-ShR*Vz*efe2*ShCz/2.0;
		}
*/
	if (qc!=0) // if chirality
		{
			memset(&Kc,0,5*5*sizeof(double));
			Kc[0][3] = -ShRx*qc*rt3*ShC/2.0-ShR*qc*rt3*ShCx/2.0;
			Kc[0][4] = ShRy*qc*rt3*ShC/2.0+ShR*qc*rt3*ShCy/2.0;
			Kc[1][2] = ShRz*qc*ShC+ShR*qc*ShCz;
			Kc[1][3] = -ShC*qc*ShRx/2.0-ShCx*qc*ShR/2.0;
			Kc[1][4] = -ShC*qc*ShRy/2.0-ShCy*qc*ShR/2.0;
			Kc[2][3] = -ShC*qc*ShRy/2.0-ShCy*qc*ShR/2.0;
			Kc[2][4] = ShC*qc*ShRx/2.0+ShCx*qc*ShR/2.0;
			Kc[3][4] = -ShRz*qc*ShC/2.0-ShR*qc*ShCz/2.0;
		}


		//L1- matrix terms
		double dot=L1*mul*(dSh[i][0]*dSh[j][0]+dSh[i][1]*dSh[j][1]+dSh[i][2]*dSh[j][2]);
		//L2-terms
		if (three_elastic_constant){//if three elastic constants used
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
				
				/////////////////////////////////////7777
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
				lK[i][j   ] += T11 +	dot	+L2_K11		+L6_K11;
				lK[i][j+4 ] += T12 +		L2_K12		+L6_K12;	
				lK[i][j+8 ] += T13 +		L2_K13		+L6_K13;
				lK[i][j+12] += T14 +		L2_K14		+L6_K14 +Kc[0][3];
				lK[i][j+16] += T15 +		L2_K15		+L6_K15 +Kc[0][4];


				lK[i+4][j+0 ] += T12 +		L2_K12		+L6_K12;
				lK[i+4][j+4 ] += T22 +	dot	+L2_K22	+L6_K22;
				lK[i+4][j+8 ] += T23 +	L2_K23 +L6_K23 +Kc[1][2];
				lK[i+4][j+12] += T24 +  L2_K24 +L6_K24 +Kc[1][3];
				lK[i+4][j+16] += T25 +	L2_K25 +L6_K25 +Kc[1][4];	

				
				lK[i+8][j+0] +=T13 +		L2_K13		+L6_K13;
				lK[i+8][j+4] +=T23 +	    L2_K23     +L6_K23;
				lK[i+8][j+8]  += T33 + dot	+L2_K33 +L6_K33 ;
				lK[i+8][j+12] += T34 + L2_K34 +L6_K34 +Kc[2][3];
				lK[i+8][j+16] += T35 + L2_K35 +L6_K35 +Kc[2][4];

			
				lK[i+12][j+0 ]+=T14 +		L2_K14		+L6_K14;
				lK[i+12][j+4 ]+=T24 +       L2_K24 +L6_K24;
				lK[i+12][j+8 ]+=T34 +       L2_K34 +L6_K34;
				lK[i+12][j+12]+=T44 +dot	+L2_K44 +L6_K44;
				lK[i+12][j+16]+=T45	+L2_K45 +L6_K45 +Kc[3][4];

				
				
				lK[i+16][j+0 ]+=T15 +		L2_K15		+L6_K15;
				lK[i+16][j+4] += T25 +	L2_K25 +L6_K25;	
				lK[i+16][j+8] += T35 + L2_K35 +L6_K35;
				lK[i+16][j+12]+=T45	+L2_K45 +L6_K45;
				lK[i+16][j+16]+=T55 +dot	+L2_K55 +L6_K55;
	
		
//Local identity matrix			
					lI[i   ][j   ]+=ShRC;
					lI[i+4 ][j+4 ]+=ShRC;
					lI[i+8 ][j+8 ]+=ShRC;
					lI[i+12][j+12]+=ShRC;
					lI[i+16][j+16]+=ShRC;

		}//end for j
	}//end for i
}//end for igp


	    
if(dt!=0)// Crank-Nicolson time stepping
	{// %0 calculates the steady state
		// the product of the mass matrix and the various q vectors
		double Mqi[20],Mqn[20],Mqd[20];
		memset(Mqi,0,20*sizeof(double));
		memset(Mqn,0,20*sizeof(double));
		memset(Mqd,0,20*sizeof(double));
		for (i=0;i<4;i++) {		//each node row
			for (j=0;j<4;j++) {	//each node column
				for (k=0;k<5;k++){//each component
					Mqi[4*k+i]+=lI[i][j]*q0[np*k+tt[j]-1];
					Mqn[4*k+i]+=lI[i][j]*qn[np*k+tt[j]-1];
					Mqd[4*k+i]+=lI[i][j]*qd[np*k+tt[j]-1];
				}
			}
		}
		for (i=0;i<20;i++) {
			for (j=0;j<20;j++) {
				 lK[i][j]=lI[i][j]*2.0*u1/dt+lK[i][j];
			}
			lL[i]=Mqn[i]*2.0*u1/dt+u1*Mqd[i]-(lL[i]+Mqi[i]*2.0*u1/dt);
		}
	}//if(dt!=0)
}//END void lokaKL
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////////Flow induced reorientation contributions to residual and stiffness matrix

void FlowReOrientation(double K[20][20],double L[20],double *q0,double *u0,int *tt, double sh[4],double dsh[4][3],double mul,int npLC, double *params)
{
	
	double mu1  =	params[9];
	double mu2	=	params[10];

	double u=0,ux=0, uy=0, uz=0;
	double v=0,vx=0, vy=0, vz=0;
	double w=0,wx=0, wy=0, wz=0;
	double q1x=0, q2x=0, q3x=0, q4x=0, q5x=0;
	double q1y=0, q2y=0, q3y=0, q4y=0, q5y=0;
	double q1z=0, q2z=0, q3z=0, q4z=0, q5z=0;
	double q1 =0, q2 =0, q3 =0, q4 =0, q5 =0;
	int i;

	for(i=0;i<4;i++)
	{
			int n1 = tt[i]-1; //shortcut to node numbers
			int n2 = tt[i]-1 +1*npLC;
			int n3 = tt[i]-1 +2*npLC;
			int n4 = tt[i]-1 +3*npLC;
			int n5 = tt[i]-1 +4*npLC;
			
			u  += sh[i]*u0[n1]; ux += dsh[i][0]*u0[n1];	uy += dsh[i][1]*u0[n1]; uz += dsh[i][2]*u0[n1];
			v  += sh[i]*u0[n2]; vx += dsh[i][0]*u0[n2];	vy += dsh[i][1]*u0[n2]; vz += dsh[i][2]*u0[n2];
			w  += sh[i]*u0[n3]; wx += dsh[i][0]*u0[n3];	wy += dsh[i][1]*u0[n3]; wz += dsh[i][2]*u0[n3];
			
			//use normal tensor basis for flow stresses 
			
			double a1 = -q0[n1]/rt6 +q0[n2]/rt2;
			double a2 = -q0[n1]/rt6 -q0[n2]/rt2;
			double a3 =  q0[n3]/rt2;
			double a4 =  q0[n4]/rt2;
			double a5 =  q0[n5]/rt2;
			q1  += a1*sh[i];	q1x +=a1*dsh[i][0];	q1y +=a1*dsh[i][1];	q1z +=a1*dsh[i][2];	
			q2  += a2*sh[i];	q2x +=a2*dsh[i][0];	q2y +=a2*dsh[i][1];	q2z +=a2*dsh[i][2];	
			q3  += a3*sh[i];	q3x +=a3*dsh[i][0];	q3y +=a3*dsh[i][1];	q3z +=a3*dsh[i][2];	
			q4  += a4*sh[i];	q4x +=a4*dsh[i][0];	q4y +=a4*dsh[i][1];	q4z +=a4*dsh[i][2];	
			q5  += a5*sh[i];	q5x +=a5*dsh[i][0];	q5y +=a5*dsh[i][1];	q5z +=a5*dsh[i][2];	
  
	
	
	}


	for(i=0;i<4;i++){
		double ShR = sh[i]*mul;
		
      L[i+0 ] += (-1/rt6*(0.5*mu2*ux+mu1*(u*q1x+v*q1y+(0.5*vx-0.5*uy)*q3-q3*(0.5*uy-0.5*vx)+w*q1z+(0.5*wx-0.5*uz)*q5-q5*(0.5*uz-0.5*wx)))-1/rt6*(0.5*mu2*vy+mu1*(u*q2x+q3*(0.5*uy-0.5*vx)-(0.5*vx-0.5*uy)*q3+v*q2y+w*q2z+(0.5*wy-0.5*vz)*q4-q4*(0.5*vz-0.5*wy)))+2.0/rt6*(0.5*mu2*wz+mu1*(u*(-q1x-q2x)+q5*(0.5*uz-0.5*wx)-(0.5*wx-0.5*uz)*q5+v*(-q1y-q2y)+q4*(0.5*vz-0.5*wy)-(0.5*wy-0.5*vz)*q4+w*(-q1z-q2z))))*ShR;
      L[i+4 ] += (1/rt2*(0.5*mu2*ux+mu1*(u*q1x+v*q1y+(0.5*vx-0.5*uy)*q3-q3*(0.5*uy-0.5*vx)+w*q1z+(0.5*wx-0.5*uz)*q5-q5*(0.5*uz-0.5*wx)))-1/rt2*(0.5*mu2*vy+mu1*(u*q2x+q3*(0.5*uy-0.5*vx)-(0.5*vx-0.5*uy)*q3+v*q2y+w*q2z+(0.5*wy-0.5*vz)*q4-q4*(0.5*vz-0.5*wy))))*ShR;
      L[i+8 ] += (1/rt2*(0.5*mu2*(0.5*vx+0.5*uy)+mu1*(u*q3x+(0.5*uy-0.5*vx)*q1+v*q3y-q2*(0.5*uy-0.5*vx)+w*q3z+(0.5*wy-0.5*vz)*q5-q4*(0.5*uz-0.5*wx)))+1/rt2*(0.5*mu2*(0.5*vx+0.5*uy)+mu1*(u*q3x-q1*(0.5*vx-0.5*uy)+v*q3y+(0.5*vx-0.5*uy)*q2+w*q3z+(0.5*wx-0.5*uz)*q4-q5*(0.5*vz-0.5*wy))))*ShR;
      L[i+12] += (1/rt2*(0.5*mu2*(0.5*wy+0.5*vz)+mu1*(u*q4x+(0.5*uz-0.5*wx)*q3-q5*(0.5*vx-0.5*uy)+v*q4y+(0.5*vz-0.5*wy)*q2+w*q4z-(-q1-q2)*(0.5*vz-0.5*wy)))+1/rt2*(0.5*mu2*(0.5*wy+0.5*vz)+mu1*(u*q4x+(0.5*uy-0.5*vx)*q5-q3*(0.5*wx-0.5*uz)+v*q4y-q2*(0.5*wy-0.5*vz)+w*q4z+(0.5*wy-0.5*vz)*(-q1-q2))))*ShR;
      L[i+16] += (1/rt2*(0.5*mu2*(0.5*wx+0.5*uz)+mu1*(u*q5x+(0.5*uz-0.5*wx)*q1+v*q5y+(0.5*vz-0.5*wy)*q3-q4*(0.5*uy-0.5*vx)+w*q5z-(-q1-q2)*(0.5*uz-0.5*wx)))+1/rt2*(0.5*mu2*(0.5*wx+0.5*uz)+mu1*(u*q5x-q1*(0.5*wx-0.5*uz)+v*q5y+(0.5*vx-0.5*uy)*q4-q3*(0.5*wy-0.5*vz)+w*q5z+(0.5*wx-0.5*uz)*(-q1-q2))))*ShR;

      



	for(int j=0;j<4;j++)
	{
		double ShC= sh[j];
		double ShCx = dsh[j][0];
		double ShCy = dsh[j][1];
		double ShCz = dsh[j][2];
		
	  K[0+i][j+0] += (-3.0/rt6*mu1*v*ShCy+(-3.0/rt6*u*ShCx-3.0/rt6*w*ShCz)*mu1)*ShR;
      K[0+i][j+4] += (-3.0/rt6*mu1*v*ShCy+(-3.0/rt6*u*ShCx-3.0/rt6*w*ShCz)*mu1)*ShR;
      K[0+i][j+8] += (-1.0/rt6*(0.1E1*vx-0.1E1*uy)-1.0/rt6*(0.1E1*uy-0.1E1*vx))*ShC*mu1*ShR;
      K[0+i][j+12] += (-1.0/rt6*(0.1E1*wy-0.1E1*vz)+2.0/rt6*(0.1E1*vz-0.1E1*wy))*ShC*mu1*ShR;
      K[0+i][j+16] += (-1.0/rt6*(0.1E1*wx-0.1E1*uz)+2.0/rt6*(0.1E1*uz-0.1E1*wx))*ShC*mu1*ShR;
      
      K[4+i][j+0] += (1/rt2*mu1*v*ShCy+(u*ShCx+w*ShCz)/rt2*mu1)*ShR;
      K[4+i][j+4] += (-1/rt2*mu1*v*ShCy+(-u*ShCx-w*ShCz)/rt2*mu1)*ShR;
      K[4+i][j+8] += (0.2E1*vx-0.2E1*uy)*ShC/rt2*mu1*ShR;
      K[4+i][j+12] += -1.0/rt2*mu1*(0.1E1*wy-0.1E1*vz)*ShR*ShC;
      K[4+i][j+16] += 1/rt2*mu1*(0.1E1*wx-0.1E1*uz)*ShR*ShC;
      
      K[8+i][j+0] += 2.0/rt2*mu1*(0.5*uy-0.5*vx)*ShR*ShC;
      K[8+i][j+4] += 2.0/rt2*mu1*(0.5*vx-0.5*uy)*ShR*ShC;
      K[8+i][j+8] += (2.0/rt2*mu1*v*ShCy+(2.0*u*ShCx+2.0*w*ShCz)/rt2*mu1)*ShR;
      K[8+i][j+12] += 2.0/rt2*mu1*(0.5*wx-0.5*uz)*ShR*ShC;
      K[8+i][j+16] += 2.0/rt2*mu1*(0.5*wy-0.5*vz)*ShR*ShC;
      
      K[12+i][j+0] += 2.0/rt2*mu1*(0.5*vz-0.5*wy)*ShR*ShC;
      K[12+i][j+4] += 2.0/rt2*mu1*(0.1E1*vz-0.1E1*wy)*ShR*ShC;
      K[12+i][j+8] += 2.0/rt2*mu1*(0.5*uz-0.5*wx)*ShR*ShC;
      K[12+i][j+12] += (2.0/rt2*mu1*v*ShCy+(2.0*u*ShCx+2.0*w*ShCz)/rt2*mu1)*ShR;
      K[12+i][j+16] += 2.0/rt2*mu1*(0.5*uy-0.5*vx)*ShR*ShC;

      K[16+i][j+0] += 2.0/rt2*mu1*(0.1E1*uz-0.1E1*wx)*ShR*ShC;
      K[16+i][j+4] += 2.0/rt2*mu1*(0.5*uz-0.5*wx)*ShR*ShC;
      K[16+i][j+8] += 2.0/rt2*mu1*(0.5*vz-0.5*wy)*ShR*ShC;
      K[16+i][j+12] += 2.0/rt2*mu1*(0.5*vx-0.5*uy)*ShR*ShC;
      K[16+i][j+16] += (2.0/rt2*mu1*v*ShCy+(2.0*u*ShCx+2.0*w*ShCz)/rt2*mu1)*ShR;









	}//end for j
	}//end for i



}

/////////////////////////////////////////////////////////////////////////////////
//*****************************************************************************
//********************pl_LOCALKL**************************************************
//*****************************************************************************

void pl_localKL(double *p, double *e,int e_num ,double *q, double *qn, double *qd, int npLC, double lK[15][15],
	double lL[15],double *params,double w, double *snorm){
	double tc=1.0;
	double lI[15][15];
	double dt,u1;
	int i,j;//,k;
	dt=params[0];
	////
	//mexPrintf("w : %f\n",w);
	double S=params[12];
		if (w<0) S = -0.5*S;   // weak homeotropic anchoring
			
	double A=1.0 / (6.0 * S);
		
//	mexPrintf("w %f A %f, S %f\n",w,A,S);
	
	//double as=(1.0/(3.0*params[12]))/2.0; //constant for isotropic part of surface energy
	double vx,vy,vz; //surface normal vector components
	
	
	///
	u1=params[9];
	memset(lK,0,15*15*sizeof(double));	//SET LOCAL MATRICES TO ZERO
	memset(lI,0,15*15*sizeof(double));
	memset(lL,0,15*sizeof(double));
	
	int gn[3]; //global node numbers of the element
	gn[0]	= (int)e[0];
	gn[1]	= (int)e[1];
	gn[2]	= (int)e[2];

	double Tii=2*A*w; //isotropic stiffness matrix contribution = simplifies to constant

	double Jdet=2*e[7];		//2*AREA OF TRIANGLE
	
	
	for (int igp=sngp;igp--;){ 
		
	double q1,q2,q3,q4,q5;
	q1=q2=q3=q4=q5=vx=vy=vz=0.0;
	for (i=0;i<3;i++){
		//mexPrintf("snorm [%f,%f,%f]\n",snorm[(gn[i]-1)*3],snorm[(gn[i]-1)*3+1],snorm[(gn[i]-1)*3+2]);
		q1+=ssh1[igp][i]*q[gn[i]-1]; //Q-tensor components with shape functions
		q2+=ssh1[igp][i]*q[gn[i]-1+npLC];
		q3+=ssh1[igp][i]*q[gn[i]-1+2*npLC];
		q4+=ssh1[igp][i]*q[gn[i]-1+3*npLC];
		q5+=ssh1[igp][i]*q[gn[i]-1+4*npLC];
				
		vx+=ssh1[igp][i]*snorm[(gn[i]-1)*3]; //surface normals with shape functions
		vy+=ssh1[igp][i]*snorm[(gn[i]-1)*3+1];
		vz+=ssh1[igp][i]*snorm[(gn[i]-1)*3+2];
		
		

	}//end for i
//mexPrintf("W=%f, v=[%f,%f,%f]\n",w,vx,vy,vz);
	

//RESIDUAL TERMS
double R1 = w*(-vx*vx*rt6 - vy*vy*rt6 + 2.0*rt6*vz*vz + 12.0*A*q1)/6.0;
double R2 = w*(vx*vx*rt2 - vy*vy*rt2  + 4.0*A*q2)/2.0;
double R3 = w*(vx*rt2*vy + 2.0*A*q3);
double R4 = w*(vy*rt2*vz + 2.0*A*q4);
double R5 = w*(vx*rt2*vz + 2.0*A*q5);

	double		mul=sw[igp]*Jdet;
  	for (i=3;i--;){
		double ShR=mul*ssh1[igp][i];
		//mexPrintf("Shr %f, mul %f, Jdet ,%f\n",1e15*ShR, 1e15*mul, 1e15*Jdet);

		for ( j=3;j--;){
		//double ShC=Sh[j];
		double Sh2=ShR*ssh1[igp][j];
				//K MATRIX
						
				lK[i][j   ]   += Sh2*Tii;
				lK[i+3][j+3 ] += Sh2*Tii;
				lK[i+6][j+6]  += Sh2*Tii;
				lK[i+9][j+9]  += Sh2*Tii;
				lK[i+12][j+12]+= Sh2*Tii;
				
				//Local Identity Matrix
				lI[i  ][j   ]+=Sh2;
				lI[i+3][j+3 ]+=Sh2;
				lI[i+6][j+6 ]+=Sh2;
				lI[i+9][j+9 ]+=Sh2;
				lI[i+12][j+12]+=Sh2;
				}//end for j


								//LOCAL RESIDUAL
				lL[i+0] +=ShR*(R1);// + S1v2);  
				lL[i+3] +=ShR*(R2);// + S2v2);
				lL[i+6] +=ShR*(R3);// + S3v2);
				lL[i+9] +=ShR*(R4);// + S4v2);
				lL[i+12]+=ShR*(R5);// + S5v2);
			}//end for i
	}//end for igp
	if(dt!=0)
	{// %0 calculates the steady state
			for ( i=15; i--;) {
			lL[i]=-lL[i];//Mqn[i]*2.0*u1/dt+u1*Mqd[i]-(lL[i]+Mqi[i]*2.0*u1/dt);
			}
	}//if(dt!=0)
	//mexPrintf("L[0] %f \n",lL[0]*1e20);
}//end void pl_localKL


/////////////////////////////////////////////////////////////////////////////////
//*****************************************************************************
//********************wk_LOCALKL**************************************************
//*****************************************************************************

void wk_localKL(double *p, double *e,int e_num ,double *q, double *qn, double *qd, int npLC, double lK[15][15],
	double lL[15],double *params, double V1[5], double V2[5]){
	double tc=1.0;
	double Sh[3];
	double lI[15][15];
	double dt,u1;
	int i,j;//,k;
	dt=params[0];
	
	double Ss=params[12];
	double W1=V1[1];//relative anchoring Strengths 
	double W2=V2[1];

	double A = (W1+W2)/(Ss*6.0); 
	//mexPrintf("A: %f \n",A);
	double w=V1[0]; //anchoring Strength
	
	u1=params[9];
	memset(lK,0,15*15*sizeof(double));	//SET LOCAL MATRICES TO ZERO
	memset(lI,0,15*15*sizeof(double));
	memset(lL,0,15*sizeof(double));
	int gn[3]; //global node numbers of the element
	
	gn[0]	= (int)e[0];
	gn[1]	= (int)e[1];
	gn[2]	= (int)e[2];
	
	double Jdet=e[7];///2.0; // determinant = half surface area		
	
	//mexPrintf("index:%i area: %g surface normal[x,y,z]=[%f,%f,%f]\n",e_num,Jdet*1e15,vx,vy,vz);
	for (int igp=0;igp<sngp;igp++){ 
	
		for( i=0;i<3;i++){
		Sh[i]=ssh1[igp][i];
		}
	
	double q1,q2,q3,q4,q5,v1x,v1y,v1z,v2x,v2y,v2z;
	q1=q2=q3=q4=q5=v1x=v1y=v1z=v2x=v2y=v2z=0;
	for (i=3;i--;){
		// Q-tensor components with shape functions
		q1+=Sh[i]*q[gn[i]-1];
		q2+=Sh[i]*q[gn[i]-1+npLC];
		q3+=Sh[i]*q[gn[i]-1+2*npLC];
		q4+=Sh[i]*q[gn[i]-1+3*npLC];
		q5+=Sh[i]*q[gn[i]-1+4*npLC];
		// vector components with shape functions
		 v1x+=Sh[i]*V1[2];// v1---- x,y and z-components
	     v1y+=Sh[i]*V1[3];
		 v1z+=Sh[i]*V1[4];

		 v2x+=Sh[i]*V2[2];// v2 ---- x,y and z-components
		 v2y+=Sh[i]*V2[3];
		 v2z+=Sh[i]*V2[4];


	}//end for i

//mexPrintf("W1 = %f, v1=[%f,%f,%f]    W2 = %f, v2=[%f,%f,%f]\n",W1,v1x,v1y,v1z,W2,v2x,v2y,v2z);

//Tii = thermotropic stiffness matrix terms
double		Tii=2*(A)*w;

//Ti = thermotropic residual terms

double	T1=2*(A*q1)*w;
double  T2=2*(A*q2)*w;
double 	T3=2*(A*q3)*w;
double 	T4=2*(A*q4)*w;
double 	T5=2*(A*q5)*w;

//vector v1 - terms
double S1v1 = (-v1x*v1x*rt6  - v1y*v1y*rt6 + 2*v1z*v1z*rt6)*w*W1/6.0;
double S2v1 = (v1x*v1x*rt2   - v1y*v1y*rt2)*w*W1/2.0;
double S3v1 = v1x * rt2 * v1y*w*W1;
double S4v1 = v1y * rt2 * v1z*w*W1;
double S5v1 = v1x * rt2 * v1z*w*W1;
//vector v2 - terms
double S1v2 = (-v2x*v2x*rt6 - v2y*v2y*rt6 +  2*v2z*v2z*rt6)*w*W2/6.0;
double S2v2 = (v2x*v2x*rt2  - v2y*v2y*rt2)*w*W2/2.0;
double S3v2 = v2x * rt2 * v2y *w*W2;
double S4v2 = v2y * rt2 * v2z *w*W2;
double S5v2 = v2x * rt2 * v2z *w*W2;




	double		mul=sw[igp]*Jdet;
  	for (i=0;i<3;i++){
		double ShR=mul*Sh[i];
		//mexPrintf("Shr %f, mul %f, Jdet ,%f\n",1e15*ShR, 1e15*mul, 1e15*Jdet);
		for ( j=0;j<3;j++){
		//double ShC=Sh[j];
		double Sh2=ShR*Sh[j];
		//K MATRIX
				lK[i][j   ]   += Sh2*Tii;
				lK[i+3][j+3 ] += Sh2*Tii;
				lK[i+6][j+6]  += Sh2*Tii;
				lK[i+9][j+9]  += Sh2*Tii;
				lK[i+12][j+12]+= Sh2*Tii;
				//Local Identity Matrix
				lI[i  ][j   ]+=Sh2;
				lI[i+3][j+3 ]+=Sh2;
				lI[i+6][j+6 ]+=Sh2;
				lI[i+9][j+9 ]+=Sh2;
				lI[i+12][j+12]+=Sh2;
				}//end for j
				//LOCAL RESIDUAL
				lL[i+0] +=ShR*(S1v1  +S1v2+T1);//  
				lL[i+3] +=ShR*(S2v1  +S2v2+T2);//
				lL[i+6] +=ShR*(S3v1  +S3v2+T3);//
				lL[i+9] +=ShR*(S4v1  +S4v2+T4);//
				lL[i+12]+=ShR*(S5v1  +S5v2+T5);//
			}//end for i
	}//end for igp
	if(dt!=0)
	{// %0 calculates the steady state
			for ( i=0; i<15; i++) {
			lL[i]=-lL[i];
		
			}
	}//if(dt!=0)
}//end void pl_localKL

void localKL_N(double *p, int t[4],double *ee,int tmat,double Jdet,double *v, int np,
			   double lK[20][20],double lL[20],double *params,int npLC,int it)
{
	int i,j;	
	double S0=params[12];
	double rt2 = sqrt(2.0);
	double rt6 = sqrt(6.0);
	double E11 = params[14];
	double E33 = params[15];
	double efe  = (2.0/S0/3.0)*(E11);//+2*e33);
	double efe2 = (4.0/S0/9.0)*(E11 - E33);
	
	memset(lK,0,20*20*sizeof(double));
	memset(lL,0,20*sizeof(double));
//separate node numbers, surface determinant and surface normal vector
	//e[0] = (int)ee[0]; e[1] = (int)ee[1]; e[2] = (int)ee[2];
	double n[3];
	n[0] = ee[4]; n[1] = ee[5]; n[2] = ee[6];				
	double eDet = ee[7];
// Jacobian
	for (int igp=0; igp<sngp; igp++) {
		double xr,xs,xt,yr,ys,yt,zr,zs,zt;
		xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;
		for (i=0; i<4; i++) 
		{
			xr+=sh1r[igp][i]*p[(t[i]-1)*3+0]*1e-6;
			xs+=sh1s[igp][i]*p[(t[i]-1)*3+0]*1e-6;
			xt+=sh1t[igp][i]*p[(t[i]-1)*3+0]*1e-6;

			yr+=sh1r[igp][i]*p[(t[i]-1)*3+1]*1e-6;
			ys+=sh1s[igp][i]*p[(t[i]-1)*3+1]*1e-6;
			yt+=sh1t[igp][i]*p[(t[i]-1)*3+1]*1e-6;

			zr+=sh1r[igp][i]*p[(t[i]-1)*3+2]*1e-6;
			zs+=sh1s[igp][i]*p[(t[i]-1)*3+2]*1e-6;
			zt+=sh1t[igp][i]*p[(t[i]-1)*3+2]*1e-6;
		}//end for i
						
		double Jinv[3][3]={(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet
						  ,(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet
						  ,(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet};

       double Sh[4],dSh[4][3];
       for(i=0;i<4;i++)
	   {
			Sh[i]=sh1[igp][i];
			dSh[i][0]=sh1r[igp][i]*Jinv[0][0]+sh1s[igp][i]*Jinv[1][0]+sh1t[igp][i]*Jinv[2][0];
            dSh[i][1]=sh1r[igp][i]*Jinv[0][1]+sh1s[igp][i]*Jinv[1][1]+sh1t[igp][i]*Jinv[2][1];
			dSh[i][2]=sh1r[igp][i]*Jinv[0][2]+sh1s[igp][i]*Jinv[1][2]+sh1t[igp][i]*Jinv[2][2];
			
		}//end for i
		
		double Vx=0, Vy=0, Vz=0;
		for(i=0;i<4;i++)
		{
			    Vx += dSh[i][0] * v[t[i]-1];
				Vy += dSh[i][1] * v[t[i]-1];
				Vz += dSh[i][2] * v[t[i]-1];
		}//end for i
		
	
	double mul=sw[igp]*eDet*2;
				
		for (i=0; i<4; i++) 
		{
		double ShR = Sh[i]*mul;	
		/*
		Lefe_1 = (rt6*(Vx*ShRx+Vy*ShRy-2.0*Vz*ShRz)*efe/6.0);
		Lefe_2 = (-rt2*(Vx*ShRx-Vy*ShRy)*efe/2.0);
		Lefe_3 = (-rt2*(Vy*ShRx+Vx*ShRy)*efe/2.0);
		Lefe_4 = (-rt2*(Vz*ShRy+Vy*ShRz)*efe/2.0);
		Lefe_5 = (-rt2*(Vz*ShRx+Vx*ShRz)*efe/2.0);
		*/
		lL[i+0] -= (rt6* (Vx*ShR*n[0]+Vy*ShR*n[1]-2.0*Vz*ShR*n[2])*efe/6.0);
		lL[i+4] -= (-rt2*(Vx*ShR*n[0]-Vy*ShR*n[1])*efe/2.0);
		lL[i+8] -= (-rt2*(Vy*ShR*n[0]+Vx*ShR*n[1])*efe/2.0);
		lL[i+12] -= (-rt2*(Vz*ShR*n[1]+Vy*ShR*n[2])*efe/2.0);
		lL[i+16] -= (-rt2*(Vz*ShR*n[0]+Vx*ShR*n[2])*efe/2.0);
		}//end for i
	
	}//end for igp
}// end void localKL


///////////////////////////////////////////////////////////////////////////////
// sparse_set
void sparse_set(double *Pr,int *Ir,int *Jc,int ri,int rj,double val)
{
// Binary search within column

	int k,k1,k2,cr=ri;
	k1=Jc[rj];
	k2=Jc[rj+1]-1;
	do{
		k=(k1+k2)>>1;
		if(cr<Ir[k])
			k2=k-1;
		else
			k1=k+1;
		if(cr==Ir[k])
			break;
	}while(k2>=k1);

	Pr[k]+=val;

	if(cr!=Ir[k])
	{
		//char str[50]="";
		//mexPrintf("\n direct_new.cpp - Error finding (%d, %d)\n",ri+1,rj+1);
		mexErrMsgTxt("direct_new.cpp - error setting sparse matrix \n ");
	}
}

////////////////////////////////////////////////////////////////////////////////
// assemble
void assemble(double *p,int *t, double *volume, int *tmat,double *q0,double *u0, double *qn, double *qd,double *v,
			  int nt, int npLC, int np, double *Pr,
			  int *Ir,int *Jc,double *L,double *params,int *per)
{
	init_shape();

//	char str[50];
	double lK[20][20];
	double lL[20];
	
	int *n;
	int ijstart;
	n=t;
	for (int it=0; it<nt; it++) {
		if( tmat[it]==4 ){
			
			localKL(p,n,&volume[it],q0,u0,qn,qd,v,npLC,np,lK,lL,params,it);
		
			ijstart=0;
	
		 for (int i=0;i<20;i++) 
		 {
			//int ri=n[i%4]+npLC*(i/4)-1;
			int ri=per[n[i%4]+npLC*(i/4)-1];
			L[ri]+=lL[i]*2e16;
			//mexPrintf("L %f, K %f \n",2e16*lL[i],2e16*lK[1][1]);
			for (int j=0;j<20;j++) 
			{
				int rj=per[n[j%4]+npLC*(j/4)-1];
				//mexPrintf("ri,rj =[%i,%i,], n=[%i,%i,%i,%i]\n",ri,rj,n[0],n[1],n[2],n[3]);
				 sparse_set(Pr,Ir,Jc,ri,rj,lK[i][j]*2e16);
			}//end for j 
		 }//end for i
		}//end for if tmat
		n+=4;
	}//end fr it
//mexPrintf("assemble OK\n");
}

////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
//************************************************************************************

void pl_assemble(double *e,int e_length,double *p,int npLC, double *q,double *qn, double *qd, int *surface_index,
			  int n_surfaces, double *W, double *params, double *L, double *Pr   , int *Ir, 
			  int *Jc,int *per, double *snorm){
	init_shape_surf();


//	char str[50];
	double lK[15][15];
	double lL[15];
	double w; // anchoring strenght
	register unsigned int ii,jj;
	int fixlc;
	register unsigned int gn[3]; //global node numbers
	//n=(int)e;
	for (int s=0; s<n_surfaces; s++){//loop over all surface indicies
	//mexPrintf("# of surfaces %i ,surface: %i index:%i weight:%f\n",n_surfaces,s,surface_index[s],W[s]);
	
	//mexPrintf("e_lenght %i\n",e_length);
	for (int it=0; it<e_length; it++) {
		fixlc=(int)e[it*8+3];
		w   =W[2*s];
		
		//mexPrintf("it:%i e(4,it)=%i\n",it,fixlc);
	
		if(((fixlc&31*2048)/2048)==surface_index[s]){ //If triangle is part of surface... 
			//mexPrintf("FixLC%i found \n",surface_index[s]);
			gn[0]=(int)e[it*8];
			gn[1]=(int)e[it*8+1];
			gn[2]=(int)e[it*8+2];

//mexPrintf("tri %i",it);
			pl_localKL(p,&e[it*8],it,q, qn, qd, npLC,lK,lL,params,w,snorm);
			//mexPrintf("surface triangle %iglobal nodes[%i,%i,%i] material %i\n",8*it,gn[0],gn[1],gn[2],fixlc);
			//ijstart=0;
//mexPrintf("OK\n");		
	
			for (int i=0;i<15;i++){ //LOCAL to GLOBAL
				int ri=per[gn[i%3]+npLC*(i/3)-1];
				//mexPrintf("local %i = global %i, L=%f\n",i,ri,lL[i]*1e20);
				L[ri]+=lL[i]*2e16;
				
			for (int j=0;j<15;j++){
				int	rj=per[gn[j%3]+npLC*(i/3)-1];
				ii=i;
				jj=j;
					
					if (i>j){//if row > column -> lower half 
					ii=j;
					jj=i;
					}
					sparse_set(Pr,Ir,Jc,ri,rj,lK[ii][jj]*2e16); //	

			}//end for j
			}//end fr i
		}//end if fixlc	
		

}//for s
}//for it
}//end void pl_assemble
//___________________________________________________________________________________________
//**********************ASSEMBLE WEAK ANCHORING SURFACES*************************************
//___________________________________________________________________________________________
void wk_assemble(double *e,int e_length,double *p,int npLC,double *q0,double *qn,double *qd,
				 int *wk_surface_index,int wk_n_surfaces,double *wk_vector1,double *wk_vector2,
				 double *params,double *L,double *Pr,int *Ir, int *Jc,int *per)
	{
	init_shape_surf();


//	char str[50];
	double lK[15][15];
	double lL[15];
	
	int ii,jj;
	int fixlc;
	int gn[3]; //global node numbers
	//n=(int)e;
	double V1[5],V2[5];
	for (int s=0; s<wk_n_surfaces; s++){//loop over all surface indicies
	//mexPrintf("# of surfaces %i ,surface: %i index:%i weight:%f\n",n_surfaces,s,surface_index[s],W[s]);
		V1[0]=V2[0]=wk_vector1[5*s+0];//anchoring strength
		V1[1]=wk_vector1[5*s+1];//relative strengths
		V2[1]=wk_vector2[5*s+1];
		V1[2]=wk_vector1[5*s+2];V1[3]=wk_vector1[5*s+3];V1[4]=wk_vector1[5*s+4]; //v1 x,y,z-components
		V2[2]=wk_vector2[5*s+2];V2[3]=wk_vector2[5*s+3];V2[4]=wk_vector2[5*s+4]; //v2 x,y,z-components

	//mexPrintf("e_lenght %i\n",e_length);
	for (int it=0; it<e_length; it++) {
		fixlc=(int)e[it*8+3];
		//mexPrintf("it:%i e(it,4)=%i\n",it,fixlc);
		//mexPrintf("%f\n",e[it*8+3]);
		
		if(((fixlc&31*2048)/2048)==wk_surface_index[s]){ //If triangle is part of surface... 
			//mexPrintf("FixLC%i found \n",surface_index[s]);
			gn[0]=(int)e[it*8];
			gn[1]=(int)e[it*8+1];
			gn[2]=(int)e[it*8+2];

			wk_localKL(p,&e[it*8],it,q0, qn, qd, npLC,lK,lL,params,V1,V2);
			//mexPrintf("global nodes[%i,%i,%i] material %i\n",gn[0],gn[1],gn[2],fixlc);
			//ijstart=0;
		
	
			for (int i=0;i<15;i++){ //LOCAL to GLOBAL
				int ri=per[gn[i%3]+npLC*(i/3)-1];
				//mexPrintf("local %i = global %i, L=%f\n",i,ri,lL[i]*1e20);
				L[ri]+=lL[i]*2e16;
				
			for (int j=0;j<15;j++){
				int	rj=per[gn[j%3]+npLC*(i/3)-1];
				ii=i;
				jj=j;
					
					if (i>j){//if row > column -> lower half 
					ii=j;
					jj=i;
					}
					sparse_set(Pr,Ir,Jc,ri,rj,lK[ii][jj]*2e16); //	

			}//end for j
			}//end fr i
		}//end if fixlc	
		

}//for s
}//for it
}//end void pl_assemble
//
// assemble Neumann boundaries
//
void assemble_surfaces(double *p, int *t,int *tmat, double *e, double *volume,double *v, int nt, int np, double *Pr, int *Ir, int *Jc, double *L, double *params, int npLC, int *e_to_t,int n_e, int *per)
{
	init_shape_N();
	double lK[20][20];
	double lL[20];
	double *ee;
	int *tt,ind,i,j;
	memset(lK,0,20*20*sizeof(double));
	memset(lL,0,20*sizeof(double));
	for (int it=0; it<n_e; it++)
	{
		// if points to a LC tet 
		//mexPrintf("it= %i, e_to_t[it]=%i\n",it,e_to_t[it]);
		if ((e_to_t[it]>=0)&&(tmat[e_to_t[it]-1]==4)) // if LC material
		{
			//mexPrintf("it = %i ",it);
			
			ee = &e[8*it]; //local surface element
			if (ee[3]==3) {mexPrintf("breaking\n");break;}
			ind = e_to_t[it]-1; // index to tet connected to current surface element
			tt = &t[4*ind]; // local tetrahedron
		
			//mexPrintf("e= [%i,%i,%i]  tt= [%i,%i,%i,%i]\n",(int)ee[0],(int)ee[1],(int)ee[2],tt[0],tt[1],tt[2],tt[3]);
			int intr=-1;//find  index to internal node
			for (i=0;i<4;i++)
			{
				if ((tt[i]!=(int)ee[0]) && (tt[i]!=(int)ee[1]) && (tt[i]!=(int)ee[2])) 
				{
					intr = i  ;
					break;
				}
			}
		
			int ti[4] = {(int)ee[0],(int)ee[1],(int)ee[2],tt[intr]}; // reordered local element
			localKL_N(p,ti,ee,4,volume[ind],v,np,lK,lL,params,npLC,it);
			
			for (i=0;i<20;i++)
			{
				int ri=per[ti[i%4]+npLC*(i/4)-1];
				L[ri]+= lL[i]*2e16;
					
				//	for (j=0; j<4; j++) {
				//	int rj = ti[j]-1;
				//	//sparse_set(Pr,Ir,Jc,ri,rj,-lK[i][j]);
				//}//end for j*/
			}//end for i
		}//end if LC element
		
	}//end for it
}//end void assemble_Neumann


void sort(int *ps,int*pe){
	int *p,*pd,*pdd;//p acending,p decending
	int val;
	// .. order the entries according to column indices
	// burble-sort is used

	for(p=ps; p<=pe; p++){
		for(pd=pe; pd>p; pd--)
		{
			pdd = pd-1;
			
			if (*pdd>*pd){
			//if (*pdd>*pd || *pdd==0){
				val=*pdd;
				*pdd = *pd;
				*pd=val;
			}
		}
	}
/*	// debug double check
	for(p=ps+1; p<=pe; p++){
		if (*(p-1)>*p) mexErrMsgTxt("Order messed up.");
	}*/

}

