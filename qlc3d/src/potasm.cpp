# include "qlc3d.h"
# include <math.h>
# include "material_numbers.h"
# include "gauss.h"
const int	npt = 4; //Number of Points per Tetrahedra

void assemble_volume(double *p,SolutionVector *v, Mesh *mesh, SparseMatrix *K, double* L);


void init_shapes()
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
void potasm(SolutionVector *v, Mesh *mesh, double *p, SparseMatrix *K, double* L)
{
	init_shapes();
	assemble_volume(p,v,mesh, K , L);




};


void localKL(double *p,int *tt, int mat, double lK[npt][npt],double lL[npt],int it, Mesh* mesh)
{
	int i,j;	
	double eps;
	
	
	
	//printf("t[%i] = [%i,%i,%i,%i]\n", it, tt[0], tt[1], tt[2], tt[3]);
	
	eps = 1;
	

	memset(lK,0,4*4*sizeof(double));
	memset(lL,0,4*sizeof(double));

	for (int igp=0; igp<ngp; igp++) {
		// Jacobian
		double xr,xs,xt,yr,ys,yt,zr,zs,zt;
		xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;
		for (i=0; i<npt; i++) {
			xr+=sh1r[igp][i]*p[(tt[i])*3+0];//*1e-6;
			xs+=sh1s[igp][i]*p[(tt[i])*3+0];//*1e-6;
			xt+=sh1t[igp][i]*p[(tt[i])*3+0];//*1e-6;

			yr+=sh1r[igp][i]*p[(tt[i])*3+1];//*1e-6;
			ys+=sh1s[igp][i]*p[(tt[i])*3+1];//*1e-6;
			yt+=sh1t[igp][i]*p[(tt[i])*3+1];//*1e-6;

			zr+=sh1r[igp][i]*p[(tt[i])*3+2];//*1e-6;
			zs+=sh1s[igp][i]*p[(tt[i])*3+2];//*1e-6;
			zt+=sh1t[igp][i]*p[(tt[i])*3+2];//*1e-6;
		}//end for i
		double Jdet= mesh->getDeterminant(it);//(xr*ys*zt-xr*zs*yt+xs*yt*zr-xs*yr*zt+xt*yr*zs-xt*ys*zr);
		if (Jdet<0) Jdet = -Jdet;// printf("negative jacobian!\n");
		
		double Jinv[3][3]={(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet
						  ,(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet
						  ,(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet};
       
       double Sh[4],dSh[4][3];
	   // x,y,z derivatives of shape functions
       for(i=0;i<4;i++)
	   {
			Sh[i]=sh1[igp][i];
            
			dSh[i][0]=sh1r[igp][i]*Jinv[0][0]+sh1s[igp][i]*Jinv[1][0]+sh1t[igp][i]*Jinv[2][0];
            dSh[i][1]=sh1r[igp][i]*Jinv[0][1]+sh1s[igp][i]*Jinv[1][1]+sh1t[igp][i]*Jinv[2][1];
			dSh[i][2]=sh1r[igp][i]*Jinv[0][2]+sh1s[igp][i]*Jinv[1][2]+sh1t[igp][i]*Jinv[2][2];
			//printf("dsh[i][0] = %f\n",sh1s[igp][i]);
		}//end for i
		
		
		//if (it == 0)
		//{
		//	printf("dshx = %f, %f, %f %f\n", dSh[0][0],dSh[1][0],dSh[2][0],dSh[3][0] );
		
		//}
		
		// Local K and L
        double mul=w[igp]*Jdet;
        //printf("mul = %f, Jdet = %f\n",mul*1e18,Jdet*1e19);		
		
		for (i=0; i<npt; i++) {
			//double ShR = sh1[igp][i];
			//lL[i]-= mul*(dSh[i][0]*P1 + dSh[i][1]*P2 + dSh[i][2]*P3);
            for (j=0; j<npt; j++) 
			{
                lK[i][j]+=eps*mul*(
                        dSh[i][0]*dSh[j][0]+
                        dSh[i][1]*dSh[j][1]+
                        dSh[i][2]*dSh[j][2]);
						
			}//end for j
		}//end for i
	}//end for igp
}// end void localKL


void assemble_volume(double *p,SolutionVector *v, Mesh *mesh, SparseMatrix *K, double* L)
{
	init_shape();
	
	double lK[npt][npt];
	double lL[npt];
	
	//int nt = mesh->getnElements();
	
	int t[4] = {0,0,0,0};
	int tmat;
	int nNodes = mesh->getnNodes();
	
	for (int it=0; it<mesh->getnElements(); it++) 
	{
		t[0] =  mesh->getNode(it,0);
		t[1] =  mesh->getNode(it,1);
		t[2] =  mesh->getNode(it,2);
		t[3] =  mesh->getNode(it,3);
		tmat = mesh->getMaterialNumber(it);
		
		
		localKL(p,t,tmat,lK,lL,it, mesh);
		
		for (int i=0; i<npt; i++)
		{
			int ri=t[i];
		
		
		if ( v->getIsFixed(ri) )
		{
			//printf("ri %i, val = %f\n",ri,v->getValue(ri));
			for (int j = 0; j<4 ; j++)
			{
				L[t[j]]+=lK[i][j]*v->getValue(ri);
				L[ri]=0;
			}
		}
			
			
			for (int j=0; j<npt; j++)
			{
				int rj=t[j];
				
				if ((v->getIsFixed(ri)) || (v->getIsFixed(rj)) )
				{
					if (ri==rj) K->sparse_set(ri,rj,1.0);
				}
				else
				{
					K->sparse_add(ri,rj,lK[i][j]);
				}
			}//end for j
			
		}//end for i
	
	}//end for it*/
	
}// end void assemble_volume








