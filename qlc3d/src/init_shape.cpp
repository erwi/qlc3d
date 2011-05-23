
#include "gauss.h"
#include "qlc3d.h"

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
	printf("sh1[%i][1] = %f\n",i,sh1[i][1]);
	}
}
