#include <energy.h>
#include <cstdio>
#include <iostream>
namespace Energy
{
	#define SIGN(x)		(x)>=0 ? 1:-1


// Gauss integration
// ---------------------------------------------------------
//     3D Gauss-Legendre weights for N = 11, D = 4
// ---------------------------------------------------------
	const int ngp=11;
	const double a=(1+sqrt(5.0/14.0))/4.0;
	const double b=(1-sqrt(5.0/14.0))/4.0;

	static double gp[ngp][4]={
		{0.25	  		, 0.25		,	0.25		,0.25},
		{11.0/14.0     	, 1.0/14.0	,	1.0/14.0	,1.0/14.0},
		{1.0/14.0      ,	11.0/14.0	,	1.0/14.0	,1.0/14.0},
		{1.0/14.0	  ,	1.0/14.0	,	11.0/14.0   ,1.0/14.0},
		{1.0/14.0	  , 1.0/14.0	,	1.0/14.0	,11.0/14.0},
		{a		  ,	a		,	b       ,b},
		{a		  , b		,   a       ,b},
		{a         , b       ,   b      ,a},
		{b		  , a       ,   a       ,b},
		{b		  , a       ,   b		,a},
		{b		  , b		,   a		,a}};

	const double w11 = -74.0/5625.0;
	const double w12 = 343.0/45000.0;
	const double w13 = 56.0/2250.0;
	static double w[ngp]={w11,w12,w12,w12,w12,w13,w13,w13,w13,w13,w13};

	static double sh1[ngp][4]; // P1 Shape functions
	static double sh1r[ngp][4]; // P1 Shape functions r-derivatives
	static double sh1s[ngp][4]; // P1 Shape functions s-derivatives
	static double sh1t[ngp][4]; //P1 shape functions t-derivative

	static double rt2 = sqrt(2.0);
	static double rt6 = sqrt(6.0);
	//static double rt3 = sqrt(3.0);

void init_shape()
{
    printf("B");
    for (int i=0; i<ngp; i++)
    {
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
        //printf("B%i\n",i);
	}

    printf("D");
}// enf void init_shape

}// end namespace Energy
void CalculateFreeEnergy(FILE* fid, Simu* simu, LC* lc, Geometry* geom, SolutionVector* v, SolutionVector* q )
{
	using std::cout;
	using std::endl;

	if (! simu->getOutputEnergy() == 1) // RETURN IF ENERGY CALCULATION NOT SET TO TRUE
		return;


	cout << " Calculating free energy ..." << endl;

	if (simu->getCurrentIteration() == 0 ){	// IF FIRST ITERATION, PRINT HEADERS
		fprintf(fid,"%% columns are:\n");
		fprintf(fid,"%% time[s],distortion,thermotropic,dielectric\n");
		fprintf(fid,"F = [ ...\n");
	}
	//
	//
	//	BULK ENERGY
	//
	//
	using namespace Energy;

        init_shape();
        printf("A");
        double e0 = 8.8541878176*1e-12;
	double S0 = lc->S0;

	double epsav = lc->eps_per / S0;
	double deleps = (lc->eps_par - lc->eps_per) / S0 ;
	double A = lc->A;
	double B = lc->B;
	double C = lc->C;
	double efe  = (2.0/S0/3.0)*(lc->e11+2*lc->e33);
	double efe2 = (4.0/S0/9.0)*(lc->e11 - lc->e33);
	efe2+= 0; //no warnings...
	double f0=(3.0* A / 4.0)*(S0*S0) + (B/4.0)*(S0*S0*S0) + (9.0*C/16.0)*(S0*S0*S0*S0);
        double* p = geom->getPtrTop();

	// energy variables
	double Fd = 0;		// distortion energy
	double Fe = 0; 		// electric energy
	double Fflx = 0;	// flexoelectric enery
	double Fth	= 0;	// thermotropic energy

	//loop over each element and calculate elastic energy contribution
        for (idx x = 0 ; x < geom->t->getnElements() ; x++ ){
		if (geom->t->getMaterialNumber(x)==MAT_DOMAIN1){//if domain1 element, change test to bitand to include more domains

		// CALCULATE ELEMENT BARYCENTRE TO CHECK WHETHER IT IS INSIDE A VALID ENERGY CALCULATION REGION
		double bary[3];
		geom->getTetBaryCentre(bary , x);
		bary[0] = ( bary[0] - fabs(simu->getEnergyRegionX() ) ) *  ( (double) SIGN(simu->getEnergyRegionX()  ) );
		bary[1] = ( bary[1] - fabs(simu->getEnergyRegionY() ) ) *  ( (double) SIGN(simu->getEnergyRegionY()  ) );
		bary[2] = ( bary[2] - fabs(simu->getEnergyRegionZ() ) ) *  ( (double) SIGN( simu->getEnergyRegionZ() ) );//*
		//cout << "bary = " << bary[0] << "," << bary[1] <<"," << bary[2] << endl;

		if ( (bary[0] > 0) && (bary[1] > 0) && (bary[2] > 0) ){ // if within EnergyRegion

		    int tt[4] = {geom->t->getNode(x,0),
						 geom->t->getNode(x,1),
						 geom->t->getNode(x,2),
						 geom->t->getNode(x,3)};

		for (int igp=0;igp<ngp;igp++) {//loop over each gauss point

		    double q1=0, q2=0, q3=0, q4=0, q5=0;
		    double q1x=0,q2x=0,q3x=0,q4x=0,q5x=0;
		    double q1y=0,q2y=0,q3y=0,q4y=0,q5y=0;
		    double q1z=0,q2z=0,q3z=0,q4z=0,q5z=0;
		    double Vx=0, Vy=0, Vz=0;

		    // Inverse Jacobian
		    double xr,xs,xt,yr,ys,yt,zr,zs,zt;
		    xr=xs=xt=yr=ys=yt=zr=zs=zt=0.0;
		    for (int i=0; i<4 ;i++) {
			xr+=sh1r[igp][i]*p[tt[i]*3+0]*1e-6;
			xs+=sh1s[igp][i]*p[tt[i]*3+0]*1e-6;
			xt+=sh1t[igp][i]*p[tt[i]*3+0]*1e-6;

			yr+=sh1r[igp][i]*p[tt[i]*3+1]*1e-6;
			ys+=sh1s[igp][i]*p[tt[i]*3+1]*1e-6;
			yt+=sh1t[igp][i]*p[tt[i]*3+1]*1e-6;

			zr+=sh1r[igp][i]*p[tt[i]*3+2]*1e-6;
			zs+=sh1s[igp][i]*p[tt[i]*3+2]*1e-6;
			zt+=sh1t[igp][i]*p[tt[i]*3+2]*1e-6;
		    }
                //----------------
                // Jacobian
		    double Jdet = geom->t->getDeterminant(x);

		    double Jinv[3][3]={
			{(zt*ys-yt*zs)/Jdet,(xt*zs-zt*xs)/Jdet,(xs*yt-ys*xt)/Jdet},
			{(yt*zr-zt*yr)/Jdet,(zt*xr-xt*zr)/Jdet,(xt*yr-yt*xr)/Jdet},
			{(yr*zs-ys*zr)/Jdet,(xs*zr-xr*zs)/Jdet,(ys*xr-xs*yr)/Jdet}};
		    // Shape function derivatives
		    double dSh[4][3];
		    for(int i = 0 ; i < 4 ; i++ ) {
			dSh[i][0]=sh1r[igp][i]*Jinv[0][0]+sh1s[igp][i]*Jinv[1][0]+sh1t[igp][i]*Jinv[2][0];
			dSh[i][1]=sh1r[igp][i]*Jinv[0][1]+sh1s[igp][i]*Jinv[1][1]+sh1t[igp][i]*Jinv[2][1];
			dSh[i][2]=sh1r[igp][i]*Jinv[0][2]+sh1s[igp][i]*Jinv[1][2]+sh1t[igp][i]*Jinv[2][2];
		    }

		    for(int i=0; i<4 ; i++){
			q1+=sh1[igp][i]*q->getValue( geom->t->getNode(x,i) , 0) ;				// q1i * Ni  = A1
			q2+=sh1[igp][i]*q->getValue( geom->t->getNode(x,i) , 1) ;
			q3+=sh1[igp][i]*q->getValue( geom->t->getNode(x,i) , 2) ;
			q4+=sh1[igp][i]*q->getValue( geom->t->getNode(x,i) , 3) ;
			q5+=sh1[igp][i]*q->getValue( geom->t->getNode(x,i) , 4) ;

			q1x+=dSh[i][0]*q->getValue( geom->t->getNode(x,i) , 0 );
			q2x+=dSh[i][0]*q->getValue( geom->t->getNode(x,i) , 1 );
			q3x+=dSh[i][0]*q->getValue( geom->t->getNode(x,i) , 2 );
			q4x+=dSh[i][0]*q->getValue( geom->t->getNode(x,i) , 3 );
			q5x+=dSh[i][0]*q->getValue( geom->t->getNode(x,i) , 4 );

			q1y+=dSh[i][1]*q->getValue( geom->t->getNode(x,i) , 0 );
			q2y+=dSh[i][1]*q->getValue( geom->t->getNode(x,i) , 1 );
			q3y+=dSh[i][1]*q->getValue( geom->t->getNode(x,i) , 2 );
			q4y+=dSh[i][1]*q->getValue( geom->t->getNode(x,i) , 3 );
			q5y+=dSh[i][1]*q->getValue( geom->t->getNode(x,i) , 4 );

			q1z+=dSh[i][2]*q->getValue( geom->t->getNode(x,i) , 0 );
			q2z+=dSh[i][2]*q->getValue( geom->t->getNode(x,i) , 1 );
			q3z+=dSh[i][2]*q->getValue( geom->t->getNode(x,i) , 2 );
			q4z+=dSh[i][2]*q->getValue( geom->t->getNode(x,i) , 3 );
			q5z+=dSh[i][2]*q->getValue( geom->t->getNode(x,i) , 4 );

                    // electric fields
			Vx+=dSh[i][0]*v->getValue( geom->t->getNode(x,i) );
			Vy+=dSh[i][1]*v->getValue( geom->t->getNode(x,i) );
			Vz+=dSh[i][2]*v->getValue( geom->t->getNode(x,i) );
		    }//end for i

			double R = q1*q1 + q2*q2 + q3*q3 + q4*q4 + q5*q5;
			double mul=w[igp]*Jdet;

			double Fd_elem =       (lc->L1/2)*(q1x*q1x +q2x*q2x +q3x*q3x +q4x*q4x +q5x*q5x
							+ q1y*q1y +q2y*q2y +q3y*q3y +q4y*q4y +q5y*q5y
							+ q1z*q1z +q2z*q2z +q3z*q3z +q4z*q4z +q5z*q5z);
				
                if (lc->L2)
				{        
                    Fd_elem+= (lc->L2/2)*(q1x*q1x/12.0+q2x*q2x/4.0+q1y*q1y/12.0+q2y*q2y/4.0+q3y*q3y/4.0+q3x*q3x/4.0+q5z*q5z/4.0+q5x*q5x/4.0+q4z*q4z/4.0+q1z*q1z/3.0+q4y*q4y/4.0+q3y*q5z/2.0
                             +q3x*q4z/2.0+q5x*q4y/2.0+q3y*q2x/2.0+q5z*q2x/2.0-q3x*q2y/2.0-q4z*q2y/2.0-q1x*
                             rt6*q2x*rt2/12.0+q1y*rt6*q2y*rt2/12.0+q5x*rt2*q1z
                             *rt6/6.0+q4y*rt2*q1z*rt6/6.0-q3y*rt2*q1x*rt6/12.0
                            -q5z*rt2*q1x*rt6/12.0-q3x*rt2*q1y*rt6/12.0-q4z*rt2*q1y*rt6/12.0);
				}
			    
				if (lc->L4 != 0){ // chiral term
			    Fd_elem+= (lc->L4/2)* ( q2*q3z - ( sqrt(3.0)*(3*q1*q4x - 3*q1x*q4 - 3*q1*q5y + 3*q1y*q5))/6 - q2z*q3 - 
				(q2*q4x)/2 + (q2x*q4)/2 - (q2*q5y)/2 + (q2y*q5)/2 - (q3*q4y)/2 + (q3y*q4)/2 + (q3*q5x)/2 - (q3x*q5)/2 - (q4*q5z)/2 + (q4z*q5)/2 );
			    		
				}
				
				if (lc->L6)
				{
							double MapleGenVar2 = -q5x*q5x*q1*rt6/12.0+q5x*q5x*q2*rt2/4.0-q3y*q3y*q1*rt6/12.0-q3y*q3y*q2*rt2/4.0+q1*rt6*q4z*q4z/6.0+q1*rt6*q1z*q1z/6.0-q4x*q4x*q1*rt6/12.0+q4x*q4x*q2*rt2/4.0-q5y*q5y*q1
                                *rt6/12.0-q3x*q3x*q1*rt6/12.0+q3x*q3x*q2*rt2/4.0-q5y*q5y*q2*rt2/4.0-q4y*q4y*q1*rt6/12.0-q4y*q4y*q2*rt2/4.0+q1*rt6*q3z*q3z/6.0+q1*rt6*q5z*q5z/6.0-q2*rt2*q1y*q1y/4.0-q2*rt2*q2y*
                                q2y/4.0-q1*rt6*q1x*q1x/12.0-q1*rt6*q2x*q2x/12.0;

                            double MapleGenVar1 = MapleGenVar2+q2*rt2*q1x*q1x/4.0+q2*rt2*q2x*q2x/4.0+q1*rt6*q2z*q2z/6.0+q5*rt2*q1x*q1z/2.0+q4*rt2*q4y*q4z/2.0+q4*rt2*q1y*q1z/2.0+q3*rt2*q2x*q2y/2.0+q4*rt2*q2y*q2z/2.0+q3*
                                rt2*q3x*q3y/2.0+q3*rt2*q5x*q5y/2.0+q3*rt2*q4x*q4y/2.0+q5*rt2*q3x*q3z/2.0+q3*rt2*q1x*q1y/2.0+q4*rt2*q5y*q5z/2.0+q5*rt2*q5x*q5z/2.0+q5*rt2*q2x*q2z/2.0+q4*rt2*q3y*q3z/2.0+q5*rt2*
                                q4x*q4z/2.0-q1*rt6*q1y*q1y/12.0-q1*rt6*q2y*q2y/12.0;

                            Fd_elem+= MapleGenVar1*(lc->L6/2);
                 }

				

			Fd += mul*Fd_elem;

			    double Fel_elem = e0*(-Vx*Vx - Vy*Vy- Vz*Vz)*epsav*0.5 +
                                e0*deleps*(Vx*Vx*q1*rt6/12.0  - Vx*Vx*q2*rt2/4.0 - Vx*Vy*q3*rt2/2.0  -  Vx*Vz*q5*rt2/2.0
                                -Vy*Vz*q4*rt2/2.0   - Vz*Vz*q1*rt6/6.0 + Vy*Vy*q2/4.0);

			Fe += mul*Fel_elem;

		    if (efe != 0){// if flexoelectric terms
			double Fflexo = efe*(2.0/3.0)*( Vx*((-1.0/6.0)*q1x*rt6 +0.5*q2x*rt2 + 0.5*q3y*rt2 + 0.5*q5z*rt2)
					  +Vy*(0.5*q3x*rt2 - (1.0/6.0)*q1y*rt6 - 0.5*q2y*rt2 + 0.5*q4z*rt2)
					  +Vz*(0.5*q1x*rt2 + 0.5*q4y*rt2 + (1.0/3)*q1z*rt6));

			Fflx += mul*Fflexo;
                }// end if flexoelectric terms

			double Fth_elem=A*(R)/2.0 +
					 B*(q5*q5*q1*rt6/4.0 - q1*rt6*q2*q2/2.0
					-q3*q3*q1*rt6/2.0 + 3.0/4.0*q5*q5*q2*rt2 + 3.0/2.0*q3*rt2*q5*q4
					+q4*q4*q1*rt6/4.0 - 3.0/4.0*q4*q4*q2*rt2 + q1*q1*q1*rt6/6.0)/3.0
					+C*(R*R)/4.0;

			Fth +=mul*( Fth_elem - f0);

		    }//end for igp, loop through gauss points
		}// end if within EnergyRegion
	    }//end if domain1 element
	}//end for loop through each element
	//
	//	END BULK ENERGY
	//
	fprintf(fid,"%e\t%e\t%e\t%e;\n",simu->getCurrentTime(), Fd , Fth , Fe);
	printf("OK\n");
}

void closeEnergyFile(FILE* fid, Simu& simu){
	if (simu.getOutputEnergy()==1){
		fprintf(fid, "];");
		fclose(fid);
	}
}
