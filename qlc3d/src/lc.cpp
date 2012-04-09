#include <lc.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
LC::LC()
{
    PhysicsFormulation = K3; // USE MORI'S 3K FORMULATION BY DEFAULT
    //printf("===========================================\n");
	A =  0.0;
	B = -0.3e6;
	C =  0.5e6;
	S0=(-B+sqrt(B*B-24*A*C))/(6*C);
	// should check for 'bad' S0 here - i.e. when it goes imaginary
	//printf("\tA = %1.2f, B = %1.2f, C = %1.2f -> S0 = %f\n", A*1e-6,B*1e-6,C*1e-6,S0);

	// elastic coefficients
	K11 = 10e-12;
	K22 = 10e-12;
	K33 = 10e-12;
	//printf("\tK11 = %2.1fpN, K22 = %2.1fpN, K33 = %2.1fpN\n",K11*1e12,K22*1e12,K33*1e12);

	L1=2.0*(K33-K11+3.0*K22)/(S0*S0*27.0);
	L2=4.0*(K11-K22)/(9.0*S0*S0);
	L3 = 0;
    L4 = 8.0*2*3.14159*K22/ ( p0 * 9.0 * S0*S0);// chiral
	L5 = 0;
	L6=4.0*(K33-K11)/(S0*S0*S0*27.0);
	//printf("\tL1 = %f, L2 = %f, L6 = %f\n",L1,L2,L6);


	// permittivity
	eps_par = 5;
	eps_per = 3;

	//flexoelectricity
	e11 = 0;
	e33 = 0;

	// viscosities
	u1 = 0.01;
	u2 = 0.01;
	gamma1 = 0.1;
	gamma2 = 0;
	alpha1 = 0;
	alpha4 = 0;
    alpha5 = 0;
    alpha6 = 0;
	//printf("\t%c%c = %f , %c|| = %f\n",238,193 , eps_per,238,eps_par);
	//printf("============================================================\n");
}

void LC::printLC()
{

	printf("\t[K11,K22,K33]\t\t=\t[%1.1f, %1.1f, %1.1f] pN\n", K11*1e12,K22*1e12,K33*1e12);
	printf("\tp0\t\t\t=\t%f microns\n",p0*1e6);
	printf("\t[A,B,C]\t\t\t=\t[%1.2f, %1.2f, %1.2f]*1e6\n", A*1e-6,B*1e-6,C*1e-6);
	printf("\tS0 \t\t\t=\t%f\n",S0);
	printf("\t[eps_par, eps_per]\t=\t[%1.1f, %1.1f]\n",eps_par,eps_per);
	printf("\t[e11,e33]\t\t= \t[%1.1f, %1.1f]*1e-12\n",e11*1e12,e33*1e12);
	printf("\t[gamma1,gamma2]\t\t=\t[%1.1f, %1.1f]\n",gamma1,gamma2);
	printf("\talpha[1,4,5,6]\t\t=\t[%1.1f, %1.1f, %1.1f, %1.1f]\n", alpha1,alpha4,alpha5,alpha6);


}// end void printLC

void LC::WriteLC(FILE* fid)
{
	if (fid!=NULL)
	{
		fprintf(fid,"#================================\n");
        fprintf(fid,"#  LC MATERIAL PARAMETERS\n");
        fprintf(fid,"#================================\n\n");
		fprintf(fid,"\tK11 = %2.4e\n",K11);//K11
		fprintf(fid,"\tK22 = %2.4e\n",K22);//K22
		fprintf(fid,"\tK33 = %2.4e\n",K33);//K33
		fprintf(fid,"\tp0  = %2.4e\n\n",p0);//p0

		fprintf(fid,"\tA   = %2.4e\n",A);//A
		fprintf(fid,"\tB   = %2.4e\n",B);//B
		fprintf(fid,"\tC   = %2.4e\n\n",C);//C
		// Dielectric and flexoelectric coefficients
		fprintf(fid,"\teps_par = %2.4f\n",eps_par);//eps_par	= 18.5; # parallel
		fprintf(fid,"\teps_per = %2.4f\n",eps_per);//eps_per = 7.0;  # perpendicular
		fprintf(fid,"\te11     = %2.4e\n",e11);//e11 = 0e-12;  # flexoelectric splay
		fprintf(fid,"\te33     = %2.4e\n\n",e33);//e33 = 0e-12;  # flexoelectric bend

		// Leslie viscous coefficients
		fprintf(fid,"\tgamma1  = %2.4f\n",gamma1);//gamma1 = 0.1;
		fprintf(fid,"\tgamma2  = %2.4f\n",gamma2);//gamma2 = 0.1;
		fprintf(fid,"\talpha1  = %2.4f\n",alpha1);//alpha1 = 1;
		fprintf(fid,"\talpha4  = %2.4f\n",alpha4);//alpha4 = 4;
		fprintf(fid,"\talpha5  = %2.4f\n",alpha5);//alpha5 = 5;
		fprintf(fid,"\talpha6  = %2.4f\n",alpha6);//alpha6 = 6;
		fprintf(fid,"\n\n");
	}
}

void LC::convert_params_n2Q()
{
	//S0 = S0=(-B+sqrt(B*B-24*A*C))/(6*C);

	S0=(-B+sqrt(B*B-24*A*C))/(6*C);

	if ( (S0<0) || (S0 > 1) ){
	    std::cout << "ERROR!, unusual equilibrium order parameter value: "<< S0 << std::endl;
	    std::cout	<< "Check Thermotropic coefficients, A:" << A
			<< " B: " << B
			<< " C:" << C << std::endl;
	    exit(1);
	}

        // WHERE DOES EXTRA FACTOR OF 2 COME IN L-TERMS?
        L1=2.0*(K33-K11+3.0*K22)/(S0*S0*27.0);
	L2=4.0*(K11-K22)/(9.0*S0*S0);
	L3 = 0;
	L4 = 0.0;
	if (p0!=0)
	{
		double q0 = 2*3.14159265/p0;
                L4 = ( 8.0 * q0 * K22 ) / ( S0 * S0 * 9.0 );
	}
	L5 = 0;
	L6=4.0*(K33-K11)/(S0*S0*S0*27.0);

// VISCOSITIES
	u1 = 4*gamma1 / (9 * S0);

}

double LC::getS0()
{
	return S0;
}

