#include <extras.h>
int FindElemByCoord(Geometry* geom , double px, double py, double pz)
{
	printf("FindElemByCoord\n");

	return 0;


}


double Interpolate(int* t, SolutionVector* sol , double D[5])
{
	
	return 	D[0]*sol->getValue(t[0]) + 
			D[1]*sol->getValue(t[1]) + 
			D[2]*sol->getValue(t[2]) +
			D[3]*sol->getValue(t[3]);

}
double Interpolate(int* t, int dim, SolutionVector* sol , double D[5])
{
	
	return 	D[0]*sol->getValue(t[0], dim) + 
			D[1]*sol->getValue(t[1], dim) + 
			D[2]*sol->getValue(t[2], dim) +
			D[3]*sol->getValue(t[3], dim);

}


void extras(Geometry* geom, SolutionVector* sol)
{
	double P[3] = {0,0,0.5};
	double D[5] = {0,0,0,0,0};
	int N = 30;
	double dx = 1.0 / N; 
	//double cumsum = 0;
	//double aver = 0;
	bool found = false;
	
	double n_dst = 1e9;
	double n_ind  = -1;
	
	double a[5] = {0,0,0,0,0};
	
	
	
	for (int i = 0 ; i < N ; i++) // loop over rows
	{
		P[1] = i*dx; // set row
		for (int j = 0 ; j < N ; j++) // loop over cols
		{
			P[0] = j*dx; // set col
			found = false;
			for ( int e = 0 ; e < geom->t->getnElements(); e++)
			{
				double dst = geom->t->ContainsCoordinate( e , geom->getPtrTop(), P, D); 
				
				
				if (dst < 1e-10)
				{
					int* t = geom->t->getPtrToElement(e);
					
					a[0] += Interpolate(t,0, sol, D);
					a[1] += Interpolate(t,1, sol, D);
					a[2] += Interpolate(t,2, sol, D);
					a[3] += Interpolate(t,3, sol, D);
					a[4] += Interpolate(t,4, sol, D);
								
					
					found = true;
					break;
				}// if found
				else if (dst < n_dst)
				{	
					n_dst = 0 ; n_ind = e;
				}
			}// end for loop over elements
			if (!found) 
			{
				printf("could not find [%f,%f,%f] in any element\n", P[0] , P[1] , P[2]); 
				printf("nearest element is %i, distance %e\n", n_ind, n_dst);
				exit(1);}
		}//end for j
	} // end for i
	
	a[0]=a[0]/ (N*N);
	a[1]=a[1]/ (N*N);
	a[2]=a[2]/ (N*N);
	a[3]=a[3]/ (N*N);
	a[4]=a[4]/ (N*N);
	
	double *n = tensortovector(a,1); // get vector data 
	
	double tlt = 180* fabs (asin(n[2])) / 3.14159;
	double twt = 180* fabs (atan(n[1] / n[0])) / 3.14159; 
	printf("n = [%f, %f, %f], S = %f, %f\n", n[0] , n[1] , n[2], n[3] , n[4]);
	printf("tlt,twt = %f , %f\n", tlt, twt);
	free(n);
	
}
