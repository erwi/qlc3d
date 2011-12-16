#include <math.h>
#include <time.h>
#include <qlc3d.h>
using namespace std;

void setNormalBox(  Box* box,
                    double* nx,
                    double* ny,
                    double* nz,
                    double* p,
                    int npLC)
{
/*! Sets the Q-tensor initial configuration (volume) within a normal Box*/
    double d_th = box->Tilt[1] *PI/180.0; // delta tilt bottom to top of box
    double d_ph = box->Twist[1] *PI/180.0;// delta twist
    double power = box->Params[0];

    for (int i = 0; i < npLC ; i++) // loop over each node
    {
	double px = p[i*3 + 0];
	double py = p[i*3 + 1];
	double pz = p[i*3 + 2];

	if ((px >= box->X[0]) && ( px <= box->X[1])) // if within X
	{
            if (( py >= box->Y[0]) &&( py <= box->Y[1])) // Y
            {
                if (( pz >= box->Z[0])&&( pz <= box->Z[1])) // Z
                {

                    // Twist, Tilt
                    double  ph = box->Twist[0] * PI /180.0 + pow(p[i*3+2]*d_ph, power);
                    double  th =  pow(( pz- box->Z[0] ) / (box->Z[1] - box->Z[0]), power ) * d_th  + box->Tilt[0]*PI/180.0;

                    nx[i] = cos(th)*cos(ph); // sets director
                    ny[i] = sin(ph)*cos(th);
                    nz[i] = sin(th);
                } // if inside z limits
            }//if inside Y limits
	} // end if inside X limits
    }// end loop over each node

}//end void setNormalBox
void setRandomBox(Box* box, double* nx, double* ny, double* nz, double* p, int npLC){
/*! Sets the Q-tensor initial configuration within a box volume. The dircor orientation is randomized */
    srand(time(NULL) );
    for (int i = 0 ; i < npLC ; i ++){
        if ((p[i*3+0] >= box->X[0]) && (p[i*3+0] <= box->X[1])) // If this node is inside the box
        if ((p[i*3+1] >= box->Y[0]) &&(p[i*3+1] <= box->Y[1]))
        if ((p[i*3+2] >= box->Z[0])&&(p[i*3+2] <= box->Z[1])){
                double th = (double) rand(); // generate random director orientations
                double ph = (double) rand();
		
                nx[i] = cos(th)*cos(ph);
                ny[i] = sin(ph)*cos(th);
                nz[i] = sin(th);
        }
    }
}

void setHedgehogBox(Box* box, double* nx, double* ny, double* nz, double* p, int npLC){
/*! Sets the Q-tensor/director initial orientation bithin a box volume. The orientation is set
  so that a single hedgehog (+1) defect is located at the centre of the box*/

    // Director componenets are set to equal vectors from box cetre to director node location
    // resulting in a hedgehog defect

    // Calculate centre coordinates of this box
    double cc[3] = { ( box->X[0] + box->X[1] ) / 2.0 ,
                     ( box->Y[0] + box->Y[1] ) / 2.0 ,
                     ( box->Z[0] + box->Z[1] ) / 2.0};

    for (unsigned int i = 0 ; i < (unsigned int) npLC ; i++){ // loop over all LC nodes
        double *pl = &p[i*3]; // shortcut to coordinates of node i

        if ( box->isInBox(pl) ){ // if node i is inside box

            // calculate vector from box centre to this node
            double v[3] = {pl[0] - cc[0] ,
                           pl[1] - cc[1] ,
                           pl[2] - cc[2] };
            // normalise vector lengt to 1
            double len = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
            len = sqrt(len);
            v[0]/= len;
            v[1]/= len;
            v[2]/= len;

            // set director to vector
            nx[i] = v[0];
            ny[i] = v[1];
            nz[i] = v[2];
        }
    }
}

void SetVolumeQ(
	SolutionVector *q,
	LC* lc,
	Boxes* boxes,
	double* p){
/*! Sets initial Q-tensor volume configuration for all boxes*/

    int npLC = q->getnDoF() ;
    // LC TILT AND TWIST IS FIRST CALCULATED AS VECTORS
    // AFTER THIS, A Q-TENSOR REPRESENTATION IS THEN CALCULATED FROM THE VECTORS

    // ALLOCATE MEMORY FOR VECTOR REPRESENTATION OF LC
    double *nx = (double*)malloc(npLC * sizeof(double));
    double *ny = (double*)malloc(npLC * sizeof(double));
    double *nz = (double*)malloc(npLC * sizeof(double));

    if (nx == NULL) {printf("error - could not allocate nx - bye"); exit(1);}
    if (ny == NULL) {printf("error - could not allocate ny - bye"); exit(1);}
    if (nz == NULL) {printf("error - could not allocate nz - bye"); exit(1);}
    // intialize all to 0 tilt, 0 twist
    for (int i = 0 ; i < npLC ; i++){
	nx[i] = 1.0;
	ny[i] = 0.0;
	nz[i] = 0.0;
    }
    // SET EACH BOX SEPARATELY
    for (int i = 0; i < boxes->n_Boxes; i++){
	if (boxes->box[i]->Type.compare("Normal") == 0)
	    setNormalBox(boxes->box[i], nx,ny,nz ,p, npLC);
	else if (boxes->box[i]->Type.compare("Random") == 0)
	    setRandomBox(boxes->box[i], nx,ny,nz, p , npLC);
        else if (boxes->box[i]->Type.compare("Hedgehog") == 0){
            setHedgehogBox(boxes->box[i], nx,ny,nz,p,npLC);
        }
        else{

            std::string type = boxes->box[i]->Type;
            printf("error in SetVolumeQ - box type %s is not supported - bye!\n", type.c_str() );
        }
    }
	
    // CONVERT VECTOR TO TENSOR
    double *a1 = (double*)malloc (npLC * sizeof (double));// memory for temporary Qxx
    double *a2 = (double*)malloc (npLC * sizeof (double));// memory for temporary Qyy
    double *a3 = (double*)malloc (npLC * sizeof (double));// memory for temporary Qxy
    double *a4 = (double*)malloc (npLC * sizeof (double));// memory for temporary Qyz
    double *a5 = (double*)malloc (npLC * sizeof (double));// memory for temporary Qxz
	
    for (int i =0; i<npLC ; i++){

        a1[i] = (lc->S0)*(3* nx[i] * nx[i] - 1)/2.0;
	a2[i] = (lc->S0)*(3* ny[i] * ny[i] - 1)/2.0;
	a3[i] = (lc->S0)*(3* nx[i] * ny[i])/2.0;
	a4[i] = (lc->S0)*(3* ny[i] * nz[i])/2.0;
	a5[i] = (lc->S0)*(3* nx[i] * nz[i])/2.0;

        //if (i == 1053)
        //{
        //    printf(" n = %f,%f,%f\n", nx[i], ny[i], nz[i]);
        //}
    }

    // convert tensor basis to symmetric traceless
    double rt6 = sqrt(6);
    double rt2 = sqrt(2);
    for (int i =0; i<npLC ; i++){
	q->Values[i+0*npLC] = (a1[i]+a2[i])*(rt6/-2);
	q->Values[i+1*npLC] = (a1[i]+(a1[i]+a2[i])/-2)*rt2;
	q->Values[i+2*npLC] = a3[i]*rt2;
	q->Values[i+3*npLC] = a4[i]*rt2;
	q->Values[i+4*npLC] = a5[i]*rt2;
    }
	
    free(a1);
    free(a2);
    free(a3);
    free(a4);
    free(a5);
    free(nx);
    free(ny);
    free(nz);
}
// end void setVolumeQ
