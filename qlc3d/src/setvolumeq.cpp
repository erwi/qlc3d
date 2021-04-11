#include <math.h>
#include <qlc3d.h>
#include <lc-representation.h>
using namespace std;

void setNormalBox(  Box &box,
                    std::vector<qlc3d::Director> &dir,
                    double* p,
                    int npLC) {
/*! Sets the Q-tensor initial configuration (volume) within a normal Box*/
    // Twist, Tilt
    double bottomTwistDegrees = box.Twist[0];
    double bottomTiltDegrees = box.Tilt[0];
    double deltaTiltDegrees = box.Tilt[1]; // delta tilt bottom to top of box
    double deltaTwistDegrees = -1 * box.Twist[1];
    double boxHeight = box.Z[1] - box.Z[0];
    double power = box.Params[0];

    for (int i = 0; i < npLC; i++) { // loop over each node
        double px = p[i * 3 + 0];
        double py = p[i * 3 + 1];
        double pz = p[i * 3 + 2];
        double S = dir[i].S();
        // TODO: add Box::contains(x, y, z)
        if ((px >= box.X[0]) && (px <= box.X[1])) { // if within X
            if ((py >= box.Y[0]) && (py <= box.Y[1])) { // Y
                if ((pz >= box.Z[0]) && (pz <= box.Z[1])) {// Z
                    // NORMALISE COORDINATE W.R.T BOX SIZE
                    double pzn = (pz - box.Z[0]) / (boxHeight);
                    double twistDegrees = bottomTwistDegrees + pow(pzn * deltaTwistDegrees, power);
                    double tiltDegrees = bottomTiltDegrees + pow(pzn * deltaTiltDegrees, power);
                    dir[i] = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, S);
                } // if inside z limits
            }//if inside Y limits
        } // end if inside X limits
    }// end loop over each node
}//end void setNormalBox

void setRandomBox(Box &box, std::vector<qlc3d::Director> &dir, double* p, int npLC){
/*! Sets the Q-tensor initial configuration within a box volume. The dircor orientation is randomized */
    srand(0); // seed with constant value so that results are repeatable. TODO: make seed user defined configuration
    for (int i = 0 ; i < npLC ; i ++) {
        double S = dir[i].S();
        if ((p[i*3+0] >= box.X[0]) && (p[i*3+0] <= box.X[1])) // If this node is inside the box
        if ((p[i*3+1] >= box.Y[0]) &&(p[i*3+1] <= box.Y[1]))
        if ((p[i*3+2] >= box.Z[0])&&(p[i*3+2] <= box.Z[1])){
                double r1 =  (double) ( rand() % 10000 ) - 5000.0;
                double r2 =  (double) ( rand() % 10000 ) - 5000.0;
                double r3 =  (double) ( rand() % 10000 ) - 5000.0;
                double len = sqrt(r1*r1 + r2*r2 + r3*r3);
                dir[i] = qlc3d::Director(r1 / len, r2 / len, r3 / len, S);
        }
    }
}

void setHedgehogBox(Box box, std::vector<qlc3d::Director> &dir, double* p, int npLC){
/*! Sets the Q-tensor/director initial orientation bithin a box volume. The orientation is set
  so that a single hedgehog (+1) defect is located at the centre of the box*/
    // Director componenets are set to equal vectors from box cetre to director node location
    // resulting in a hedgehog defect
    // Calculate centre coordinates of this box
    double cc[3] = { ( box.X[0] + box.X[1] ) / 2.0 ,
                     ( box.Y[0] + box.Y[1] ) / 2.0 ,
                     ( box.Z[0] + box.Z[1] ) / 2.0};

    for (unsigned int i = 0 ; i < (unsigned int) npLC ; i++){ // loop over all LC nodes
        double S = dir[i].S();
        double *pl = &p[i*3]; // shortcut to coordinates of node i

        if ( box.isInBox(pl) ) { // if node i is inside box

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
            dir[i] = qlc3d::Director(v[0], v[1], v[2], S);
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
    // AFTER THIS, A "T-TENSOR" REPRESENTATION IS THEN CALCULATED FROM THE VECTORS

    // By default, when no boxes exist, the director is initialised along (1, 0, 0) at equilibrium order.
    auto defaultDirector = qlc3d::Director(1, 0, 0, lc->S0);
    std::vector<qlc3d::Director> dir(npLC, defaultDirector);

    // override the director within each box
    for (int i = 0; i < boxes->n_Boxes; i++) {
        Box &b = *boxes->box[i];
        switch (b.Type) {
            case Box::Normal:
                setNormalBox(b, dir, p, npLC);
                break;
            case Box::Random:
                setRandomBox(b, dir, p, npLC);
                break;
            case Box::Hedgehog:
                setHedgehogBox(b, dir, p, npLC);
                break;
            default:
                throw std::invalid_argument("unsupported box type " + b.TypeString);
        }
    }

    // convert the director to tensor
    for (int i = 0; i < npLC; i++) {
        q->setValue(i, qlc3d::TTensor::fromDirector(dir[i]));
    }
}
// end void setVolumeQ
