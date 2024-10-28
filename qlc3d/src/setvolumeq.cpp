#include <math.h>
#include <qlc3d.h>
#include <lc-representation.h>
#include <util/logging.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
using namespace std;

void setNormalBox(  const Box &box,
                    std::vector<qlc3d::Director> &dir,
                    const Coordinates& coordinates,
                    int npLC) {
/*! Sets the Q-tensor initial configuration (volume) within a normal Box*/
    // Twist, Tilt
    double bottomTwistDegrees = box.Twist[0];
    double bottomTiltDegrees = box.Tilt[0];
    double deltaTiltDegrees = box.Tilt[1]; // delta tilt bottom to top of box
    double deltaTwistDegrees = box.Twist[1];

    const AABox &bbox = box.getBoundingBox();
    double boxHeight = bbox.getZMax() - bbox.getZMin();
    double power = box.getParam(0, 1.0);

    for (int i = 0; i < npLC; i++) { // loop over each node
        Vec3 p = coordinates.getPoint(i);
        double S = dir[i].S();

        if (box.contains(p)) {
            double pzn = (p.z() - bbox.getZMin()) / (boxHeight);
            double twistDegrees = bottomTwistDegrees + pow(pzn * deltaTwistDegrees, power);
            double tiltDegrees = bottomTiltDegrees + pow(pzn * deltaTiltDegrees, power);
            dir[i] = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, S);
        }
    } // end loop over each node
} //end void setNormalBox

void setRandomBox(const Box &box, std::vector<qlc3d::Director> &dir, const Coordinates &coordinates, int npLC){
/*! Sets the Q-tensor initial configuration within a box volume. The director orientation is randomized */
    srand(0); // seed with constant value so that results are repeatable. TODO: make seed user defined configuration
    for (int i = 0 ; i < npLC ; i ++) {
        double S = dir[i].S();
        Vec3 p = coordinates.getPoint(i);
        if (box.contains(p)) {
                double r1 =  (double) ( rand() % 10000 ) - 5000.0;
                double r2 =  (double) ( rand() % 10000 ) - 5000.0;
                double r3 =  (double) ( rand() % 10000 ) - 5000.0;
                double len = sqrt(r1*r1 + r2*r2 + r3*r3);
                dir[i] = qlc3d::Director(r1 / len, r2 / len, r3 / len, S);
        }
    }
}
/*!
 * Sets the Q-tensor/director initial orientation within a box volume. The orientation is set
 * so that a single hedgehog (+1) defect is located at the centre of the box
*/
void setHedgehogBox(Box box, std::vector<qlc3d::Director> &dir, const Coordinates &coordinates, int npLC) {
    // Director componenets are set to equal vectors from box cetre to director node location
    // resulting in a hedgehog defect
    // Calculate centre coordinates of this box
    Vec3 centre = box.centroid();
    for (unsigned int i = 0 ; i < (unsigned int) npLC ; i++){ // loop over all LC nodes
        double S = dir[i].S();
        Vec3 pl = coordinates.getPoint(i);

        if (box.contains(pl) ) { // if node i is inside box
            // calculate vector from box centre to this node
            Vec3 v = pl - centre;
            dir[i] = qlc3d::Director(v.normalized(), S);
        }
    }
}
/*!
 * Sets initial Q-tensor volume configuration for all boxes
 */
void setVolumeQ(
        SolutionVector &q,
        double S0,
        const Boxes &boxes,
        const Coordinates & coordinates) {
    Log::info("Setting initial LC configuration for {} boxes.", boxes.getBoxCount());

    assert(q.getnDimensions() == 5);
    assert(q.getnDoF() > 0);

    int npLC = q.getnDoF() ;
    // LC TILT AND TWIST IS FIRST CALCULATED AS VECTORS
    // AFTER THIS, A "T-TENSOR" REPRESENTATION IS THEN CALCULATED FROM THE VECTORS

    // By default, when no boxes exist, the director is initialised along (1, 0, 0) at equilibrium order.
    auto defaultDirector = qlc3d::Director(1, 0, 0, S0);
    std::vector<qlc3d::Director> dir(npLC, defaultDirector);

    // override the director within each box
    for (size_t i = 0; i < boxes.getBoxCount(); i++) {
        const Box &b = boxes.getBox(i);
        Log::info("BOX{}:{}.", i + 1, b.toString());
        switch (b.getType()) {
            case Normal:
                setNormalBox(b, dir, coordinates, npLC);
                break;
            case Random:
                setRandomBox(b, dir, coordinates, npLC);
                break;
            case Hedgehog:
                setHedgehogBox(b, dir, coordinates, npLC);
                break;
            default:
                throw std::invalid_argument("unsupported box type " + b.getTypeString());
        }
    }

    // convert the director to tensor
    for (int i = 0; i < npLC; i++) {
        q.setValue(i, qlc3d::TTensor::fromDirector(dir[i]));
    }
}
// end void setVolumeQ
