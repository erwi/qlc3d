#include <box.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <stringenum.h>
#include <settings_file_keys.h>
#include <cassert>
#include <geom/vec3.h>
#include <lc-representation.h>
#include <util/exception.h>
#include <util/logging.h>
#include <solutionvector.h>
#include <geom/coordinates.h>

const std::vector<std::string> Box::VALID_TYPES = {"Normal", "Random", "Hedgehog"};
const std::string Box::DEFAULT_TYPE = Box::VALID_TYPES[0];
const std::vector<double> Box::DEFAULT_PARAMS = {};
const std::vector<double> Box::DEFAULT_X_Y_Z = {0,0};
const std::vector<double> Box::DEFAULT_TILT_TWIST = {0,0};

using std::vector;
using std::cerr;
Box::Box(int boxnum) {
    BoxNumber = boxnum;
    Params = Box::DEFAULT_PARAMS;
    Tilt = Box::DEFAULT_TILT_TWIST;
    Twist = Box::DEFAULT_TILT_TWIST;
    setBoxType(Box::DEFAULT_TYPE);
    boundingBox = AABox(DEFAULT_X_Y_Z[0], DEFAULT_X_Y_Z[1], DEFAULT_X_Y_Z[0], DEFAULT_X_Y_Z[1], DEFAULT_X_Y_Z[0], DEFAULT_X_Y_Z[1]);
}

void Box::setX(std:: vector<double> x) {
  if (x.size() != 2) {
    throw new std::invalid_argument("Invalid X length " + std::to_string(x.size()) + " expected 2");
  }
  boundingBox.setBoundsX(x[0], x[1]);
}

void Box::setY(std:: vector<double> y) {
  if (y.size() != 2) {
    throw new std::invalid_argument("Invalid Y length " + std::to_string(y.size()) + " expected 2");
  }
  boundingBox.setBoundsY(y[0], y[1]);
}

void Box::setZ(std:: vector<double> z) {
  if (z.size() != 2) {
    throw new std::invalid_argument("Invalid Z length " + std::to_string(z.size()) + " expected 2");
  }
  boundingBox.setBoundsZ(z[0], z[1]);
}

void Box::setTilt(std:: vector<double> tlt) {
    if (tlt.size() == 2) {
        Tilt[0] = tlt[0];
        Tilt[1] = tlt[1];
    } else {
        std::cerr << "error, Box::setTilt, invalid Tilt length - bye!" << std::endl;
        exit(1);
    }
}
void Box::setTwist(std:: vector<double> twt) {
    if (twt.size() == 2) {
        Twist[0] = twt[0];
        Twist[1] = twt[1];
    } else {
        std::cerr << "error, Box::setTwist, invalid Twist length - bye!" << std::endl;
        exit(1);
    }
}

void Box::setRandomGenerator(RandomGenerator &rg) {
  randomGenerator = rg;
}

bool Box::contains(double *coords) const {
    return boundingBox.contains(coords[0], coords[1], coords[2]);
}

Vec3 Box::centroid() const {
  return boundingBox.center();
}

bool Box::contains(const Vec3 &p) const {
  double coords[3] = {p.x(), p.y(), p.z()};
  return this->contains(coords);
}

void Box::setBoxType(const std::string &bt) {
/*!Sets the current box type from type name string.*/

    std::string typeKey = wildcardToNum(SFK_BOX_TYPE, this->BoxNumber);
    StringEnum<BoxType> validator(typeKey, Box::VALID_TYPES);
    try {
        type = validator.getEnumValue(bt);
    } catch (...) {
        validator.printErrorMessage(bt);
        RUNTIME_ERROR("Error setting box type to " + bt);
    }
}

void Box::setParams(const std::vector<double> &params) {
    assert(Params.empty());
    std::copy(params.begin(), params.end(), std::back_inserter(Params));
}

double Box::getParam(unsigned int i, double defaultValue) const {
    if (i < this->Params.size()) {
        return this->Params[i];
    } else {
        return defaultValue;
    }
}

std::string Box::toString() const {
    return "Type=" + getTypeString() + ", bounds=["
    + std::to_string(boundingBox.getXMin()) + ", " + std::to_string(boundingBox.getXMax()) + ", "
    + std::to_string(boundingBox.getYMin()) + ", " + std::to_string(boundingBox.getYMax()) + ", "
    + std::to_string(boundingBox.getZMin()) + ", " + std::to_string(boundingBox.getZMax()) + "]";
}

qlc3d::Director Box::getDirectorForNormalBox(const Vec3 &p) const {
  double bottomTwistDegrees = Twist[0];
  double bottomTiltDegrees = Tilt[0];
  double deltaTiltDegrees = Tilt[1]; // delta tilt bottom to top of box
  double deltaTwistDegrees = Twist[1];

  double boxHeight = boundingBox.getZMax() - boundingBox.getZMin();
  double power = getParam(0, 1.0);

  double pzn = (p.z() - boundingBox.getZMin()) / (boxHeight);
  double twistDegrees = bottomTwistDegrees + pow(pzn * deltaTwistDegrees, power);
  double tiltDegrees = bottomTiltDegrees + pow(pzn * deltaTiltDegrees, power);
  return qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, 1.0);
}

qlc3d::Director Box::getDirectorForRandomBox(const Vec3 &p) const {
  Vec3 n( (double) ( rand() % 10000 ) - 5000.0,
          (double) ( rand() % 10000 ) - 5000.0,
          (double) ( rand() % 10000 ) - 5000.0);

  return {n, 1.0};
}

qlc3d::Director Box::getDirectorForHedgehogBox(const Vec3 &p) const {
  Vec3 centre = centroid();
  Vec3 n = p - centre;
  n.normalize();

  return {n, 1.0};
}

qlc3d::Director Box::getDirectorAt(const Vec3 &p) const {
  if (!contains(p)) {
    RUNTIME_ERROR("Error, Box::getDirectorAt, point not in box" + p.toString());
  }

  switch(getType()) {
    case Normal:
      return getDirectorForNormalBox(p);
    case Random:
      return getDirectorForRandomBox(p);
    case Hedgehog:
      return getDirectorForHedgehogBox(p);
    default:
      RUNTIME_ERROR("Error, Box::getDirectorAt, unsupported box type " + getTypeString());
  }
}

//===================================================================
InitialVolumeOrientation::InitialVolumeOrientation() {
  randomGenerator = RandomGenerator();
}

InitialVolumeOrientation::InitialVolumeOrientation(RandomGenerator &rg) {
    randomGenerator = rg;
}

void InitialVolumeOrientation::addBox(Box *b) {
    boxRegions.push_back(std::unique_ptr<Box>(b));
}

void InitialVolumeOrientation::addBox(const int &boxNum,
                                      const std::string &boxType,
                                      const std::vector<double> &params,
                                      const std::vector<double> &x,
                                      const std::vector<double> &y,
                                      const std::vector<double> &z,
                                      const std::vector<double> &tilt,
                                      const std::vector<double> &twist) {
    Box *b = new Box(boxNum);
    b->setBoxType(boxType);
    b->setParams(params);
    b->setX(x);
    b->setY(y);
    b->setZ(z);
    b->setTilt(tilt);
    b->setTwist(twist);
    b->setRandomGenerator(randomGenerator);
    this->addBox(b);
}

void InitialVolumeOrientation::setVolumeQ(SolutionVector &q, double S0, const Coordinates &coordinates) const {
  Log::info("Setting initial LC configuration for {} boxes.", boxRegions.size());

  if (q.getnDimensions() != 5) {
    RUNTIME_ERROR("Error, Boxes::setVolumeQ, invalid q dimensions " + std::to_string(q.getnDimensions()));
  }
  if (q.getnDoF() == 0) {
    RUNTIME_ERROR("Error, Boxes::setVolumeQ, invalid q DoF " + std::to_string(q.getnDoF()));
  }

  const unsigned int npLC = q.getnDoF();

  // By default, when no boxes exist, the director is initialised along (1, 0, 0) at equilibrium order.
  // TODO: this could be a user defined configuration: a "background director orientation"
  auto defaultDirector = qlc3d::Director(1, 0, 0, S0);
  std::vector<qlc3d::Director> dir(npLC, defaultDirector);

  for (unsigned int i = 0; i < npLC; i++) {
    auto p = coordinates.getPoint(i);

    for (auto &b : boxRegions) {
      if (b->contains(p)) {
        auto d = b->getDirectorAt(p);
        dir[i] = qlc3d::Director(d.nx(), d.ny(), d.nz(), S0);
      }
    }
  }

  // convert the director to tensor
  for (unsigned int i = 0; i < npLC; i++) {
    q.setValue(i, qlc3d::TTensor::fromDirector(dir[i]));
  }
}


