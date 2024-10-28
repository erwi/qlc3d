#include <box.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <stringenum.h>
#include <settings_file_keys.h>
#include <cassert>
#include <geom/vec3.h>
#include "util/exception.h"

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

//===================================================================
Boxes::Boxes() {
}

void Boxes::addBox(Box *b) {
    box.push_back(std::unique_ptr<Box>(b));
}

void Boxes::addBox(const int &boxNum,
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
    this->addBox(b);
}


