#include <alignment.h>
#include <algorithm>
#include <iostream>
#include <stringenum.h>
#include <settings_file_keys.h>
#include <util/exception.h>
#include <geom/vec3.h>

const std::vector<std::string> Surface::VALID_ANCHORING_TYPES = {"Strong", "Weak", "Homeotropic",
                                                       "Degenerate", "Freeze", "Polymerise",
                                                       "ManualNodes"};
const std::string Surface::DEFAULT_ANCHORING_TYPE = Surface::VALID_ANCHORING_TYPES[0];
const double Surface::DEFAULT_ANCHORING_STRENGTH = 1e-4;
const double Surface::DEFAULT_ANCHORING_K1 = 1;
const double Surface::DEFAULT_ANCHORING_K2 = 1;
const bool Surface::DEFAULT_ANCHORING_OVERRIDE_VOLUME = true;
const std::vector<double> Surface::DEFAULT_ANCHORING_EASY = {0,0,0};
const std::vector<double> Surface::DEFAULT_ANCHORING_PARAMS = {};

Surface::Surface(AnchoringType anchoringType,
                 double strength, double k1, double k2,
                 const double easyAnglesDegrees[3],
                 //const Vec3 &v1, const Vec3 &v2,
                 bool usesSurfaceNormal,
                 bool isFixed,
                 bool overrideVolume,
                 unsigned int fixLcNumber)
    : Type(anchoringType), Strength(strength), K1(k1), K2(k2),
    easyAnglesDegrees{easyAnglesDegrees[0], easyAnglesDegrees[1], easyAnglesDegrees[2]},

    v1(calculateV1(easyAnglesDegrees[0], easyAnglesDegrees[1], easyAnglesDegrees[2])),
    v2(calculateV2(easyAnglesDegrees[0], easyAnglesDegrees[1], easyAnglesDegrees[2])),
    e(v1.cross(v2)),
    usesSurfaceNormal(usesSurfaceNormal), isFixed(isFixed), overrideVolume(overrideVolume), fixLcNumber(fixLcNumber) {

}

Vec3 Surface::calculateV1(double tiltDegrees, double twistDecrees, double rotDegrees) {
  double a = twistDecrees * PI / 180.0; // twist
  double b = tiltDegrees * PI / 180.0; // tilt
  double g = rotDegrees * PI / 180.0; // rotation around

  return {
          sin(a) * sin(g) - cos(a) * sin(b) * cos(g),
          cos(a)*sin(g) + sin(a)*sin(b)*cos(g),
          cos(b)*cos(g)
  };
}

Vec3 Surface::calculateV2(double tiltDegrees, double twistDegrees, double rotDegrees) {
  double a = twistDegrees * PI / 180.0; // twist
  double b = tiltDegrees * PI / 180.0; // tilt
  double g = rotDegrees * PI / 180.0; // rotation around

  return {
          sin(a)*cos(g) + cos(a)*sin(b)*sin(g),
          cos(a)*cos(g) - sin(a)*sin(b)*sin(g),
          -sin(b)*sin(g),
  };
}


std::string Surface::getAnchoringTypeName() const {
    return Surface::VALID_ANCHORING_TYPES[this->Type];
}

std::string Surface::toString() const {
  return fmt::format("FIXLC{}, Anchoring:{}, Strength:{}, K1:{}, K2:{},"
                     " EasyAngles:[{},{},{}], v1:[{}], v2:[{}], e:[{}]",
                     getFixLCNumber(), getAnchoringTypeName(), getStrength(), getK1(), getK2(),
                     easyAnglesDegrees[0], easyAnglesDegrees[1], easyAnglesDegrees[2],
                     v1, v2, e);
}

AnchoringType Surface::getAnchoringType() const {
    return this->Type;
}

double Surface::getStrength() const		{		return Strength;}
double Surface::getK1() const				{		return K1;}
double Surface::getK2() const				{		return K2;}
double Surface::getEasyTilt() const{			return easyAnglesDegrees[0];}
double Surface::getEasyTwist() const{			return easyAnglesDegrees[1];}

bool	Surface::getUsesSurfaceNormal() const {
    return usesSurfaceNormal;
}
bool    Surface::isStrong() const {
    return isFixed;
}

Surface Surface::ofStrongAnchoring(unsigned int fixLcNumber_, double tiltDegrees, double twistDegrees) {
  if (fixLcNumber_ == 0 || fixLcNumber_ > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber_));
  }
  double easyAnglesDegrees[3] = {tiltDegrees, twistDegrees, 0};
  return {AnchoringType::Strong,
          DEFAULT_ANCHORING_STRENGTH, 1, 1,
          easyAnglesDegrees,
          false, true, true,
          fixLcNumber_};
}

Surface Surface::ofPlanarDegenerate(unsigned int fixLcNumber, double strength) {
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  double easyAnglesDegrees[3] = {0, 0, 0};

  return {AnchoringType::Degenerate,
          strength, 1, 1,
          easyAnglesDegrees, // actually calculated for each node from surface normal
          true, false, true, fixLcNumber};
}

Surface Surface::ofHomeotropic(unsigned int fixLcNumber) {
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  double easyAnglesDegrees[3] = {0, 0, 0};

  return {AnchoringType::Homeotropic,
          DEFAULT_ANCHORING_STRENGTH, 1, 1,
          easyAnglesDegrees, // actually calculated for each node from surface normal
          true,
          true,
          true,
          fixLcNumber};
}

Surface Surface::ofWeakAnchoring(unsigned int fixLcNumber, double tiltDegrees, double twistDegrees, double strength,
                                 double k1, double k2) {
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  double easyAnglesDegrees[3] = {tiltDegrees, twistDegrees, 0};
  return {AnchoringType::Weak,
          strength, k1, k2,
          easyAnglesDegrees,
          false, false, true, fixLcNumber

  };
}

//====================================================
//
//		Alignment
//
//=====================================================

Alignment::Alignment(){	setnSurfaces(0);}//end constructor

void Alignment::setnSurfaces(int n){	n_surfaces = n;}

void Alignment::addSurface(Surface s){
  // check against adding multiple surfaces with the same FixLC number
  if (hasSurface(s.getFixLCNumber())) {
    throw std::invalid_argument("Surface with FixLC number " + std::to_string(s.getFixLCNumber()) + " already exists");
  }

  surface.push_back(s);
  n_surfaces++;
}

const Surface& Alignment::getSurface(const idx &i) const {
    if (i >= (idx) surface.size()) {
        throw std::invalid_argument("no such surface " + std::to_string(i) +
        " in " + __PRETTY_FUNCTION__ );
    }

    return surface[i];
}


int Alignment::getnSurfaces() const {	return n_surfaces;}

bool Alignment::IsStrong(int i) const {
    if (i >= (int) surface.size() ) {
        RUNTIME_ERROR("Invalid alignment surface FIXLC" + std::to_string(i + 1) + ", number of alignment surfaces is " + std::to_string(surface.size()));
    }
    return getSurface(i).isStrong();
}

bool Alignment::weakSurfacesExist() const {
  for (const auto &s : surface) {
    if (!s.isStrong()) {
      return true;
    }
  }
  return false;
}

double Alignment::getStrength(int n) { return surface[n-1].getStrength(); }
double Alignment::getK1(int n) { return surface[n-1].getK1(); }
double Alignment::getK2(int n) { return surface[n-1].getK2(); }

bool Alignment::getUsesSurfaceNormal(int n) const {
  return surface[n-1].getUsesSurfaceNormal();
}

bool Alignment::hasSurface(unsigned int fixLcNumber) const {
    for (const auto &s : surface) {
        if (s.getFixLCNumber() == fixLcNumber) {
            return true;
        }
    }
    return false;
}
