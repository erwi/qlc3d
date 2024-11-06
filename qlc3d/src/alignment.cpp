#include <alignment.h>
#include <util/exception.h>
#include <geometry.h>
#include <geom/vec3.h>
#include <geom/coordinates.h>
#include <lc-representation.h>
#include <solutionvector.h>
#include <util/logging.h>
#include <expression.h>

const std::vector<std::string> Surface::VALID_ANCHORING_TYPES = {"Strong", "Weak", "Homeotropic",
                                                       "Degenerate", "Freeze", "Polymerise",
                                                       "ManualNodes", "WeakHomeotropic"};
const std::string Surface::DEFAULT_ANCHORING_TYPE = Surface::VALID_ANCHORING_TYPES[0];
const double Surface::DEFAULT_ANCHORING_STRENGTH = 1e-4;
const double Surface::DEFAULT_ANCHORING_K1 = 1;
const double Surface::DEFAULT_ANCHORING_K2 = 1;
const bool Surface::DEFAULT_ANCHORING_OVERRIDE_VOLUME = true;
const std::vector<double> Surface::DEFAULT_ANCHORING_EASY = {0,0,0};
const std::vector<double> Surface::DEFAULT_ANCHORING_PARAMS = {};

// <editor-fold desc="Surface">
Surface::Surface(AnchoringType anchoringType,
                 double strength, double k1, double k2,
                 const double easyAnglesDegrees[3],
                 bool overrideVolume,
                 unsigned int fixLcNumber)
    : Type(anchoringType), Strength(strength), K1(k1), K2(k2),
    easyAnglesDegrees{easyAnglesDegrees[0], easyAnglesDegrees[1], easyAnglesDegrees[2]},

    v1(calculateV1(easyAnglesDegrees[0], easyAnglesDegrees[1], easyAnglesDegrees[2])),
    v2(calculateV2(easyAnglesDegrees[0], easyAnglesDegrees[1], easyAnglesDegrees[2])),
    e(v1.cross(v2)), overrideVolume(overrideVolume), fixLcNumber(fixLcNumber) {
}

Vec3 Surface::calculateV1(double tiltDegrees, double twistDegrees, double rotDegrees) {
  double a = -twistDegrees * PI / 180.0; // twist
  double b = tiltDegrees * PI / 180.0; // tilt
  double g = rotDegrees * PI / 180.0; // rotation around

  return {
    sin(a)*cos(g) + cos(a)*sin(b)*sin(g),
    cos(a)*cos(g) - sin(a)*sin(b)*sin(g),
    -cos(b)*sin(g)
  };
}

Vec3 Surface::calculateV2(double tiltDegrees, double twistDegrees, double rotDegrees) {
  double a = -twistDegrees * PI / 180.0; // twist
  double b = tiltDegrees * PI / 180.0; // tilt
  double g = rotDegrees * PI / 180.0; // rotation around

  return {
    sin(a) * sin(g) - cos(a) * sin(b) * cos(g),
    cos(a) * sin(g) + sin(a) * sin(b) * cos(g),
    cos(b) * cos(g)
  };
}


std::string Surface::getAnchoringTypeName() const {
  switch (getAnchoringType()) {
    case AnchoringType::Strong:
      return "Strong";
    case AnchoringType::Weak:
      return "Weak";
    case AnchoringType::Homeotropic:
      return "Homeotropic";
    case AnchoringType::Degenerate:
      return "Degenerate";
    case AnchoringType::Freeze:
      return "Freeze";
    case AnchoringType::Polymerise:
      return "Polymerise";
    case AnchoringType::ManualNodes:
      return "ManualNodes";
    case AnchoringType::WeakHomeotropic:
      return "WeakHomeotropic";
    default:
      RUNTIME_ERROR("Unhandled anchoring type " + std::to_string((int) getAnchoringType()) + ".");
  }

    return Surface::VALID_ANCHORING_TYPES[this->Type];
}

std::string Surface::toString() const {

  std::string easyTilt = fmt::format("{}", easyAnglesDegrees[0]);
  if (tiltExpression_.has_value()) {
    easyTilt = tiltExpression_.value().getExpression();
  }

  std::string easyTwist = fmt::format("{}", easyAnglesDegrees[1]);
  if (twistExpression_.has_value()) {
    easyTwist = twistExpression_.value().getExpression();
  }

  return fmt::format("FIXLC{}, Anchoring:{}, Strength:{}, K1:{}, K2:{},"
                     " EasyAngles:[{}, {}, {}], v1:[{}], v2:[{}], e:[{}]",
                     getFixLCNumber(), getAnchoringTypeName(), getStrength(), getK1(), getK2(),
                     easyTilt, easyTwist, easyAnglesDegrees[2],
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

bool	Surface::usesSurfaceNormal() const {
  AnchoringType type = getAnchoringType();
  return type == AnchoringType::Degenerate || type == AnchoringType::Homeotropic || type == AnchoringType::WeakHomeotropic;
}
bool    Surface::isStrong() const {
  AnchoringType type = getAnchoringType();
  return type != AnchoringType::Weak && type != AnchoringType::Degenerate && type != AnchoringType::WeakHomeotropic;
}


void Surface::setManualNodesAnchoring(SolutionVector &q, double S0) const {
  double tilt = getEasyTilt();
  double twist = getEasyTwist();

  // node indices are stored in Params as doubles, convert to idx
  std::set<idx> nodes_idx;
  for (double i : Params) {
    if ((i < 0 )){
      RUNTIME_ERROR("Negative node index when setting ManualNodesAnchoring.");
    }
    idx nodeIdx = static_cast<idx>(i);
    nodes_idx.insert(nodeIdx);
  }
  auto director = qlc3d::Director::fromDegreeAngles(tilt, twist, S0);

  for (auto &i : nodes_idx) {
    q.setValue(i, director);
  }
}

double Surface::getEasyTiltAngleAt(const Vec3 &p) const {
  if (tiltExpression_.has_value()) {
    return tiltExpression_.value().evaluate(p);
  } else {
    return getEasyTilt();
  }
}

double Surface::getEasyTwistAngleAt(const Vec3 &p) const {
  if (twistExpression_.has_value()) {
    return twistExpression_->evaluate(p);
  } else {
    return getEasyTwist();
  }
}

Vec3 Surface::getEasyDirectionAt(const Vec3 &p) const {
  double tilt = getEasyTiltAngleAt(p);
  double twist = getEasyTwistAngleAt(p);
  return Vec3::fromDegreeAngles(tilt, twist);
}

void Surface::setFromTiltAndTwistAngles(SolutionVector &q, double S0, const Geometry &geom) const {
  std::set<idx> indSurfaceNodes = geom.getTriangles().listFixLCSurfaceNodes(getFixLCNumber());
  const Coordinates &coordinates = geom.getCoordinates();

  for (auto &i : indSurfaceNodes) {
    auto &p = coordinates.getPoint(i);
    double tilt = getEasyTiltAngleAt(p);
    double twist = getEasyTwistAngleAt(p);
    auto director = qlc3d::Director::fromDegreeAngles(tilt, twist, S0);
    q.setValue(i, director);
  }
}

void Surface::setHomeotropicOrientation(SolutionVector &q, double S0, const Geometry &geom) const {
  auto indSurface = geom.getTriangles().listFixLCSurfaceNodes(getFixLCNumber());
  for (auto &i : indSurface) {
    Vec3 normal = geom.getNodeNormal(i);
    qlc3d::Director n(normal, S0);
    q.setValue(i, n);
  }
}

void Surface::setTiltAngleExpression(const std::string &tiltExpression) {
  tiltExpression_.emplace(tiltExpression);
}

void Surface::setTwistAngleExpression(const std::string &twistExpression) {
  twistExpression_.emplace(twistExpression);
}

void Surface::setAlignmentOrientation(SolutionVector &q, double S0, const Geometry &geom) const {
  Log::info("Setting initial LC configuration for FIXLC{} = {}.", getFixLCNumber(), toString());
  if (!getOverrideVolume()) {
    Log::info("not setting initial orientation for FixLC{} because overrideVolume is false", getFixLCNumber());
    return;
  }

  if (getAnchoringType() == ManualNodes) {
    setManualNodesAnchoring(q, S0);
  } else if (getAnchoringType() == Strong || getAnchoringType() == Weak ||
             (getAnchoringType() == Degenerate && getStrength() >= 0)) {
    setFromTiltAndTwistAngles(q, S0, geom);
  } else if (getAnchoringType() == Homeotropic ||
             (getAnchoringType() == Degenerate && getStrength() < 0) ||
             (getAnchoringType() == WeakHomeotropic)) {
    setHomeotropicOrientation(q, S0, geom);
    return;
  } else if (getAnchoringType() == Freeze || getAnchoringType() == Polymerise) {
    Log::info("Using current bulk q-tensor values for FIXLC{}.", getFixLCNumber());
  } else {
    RUNTIME_ERROR("Unhandled anchoring type " + getAnchoringTypeName() + ".");
  }
}

Surface Surface::ofStrongAnchoring(unsigned int fixLcNumber_, double tiltDegrees, double twistDegrees) {
  if (fixLcNumber_ == 0 || fixLcNumber_ > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber_));
  }
  double easyAnglesDegrees[3] = {tiltDegrees, twistDegrees, 0};
  return {AnchoringType::Strong,
          std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(),
          easyAnglesDegrees,
          true,
          fixLcNumber_};
}

Surface Surface::ofStrongAnchoring(unsigned int fixLcNumber, const std::string &tiltExpression,
                                   const std::string &twistExpression) {
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  double easyAnglesDegrees[3] = {0, 0, 0};
  Surface surf = {AnchoringType::Strong,
          std::numeric_limits<double>::infinity(), 1, 1,
          easyAnglesDegrees,
          true,
          fixLcNumber};
  surf.setTiltAngleExpression(tiltExpression);
  surf.setTwistAngleExpression(twistExpression);
  // set v1, v2, e to all NaN, since they may vary with position so are invalid anyway
  surf.v1.set(std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN());
  surf.v2.set(std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN());
  surf.e.set(std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN());
  return surf;
}

Surface Surface::ofPlanarDegenerate(unsigned int fixLcNumber, double strength) {
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  double easyAnglesDegrees[3] = {0,
                                 std::numeric_limits<double>::signaling_NaN(),
                                 std::numeric_limits<double>::signaling_NaN()};
  // For degenerate anchoring, the surface normal (v2) repels the LC director
  return {AnchoringType::Degenerate,
          strength, 0, 1,
          easyAnglesDegrees, // actually calculated for each node from surface normal
          false,
          fixLcNumber};
}

Surface Surface::ofStrongHomeotropic(unsigned int fixLcNumber) {
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  double easyAnglesDegrees[3] = {90,
                                 std::numeric_limits<double>::signaling_NaN(),
                                 std::numeric_limits<double>::signaling_NaN()};

  return {AnchoringType::Homeotropic,
          std::numeric_limits<double>::infinity(), 0, 1,
          easyAnglesDegrees, // actually calculated for each node from surface normal
          true,
          fixLcNumber};
}

Surface Surface::ofWeakHomeotropic(unsigned int fixLCNumber, double strength) {
  if (fixLCNumber == 0 || fixLCNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLCNumber));
  }
  double easyAnglesDegrees[3] = {90,
                                 std::numeric_limits<double>::signaling_NaN(),
                                 std::numeric_limits<double>::signaling_NaN()};

  return {AnchoringType::WeakHomeotropic,
          strength, 0, 1,
          easyAnglesDegrees,
          true,
          fixLCNumber};
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
          true,
          fixLcNumber
  };
}

Surface Surface::ofFreeze(unsigned int fixLcNumber) {
  // essentially this is a strong anchoring surface, but it does not override any previous LC orientation
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  const auto nan = std::numeric_limits<double>::signaling_NaN();
  double easyAnglesDegrees[3] = {nan, nan, nan};
  Surface surf = {AnchoringType::Freeze,
                  std::numeric_limits<double>::infinity(),
                  nan, nan,
                  easyAnglesDegrees,
                  false,
                  fixLcNumber};
  // set v1, v2, e to all NaN, since they may vary with position so are invalid anyway
  surf.v1.set(nan, nan, nan);
  surf.v2.set(nan, nan, nan);
  surf.e.set(nan, nan, nan);
  return surf;
}

// </editor-fold>

// <editor-fold desc="Alignment">
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
    RUNTIME_ERROR(fmt::format("No such surface {}. Surfaces size={}", i, surface.size()));
  }
  return surface[i];
}

//Surface& Alignment::getSurface(const idx &i) {
//  return const_cast<Surface&>(static_cast<const Alignment&>(*this).getSurface(i));
//}

int Alignment::getnSurfaces() const {	return n_surfaces;}

bool Alignment::isStrong(int n) const {
    if (n >= (int) surface.size() ) {
        RUNTIME_ERROR("Invalid alignment surface FIXLC" + std::to_string(n + 1) + ", number of alignment surfaces is " + std::to_string(surface.size()));
    }
    return getSurface(n).isStrong();
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
  return surface[n - 1].usesSurfaceNormal();
}

bool Alignment::hasSurface(unsigned int fixLcNumber) const {
    for (const auto &s : surface) {
        if (s.getFixLCNumber() == fixLcNumber) {
            return true;
        }
    }
    return false;
}

std::unordered_map<unsigned int, Surface> Alignment::getWeakSurfacesByFixLcNumber() const {
  std::unordered_map<unsigned int, Surface> weakSurfaces;
  for (const auto &s : surface) {
    if (!s.isStrong()) {
      unsigned int fixLcNumber = s.getFixLCNumber();
      weakSurfaces.insert({fixLcNumber, s});
    }
  }
  return weakSurfaces;
}

// </editor-fold>