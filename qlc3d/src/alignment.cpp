#include <alignment.h>
#include <algorithm>
#include <iostream>
#include <stringenum.h>
#include <settings_file_keys.h>
#include <util/exception.h>

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

Surface::Surface(int fxlcnum, const std::string &type) {
  fixLcNumber = fxlcnum;
  Strength = DEFAULT_ANCHORING_STRENGTH;
  K1 = DEFAULT_ANCHORING_K1;
  K2 = DEFAULT_ANCHORING_K2;
  setAnchoringType(type);
  this->setEasyAngles(DEFAULT_ANCHORING_EASY);
  UsesSurfaceNormal = false;
  isFixed = true;
  overrideVolume = DEFAULT_ANCHORING_OVERRIDE_VOLUME;
}

void Surface::setAnchoringType(const std::string &atype) {
    std::string typeKey = wildcardToNum(SFK_FIXLC_ANCHORING, getFixLCNumber());
    StringEnum<AnchoringType> validator(typeKey, Surface::VALID_ANCHORING_TYPES);
    this->Type = validator.getEnumValue(atype);
    switch (Type) {
        case Strong :
            isFixed = true;
            UsesSurfaceNormal = false;
            break;
        case Homeotropic :
            isFixed = true;
            UsesSurfaceNormal = true;
            break;
        case Weak :
            isFixed = false;
            UsesSurfaceNormal = false;
            break;
        case Degenerate :
            isFixed = false;
            UsesSurfaceNormal = true;
            break;
        case Freeze :
            isFixed = true;
            UsesSurfaceNormal = false;
            break;
        case ManualNodes :
            isFixed = true;
            UsesSurfaceNormal = false;
            break;
        case Polymerise:
            isFixed = true;
            UsesSurfaceNormal = false;
            break;
        default:
            throw std::runtime_error("unhandled anchoring type in " + std::string(__PRETTY_FUNCTION__));
    }
}// end setAnchoringType
void Surface::setStrength(double str) {	Strength = str; }
void Surface::setK1(double k1) { K1 = k1; }
void Surface::setK2(double k2) { K2 = k2; }

/**
* set easy direction, tilt, twist & rotation angles in degrees
*/
void Surface::setEasyAngles(const std::vector<double> &e) {
    if (e.size() > 3) {
        std::string errorMsg = std::string("error in Surface::setEasyAngles, too many easy angles,") +
            " got " + std::to_string(e.size()) + ", but can only handle up to 3 (tilt, twist, rotation)";
        throw std::invalid_argument(errorMsg);
    }
  easyAnglesDegrees[0] = 0.0; easyAnglesDegrees[1] = 0.0; easyAnglesDegrees[2] = 0.0;
    for (size_t i = 0; i < e.size() ; i++) {
      easyAnglesDegrees[i] = e[i];
    }

    this->calcV1V2(); // calculates primary anchoring axes from easy angles
}
void Surface::setv1( double v[3] ) { v1[0] = v[0]; v1[1] = v[1] ; v1[2] = v[2]; }
void Surface::setv2( double v[3] ) { v2[0] = v[0]; v2[1] = v[1] ; v2[2] = v[2]; }
void Surface::setEasyVector( double v[3]){	e[0] = v[0]; e[1] = v[1]; e[2] = v[2]; }

void Surface::calcV1V2(){
    /*! Calculates v1 and v2 vectors given tilt and twist angles
 *  Rotation matrices are given in Willman, IEEE Trans. Electron Dev. 54, 10, 2007 
 * */

    if (this->Type != Weak) // vectors are only set for 'Weak' anchoring type (why?)
        return;
    double a = easyAnglesDegrees[1] * PI / 180.0; // twist
    double b = easyAnglesDegrees[0] * PI / 180.0; // tilt
    double g = easyAnglesDegrees[2] * PI / 180.0; // rotation around
    double k[3] = {0.0, 0.0, 0.0};
    double l[3] = {0.0, 0.0, 0.0};
    // apply rotation matrices
    k[0] = sin(a)*cos(g) + cos(a)*sin(b)*sin(g);
    k[1] = cos(a)*cos(g) - sin(a)*sin(b)*sin(g);
    k[2] = -sin(b)*sin(g);

    l[0] = sin(a)*sin(g) - cos(a)*sin(b)*cos(g);
    l[1] = cos(a)*sin(g) + sin(a)*sin(b)*cos(g);
    l[2] = cos(b)*cos(g);
    // v1 is surface normal = (0,0,1) for flat bottom surface
    this->setv1( l );
    this->setv2( k );
    // calculate easy direction vector e = v1 x v2
    double e[3] = {0,0,0};
    e[0] =  v1[1]*v2[2] - v2[1]*v1[2];  // OUT OF BOUNDS!!!!
    e[1] = -v1[0]*v2[2] + v2[0]*v1[2];
    e[2] =  v1[0]*v2[1] - v2[0]*v1[1];
    this->setEasyVector( e );
}

std::string Surface::getAnchoringTypeName() const {
    return Surface::VALID_ANCHORING_TYPES[this->Type];
}

std::string Surface::toString() const {
  return fmt::format("FIXLC{}, Anchoring:{}, Strength:{}, K1:{}, K2:{},"
                     " EasyAngles:[{},{},{}], v1:[{},{},{}], v2:[{},{},{}]",
                     getFixLCNumber(), getAnchoringTypeName(), getStrength(), getK1(), getK2(),
                     easyAnglesDegrees[0], easyAnglesDegrees[1], easyAnglesDegrees[2],
                     v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
}

AnchoringType Surface::getAnchoringType() const {
    return this->Type;
}

double Surface::getStrength() const		{		return Strength;}
double Surface::getK1() const				{		return K1;}
double Surface::getK2() const				{		return K2;}
double Surface::getEasyTilt() const{			return easyAnglesDegrees[0];}
double Surface::getEasyTwist() const{			return easyAnglesDegrees[1];}
double* Surface::getPtrTov1(){			return &v1[0];}
double* Surface::getPtrTov2(){			return &v2[0];}
bool	Surface::getUsesSurfaceNormal() const {
    return UsesSurfaceNormal;
}
bool    Surface::isStrong() const {
    return isFixed;
}

Surface* Surface::ofStrongAnchoring(unsigned int fixLcNumber_, double tiltDegrees, double twistDegrees) {
  if (fixLcNumber_ == 0 || fixLcNumber_ > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber_));
  }
  auto *s = new Surface((int) fixLcNumber_, "Strong");
  s->setAnchoringType("Strong");
  s->setStrength(0);
  s->setEasyAngles({tiltDegrees, twistDegrees, 0});
  s->setK1(1.);
  s->setK2(1.);
  return s;
}

Surface* Surface::ofPlanarDegenerate(unsigned int fixLcNumber, double strength) {
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  auto *s = new Surface((int) fixLcNumber, "Degenerate");
  s->setStrength(strength);
  s->setEasyAngles({0, 0, 0});
  s->setK1(1.);
  s->setK2(1.);
  return s;
}

Surface* Surface::ofHomeotropic(unsigned int fixLcNumber) {
  if (fixLcNumber == 0 || fixLcNumber > 9) {
    RUNTIME_ERROR(fmt::format("FixLC number must be in range 1 - 9, got {}", fixLcNumber));
  }
  auto *s = new Surface((int) fixLcNumber, "Homeotropic");

  s->setStrength(0);
  s->setEasyAngles({0, 0, 0});
  s->setK1(1.);
  s->setK2(1.);
  return s;
}

//====================================================
//
//		Alignment
//
//=====================================================

Alignment::Alignment(){	setnSurfaces(0);}//end constructor
Alignment::~Alignment(){
    std::vector<Surface*>::iterator itr;
    for(itr = surface.begin() ; itr!= surface.end() ; ++itr){
        delete (*itr);
    }
}

void Alignment::setnSurfaces(int n){	n_surfaces = n;}

void Alignment::addSurface(Surface* s){
  // check against adding multiple surfaces with the same FixLC number
  if (hasSurface(s->getFixLCNumber())) {
    throw std::invalid_argument("Surface with FixLC number " + std::to_string(s->getFixLCNumber()) + " already exists");
  }

  surface.push_back(s);
  n_surfaces++;
}

void Alignment::addSurface(const int fixLcNumber,
                           const std::string &anchoring,
                           const double &strength,
                           const std::vector<double> &easyAnglesDegrees,
                           const double &k1,
                           const double &k2,
                           const std::vector<double> &params,
                           const bool overrideVolume) {

    auto s = new Surface(fixLcNumber, anchoring);
    s->setStrength(strength);
    s->setEasyAngles(easyAnglesDegrees);
    s->setK1(k1);
    s->setK2(k2);
    s->setEnforce(overrideVolume);
    s->Params = params;
    this->addSurface(s);
}

const Surface& Alignment::getSurface(const idx &i) const {
    if (i >= (idx) this->surface.size()) {
        throw std::invalid_argument("no such surface " + std::to_string(i) +
        " in " + __PRETTY_FUNCTION__ );
    }

    return *this->surface[i];
}


int Alignment::getnSurfaces() const {	return n_surfaces;}

bool Alignment::IsStrong(int i) const {
    if (i >= (int) surface.size() ) {
        RUNTIME_ERROR("Invalid alignment surface FIXLC" + std::to_string(i + 1) + ", number of alignment surfaces is " + std::to_string(surface.size()));
    }
    return getSurface(i).isStrong();
}

bool Alignment::WeakSurfacesExist(){
/*!
returns true if any non-fixed surfaces have been defined
*/
    for(int n = 0 ; n < getnSurfaces() ; n++ )
        if ( !IsStrong(n) )
            return true;
    return false;
}

double Alignment::getStrength(int n) 	{return surface[n-1]->getStrength();}
double Alignment::getK1(int n)			{return surface[n-1]->getK1();}
double Alignment::getK2(int n)			{return surface[n-1]->getK2();}
double* Alignment::getPtrTov1(int n)		{return surface[n-1]->getPtrTov1();}
double* Alignment::getPtrTov2(int n)		{return surface[n-1]->getPtrTov2();}
bool Alignment::getUsesSurfaceNormal(int n)  {return surface[n-1]->getUsesSurfaceNormal();}

bool Alignment::hasSurface(unsigned int fixLcNumber) const {
    for (const auto &s : surface) {
        if (s->getFixLCNumber() == fixLcNumber) {
            return true;
        }
    }
    return false;
}
