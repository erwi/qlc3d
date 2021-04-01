#include <alignment.h>
#include <algorithm>
#include <iostream>
#include <reader.h>
#include <stringenum.h>
#include <settings_file_keys.h>
const vector<string> Surface::VALID_ANCHORING_TYPES = {"Strong", "Weak", "Homeotropic",
                                                       "Degenerate", "Freeze", "Polymerise",
                                                       "ManualNodes"};
const string Surface::DEFAULT_ANCHORING_TYPE = Surface::VALID_ANCHORING_TYPES[0];
const double Surface::DEFAULT_ANCHORING_STRENGTH = 1e-4;
const double Surface::DEFAULT_ANCHORING_K1 = 1;
const double Surface::DEFAULT_ANCHORING_K2 = 1;
const vector<double> Surface::DEFAULT_ANCHORING_EASY = {0,0,0};
const vector<double> Surface::DEFAULT_ANCHORING_PARAMS = {};
Surface::Surface(int fxlcnum) {
    FixLCNumber = fxlcnum;
    Strength = DEFAULT_ANCHORING_STRENGTH;
    K1 = 1;
    K2 = 1;
    this->setEasyAngles(DEFAULT_ANCHORING_EASY);
    UsesSurfaceNormal = false;
    isFixed = true;
}

void Surface::setAnchoringType(const std::string &atype){

    std::string typeKey = wildcardToNum(SFK_FIXLC_ANCHORING, this->FixLCNumber);
    StringEnum<AnchoringType> validator(typeKey, Surface::VALID_ANCHORING_TYPES);
    try {
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
            std::cerr << "unhandled anchoring type in " << __PRETTY_FUNCTION__ << std::endl;
            std::exit(1);
        }
    } catch (...) {
        validator.printErrorMessage(atype);
        std::exit(1);
    }
}// end setAnchoringType
void Surface::setStrength(double str){	Strength = str;}
void Surface::setK1(double k1){			K1 = k1;}
void Surface::setK2(double k2){			K2 = k2;}

void Surface::setEasyAngles(const std::vector<double> &e){
    /*!
     * set easy direction, tilt, twist & rotation
     */
    if (e.size()>3){
        std::cerr << "error in Surface::setEasyAngles, too many easy angles" << std::endl;
        std::cerr << "got " << e.size() << " but can only handle up to 3 (tilt, twist, rotation) - bye!" << std::endl;
        exit(1);
    }
    Easy[0] = 0.0; Easy[1] = 0.0; Easy[2] = 0.0;
    for (size_t i = 0; i < e.size() ; i++)
        Easy[i] = e[i];
    //
    this->calcV1V2(); // calculates primary anchoring axes from easy angles
}



void Surface::setv1( double v[3] ){		v1[0] = v[0]; v1[1] = v[1] ; v1[2] = v[2];}
void Surface::setv2( double v[3] ){		v2[0] = v[0]; v2[1] = v[1] ; v2[2] = v[2];}
void Surface::setEasyVector( double v[3]){	e[0] = v[0]; e[1] = v[1]; e[2] = v[2];}

void Surface::calcEasyVector(){
    /*! Calculates easy vector e from easy angles, tilt and twist*/
    std::cerr << "unimplemented method: " << __PRETTY_FUNCTION__ << std::endl;
    std::exit(1);
}
void Surface::calcV1V2(){
    /*! Calculates v1 and v2 vectors given tilt and twist angles
 *  Rotation matrices are given in Willman, IEEE Trans. Electron Dev. 54, 10, 2007 
 * */

    if (this->Type != Weak) // vectors are only set for 'Weak' anchoring type (why?)
        return;
    double a = Easy[1] * PI / 180.0; // twist
    double b = Easy[0] * PI / 180.0; // tilt
    double g = Easy[2] * PI / 180.0; // rotation around
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

AnchoringType Surface::getAnchoringType() const {
    return this->Type;
}

//unsigned int Surface::getAnchoringNum() const {
//    return static_cast<unsigned int> (this->Type);
//}
double Surface::getStrength() const		{		return Strength;}
double Surface::getK1() const				{		return K1;}
double Surface::getK2() const				{		return K2;}
double Surface::getEasyTilt() const{			return Easy[0];}
double Surface::getEasyTwist() const{			return Easy[1];}
double Surface::getEasyRot() const{			return Easy[2];}
double* Surface::getPtrTov1(){			return &v1[0];}
double* Surface::getPtrTov2(){			return &v2[0];}
bool	Surface::getUsesSurfaceNormal() const {
    return UsesSurfaceNormal;
}
bool    Surface::isStrong() const {
    return isFixed;
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
    surface.push_back(s);
    n_surfaces++;
}

void Alignment::addSurface(const int fixLcNumber,
                           const string &anchoring,
                           const double &strength,
                           const vector<double> &easy,
                           const double &k1,
                           const double &k2,
                           const vector<double> &params) {

    Surface *s = new Surface(fixLcNumber);
    s->setAnchoringType(anchoring);
    s->setStrength(strength);
    s->setEasyAngles(easy);
    s->setK1(k1);
    s->setK2(k2);
    s->Params = params;
    this->addSurface(s);
}


const Surface& Alignment::getSurface(const idx &i) const {
    if (i >= (idx) this->surface.size()) {
        std::cerr << "no such surface " << i << " in " << __PRETTY_FUNCTION__ << std::endl;
        std::exit(1);
    }

    return *this->surface[i];
}


int Alignment::getnSurfaces(){	return n_surfaces;}

bool Alignment::IsStrong(int i) {
    if (i >= (int) surface.size() ){
        printf("error - Alignment::IsStrong - surface %i does not exist\n", i+1 );
        printf("number of surfaces = %i\n", (int) surface.size() );
        exit(1);
    }
    return getSurface(i).isStrong();
    //return surface[i]->getisFixed();
}

AnchoringType Alignment::getTypeOfSurface(const idx &n) const {
    return getSurface(n).getAnchoringType();
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
