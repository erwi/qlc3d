#include <alignment.h>
#include <algorithm>
#include <iostream>
#include <reader.h>

const double Surface::DEFAULT_ANCHORING_STRENGTH = 1e-4;


Surface::Surface(int fxlcnum)
{

    Anchoring = "Strong";
    AnchoringNum = ANCHORING_STRONG;
    FixLCNumber = fxlcnum;
    Strength = DEFAULT_ANCHORING_STRENGTH;
    K1 = 1;
    K2 = 1;
    Easy[0] = 0;
    Easy[1] = 0;
    Easy[2] = 0;

    v1[0] = 0;
    v1[1] = 0;
    v1[2] = 1;

    v2[0] = 0;
    v2[1] = 1;
    v2[2] = 0;

    e[0] = 1;
    e[1] = 0;
    e[2] = 0;
    UsesSurfaceNormal = false;
    isFixed = true;
}

void Surface::setAnchoringType(std::string &atype){
    Anchoring = atype;

    // MAKE SURE LOWERCASE
    std::transform(atype.begin(), atype.end() , atype.begin(), ::tolower);
    if (atype.compare("strong") == 0) {
        Anchoring = atype;
        AnchoringNum = ANCHORING_STRONG;
        isFixed = true;
    }
    else if (atype.compare("homeotropic") == 0) {
        Anchoring = atype;
        AnchoringNum = ANCHORING_HOMEOTROPIC;
        isFixed = true;
    }
    else if (atype.compare("weak") == 0) {
        Anchoring = atype;
        AnchoringNum = ANCHORING_WEAK;
        isFixed = false;
    }
    else if (atype.compare("degenerate") == 0) {
        Anchoring = atype;
        AnchoringNum = ANCHORING_DEGENERATE;
        setUsesSurfaceNormal(true);
        isFixed = false;
    }
    else if (atype.compare("freeze") == 0 ) {
        Anchoring = atype;
        AnchoringNum = ANCHORING_FREEZE;
        isFixed = true;
        setUsesSurfaceNormal(false);
    }
    else if (atype.compare("polymerise") == 0 ) {
        Anchoring = atype;
        AnchoringNum = ANCHORING_POLYMERISE;
        isFixed = true;
        setUsesSurfaceNormal(false);
    }
    else if (atype.compare("manualnodes") == 0) {
        Anchoring = atype;
        AnchoringNum = ANCHORING_MANUAL_NODES;
        isFixed = true;
        setUsesSurfaceNormal(false);
    }
    else{
        using std::cerr;
        using std::endl;
        cerr << "error specifying FIXLC"<< FixLCNumber<<".Anchoring, ";
        cerr << "\""<<Anchoring << "\" is not a valid anchoring type." << endl;
        cerr << "valid types are:\n" <<
                "\tstrong\n"<<
                "\tweak\n"<<
                "\tdegenerate\n"<<
                "\tfreeze\n"<<
                "\tpolymerise\n" <<
                "Check your settings file for typos - bye!\n " <<endl;
        exit(1);
    }

}// end setAnchoringType
void Surface::setStrength(double str){	Strength = str;}
void Surface::setK1(double k1){			K1 = k1;}
void Surface::setK2(double k2){			K2 = k2;}


//void Surface::setEasyAngles(double ttr[3]){		Easy[0] = ttr[0];	Easy[1] = ttr[1];	Easy[2] = ttr[2];}
void Surface::setEasyAngles(std::vector<double> &e){
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

    if (this->AnchoringNum != ANCHORING_WEAK ){ // vectors are only set for 'Weak' anchoring type
        return;
    }

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
void Surface::setUsesSurfaceNormal(bool sn) { UsesSurfaceNormal = sn;}

string Surface::getAnchoringType()	{		return Anchoring;}
unsigned int Surface::getAnchoringNum() {	return AnchoringNum;}
double Surface::getStrength()		{		return Strength;}
double Surface::getK1()				{		return K1;}
double Surface::getK2()				{		return K2;}
double Surface::getEasyTilt(){			return Easy[0];}
double Surface::getEasyTwist(){			return Easy[1];}
double Surface::getEasyRot(){			return Easy[2];}
double* Surface::getPtrTov1(){			return &v1[0];}
double* Surface::getPtrTov2(){			return &v2[0];}
bool	Surface::getUsesSurfaceNormal(){return UsesSurfaceNormal;}
bool    Surface::getisFixed()           {return isFixed;}
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


int Alignment::getnSurfaces(){	return n_surfaces;}

bool Alignment::IsStrong(int i){
    if (i >= (int) surface.size() ){
        printf("error - Alignment::IsStrong - surface %i does not exist\n", i+1 );
        printf("number of surfaces = %i\n", (int) surface.size() );
        exit(1);
    }
    return surface[i]->getisFixed();
}

unsigned int Alignment::getAnchoringNum(const int &n){
    // returns anchoring number of alignment surface n
    // Anchoring number is a descriptor that maps to a type.
    // It should be made into an enumerator
    vector <Surface*>::iterator itr;
    itr = surface.begin();
    if ( (n < getnSurfaces() )  && (n>=0) ){
        return (*(itr+n))->getAnchoringNum();
    }
    else{
        printf("error - Alignment::getAnchoringNum(%i) - %i is an invalid surface number, bye!\n", n ,n);
        exit(1);
    }

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
