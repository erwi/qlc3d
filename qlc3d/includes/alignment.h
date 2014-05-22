#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include<stdio.h>
#include<vector>
#include<string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//using namespace std;

#define ANCHORING_STRONG    1
#define ANCHORING_WEAK      2
#define ANCHORING_HOMEOTROPIC   3
#define ANCHORING_DEGENERATE    4
#define ANCHORING_FREEZE    5   // anchoring type where initial Q-Tensor orientation is frozen
#define ANCHORING_POLYMERISE    6
#define ANCHORING_MANUAL_NODES 7
#ifndef PI
#define PI 3.14159265358979323846264338327950288419716939937510
#endif
using std::string;
using std::vector;
class Surface {
    /*!The surface class represents a single FIXLC anchoring surface.*/
private:


    string Anchoring;                   // name of anchoring type
    unsigned int AnchoringNum;          // number 1, 2, 3 or 4 corresponding to anchoring type, as defined above
    double Strength;
    double K1;
    double K2;
    double Easy[3];                     // Easy direction angles
    double v1[3];                       // First principal axis of anchoring
    double v2[3];                       // Second principal axis of anchoring
    double e[3];                        // Easy direction vector
    bool UsesSurfaceNormal;             // whether to use local surface normal vector or v1 and v2
    bool isFixed;                       // whether this surface is fixed or not
    void calcV1V2();                    // calculates v1 and v2 values from easy angles
public:
    static const string DEFAULT_ANCHORING_TYPE;
    static const double DEFAULT_ANCHORING_STRENGTH;
    static const double DEFAULT_ANCHORING_K1;
    static const double DEFAULT_ANCHORING_K2;
    static const vector<double> DEFAULT_ANCHORING_EASY;
    static const vector<double> DEFAULT_ANCHORING_PARAMS;
    int FixLCNumber;
    std::vector<double> Params;         // holds optional parameters, but is mostly empty
    Surface(int fxlcnum);
    void setAnchoringType(const std::string &atype);
    void setStrength(double str);
    void setK1(double k1);
    void setK2(double k2);
    void setEasyAngles(const std::vector<double> &e);
    void setv1(double v[3]);
    void setv2(double v[3]);
    void setEasyVector(double v[3]);
    void calcEasyVector();              // calculates Easy Vector from easy angles

    void setUsesSurfaceNormal(bool sn);
    string getAnchoringType();
    unsigned int getAnchoringNum();     //
    double getStrength();
    double getK1();
    double getK2();
    double getEasyTilt();
    double getEasyTwist();
    double getEasyRot();
    double *getPtrTov1();
    double *getPtrTov2();
    bool    getUsesSurfaceNormal();
    bool    getisFixed();

    friend class Alignment;
};

class Reader; // forward declaration of reader class

class Alignment {
    /*! A collection of Surface objects, each representing a FIXLC surface*/
private:
    int n_surfaces;
    void addSurface(Surface *s);
public:
    vector<Surface *> surface;
    Alignment();
    ~Alignment();

    void addSurface(const int fixLcNumber,
                    const string &anchoring,
                    const double &strength,
                    const vector<double> &easy,
                    const double &k1,
                    const double &k2,
                    const vector<double> &params);

    void setnSurfaces(int n);
    double getStrength(int n);  // get strength of FixLCn
    double getK1(int n);        // get K1 of FixLCn
    double getK2(int n);        // get K2 of FixLCn
    double *getPtrTov1(int n);  // get pointer to v1 of FixLCn
    double *getPtrTov2(int n);  // get pointer to v2 of FicLCn
    int getnSurfaces();
    bool IsStrong(int n);       // is surface n strong?
    unsigned int getAnchoringNum(const int &n); // // returns anchoring number of alignment surface n
    bool getUsesSurfaceNormal(int n); // if n uses surface normal instead of v1 and v2 ?
    bool WeakSurfacesExist();   // checks if weak surfaces are defined
};


#endif
