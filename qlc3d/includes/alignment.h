#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include<stdio.h>
#include<vector>
#include<string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <globals.h>


enum AnchoringType {Strong = 0, Weak = 1, Homeotropic = 2,
                     Degenerate = 3, Freeze = 5, Polymerise = 6,
                     ManualNodes = 7, AnchoringTypesCount};

class Surface {
    /*!The surface class represents a single FIXLC anchoring surface.*/
private:


    //string Anchoring;                   // name of anchoring type TODO recover this using index to VALID_TYPES instead
    //unsigned int AnchoringNum;          // number 1, 2, 3 or 4 corresponding to anchoring type, as defined above
    AnchoringType Type;
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
    static const std::vector<std::string> VALID_ANCHORING_TYPES;
    static const std::string DEFAULT_ANCHORING_TYPE;
    static const double DEFAULT_ANCHORING_STRENGTH;
    static const double DEFAULT_ANCHORING_K1;
    static const double DEFAULT_ANCHORING_K2;
    static const std::vector<double> DEFAULT_ANCHORING_EASY;
    static const std::vector<double> DEFAULT_ANCHORING_PARAMS;
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

    std::string getAnchoringTypeName() const;
    AnchoringType getAnchoringType() const;
    double getStrength() const;
    double getK1() const;
    double getK2() const;
    double getEasyTilt() const;
    double getEasyTwist() const;
    double getEasyRot() const;
    double *getPtrTov1();
    double *getPtrTov2();
    bool    getUsesSurfaceNormal() const;
    bool    isStrong() const;

    friend class Alignment;
};

class Reader; // forward declaration of reader class

class Alignment {
    /*! A collection of Surface objects, each representing a FIXLC surface*/
private:
    int n_surfaces; // TODO: why not just return surface.size() ??
    void addSurface(Surface *s);
    void setnSurfaces(int n);
public:
    std::vector<Surface *> surface;
    Alignment();
    ~Alignment();

    void addSurface(const int fixLcNumber,
                    const std::string &anchoring,
                    const double &strength,
                    const std::vector<double> &easy,
                    const double &k1,
                    const double &k2,
                    const std::vector<double> &params);
    const Surface & getSurface(const idx& i) const; // returns read-only reference to i'th surface

    double getStrength(int n);  // get strength of FixLCn
    double getK1(int n);        // get K1 of FixLCn
    double getK2(int n);        // get K2 of FixLCn
    double *getPtrTov1(int n);  // get pointer to v1 of FixLCn
    double *getPtrTov2(int n);  // get pointer to v2 of FicLCn
    int getnSurfaces();
    bool IsStrong(int n) const;       // is surface n strong?
    //unsigned int getAnchoringNum(const int &n); // // returns anchoring number of alignment surface n
    AnchoringType getTypeOfSurface(const idx &n) const; //returns type of n'th surface
    bool getUsesSurfaceNormal(int n); // if n uses surface normal instead of v1 and v2 ?
    bool WeakSurfacesExist();   // checks if weak surfaces are defined
};


#endif
