#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include<stdio.h>
#include<vector>
#include<string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <globals.h>

class Vec3;
enum AnchoringType {Strong = 0, Weak = 1, Homeotropic = 2,
                     Degenerate = 3, Freeze = 5, Polymerise = 6,
                     ManualNodes = 7, AnchoringTypesCount};

class Surface {
    /*!The surface class represents a single FIXLC anchoring surface.*/
private:
    AnchoringType Type;
    double Strength;
    double K1;
    double K2;
    double easyAnglesDegrees[3] = {0, 0, 0};                     // Easy direction angles in degrees
    double v1[3] = {0, 0, 0};                       // First principal axis of anchoring
    double v2[3] = {0, 0, 0};                       // Second principal axis of anchoring
    double e[3] = {0, 0, 0};                        // Easy direction vector
    bool UsesSurfaceNormal;             // whether to use local surface normal vector or v1 and v2
    bool isFixed;                       // whether this surface is fixed or not
    bool overrideVolume;               // whether to override volumes at startup. This is set to true by default
    void calcV1V2();                    // calculates v1 and v2 values from easy angles
    unsigned int fixLcNumber = 0;
    void setAnchoringType(const std::string &atype);
public:
    static const std::vector<std::string> VALID_ANCHORING_TYPES;
    static const std::string DEFAULT_ANCHORING_TYPE;
    static const double DEFAULT_ANCHORING_STRENGTH;
    static const double DEFAULT_ANCHORING_K1;
    static const double DEFAULT_ANCHORING_K2;
    static const std::vector<double> DEFAULT_ANCHORING_EASY;
    static const std::vector<double> DEFAULT_ANCHORING_PARAMS;
    static const bool DEFAULT_ANCHORING_OVERRIDE_VOLUME;
    std::vector<double> Params;         // holds optional parameters, but is mostly empty

    Surface(int fxlcnum, const std::string &type);
    void setStrength(double str);
    void setK1(double k1);
    void setK2(double k2);
    void setEasyAngles(const std::vector<double> &e);
    void setv1(double v[3]);
    void setv2(double v[3]);
    void setEasyVector(double v[3]);
    void setEnforce(bool overrideVolume) { this->overrideVolume = overrideVolume; }

  [[nodiscard]] std::string getAnchoringTypeName() const;
  [[nodiscard]] AnchoringType getAnchoringType() const;
  [[nodiscard]] double getStrength() const;
  [[nodiscard]] double getK1() const;
  [[nodiscard]] double getK2() const;
  [[nodiscard]] double getEasyTilt() const;
  [[nodiscard]] double getEasyTwist() const;
  double *getPtrTov1();
  double *getPtrTov2();
  [[nodiscard]] bool getUsesSurfaceNormal() const;
  [[nodiscard]] bool isStrong() const;
  [[nodiscard]] bool getOverrideVolume() const { return overrideVolume; }

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] unsigned int getFixLCNumber() const { return fixLcNumber; }
  void setFixLCNumber(unsigned int fixLcNumber);
  friend class Alignment;

  [[nodiscard]] static Surface* ofStrongAnchoring(unsigned int fixLcNumber, double tiltDegrees, double twistDegrees);
  [[nodiscard]] static Surface* ofPlanarDegenerate(unsigned int fixLcNumber, double strength);
  [[nodiscard]] static Surface* ofHomeotropic(unsigned int fixLcNumber);
};

class Reader; // forward declaration of reader class

class Alignment {
    /*! A collection of Surface objects, each representing a FIXLC surface*/
private:
    int n_surfaces; // TODO: why not just return surface.size() ??
    void setnSurfaces(int n);
public:
    std::vector<Surface *> surface;
    Alignment();
    ~Alignment();

    void addSurface(const int fixLcNumber,
                    const std::string &anchoring,
                    const double &strength,
                    const std::vector<double> &easyAnglesDegrees,
                    const double &k1,
                    const double &k2,
                    const std::vector<double> &params,
                    const bool overrideVolume = true);

    void addSurface(Surface *s);

    [[nodiscard]] const Surface & getSurface(const idx& i) const; // returns read-only reference to i'th surface

    double getStrength(int n);  // get strength of FixLCn
    double getK1(int n);        // get K1 of FixLCn
    double getK2(int n);        // get K2 of FixLCn
    double *getPtrTov1(int n);  // get pointer to v1 of FixLCn
    double *getPtrTov2(int n);  // get pointer to v2 of FicLCn
    [[nodiscard]] int getnSurfaces() const;
    [[nodiscard]] bool IsStrong(int n) const;       // is surface n strong?
    bool getUsesSurfaceNormal(int n); // if n uses surface normal instead of v1 and v2 ?
    bool WeakSurfacesExist();   // checks if weak surfaces are defined
    [[nodiscard]] bool hasSurface(unsigned int fixLcNumber) const;
};
#endif
