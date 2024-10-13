#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include <stdio.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <globals.h>
#include <geom/vec3.h>

class Vec3;
enum AnchoringType {Strong = 0, Weak = 1, Homeotropic = 2,
                     Degenerate = 3, Freeze = 5, Polymerise = 6,
                     ManualNodes = 7, WeakHomeotropic = 8};

/**The surface class represents a single FIXLC anchoring surface.*/
class Surface {

private:
    AnchoringType Type;
    double Strength;
    double K1;
    double K2;
    double easyAnglesDegrees[3]{};                     // Easy direction angles in degrees
    Vec3 v1 = {0, 0, 0};                       // First principal axis of anchoring
    Vec3 v2 = {0, 0, 0};                       // Second principal axis of anchoring
    Vec3 e = {0, 0, 0};                        // Easy direction vector
    bool overrideVolume;               // whether to override volumes at startup. This is set to true by default
    unsigned int fixLcNumber = 0;

    /**
    * Calculates v1 and v2 vectors given tilt and twist angles Rotation matrices are given in
    * Willman, IEEE Trans. Electron Dev. 54, 10, 2007
    **/
    [[nodiscard]] static Vec3 calculateV1(double tiltDegrees, double twistDegrees, double rotDegrees = 0);
    [[nodiscard]] static Vec3 calculateV2(double tiltDegrees, double twistDegrees, double rotDegrees = 0);
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
    Surface(AnchoringType anchoringType, double strength, double k1, double k2,
            const double easyAnglesDegrees[3],
            bool overrideVolume, unsigned int fixLcNumber);

    //Surface(const Surface &s);

    [[nodiscard]] static Surface ofStrongAnchoring(unsigned int fixLcNumber, double tiltDegrees, double twistDegrees);
    [[nodiscard]] static Surface ofPlanarDegenerate(unsigned int fixLcNumber, double strength);
    [[nodiscard]] static Surface ofStrongHomeotropic(unsigned int fixLcNumber);
    [[nodiscard]] static Surface ofWeakHomeotropic(unsigned int fixLCNumber, double strength);
    [[nodiscard]] static Surface ofWeakAnchoring(unsigned int fixLcNumber,
                                               double tiltDegrees, double twistDegrees, double strength, double k1, double k2);

  [[nodiscard]] std::string getAnchoringTypeName() const;
  [[nodiscard]] AnchoringType getAnchoringType() const;
  [[nodiscard]] double getStrength() const;
  [[nodiscard]] double getK1() const;
  [[nodiscard]] double getK2() const;
  [[nodiscard]] double getEasyTilt() const;
  [[nodiscard]] double getEasyTwist() const;
  [[nodiscard]] const Vec3& getV1() const { return v1; }
  [[nodiscard]] const Vec3& getV2() const { return v2; }
  [[nodiscard]] const Vec3& getEasyVector() const { return e; }

  [[nodiscard]] bool usesSurfaceNormal() const;
  [[nodiscard]] bool isStrong() const;
  [[nodiscard]] bool getOverrideVolume() const { return overrideVolume; }

  [[nodiscard]] std::string toString() const;
  [[nodiscard]] unsigned int getFixLCNumber() const { return fixLcNumber; }
  friend class Alignment;
};

class Reader; // forward declaration of reader class

class Alignment {
    /*! A collection of Surface objects, each representing a FIXLC surface*/
private:
    int n_surfaces; // TODO: why not just return surface.size() ??
    void setnSurfaces(int n);
public:
    std::vector<Surface> surface;
    Alignment();

    void addSurface(Surface s);

    [[nodiscard]] const Surface & getSurface(const idx& i) const;

    double getStrength(int n);  // get strength of FixLCn
    double getK1(int n);        // get K1 of FixLCn
    double getK2(int n);        // get K2 of FixLCn

    [[nodiscard]] int getnSurfaces() const;
    [[nodiscard]] bool isStrong(int n) const;       // is surface n strong?
    [[nodiscard]] bool getUsesSurfaceNormal(int n) const; // if n uses surface normal instead of v1 and v2 ?
    [[nodiscard]] bool weakSurfacesExist() const;
    [[nodiscard]] bool hasSurface(unsigned int fixLcNumber) const;
    [[nodiscard]] std::unordered_map<unsigned int, Surface> getWeakSurfacesByFixLcNumber() const;
};
#endif
