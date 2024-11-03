#ifndef BOX_H
#define BOX_H
#include <geom/aabox.h>
#include <util/randomgenerator.h>
#include <expression.h>
#include <string>
#include <vector>
#include <memory>

class Vec3;
class SolutionVector;
class Coordinates;

enum BoxType {Normal, Random, Hedgehog, InvalidType};

namespace qlc3d {
  class Director;
}


/*!
 * The Box class is used when defining initial LC orientations. It represents a cuboidal sub region within the modelled
 * volume, where the LC orientation properties are set prior to starting the simulation proper.
 */
class Box {
  std::vector<double> Params;
  std::unique_ptr<CartesianExpression> tiltExpression_;
  std::unique_ptr<CartesianExpression> twistExpression_;
  BoxType type;
  AABox boundingBox;
  RandomGenerator randomGenerator;

  [[nodiscard]] double getTiltAt(const Vec3 &p) const;
  [[nodiscard]] double getTwistAt(const Vec3 &p) const;

  [[nodiscard]] qlc3d::Director getDirectorForNormalBox(const Vec3 &p) const;
  [[nodiscard]] qlc3d::Director getDirectorForRandomBox(const Vec3 &p) const;
  [[nodiscard]] qlc3d::Director getDirectorForHedgehogBox(const Vec3 &p) const;

public:
    // Declare default box types
    const static std::vector<std::string> VALID_TYPES;   // list of valid type strings
    const static std::string DEFAULT_TYPE;
    const static std::string DEFAULT_TILT_ANGLE_DEGREES;
    const static std::string DEFAULT_TWIST_ANGLE_DEGREES;
    const static std::vector<double> DEFAULT_PARAMS;
    const static std::vector<double> DEFAULT_X_Y_Z;      // default is same in all dimensions
    const static std::vector<double> DEFAULT_TILT_TWIST; // same defaults bot tilt and twist

    int BoxNumber;

    std::vector<double> Tilt;
    std::vector<double> Twist;
    Box(int boxnum);
    void setX(std::vector<double> x);
    void setY(std::vector<double> y);
    void setZ(std::vector<double> z);
    void setTilt(std::vector<double> tlt);
    void setTwist(std::vector<double> twt);
    void setRandomGenerator(RandomGenerator &rg);
    void setTiltExpression(const std::string &expression);
    void setTwistExpression(const std::string &expression);

    // checks whether [x,y,z] coordinates in array of size 3 are inside the box
    bool contains(double *coords) const;
    bool contains(const Vec3 &coords) const;
    void setBoxType(const std::string &bt);
    [[nodiscard]] qlc3d::Director getDirectorAt(const Vec3 &p) const;

    [[nodiscard]] Vec3 centroid() const;

    [[nodiscard]] BoxType getType() const { return type; }
    [[nodiscard]] std::string getTypeString() const { return VALID_TYPES[type]; }
    [[nodiscard]] const AABox & getBoundingBox() const { return boundingBox; }
    /** Set the params vector. Note: this can only be done once */
    void setParams(const std::vector<double> &p);
    /** return the i'th parameted or default if out of range */
    [[nodiscard]] double getParam(unsigned int i, double defaultValue) const;
    [[nodiscard]] std::string toString() const;

};

class BoxBuilder {
private:
  int boxNum_;
  std::string boxType_;
  std::vector<double> params_;
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> z_;
  std::vector<double> tilt_;
  std::vector<double> twist_;
  std::string tiltExpression_;
  std::string twistExpression_;

public:
  BoxBuilder(int boxNum) : boxNum_(boxNum), boxType_(Box::DEFAULT_TYPE), params_(Box::DEFAULT_PARAMS),
                           x_(Box::DEFAULT_X_Y_Z), y_(Box::DEFAULT_X_Y_Z), z_(Box::DEFAULT_X_Y_Z),
                           tilt_(Box::DEFAULT_TILT_TWIST), twist_(Box::DEFAULT_TILT_TWIST),
                           tiltExpression_(""), twistExpression_("") {}

  BoxBuilder &setBoxType(const std::string &type) {
    boxType_ = type;
    return *this;
  }

  BoxBuilder &setParams(const std::vector<double> &params) {
    params_ = params;
    return *this;
  }

  BoxBuilder &setX(const std::vector<double> &x) {
    if (x.size() != 2) {
      throw std::invalid_argument("Invalid X length " + std::to_string(x.size()) + " expected 2");
    }
    x_ = x;
    return *this;
  }

  BoxBuilder &setY(const std::vector<double> &y) {
    if (y.size() != 2) {
      throw std::invalid_argument("Invalid Y length " + std::to_string(y.size()) + " expected 2");
    }
    y_ = y;
    return *this;
  }

  BoxBuilder &setZ(const std::vector<double> &z) {
    if (z.size() != 2) {
      throw std::invalid_argument("Invalid Z length " + std::to_string(z.size()) + " expected 2");
    }
    z_ = z;
    return *this;
  }

  BoxBuilder &setTilt(const std::vector<double> &tilt) {
    if (tilt.size() != 2) {
      throw std::invalid_argument("Invalid Tilt length " + std::to_string(tilt.size()) + " expected 2");
    }
    tilt_ = tilt;
    return *this;
  }

  BoxBuilder &setTwist(const std::vector<double> &twist) {
    if (twist.size() != 2) {
      throw std::invalid_argument("Invalid Twist length " + std::to_string(twist.size()) + " expected 2");
    }
    twist_ = twist;
    return *this;
  }

  BoxBuilder &setTiltExpression(const std::string &expression) {
    tiltExpression_ = expression;
    return *this;
  }

  BoxBuilder &setTwistExpression(const std::string &expression) {
    twistExpression_ = expression;
    return *this;
  }

  [[nodiscard]] std::unique_ptr<Box> build() const;
};
/*
 * A class that defines the LC initial orientation in the volume.
 */
class InitialVolumeOrientation {
  std::vector<std::unique_ptr<Box>> boxRegions;
  RandomGenerator randomGenerator;
public:
    InitialVolumeOrientation();
    InitialVolumeOrientation(RandomGenerator &rg);
    void addBox(Box *b);
    void addBox(std::unique_ptr<Box> b);
    void addBox(const int &boxNum,
                const std::string &boxType,
                const std::vector<double> &params,
                const std::vector<double> &x,
                const std::vector<double> &y,
                const std::vector<double> &z,
                const std::vector<double> &tilt,
                const std::vector<double> &twist);
  [[nodiscard]] size_t getBoxCount() const { return boxRegions.size(); }
  [[nodiscard]] const Box & getBox(unsigned int i) const {return *boxRegions[i]; };

  /** Initialises the Q-tensor for all boxes */
  void setVolumeQ(SolutionVector &q, double S0, const Coordinates &coordinates) const;
};
#endif
