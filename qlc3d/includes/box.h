#ifndef BOX_H
#define BOX_H
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>

class Vec3;
/*!
 * The Box class is used when defining initial LC orientations. It represents a cuboidal sub region within the modelled
 * volume, where the LC orientation properties are set prior to starting the simulation proper.
 */
class Box {
    std::vector<double> Params;
public:
    // Declare default box types
    const static std::vector<std::string> VALID_TYPES;   // list of valid type strings
    const static std::string DEFAULT_TYPE;
    const static std::vector<double> DEFAULT_PARAMS;
    const static std::vector<double> DEFAULT_X_Y_Z;      // default is same in all dimensions
    const static std::vector<double> DEFAULT_TILT_TWIST; // same defaults bot tilt and twist

    enum BoxTypes {Normal, Random, Hedgehog, InvalidType};
    BoxTypes Type;
    std::string TypeString;
    int BoxNumber;

    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    std::vector<double> Tilt;
    std::vector<double> Twist;
    Box(int boxnum);
    void setX(std::vector<double> x);
    void setY(std::vector<double> y);
    void setZ(std::vector<double> z);
    void setTilt(std::vector<double> tlt);
    void setTwist(std::vector<double> twt);

    // checks whether [x,y,z] coordinates in array of size 3 are inside the box
    bool contains(double *coords) const;
    bool contains(const Vec3 &coords) const;
    void setBoxType(const std::string &bt);
    [[nodiscard]] Vec3 centroid() const;
    /** Set the params vector. Note: this can only be done once */
    void setParams(const std::vector<double> &p);
    /** return the i'th parameted or default if out of range */
    [[nodiscard]] double getParam(int i, double defaultValue) const;

    [[nodiscard]] std::string toString() const;
};


class Boxes {
    /*! A collection of multiple Box regions*/
public:
    std::vector<Box *> box;
    int n_Boxes;
    Boxes();
    ~Boxes();
    void addBox(Box *b);
    void addBox(const int &boxNum,
                const std::string &boxType,
                const std::vector<double> &params,
                const std::vector<double> &x,
                const std::vector<double> &y,
                const std::vector<double> &z,
                const std::vector<double> &tilt,
                const std::vector<double> &twist);


};
#endif
