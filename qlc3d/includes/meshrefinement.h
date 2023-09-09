#ifndef MESHREFINEMENT_H
#define MESHREFINEMENT_H
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <material_numbers.h>

/**
 * This class contains information about a mesh region that may need refinement. It is used by the
 * mesh refinement algorithm. This is also the the void pointer "eventData", in the Event class, for mesh
 * refinement events.
 */
class RefInfo
{
public:
    enum Type {None = 0, Change = 1, Sphere=2, Box = 3};
private:
    Type _type;
    long long int _iter; // WHEN TO REFINE, IF BOTH ZERO ASSUME REPEATING
    double _time;   //
    unsigned int _refIter;   // HOW MANY REFINMENT STEPS THIS OBJECT DESCRIBES
    std::vector<double> _values; // THIS WILL HOLD DIFFERENT VALUES DEPENDING ON type_
    std::vector<double> _x;      // LISTS OF COORDINATES FOR A SPECIFIED REGION
    std::vector<double> _y;
    std::vector<double> _z;
    // Private constructors, use the factory method make() instead
    RefInfo() { }
    RefInfo(const std::string& Type);
    void setType( std::string Type );
    void setRefIter();
    void setIteration(const long int i) { _iter = i; }
    void setTime(const double t) {_time = t; }
    /** set values as deep copy */
    void setValues(const std::vector<double>& values);
    /** set coordinates as deep copy */
    void setCoords(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& z);
    static void validate(const RefInfo& refinfo);  // TRIES TO VALIDATE TO MAKE SURE INFO PROVIDED MAKES SENSE
public:
    RefInfo(const RefInfo& other);      // COPY CONSTRUCTOR

    [[nodiscard]] unsigned int getRefIter() const { return _refIter; }
    [[nodiscard]] double getRefTime() const { return _time; }
    [[nodiscard]] double getValue(const size_t i) const;
    [[nodiscard]] const std::vector<double> & getX() const { return _x; }
    [[nodiscard]] const std::vector<double> & getY() const { return _y; }
    [[nodiscard]] const std::vector<double> & getZ() const { return _z; }
    [[nodiscard]] Type getType() const { return _type; }

    static RefInfo *make(const std::string& type,       // A Factory method. simplified creation
                         long int iteration,            // sets values and calls validate()
                         double time,
                         const std::vector<double> &values,
                         const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &z);

    /** Factory method that creates raw pointer to RefInfo configured to occur periodically */
    static RefInfo* ofPeriodicMeshRefinement(const std::string &type,
                                             const std::vector<double> &values,
                                             const std::vector<double> &x,
                                             const std::vector<double> &y,
                                             const std::vector<double> &z);
};

/**
 * RefinementConfig represents the data in a settings file for a single REFINEMENT object. These
 * are later translated into refinement events.
 *
 * TODO: this seems to essentially duplicate RefInfo. Are both really required?
 */
class RefinementConfig {
    /** throws exception if an invalid RefinementConfig is detected */
    void validate();
public:
    const std::string type_;
    const std::vector<int> iterations_;
    const std::vector<double> times_;
    const std::vector<double> values_;
    const std::vector<double> x_;
    const std::vector<double> y_;
    const std::vector<double> z_;

    RefinementConfig(std::string type, std::vector<int> iterations, std::vector<double> times,
                     std::vector<double> values,
                     std::vector<double> x, std::vector<double> y, std::vector<double> z) :
            type_ { type },
            iterations_ { iterations },
            times_ { times },
            values_ { values },
            x_ { x }, y_ { y }, z_ { z } {
        validate();
    }

    /**
     * if no iteration or time is defined, then this is assumed to occur periodically, with a period defined as
     * repRefIter and/or repRefTime
     */
    bool occursPeriodically() const { return iterations_.empty() && times_.empty(); }
};

class MeshRefinement {
    unsigned int repRefIter_ = 0;    // iteration period for repeating mesh refinement
    double repRefTime_ = 0;          // time period for repeating mesh refinement
    std::vector<RefinementConfig> refinementConfigs_;

public:
    void setRepRefIter(unsigned int iter) { repRefIter_ = iter; }
    [[nodiscard]] unsigned int getRepRefIter() const { return repRefIter_; }

    void setRepRefTime(double time) { repRefTime_ = time; };
    [[nodiscard]] double getRepRefTime() const { return repRefTime_; }

    void setRefinementConfig(std::vector<RefinementConfig> &&ref);
    [[nodiscard]] const std::vector<RefinementConfig>& getRefinementConfig() const { return refinementConfigs_; }
};
#endif
