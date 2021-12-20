#ifndef MESHREFINEMENT_H
#define MESHREFINEMENT_H
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <material_numbers.h>

enum RefReg_Type { RefReg_Sphere , RefReg_Line, RefReg_Box, RefReg_Surface};

/*!
 * RefReg is a refinement region. It can be of dfferent types. The different types
 * provide different shapes and functionalities
*/
/*
class RefReg {

public:
    RefReg_Type     Type;
    unsigned int material_num;  // material number to be refined when Type is Refreg_Surface
    vector <double>     Params;
    vector <double>     Distance;
    vector <double>     X;
    vector <double>     Y;
    vector <double>     Z;
    RefReg();
    ~RefReg();
    void setType(const RefReg_Type &type) {
        Type = type;
    }
    bool setType(const std::string &type);
    void addParams(const double &P) {
        Params.push_back(P);
    }
    void addDistance(const double &D) {
        Distance.push_back(D);
    }
    void addDistance(const std::vector<double> D);
    void addCoordinate(const double &x, const double &y , const double &z) {
        X.push_back(x);
        Y.push_back(y);
        Z.push_back(z);
    }
    void addX(const double &x) {
        X.push_back(x);
    }
    void addY(const double &y) {
        Y.push_back(y);
    }
    void addZ(const double &z) {
        Z.push_back(z);
    }
    void addX(const std::vector <double> &x);
    void addY(const std::vector <double> &y);
    void addZ(const std::vector <double> &z);
    void PrintRefinementRegion();

    int getNumIterations(); // number of refinement iterations applied to same region

    double getDistance(const int &iter = 0);

    RefReg_Type getType() {
        return Type;
    }
};
*/
/*
enum AutoRef_Type {Change};
class AutoRef {
    std::vector <double> MaxValue;
    double  MinSize; // minimum element size that can be selected as a red tet
    int RefIter;
public:
    AutoRef_Type Type;
    AutoRef();
    void addMaxValue(const double &val) {
        MaxValue.push_back(val);
    }
    void addMaxValue(const std::vector<double> &vals);
    double getMaxValue(const unsigned int &iter);
    double getMinSize() {
        return MinSize;
    }
    void setMinSize(const double &val) {
        MinSize = val;
    }
    bool setType(const std::string &type);
    unsigned int getNumIterations();
    void printAutoref();

    inline int getRefIter() {
        return RefIter;
    }
    inline void setRefIter(const int &iteration) {
        RefIter = iteration;
    }
    inline bool isRefIter(const int &iteration) {
        if (iteration == 0) {
            return false;   // zero RefIter = refine only at end of simu
        }
        if (iteration % RefIter == 0) {
            return true;   // Returns true if this is a refinement iteration
        }
        return false; // default is no
    }
};

// EndRef is performed when simulation has reached final solution.
class EndRef: public AutoRef {
    int EndRefIteration;  // counter of how many times end-refinement has been performed
public:
    EndRef();
    ~EndRef();
    int getEndRefIteration() {
        return EndRefIteration;
    }
    void incrementEndRefIteration() {
        EndRefIteration ++;
    }
};
*/

class RefInfo;
/**
 * RefinementConfig represents the data in a settings file for a single REFINEMENT object. These
 * are later translated into refinement events.
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
    unsigned int    refiter;        // mesh refinement counter
    bool            needs_new_mesh; // true if mesh file has been modified
    unsigned int repRefIter_ = 0;    // iteration period for repeating mesh refinement
    double repRefTime_ = 0;          // time period for repeating mesh refinement

    std::vector<RefinementConfig> refinementConfigs_;

public:
    //std::vector <RefReg> RefinementRegion;
    //std::vector <AutoRef> AutoRefinement;
    //std::vector <EndRef> EndRefinement;
    MeshRefinement();
    //~MeshRefinement() = default;
    //void addRefinementRegion(const RefReg reg);
    //void addAutorefinement(AutoRef reg);
    //void addEndrefinement(EndRef reg);
    void PrintRefinementRegions();
    int getMaxNumRefIterations();
    int getMaxNumAutoRefIterations();
    int getMaxNumEndRefIterations();
    void incRefiter();
    unsigned int getRefiter() {
        return refiter;   // is this global nth time mesh is refined??
    }
    bool isNeedsNewMesh() {
        return needs_new_mesh;
    }
    void setNeedsNewMesh(bool n) {
        needs_new_mesh = n;
    }

    void setRepRefIter(unsigned int iter) { repRefIter_ = iter; }
    [[nodiscard]] unsigned int getRepRefIter() const { return repRefIter_; }

    void setRepRefTime(double time) { repRefTime_ = time; };
    [[nodiscard]] double getRepRefTime() const { return repRefTime_; }

    void setRefinementConfig(std::vector<RefinementConfig> &&ref);
    [[nodiscard]] const std::vector<RefinementConfig>& getRefinementConfig() const { return refinementConfigs_; }
    //bool isRefinementIteration(const int &iteration);  // checks if this iteration is an refinement iteration
};
#endif
