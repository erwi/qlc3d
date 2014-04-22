#ifndef MESHREFINEMENT_H
#define MESHREFINEMENT_H
#include <list>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <material_numbers.h>
using std::vector;
using std::list;
using std::cout;
using std::endl;

enum RefReg_Type { RefReg_Sphere , RefReg_Line, RefReg_Box, RefReg_Surface};

class RefReg {
    /*! RefReg is a refinement region. It can be of dfferent types. The different types
    provide different shapes and functionalities
    */

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
    void addX(const vector <double> &x);
    void addY(const vector <double> &y);
    void addZ(const vector <double> &z);
    void PrintRefinementRegion();

    int getNumIterations(); // number of refinement iterations applied to same region

    double getDistance(const int &iter = 0);

    RefReg_Type getType() {
        return Type;
    }
};

enum AutoRef_Type {Change};
class AutoRef {
    vector <double> MaxValue;
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

class MeshRefinement {
    unsigned int    refiter;        // mesh refinement counter
    bool            needs_new_mesh; // true if mesh file has been modified
public:
    vector <RefReg> RefinementRegion;
    vector <AutoRef> AutoRefinement;
    vector <EndRef> EndRefinement;
    MeshRefinement();
    ~MeshRefinement();
    void addRefinementRegion(const RefReg reg);
    void addAutorefinement(AutoRef reg);
    void addEndrefinement(EndRef reg);
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
    bool isRefinementIteration(const int &iteration);  // checks if this iteration is an refinement iteration
};
#endif
