#ifndef REFINFO_H
#define REFINFO_H
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <iostream>

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
    // Privatre constructors, use the factory method make() instead
    RefInfo(){}
    RefInfo(const std::string& Type);
    void setType( std::string Type );
    void setRefIter();
    void setIteration(const long int i){_iter = i;}
    void setTime(const double t){_time = t;}
    /** set values as deep copy */
    void setValues(const std::vector<double>& values);
    /** set coordinates as deep copy */
    void setCoords(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& z);
    static void validate(const RefInfo& refinfo );  // TRIES TO VALIDATE TO MAKE SURE INFO PROVIDED MAKES SENSE
public:
    RefInfo(const RefInfo& other);      // COPY CONSTRUCTOR

    unsigned int    getRefIter()const {return _refIter;}
    double getValue(const size_t i) const;//
    void   getCoords(std::vector<double>& x,    // RETURNS ALL COORDINATES IN A VECTOR
                     std::vector<double>& y,
                     std::vector<double>& z) const ;
    void   getCoord( double& x, double &y, double& z)const; // RETURNS ONLY FIRST COORD
    Type   getType() const {return _type;}
    void   printRefInfo(FILE* fid = stdout) const;

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
#endif // REFINFO_H
