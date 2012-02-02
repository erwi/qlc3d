#ifndef REFINFO_H
#define REFINFO_H
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <vector>
class RefInfo
{
public:
    enum Type {None = 0, Change = 1, Sphere=2};
private:
    Type type_;
    long int iter_; // WHEN TO REFINE, IF BOTH ZERO ASSUME REPEATING
    double time_;   //
    int refIter_;    // HOW MANY REFINMENT STEPS THIS OBJECT DESCRIBES

    std::vector<double> values_; // THIS WILL HOLD DIFFERENT VALUES DEPENDING ON type_
    std::vector<double> X_;      // LISTS OF COORDINATES FOR A SPECIFIED REGION
    std::vector<double> Y_;
    std::vector<double> Z_;



    RefInfo(){}
    void setType( std::string Type );
    void setRefIter();

    //void error(const char* var);    // PRINTS ERROR MESSAGE FOR VARIABLE var AND EXITS
public:
    RefInfo(const std::string& Type);
    RefInfo(const RefInfo& other);      // COPY CONSTRUCTOR

    void setIteration(const long int i){iter_ = i;}
    void setTime(const double t){time_ = t;}
    void setValues(std::vector<double>& values);
    void setCoords(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& z);



    int    getRefIter()const {return refIter_;}
    double getValue(const size_t i) const;//
    void   getCoords(std::vector<double>& x,    // RETURNS ALL COORDINATES IN A VECTOR
                     std::vector<double>& y,
                     std::vector<double>& z) const ;
    void   getCoord( double& x, double &y, double& z)const; // RETURNS ONLY FIRST COORD
    Type   getType() const {return type_;}
    void   printRefInfo(FILE* fid = stdout) const;
    static void validate(const RefInfo& refinfo );  // TRIES TO VALIDATE TO MAKE SURE INFO PROVIDED MAKES SENSE
};

#endif // REFINFO_H
