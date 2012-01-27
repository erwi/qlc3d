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
    enum Type {None = 0, Change};
private:
    Type type_;
    long int iter_; // WHEN TO REFINE, IF BOTH ZERO ASSUME REPEATING
    double time_;   //
    int refIter_;    // HOW MANY REFINMENT STEPS THIS OBJECT DESCRIBES

    std::vector<double> values_; // THIS WILL HOLD DIFFERENT VALUES DEPENDING ON type_

    RefInfo(){}
    void setType( std::string Type );
    void setRefIter();
public:
    RefInfo(const std::string& Type);
    RefInfo(const RefInfo& other);      // COPY CONSTRUCTOR

    void setIteration(const long int i){iter_ = i;}
    void setTime(const double t){time_ = t;}
    void setValues(std::vector<double>& values);
    int getRefIter()const {return refIter_;}
    void printRefInfo(FILE* fid = stdout) const;
};

#endif // REFINFO_H
