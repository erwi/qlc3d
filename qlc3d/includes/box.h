#ifndef BOX_H
#define BOX_H
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
class Box
{
	public:

    enum BoxTypes {Normal, Random, Hedgehog};
    BoxTypes Type;
    std::string TypeString;
    int BoxNumber;

    double Params[2];
	double X[2];
	double Y[2];
	double Z[2];

	double Tilt[2];
	double Twist[2];

    Box(int boxnum);
	void printBox();

    void setParams(std::vector<double> p);
    void setX(std::vector<double> x);
    void setY(std::vector<double> y);
    void setZ(std::vector<double> z);
    void setTilt(std::vector<double> tlt);
    void setTwist(std::vector<double> twt);
    bool isInBox(double* coords);           // checks whether [x,y,z] coordinates in array of size 3 are inside the box
    void setBoxType(std::string &bt);
};



class Boxes
{
	public:
		std::vector<Box*> box;
		int n_Boxes;
	Boxes();
	~Boxes();

	void addBox(Box* b);
	void printBoxes();

};




#endif
