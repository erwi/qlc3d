#ifndef ALIGNMENT_H
#define ALIGNMENT_H
#include<stdio.h>
#include<vector>
#include<string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
using namespace std;

#define ANCHORING_STRONG		1
#define ANCHORING_WEAK			2
#define ANCHORING_HOMEOTROPIC	3
#define ANCHORING_DEGENERATE	4
#define PI 3.14159265358979323846264338327950288419716939937510 

class Surface
{
	private:
		string Anchoring;					// name of anchoring type
		unsigned int AnchoringNum;			// number 1, 2, 3 or 4 corresponding to anchoring type, as defined above
		double Strength;
		double K1;
		double K2;
		double Easy[3];
		double v1[3];
		double v2[3];
		double e[3];
		bool 	UsesSurfaceNormal; // whether to use local surface normal vector or v1 and v2
        bool isFixed;   // whether this surface is fixed or not
	public:
		//string Anchoring;



		Surface();
		void printSurface();
		void WriteSurface(FILE* fid, int surfnum);


		void setAnchoringType(string atype);
		void setStrength(double str);
		void setK1(double k1);
		void setK2(double k2);
		void setEasyAngles(double ttr[3]);
		void setv1( double v[3] );
		void setv2( double v[3] );
		void setEasyVector( double v[3]);
		void calcEasyVector(); // calculates Easy Vector from easy angles
		void calcV1V2();// calculates v1 and v2 values from easy angles
		void setUsesSurfaceNormal(bool sn);

		string getAnchoringType();
		unsigned int getAnchoringNum();		//
		double getStrength();
		double getK1();
		double getK2();
		double getEasyTilt();
		double getEasyTwist();
		double getEasyRot();
		double* getPtrTov1();
		double* getPtrTov2();
		bool	getUsesSurfaceNormal();
		bool    getisFixed();
};

class Alignment
{
	private:
		int n_surfaces;

	public:
		//int n_Surfaces;
		vector<Surface*> surface;

	Alignment();
	~Alignment();

	void addSurface();
	void addSurface(Surface* s);
	void printAlignment();
	void WriteAlignment(FILE* fid); // writes alignment settings to textfile
	void setnSurfaces(int n);

	double getStrength(int n);	// get strength of FixLCn
	double getK1(int n);		// get K1 of FixLCn
	double getK2(int n);		// get K2 of FixLCn
	double* getPtrTov1(int n);	// get pointer to v1 of FixLCn
	double* getPtrTov2(int n);  // get pointer to v2 of FicLCn

	int getnSurfaces();
	bool IsStrong(int n); 		// is surface n strong?
	unsigned int getAnchoringNum(const int &n); // // returns anchoring number of alignment surface n
	bool getUsesSurfaceNormal(int n); // if n uses surface normal instead of v1 and v2 ?
	bool WeakSurfacesExist(); 	// checks if weak surfaces are defined

};


#endif
