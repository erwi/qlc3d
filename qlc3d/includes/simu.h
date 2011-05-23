#ifndef SIMU_H
#define SIMU_H
#include <string>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define SIMU_N  0
#define SIMU_END_SIMULATION -20000000 // magig minus bignum

#define SIMU_OUTPUT_FORMAT_BINARY 	0
#define SIMU_OUTPUT_FORMAT_TEXT		1

enum PotentialConsistency {Off, Loop};

using namespace std;
class Simu
{
	//using namespace std;
	private:
		PotentialConsistency PotCons;
		double  TargetPotCons; // minimum potential consistency when PotCons == Loop

		double	MaxError;
		double 	CurrentTime;
		double  CurrentChange; // change in Q, for determining convergence

		// ADAPTIVE TIME-STEPPING VARIABLES
		double  TargetdQ;
		double  dtLimits[2];
		double  dtFunction[4];
		double 	Maxdt;




		string 	EndCriterion;
		string  LoadQ;
		string  SaveDir;
		string 	LoadDir;

		double 	EndValue;
		double  StretchVector[3];
		double  EnergyRegion[3]; // x,y,z coordinates for determining regions above/below which to calculate energy
		int 	CurrentIteration;
		bool	AssembleMatrix;
                bool    MeshModified;
                int     MeshNumber;     // counts number of modifications. This number is appended to the end of mesh name
                int	OutputEnergy;	// boolean whether or not to calculate energy

                int	OutputFormat;
                int	SaveIter; // determines frequency of saving intermediate result files

	public:
		string MeshName;
		double dt;


//	int OutputEnergy;

	Simu();
	void PrintSimu();
	void WriteSimu(FILE* fid); // writes simu settings structure to file

	void setMeshName(string meshname);
	void setSaveDir(string savedir);
	void setLoadDir(string loaddir);
	void setMaxError(double me);
	void setEndCriterion( string ec);
	void setEndValue(double ev);
    void resetEndCriterion();       // starts simulation from beginning
    void setdt(double td);
	void setdtLimits(const double& min, const double& max);
	void setdtFunction(double* f4); // array of length 4
	void getdtFunction(double* f );
	void setCurrentTime(double ct);
	void setCurrentIteration( int i);
	void setAssembleMatrix(bool yn);
	void setCurrentChange(double ch);
	void setLoadQ(string qbackup);
	void setMaxdt(double maxdt);
	void setTargetdQ(const double& dq);


	void setOutputEnergy(int ope);
	void setOutputFormat(int opf);
	void setStretchVectorX(double sx);
	void setStretchVectorY(double sy);
	void setStretchVectorZ(double sz);
	inline void setEnergyRegionX(const double& x) { EnergyRegion[0] = x;}
	inline void setEnergyRegionY(const double& y) { EnergyRegion[1] = y;}
	inline void setEnergyRegionZ(const double& z) { EnergyRegion[2] = z;}
	inline void setSaveIter(const int& saveIter) {this->SaveIter = saveIter;}

	inline void setPotCons(PotentialConsistency pc) {PotCons = pc;}
	inline void setTargetPotCons(const double& tpc) {TargetPotCons = tpc;}

	PotentialConsistency getPotCons() {return PotCons;}
	double getTargetPotCons() {return TargetPotCons;}

	string 	getLoadQ();
	string  getLoadDir();
	string  getSaveDir();
	string  getEndCriterion();
    string  getMeshName();  // returns mesh filename, with MeshNumber appended
	string  getMeshFileNameOnly(); // returns mesh filename, without directory

	double getMaxError();
	inline double getdt() {return dt;}
	inline double getTargetdQ(){return TargetdQ;}
	inline double getMaxdt(){return dtLimits[1];}
	inline double getMindt(){return dtLimits[0];}
	double getCurrentTime();
	double getCurrentChange();
	double getEndValue();

	double getStretchVectorX();
	double getStretchVectorY();
	double getStretchVectorZ();

	inline double getEnergyRegionX(){return EnergyRegion[0];}
	inline double getEnergyRegionY(){return EnergyRegion[1];}
	inline double getEnergyRegionZ(){return EnergyRegion[2];}
	int getSaveIter(){ return SaveIter;}

	int	getCurrentIteration();
	int getOutputEnergy();
	int getOutputFormat();

	void IncrementCurrentTime();
	void IncrementCurrentIteration();
	bool IsRunning();
    inline bool IsAssembleMatrix() {return AssembleMatrix;}

    inline bool IsMeshModified(){return MeshModified;}
    inline void setMeshModified(const bool& mm){MeshModified = mm;}
    inline int  getMeshNumber(){ return MeshNumber;}
    inline void setMeshNumber(const int& n) {MeshNumber = n;}
    inline void IncrementMeshNumber(){ MeshNumber++;}



};
#endif
