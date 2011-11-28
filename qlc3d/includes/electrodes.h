#ifndef ELECTRODES_H
#define ELECTRODES_H
#include <stdlib.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using std::vector;

class Electrode
{
private:
	int nTimes;
	double currentPotential;
public:

	vector<double> Potential;	// potential in Volts
	vector<double> Time;		// time in seconds

	Electrode();
	void setPotential(std::vector<double> pot);
	void setTime(std::vector<double> tme);
	void PrintElectrode();
	int getnTimes();
	double getCurrentPotential() {return currentPotential; }
	void setCurrentPotential(const double& cp) { currentPotential = cp ; }
	//void removeSwitching( const int& idx); // removes idx'th swithing event
};

class Electrodes
{
	private:
		bool CalcPot; // flag for checking whether it is necessary to calculate the potential


	public:
		int nElectrodes;
		vector<Electrode*> E;
		vector <double> eps_dielectric;
                double EField[3];       // Contains the x,y,z components of a uniform E-field


		//void Electroodes();
		Electrodes();
		~Electrodes();
		void printElectrodes();
		void AddElectrode();
		void AddElectrode(Electrode* El);

		//void AddDielectricPermittivity(double eps); // adds relative dielectric permittivity
		double getDielectricPermittivity(int i); 	// gets realtive dielectric permittivity of dielectric#i

		void setCalcPot(bool yn);
		bool getCalcPot();
                bool isEField();    // returns true if uniform E-field has been defined
		int getnElectrodes();
		void WriteElectrodes(FILE* fid);			// writes electrode settings to file fid

};
#endif
