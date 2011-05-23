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

public:

	vector<double> Potential;
	vector<double> Time;

	Electrode();
	void setPotential(std::vector<double> pot);
	void setTime(std::vector<double> tme);
	void PrintElectrode();
	int getnTimes();
};

class Electrodes
{
	private:
		bool CalcPot; // flag for checking whether it is necessary to calculate the potential


	public:
		int nElectrodes;
		vector<Electrode*> E;
		vector <double> eps_dielectric;
		//Electrode E[9];

		//void Electroodes();
		Electrodes();
		~Electrodes();
		void PrintElectrodes();
		void AddElectrode();
		void AddElectrode(Electrode* El);

		//void AddDielectricPermittivity(double eps); // adds relative dielectric permittivity
		double getDielectricPermittivity(int i); 	// gets realtive dielectric permittivity of dielectric#i

		void setCalcPot(bool yn);
		bool getCalcPot();

		int getnElectrodes();
		void WriteElectrodes(FILE* fid);			// writes electrode settings to file fid

};
#endif
