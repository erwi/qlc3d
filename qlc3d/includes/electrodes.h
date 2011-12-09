#ifndef ELECTRODES_H
#define ELECTRODES_H
#include <stdlib.h>
#include <vector>
#include <list>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <eventlist.h>







class Electrodes
{
    private:
        bool CalcPot; // flag for checking whether it is necessary to calculate the potential
        std::vector<double> potentials_;    // keeps current potential values for each electrode

    public:
        int nElectrodes;
  //      vector<Electrode> E;       // THIS USED TO BE VECTOR OF POINTERS - CHECK ERRORS
        std::vector <double> eps_dielectric;
        double EField[3];       // Contains the x,y,z components of a uniform E-field

        Electrodes();
        ~Electrodes();
        void printElectrodes()const;
        void AddElectrode();
        //void AddElectrode(Electrode* El);
        double getDielectricPermittivity(int i) const; 	// gets realtive dielectric permittivity of dielectric#i
        void setCalcPot(bool yn);
        bool getCalcPot()const ;
        bool isEField()const ;    // returns true if uniform E-field has been defined
        //int getnElectrodes() const;
        void WriteElectrodes(FILE* fid)const;			// writes electrode settings to file fid

// NEW METHODS
        double getCurrentElectrodePotential(const size_t& eln) const; // get current potential for electrode eln
        size_t getnElectrodes()const {return nElectrodes;}
        void setnElectrodes(const size_t& numE)
        {
            potentials_.resize(numE,0.0);
            nElectrodes = numE;
        }
        void setElectrodePotential( const size_t& eNum, const double& pot );
};
#endif
