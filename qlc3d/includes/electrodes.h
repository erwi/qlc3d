#ifndef ELECTRODES_H
#define ELECTRODES_H
#include <stdlib.h>
#include <vector>
#include <list>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <eventlist.h>
#include <limits>

class SwitchingInstance {
    /*!A SwitchingInstance contains information needed to describe a single
    electrode switching event*/
public:
    SwitchingInstance(const double t, const double pot, const size_t en):
        time(t), potential(pot), electrodeNumber(en) {
    }
    static const size_t UNIFORM_E_FIELD;
    double time;
    double potential;
    size_t electrodeNumber;
};

class Electrodes {
    /*!A class with parameters related to potential calculations and electrodes*/
private:
    std::vector<double> currentElectrodePotentials;    // keeps current potential values for each electrode
    int nElectrodes;
public:

    std::vector <double> eps_dielectric;
    double EField[3];                   // Contains the x,y,z components of a uniform E-field

    Electrodes();
    ~Electrodes();
    void printElectrodes()const;
    void AddElectrode();
    double getDielectricPermittivity(int i) const;  // gets relative dielectric permittivity of dielectric#i
    bool getCalcPot() const;
    bool isEField()const ;              // returns true if uniform E-field has been defined
    void WriteElectrodes(FILE *fid)const;           // writes electrode settings to file fid


    double getCurrentElectrodePotential(const size_t &eln) const; // get current potential for electrode eln
    size_t getnElectrodes() const { return nElectrodes; }
    void setnElectrodes(const size_t &numE);
    void setEField(const std::vector<double> &vec3);
    void setElectrodePotential(const size_t &electrodeNumber, const double &potentialValue);
    //void setImplicitVariables(); // SETS FLAGS THAT DEPEND ON EXPLICITLY DEFINED VALUES

};
#endif
