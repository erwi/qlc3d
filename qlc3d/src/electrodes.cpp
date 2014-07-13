#include <electrodes.h>
#include <algorithm>
#include <iostream>
using std::cerr;
using std::endl;

// SET MAGIC CONSTANT FOR INDICATING UNIFORM ELECTRIC FIELD
const size_t SwitchingInstance::UNIFORM_E_FIELD = numeric_limits<size_t>::max();

//===========================================================
//
//		ELECTRODES
//
//

Electrodes::Electrodes():nElectrodes(0) {
    eps_dielectric.push_back(1.0);
    EField[0] = 0.0;
    EField[1] = 0.0;
    EField[2] = 0.0;
}
Electrodes::~Electrodes() { }

double Electrodes:: getDielectricPermittivity(int i) const {
	if (i < (int) eps_dielectric.size() )
		return eps_dielectric[i];
    else {
        cerr << "error - Electrodes::getDielectricPermittivity - trying to access dielectric "
             << i << "  when only "<< eps_dielectric.size() << " exist " << endl;
        std::exit(1);
    }
}


bool Electrodes::getCalcPot() const {
/*! Returns whether potential calculations are necessary */
    if (isEField()) { // if have fixed E-field
        return false;
    }
    if (0 == getnElectrodes()) { // if have no electrodes
        return false;
    }
    // Otherwise
    return true;
}

bool Electrodes::isEField() const {
    // CHECKS WHETHER UNIFORM E-FIELD HAS BEEN DEFINED
    if ( ( this->EField[0]!= 0.0) ||
         ( this->EField[1]!= 0.0) ||
         ( this->EField[2]!= 0.0) )
        return true;
    else
        return false;
}
void Electrodes::setEField(std::vector<double> const &vec3) {
    if (vec3.size() != 3) {
        cerr << "error - specified electric filed must havse 3 components, got " << vec3.size() << endl;
        std::exit(1);
    }
    this->EField[0] = vec3[0];
    this->EField[1] = vec3[1];
    this->EField[2] = vec3[2];
}

void Electrodes::setElectrodePotential(const size_t &electrodeNumber, const double &potentialValue) {
/*! SETS THE CURRENT POTENTIAL OF ELECTRODE electrodeNumbe TO VALUE potentialValue */
    if ( electrodeNumber >= getnElectrodes() ) {
        std::cerr << "error in " << __func__ << " Electrode " << electrodeNumber + 1 << " does not exist\nbye!" << std::endl;
        std::exit(1);
    }
    currentElectrodePotentials[electrodeNumber] = potentialValue;
}

double Electrodes::getCurrentElectrodePotential(const size_t &eNum) const {
// RETURNS THE CURRENT POTENTIAL OF ELECTRODE eln
    if ( eNum >= getnElectrodes() ) {
        printf("\nerror in %s, can't find potential for electrode %u.\n",__func__, (unsigned int) eNum + 1 );
        printf("Potentials have been specified for only %u electrodes.\n", (unsigned int) getnElectrodes());
        printf("Check your settings file - bye!\n");
        exit(1);
    }
    return currentElectrodePotentials[eNum];
}
