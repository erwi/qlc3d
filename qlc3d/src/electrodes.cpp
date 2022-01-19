#include <electrodes.h>
#include <algorithm>
#include <util/exception.h>
#include <fmt/format.h>

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
        RUNTIME_ERROR(fmt::format("Trying to access dielectic {} when only {} dielected materials exist.", i, eps_dielectric.size()));
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
        RUNTIME_ERROR(fmt::format("Specified electric field must have 3 components, got {}."))
    }
    this->EField[0] = vec3[0];
    this->EField[1] = vec3[1];
    this->EField[2] = vec3[2];
}

void Electrodes::setElectrodePotential(const size_t &electrodeNumber, const double &potentialValue) {
/*! SETS THE CURRENT POTENTIAL OF ELECTRODE electrodeNumbe TO VALUE potentialValue */
    if (electrodeNumber >= getnElectrodes()) {
        RUNTIME_ERROR(fmt::format("Electrode {} does not exist.", electrodeNumber + 1));
    }
    currentElectrodePotentials[electrodeNumber] = potentialValue;
}

double Electrodes::getCurrentElectrodePotential(const size_t &eNum) const {
// RETURNS THE CURRENT POTENTIAL OF ELECTRODE eln
    if ( eNum >= getnElectrodes() ) {
        RUNTIME_ERROR(fmt::format("Can't find potential for electrode {}. Potentials have been specified only "
                                  "for {} electrodes.", eNum + 1, getnElectrodes()));
    }
    return currentElectrodePotentials[eNum];
}
