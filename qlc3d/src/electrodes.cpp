#include <electrodes.h>
#include <algorithm>
#include <iostream>
using std::cerr;
using std::endl;
/*
Electrode::Electrode()
{
        //nTimes = 0;
        //Potential.clear();
        //Time.clear();
        //currentPotential = 0.;

}

void Electrode::addSwitching(const double &Time, const double &Potential)
{
    // ADDS NEW SWITCHING INSTANCE TO THIS ELECTRODE
    switching_.push_back( Switching(Time, Potential ) );
    switching_.sort();
}

void Electrode::printElectrode(const int& num) const
{

// DEBUG PRINTOUT OF THIS ELECTRODES SWITCHING INSTANCES
// NUM IS AN OPTIONAL PARAMETER WITH DEFAULT VALUE -1

    if ( num>=0 )
        printf("Electrode %i:\n", num);

    std::list<Switching> ::const_iterator itr = switching_.begin();
    for ( ; itr != switching_.end() ; itr++ )
    {
        printf("Time = %e, Pot = %e\n", itr->getTime(), itr->getPotential() );
    }

}
*/
// SET MAGIC CONSTANT FOR INDICATING UNIFORM ELECTRIC FIELD
const size_t SwitchingInstance::UNIFORM_E_FIELD = numeric_limits<size_t>::max();

//===========================================================
//
//		ELECTRODES
//
//

Electrodes::Electrodes():
    CalcPot(false),
    nElectrodes(0) {
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

//void Electrodes::setCalcPot(bool yn)	{	CalcPot = yn;}
bool Electrodes::getCalcPot()const		{	return CalcPot;}
//int Electrodes::getnElectrodes()const	{   return nElectrodes;}
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
        cerr << "error - specified electric filed must have 3 components, got " << vec3.size() << endl;
        std::exit(1);
    }
    this->EField[0] = vec3[0];
    this->EField[1] = vec3[1];
    this->EField[2] = vec3[2];
}

void Electrodes::setElectrodePotential(const size_t &eNum, const double &pot)
{
// SETS THE CURRENT POTENTIAL OF ELECTRODE eNum TO VALUE pot

    if ( eNum >= getnElectrodes() )
    {
        printf("error in %s, cant set potential of electrode %u.\n" , __func__ , (unsigned int) eNum + 1);
        printf("Total number of electrodes is %u\n", (unsigned int) getnElectrodes() );
        exit(1);
    }
    potentials_[eNum] = pot;

// CHECKS WHETHER POTENTIAL CALCULATION IS NEEDED
    CalcPot = false;
    double minPot = *std::min_element(potentials_.begin() , potentials_.end() );
    double maxPot = *std::max_element(potentials_.begin() , potentials_.end() );
    if (minPot != maxPot)
        CalcPot = true;


}

double Electrodes::getCurrentElectrodePotential(const size_t &eNum) const
{
// RETURNS THE CURRENT POTENTIAL OF ELECTRODE eln
    if ( eNum >= getnElectrodes() )
    {
        printf("\nerror in %s, can't find potential for electrode %u.\n",__func__, (unsigned int) eNum + 1 );
        printf("Potentials have been specified for only %u electrodes.\n", (unsigned int) getnElectrodes());
        printf("Check your settings file - bye!\n");
        exit(1);
    }
    return potentials_[eNum];
}

void Electrodes::setImplicitVariables()
{
    /*!
      SETS FLAGS THAT DEPEND ON EXPLICITLY DEFINED VALUES.
      THIS SHOULD BE CALLED AFTER OTHER VALUES HAVE BEEN SET.
    */

    // SET CalcPot FLAG THAT INDICATES WHETHER POTENTIAL
    // NEEDS TO BE CALCULATED AND REQUIRES A MATRIX
    if ( getnElectrodes() && !isEField() )
        CalcPot = true;
    else
        CalcPot = false;



}
