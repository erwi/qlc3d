#include <electrodes.h>
#include <algorithm>
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

//===========================================================
//
//		ELECTRODES
//
//

Electrodes::Electrodes():
    CalcPot(false),
    nElectrodes(0)
{

    eps_dielectric.push_back(1.0);
    EField[0] = 0.0;
    EField[1] = 0.0;
    EField[2] = 0.0;

}
Electrodes::~Electrodes(){

    //std::vector<Electrode*>::iterator itr;

//	for(itr = E.begin() ; itr!=E.end() ; itr++)
//		delete (*itr);

}

/*
void Electrodes::AddElectrode()
{
	// Adds new empty electrode
	nElectrodes++;
	E.push_back(new Electrode);
}
*/
//void Electrodes::AddElectrode(Electrode* El)
//{
    // adds non-empty electrode
//	nElectrodes++;
//    E.push_back( El );
//}
 
/*
void Electrodes::AddDielectricPermittivity(double eps){	eps_dielectric.push_back(eps); }
*/

double Electrodes:: getDielectricPermittivity(int i)const
{
	if (i < (int) eps_dielectric.size() )
		return eps_dielectric[i];
	else
		{
			printf("error - Electrodes::getDielectricPermittivity - trying to access dielectric %i, when only %i dielectric are defined - bye !\n", i , (int) eps_dielectric.size() );
			exit(1);
		}
}

void Electrodes::printElectrodes()const
{
/*
    for (int i = 0; i < nElectrodes; i++)
    {
        printf("E%i:\n",i+1);
        E[i]->PrintElectrode();
    }

    std::cout << "eps_dielectric = " << std::endl;
    for (int i = 0 ; i < (int) eps_dielectric.size() ; i++)
        std::cout <<" " << eps_dielectric[i];

    std::cout << std::endl;

    std::cout << "EField = [" << EField[0] <<"," <<EField[1] << "," << EField[2] <<"] V/um" << std::endl;
    */
}
void Electrodes::setCalcPot(bool yn)	{	CalcPot = yn;}
bool Electrodes::getCalcPot()const		{	return CalcPot;}
//int Electrodes::getnElectrodes()const	{   return nElectrodes;}
bool Electrodes::isEField()const
{
    // CHECKS WHETHER UNIFORM E-FIELD HAS BEEN DEFINED

    if ( ( this->EField[0]!= 0.0) ||
         ( this->EField[1]!= 0.0) ||
         ( this->EField[2]!= 0.0) )
        return true;
    else
        return false;


}

void Electrodes::WriteElectrodes(FILE* fid) const
{

}

void Electrodes::setElectrodePotential(const size_t &eNum, const double &pot)
{
// SETS THE CURRENT POTENTIAL OF ELECTRODE eNum TO VALUE pot

    if ( eNum >= getnElectrodes() )
    {
        printf("error in %s, cant set potential of electrode %u.\n" , __func__ , eNum + 1);
        printf("Total number of electrodes is %u\n", getnElectrodes() );
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
        printf("error in %s, can't get potential of electrode %u\n",__func__, eNum + 1 );
        printf("Total number of electrodes is %u\n", getnElectrodes() );
        exit(1);
    }
    return potentials_[eNum];
}
