#include "../includes/electrodes.h"

Electrode::Electrode()
{
	nTimes = 0;
	Potential.clear();
	Time.clear();


}


void Electrode::PrintElectrode()
{

	int npots = Potential.size();
	int ntimes = Time.size();
	printf("     number of potentials, times : %i, %i\n",npots,ntimes);

        if (npots!=ntimes){
		printf("Problem with configuration file - number of potentials must equal number of times, bye!\n");
		exit(1);
	}

        if (npots>0){
            printf("     Pots:");
            for (int i = 0; i<npots ; i++)
		printf(" %f", Potential[i]);

            printf("\n     Time:");
            for (int i = 0; i<npots ; i++)
		printf(" %f", Time[i]);

		printf("\n");
	} // end if npots>0

}
int Electrode::getnTimes()
{
	return Time.size();

}
void Electrode::setPotential(std::vector<double> pot){

    Potential = pot;

}

void Electrode::setTime(std::vector<double> tme){
    Time = tme;
    nTimes = Time.size();
    if (getnTimes() != (int) Potential.size() ){
        std::cout << "error, Electrode times don't seem to match potentials - bye!" << std::endl;
        exit(1);
    }

}



//===========================================================
//
//		ELECTRODES
//
//

Electrodes::Electrodes()
{
	nElectrodes =0;
	CalcPot = false;
        eps_dielectric.push_back(1.0);
}
Electrodes::~Electrodes(){
	std::vector<Electrode*>::iterator itr;

	for(itr = E.begin() ; itr!=E.end() ; itr++)
		delete (*itr);

}

void Electrodes::AddElectrode()
{
	//printf("adding Electrode to electrodes");
	nElectrodes++;
	E.push_back(new Electrode);
	//printf(" OK\n");
}
void Electrodes::AddElectrode(Electrode* El)
{
    nElectrodes++;
    E.push_back( El );
    //std::cout << "added";

}

/*
void Electrodes::AddDielectricPermittivity(double eps){	eps_dielectric.push_back(eps); }
*/
double Electrodes:: getDielectricPermittivity(int i)
{
	if (i < (int) eps_dielectric.size() )
		return eps_dielectric[i];
	else
		{
			printf("error - Electrodes::getDielectricPermittivity - trying to access dielectric %i, when only %i dielectric are defined - bye !\n", i , (int) eps_dielectric.size() );
			exit(1);
		}
}

void Electrodes::PrintElectrodes()
{
	for (int i = 0; i < nElectrodes; i++)
	{
		printf("E%i:\n",i+1);
		E[i]->PrintElectrode();
	}

    std::cout << "eps_dielectric = " << std::endl;
    for (int i = 0 ; i < (int) eps_dielectric.size() ; i++)
        std::cout <<" " << eps_dielectric[i];

    std::cout << std::endl;
}
void Electrodes::setCalcPot(bool yn)	{	CalcPot = yn;}
bool Electrodes::getCalcPot()			{	return CalcPot;}
int Electrodes::getnElectrodes()		{   return nElectrodes;}

void Electrodes::WriteElectrodes(FILE* fid)
{

	if (fid!=NULL)
	{
		fprintf(fid,"#======================\n");
		fprintf(fid,"#  ELECTRIC STUFF\n" );
		fprintf(fid,"#======================\n\n");
			for (int i = 0 ; i < getnElectrodes() ; i ++ ) // write electrode #i
			{

				char str[10];
				sprintf(str, "E%i", i+1);


				fprintf(fid,"\t%s.Time = [%2.4f", str, E[i]->Time[0]);
				for (int t = 1 ; t < (int) E[i]->Time.size() ; t++ )
				{	fprintf(fid,",%2.4f",E[i]->Time[t]);	}
				fprintf(fid,"]\n\t%s.Pot  = [%f",str , E[i]->Potential[0]);
				for (int p = 1 ; p < (int) E[i]->Potential.size() ; p++)
				{	fprintf(fid,",%2.4f",E[i]->Potential[p]);	}
				fprintf(fid,"]\n\n");

			}

		if (eps_dielectric.size()>0) // if dielectric permittivities are defined
		{
			fprintf(fid, "\teps_dielectric = [%2.4f",eps_dielectric[0]);
			for (int e = 1 ; e < (int) eps_dielectric.size() ; e++)
				fprintf(fid, ",%2.4f",eps_dielectric[e]);

			fprintf(fid, "]\n");

		}


	}
}
