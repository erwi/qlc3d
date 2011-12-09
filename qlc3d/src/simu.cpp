#include <simu.h>



Simu::Simu()
{
    PotCons				= Loop;
    TargetPotCons		= 1e-3;
    EndValue 			= 0;
    EndValue_orig       = 0;
    MeshName 		= "";
    LoadQ		= "";
    CurrentDir          = "";
    SaveDir		= "res";
    SaveIter		= 1;
    MeshNumber          = 0;
    LoadDir		= "";
    EndCriterion        = "Time";
		
    dt                  = 1e-6;
    dtLimits[0]         = 1e-9;
    dtLimits[1]         = 1e-3;
    TargetdQ            = 1e-3;

    dtFunction[0]       = 0.5;	// value of R where S = 0
    dtFunction[1]       = 0.8; // S = R (min)
    dtFunction[2]       = 1.2; // S = R (max)
    dtFunction[3]       = 10.0;  // S = 2

    Maxdt               = 1e-4; // deprecated, use dtLimits[1]
    restrictedTimeStep     = false;

    CurrentTime 		= 0;
    CurrentIteration 	= 0;
    CurrentChange		= 1;
    setMaxError(1e-7);
    AssembleMatrix 		= true;
    MeshModified        = true; // Always assume mesh has changed -> output mesh file
    OutputEnergy		= 0;
    StretchVector[0]    = 1.0;	StretchVector[1]	= 1.0;	StretchVector[2]	= 1.0;
    EnergyRegion [0]	= 0.0;	EnergyRegion [1]	= 0.0;	EnergyRegion [2]	= 0.0;
    OutputFormat		= SIMU_OUTPUT_FORMAT_BINARY;
		//OutputEnergy_fid	= NULL;
}
void Simu::PrintSimu(){}

void Simu::WriteSimu(FILE* fid)
{
    if (fid != NULL)
    {
        fprintf(fid,"\tMeshName =  %s\n",MeshName.c_str());				// MESHNAME
        fprintf(fid,"\tEndCriterion = %s\n",EndCriterion.c_str());		// ENDCRITERION
        fprintf(fid,"\tEndValue = %2.4e\n",EndValue);							// ENDVALUE

        fprintf(fid,"\tdt = %2.4e\n",dt);					// DT
        fprintf(fid,"\tdtLimits = [%2.4e,%2.4e]\n" , dtLimits[0], dtLimits[1] );
        fprintf(fid,"\tTargetdQ = %2.4e\n",TargetdQ);
        fprintf(fid,"\tMaxdt = %2.4e\n",Maxdt);					// MAXDT
        fprintf(fid,"\tMaxError = %2.4e\n",MaxError);				// MAXERROR
        printf("\tdtFunction = [%f,%f,%f,%f]\n", dtFunction[0], dtFunction[1], dtFunction[2], dtFunction[3]);

       	if ( !LoadQ.empty() )
	    fprintf(fid,"\tLoadQ = %s;\n",LoadQ.c_str());				// LOADQ
        if ( !getSaveDir().empty() )
            fprintf(fid,"\tSaveDir = %s\n", getSaveDir().c_str() );		// SAVEDIR

        fprintf(fid,"\tSaveIter = %i\n", getSaveIter() );
        fprintf(fid,"\tOuputFormat = %i\n", getOutputFormat() );				// OUTPUT FORMAT
        fprintf(fid,"\tOutputEnergy = %i\n",getOutputEnergy());			// OUTPUTENERGY
        fprintf(fid,"\tEnergyRegion = [%f, %f, %f]\n", EnergyRegion[0] , EnergyRegion[1] , EnergyRegion[2]); // ENERGY REGION VECTOR
        fprintf(fid,"\tStretchVector = [%f, %f, %f]\n", StretchVector[0] , StretchVector[1] , StretchVector[2] ); // MESH SCALING

        // POTENTIAL CONSISTENCY
        fprintf(fid,"\tPotCons = ");
        switch(PotCons){
            case(Off):
                fprintf(fid,"Off\n");
                break;
            case(Loop):
                fprintf(fid,"Loop\n");
                break;
            default:
                fprintf(fid,"error - Simu::PrintSimu, unknown PotCons\n");
        }

        fprintf(fid, "\tTargetPotCons = %e\n",TargetPotCons);

	} // end if fid != NULL
}
void Simu::setSaveDir(string savedir) {SaveDir = savedir;}
void Simu::setLoadDir(string loaddir) {LoadDir = loaddir;}
void Simu::setMaxError(double me){
	if (me>0)
		MaxError = me;
	else{
		printf("error - Simu:setMaxError - negative MaxError - bye!");
		exit(1);
	}
}
void Simu::setdt(double td){
// clamp time step between min max values if not goinf for steady state (dt=0)
    if (td > 0){
    	dt = td;

	if ( td < dtLimits[0])
	    dt = dtLimits[0];
	if ( td > dtLimits[1])
	    dt = dtLimits[1];
    }
    else{
	dt = 0;
    }
}
void Simu::setdtLimits(const double &min, const double &max){
    if ((min>max) || (min<=0)){
	cout << "error - Simu::setdtLimits - invalid dtLimits - bye! \n" << endl;
	exit(1);
    }

    dtLimits[0] = min;
    dtLimits[1] = max;

}

void Simu::setTargetdQ(const double &dq){
    if (dq <= 0 ){
	cout << "error - Simu::setTargetdQ - TargetdQ must be larget than 0 - bye!" << endl;
	exit(1);
    }
    TargetdQ = dq;
}
void Simu::setdtFunction(double *f){
    dtFunction[0] = f[0];
    dtFunction[1] = f[1];
    dtFunction[2] = f[2];
    dtFunction[3] = f[3];
}

void Simu::getdtFunction(double* f){
    f[0] = dtFunction[0];
    f[1] = dtFunction[1];
    f[2] = dtFunction[2];
    f[3] = dtFunction[3];
}

void Simu::setCurrentIteration(int i)
{
	CurrentIteration = i;
}
void Simu::setMaxdt(double maxdt)	{		Maxdt=maxdt;		}
void Simu::setEndCriterion(string ec)
{
	//std::cout << "ec :" << ec << std::endl;
	if (ec.compare("Iterations") == 0){

		EndCriterion = ec;
	}
	else if (ec.compare("Time") == 0){
		EndCriterion = ec;
	}
	else if (ec.compare("Change") == 0){
		EndCriterion = ec;
	}
	else{
		printf("error - Simu::setEndCriterion - %s is an unknown EndCriterion, bye!\n",ec.c_str());
		exit(1);
	}
}
void Simu::setEndValue(double ev)
{
	// this should prevent end-refinement from overwriting original value of endvalue
	// when resetEndCriterion is called
	if (EndValue == 0 ) 
		EndValue_orig = ev;
		
	EndValue = ev;
}
void Simu::resetEndCriterion()
{
/*! Resets simu varibales so that simulation is started from beginning, while continuing
with the current result. This is used after an end-refinement is performed to force simulation to do
at least one more step.

For example, is end criterion is MaxChange, Current chancge is
set to something large so that simlation is effectively restarted
*/
    if ( this->EndCriterion.compare("Change") == 0 )
    {
        // Force more iterations by multiplying current change by 1000.
        this->setCurrentChange( this->getEndValue()*1000.f );
    }
    else if (this->EndCriterion.compare("Iterations") == 0 )
    {
        this->setEndValue( this->getEndValue() + this->EndValue_orig ); // end counter by original value defined i nsettings file
    }
	else if (this->EndCriterion.compare("Time") == 0 )
	{
		this->setEndValue( this->EndValue + this->EndValue_orig );
	}
    

}

void Simu::setCurrentTime(double ct)		{   CurrentTime = ct;}
void Simu::setAssembleMatrix(bool yn)       {	AssembleMatrix = yn;}
void Simu::setLoadQ(string qbackup)			{	LoadQ  = qbackup;}
void Simu::setCurrentChange(double ch) 		{	CurrentChange = ch;}
void Simu::setMeshName(string meshname)     {	MeshName = meshname;}
void Simu::setOutputEnergy(int ope)
{
	OutputEnergy = ope;
	if (ope == 1) // if output energy == true, open file pointer
	{
		//string filename = "res/energy.txt";
		//OutputEnergy_fid = fopen(filename.c_str(),"wt");
		//if (OutputEnergy_fid == NULL)
		//{
		//	printf("error - Simu:setOutputEnergy(int) - could not open %s - bye!\n",filename.c_str());
		//	exit(1);
		//}
	}
}
void Simu::setOutputFormat(int opf){
	if ( 	(opf == SIMU_OUTPUT_FORMAT_BINARY) ||
			(opf == SIMU_OUTPUT_FORMAT_TEXT) )
			{
				OutputFormat = opf;
			}
			else
			{
				printf("error - Simu::setOutputFormat, %i is not a valid output format - bye!\n", opf);
				fflush(stdout);
				exit(1);
			}
}// end setOutputFormat
void Simu::setStretchVectorX(double sx)
{
	if (sx>0)
		StretchVector[0] = sx;
	else
	{
		printf("error - invalid StretchVector x component - bye!\n");
		exit(1);
	}
}
void Simu::setStretchVectorY(double sy)
{
	if (sy>0)
		StretchVector[1] = sy;
	else
	{
		printf("error - invalid StretchVector y component - bye!\n");
		exit(1);
	}
}
void Simu::setStretchVectorZ(double sz)
{
	if (sz>0)
		StretchVector[2] = sz;
	else
	{
		printf("error - invalid StretchVector z component - bye!\n");
		exit(1);
	}
}


string Simu::getMeshName() const
{
    // returns mesh name with mesh number appended
    std::string meshname(this->MeshName );
    std::stringstream ss;
    ss << MeshNumber;

    std::string num;
    ss >> num;

    size_t pos = meshname.find_last_of("."); // position of separator point

    meshname.insert(pos, num );

    return meshname;
}

string Simu::getMeshFileNameOnly(){
    // GET FILE NAME OF CURRENT MESH, REMOVE DIRECTORY INFORMATION
    std::string meshname( this->getMeshName() );
    size_t pos = meshname.find_last_of("/");
    if (pos != std::string::npos)
	meshname = meshname.substr(pos+1 , meshname.length()+1 - pos);

    return meshname;
}

void Simu::IncrementCurrentTime(){	CurrentTime+=dt;}
void Simu::IncrementCurrentIteration()		{	CurrentIteration ++;}
bool Simu::IsRunning()const{
    if ( EndCriterion.compare("Iterations") == 0 )
    {
        if (CurrentIteration > (int) EndValue) return false;
    }
    else if ( EndCriterion.compare("Time") == 0)
    {
        if (CurrentTime > EndValue) return false;
    }
    else if (EndCriterion.compare("Change") == 0)
    {
        if (dt >0)
        {
            printf("\tdQ / dt = %1.3e , EndValue = %1.3e\n", fabs( getCurrentChange() ), EndValue );
            if ( fabs(getCurrentChange()  ) <= EndValue) // if dQ / dt < maxchange
                return false;
        }

        if (fabs(getCurrentChange()) <= EndValue )
            return false;
	}
    else
    {
        printf("error - Simu::IsRunning() - unknowns EndCriterion %s , bye!\n",EndCriterion.c_str());
        exit(1);
    }
    return true;
}


