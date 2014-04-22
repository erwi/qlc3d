#include <simu.h>
#include <algorithm>
#include <omp.h>
#include "stringenum.h"
#include <reader.h>
#include <globals.h>
Simu::Simu():
PotCons(Off),
    TargetPotCons(1e-3),
    MaxError(1e-7),
    CurrentTime(0),
    CurrentChange(0),
    TargetdQ(1e-3),
    Maxdt(1e-4),
    EndCriterion(Time),
    LoadQ(""),
    CurrentDir(""),
    SaveDir("res"),
    LoadDir(""),
    EndValue(1e-3),
    EndValue_orig(0),
    CurrentIteration(0),

    AssembleMatrix(true),
    MeshModified(true),
    MeshNumber(0),
    OutputEnergy(0),
    OutputFormat(SIMU_OUTPUT_FORMAT_BINARY),
    SaveIter(0),
    SaveFormat( LCview ),
    numAsseblyThreads(0),       // 0 MEANS USE ALL AVAILABLE, AND IS DEFAULT
    numMatrixSolverThreads(0),
    MeshName(""),
    dt(1e-9)
{

    dtLimits[0]         = 1e-9;
    dtLimits[1]         = 1e-3;

    dtFunction[0]       = 0.5;	// value of R where S = 0
    dtFunction[1]       = 0.8; // S = R (min)
    dtFunction[2]       = 1.2; // S = R (max)
    dtFunction[3]       = 10.0;  // S = 2

    restrictedTimeStep     = false;

    StretchVector[0]    = 1.0;	StretchVector[1]	= 1.0;	StretchVector[2]	= 1.0;
    EnergyRegion [0]	= 0.0;	EnergyRegion [1]	= 0.0;	EnergyRegion [2]	= 0.0;
    OutputFormat		= SIMU_OUTPUT_FORMAT_BINARY;

    RegularGridSize[0] = 0; RegularGridSize[1] = 0; RegularGridSize[2] = 0;

    // SET THREAD COUNTS TO MAXIMUM DETECTED BY OPENMP
#ifndef DEBUG
    numAsseblyThreads = omp_get_max_threads();
    numMatrixSolverThreads = omp_get_max_threads();
#endif

    // SET MATRIX SOLVER TO AUTO. THIS WILL BE CHANGED LATER
    // IN THE PROGRAM, DEPENDING ON THE PROPERTIES OF THE MATRIX
    // OR AS SPECIFIED IN THE SETTINGS FILE
    QMatrixSolver = Auto;

    //std::vector<int> temp;
    //temp.push_back(Simu::None);
    //temp.push_back(Simu::LCview);
    //temp.push_back(Simu::RegularVTK);
    //temp.push_back(Simu::RegularVecMat);
    //temp.push_back(Simu::DirStackZ);
    //temp.push_back(Simu::LCviewTXT);
    //validSaveFormatStrings.setIntValues(temp);

}
void Simu::PrintSimu(){}

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

void Simu::setdtForced(const double &dt)
{
// FORCES VALUE OF DT
#ifdef DEBUG
    if (dt<=0)
    {
        printf("warning in %s, trying to set dt to %e\n", __func__, dt);
        exit(1);
    }
#endif

    this->dt = dt;
}


void Simu::setdtLimits(const vector<double> &vec2) {
    double min = vec2.at(0);
    double max = vec2.at(1);
    if ((min>max) || (min<=0)) {
        cerr << "error - Simu::setdtLimits - invalid dtLimits - bye! \n" << endl;
        exit(qlc3d_GLOBALS::ERROR_CODE_BAD_SETTINGS_FILE);
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
void Simu::setdtFunction(const vector<double> &vec4) {
    dtFunction[0] = vec4.at(0);
    dtFunction[1] = vec4.at(1);
    dtFunction[2] = vec4.at(2);
    dtFunction[3] = vec4.at(3);
}
void Simu::setRegularGridSize(const vector<unsigned int> &vec3) {
    this->RegularGridSize[0] = vec3.at(0);
    this->RegularGridSize[1] = vec3.at(1);
    this->RegularGridSize[2] = vec3.at(2);
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

    StringEnum<Simu::EndCriteria> validator("EndCriterion","Iterations,Time,Change" );
   try {
       EndCriterion = validator.getEnumValue(ec);
   }
   catch (...){
       validator.printErrorMessage(ec);
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
    if ( EndCriterion == Change )
    {
        // Force more iterations by multiplying current change by 1000.
        this->setCurrentChange( this->getEndValue()*1000.f );
    }
    else if (EndCriterion == Iterations )
    {
        this->setEndValue( this->getEndValue() + this->EndValue_orig ); // end counter by original value defined i nsettings file
    }
    else if (EndCriterion == Time )
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

void Simu::setAsseblyThreadCount(unsigned int numT)
{
    numAsseblyThreads = numT;
#ifndef DEBUG
    if (numT == 0 )
        numAsseblyThreads = omp_get_max_threads();
#endif
}

void Simu::setMatrixSolverThreadCount(unsigned int numT)
{
    numMatrixSolverThreads = numT;
#ifndef DEBUG
    if (numT == 0)
        numMatrixSolverThreads = omp_get_max_threads();
#endif
}

void Simu::setQMatrixSolver(const string &solver)
{
    StringEnum<Simu::QMatrixSolvers> validator("QMatrixSolver","Auto,PCG,GMRES");

    try{
        QMatrixSolver = validator.getEnumValue(solver);
    }
    catch(...)
    {
        validator.printErrorMessage(solver);
        exit(1);
    }
}

void Simu::setStretchVector(const std::vector<double> vec3) {
    for (int i = 0; i < 3; i++)
        StretchVector[i] = vec3.at(i);
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
    if ( EndCriterion == Iterations)
    {
        if (CurrentIteration > (int) EndValue) return false;
    }
    else if ( EndCriterion == Time )
    {
        if (CurrentTime > EndValue) return false;
    }
    else if (EndCriterion == Change)
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
        printf("error - Simu::IsRunning() - unknowns EndCriterion - bye!\n");
        exit(1);
    }
    return true;
}

void Simu::addSaveFormat(std::string format)
{
    /*! Add output file firmat to write*/
    int SaveFormatValues[6] = {None,LCview,RegularVTK,RegularVecMat,DirStackZ, LCviewTXT};
    StringEnum<Simu::SaveFormats> validator("SaveFormat",
                                            "None,LCView,RegularVTK,RegularVecMat,DirStackZ,LCviewTXT",
                                            SaveFormatValues);
    try{
        SaveFormat = SaveFormat | validator.getEnumValue(format);
        std::cout << "added SaveFormat : " << format << std::endl;
    }
    catch(...){
        validator.printErrorMessage(format);
        exit(1);
    }
}

void Simu::setSaveFormats(const std::vector<string> saveFormats) {
    for (const auto & s :saveFormats)
        this->addSaveFormat(s);
}


void Simu::readSettingsFile(Reader &reader) {
/*! loads values from settings file reader */
    try {
        MeshName = reader.getValueByKey<string>("MeshName");
        std::string key;
        // Read optional string values
        if (reader.containsKey("LoadQ"))
            this->LoadQ = reader.get<string>();
        if (reader.containsKey("SaveDir"))
            this->SaveDir = reader.get<string>();
        if (reader.containsKey("QMatrixSolver"))
            this->setQMatrixSolver(reader.get<string>());
        // string arrays
        if (reader.containsKey("SaveFormat"))
            this->setSaveFormats(reader.get<vector<string> >());


        // Read optional scalar values
        if (reader.containsKey("EndValue"))
            this->setEndValue(reader.get<double>());
        if (reader.containsKey("dt"))
            this->setdt(reader.get<double>());
        if (reader.containsKey("TargetdQ"))
            this->setTargetdQ(reader.get<double>());
        if (reader.containsKey("Maxdt"))
            this->setMaxdt(reader.get<double>());
        if (reader.containsKey("MaxError"))
            this->setMaxError(reader.get<double>());

        if (reader.containsKey("OutputEnergy"))
            this->setOutputEnergy(reader.get<int>());
        if (reader.containsKey("OutputFormat"))
            this->setOutputFormat(reader.get<int>());
        if (reader.containsKey("SaveIter"))
            this->setSaveIter(reader.get<int>());
        if (reader.containsKey("NumAssemblyThreads"))
            this->setAsseblyThreadCount(reader.get<unsigned int>());
        if (reader.containsKey("NumMatrixSolverThreads"))
            this->setMatrixSolverThreadCount(reader.get<unsigned int>());

        // Number arrays
        if (reader.containsKey("StretchVector"))
            this->setStretchVector(reader.get<vector<double> >());
        if (reader.containsKey("dtLimits"))
            this->setdtLimits(reader.get<vector<double> > ());
        if (reader.containsKey("dtFunction"))
            this->setdtFunction(reader.get<vector<double> >());
        if (reader.containsKey("RegularGridSize"))
            this->setRegularGridSize(reader.get<vector<idx> >());


    } catch (ReaderError e) {
       e.printError();
   }
}
