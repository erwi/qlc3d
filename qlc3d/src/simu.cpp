#include <simu.h>
#include <algorithm>
#include <omp.h>
#include "stringenum.h"
#include <reader.h>
#include <globals.h>
#include <settings_file_keys.h>

// Define valid enum string keys
const vector<string> Simu::VALID_END_CRITERIA = {"Iterations", "Time", "Change"};
const vector<string> Simu::VALID_SAVE_FORMATS = {"None","LCView","RegularVTK","RegularVecMat","DirStackZ","LCviewTXT"};
const vector<string> Simu::VALID_Q_MATRIX_SOLVERS = {"Auto", "PCG", "GMRES"};
// Define values of default Simu parameters
const string Simu::DEFAULT_LOAD_Q = "";
const string Simu::DEFAULT_SAVE_DIR = "res";
const string Simu::DEFAULT_Q_MATRIX_SOLVER = "Auto";
const vector<string> Simu::DEFAULT_SAVE_FORMATS = {"LCView"};
const string Simu::DEFAULT_END_CRITERION = Simu::VALID_END_CRITERIA[1]; // "Time"
const double Simu::DEFAULT_END_VALUE = 1e-3;
const double Simu::DEFAULT_DT = 1e-9;
const double Simu::DEFAULT_TARGET_DQ = 1e-3;
const double Simu::DEFAULT_MAX_DT = 1e-4;
const double Simu::DEFAULT_MAX_ERROR = 1e-7;
const vector<double> Simu::DEFAULT_STRETCH_VECTOR = {1,1,1};
const vector<double> Simu::DEFAULT_DT_LIMITS = {Simu::DEFAULT_DT, Simu::DEFAULT_MAX_DT};
const vector<double> Simu::DEFAULT_DT_FUNCTION = {0.5, 0.8, 1.2, 10};
const vector<idx> Simu::DEFAULT_REGULAR_GRID_SIZE = {0,0,0};
const int Simu::DEFAULT_OUTPUT_ENERGY = 0;
const int Simu::DEFAULT_OUTPUT_FORMAT = SIMU_OUTPUT_FORMAT_BINARY; // TODO: Get rid of this
const int Simu::DEFAULT_SAVE_ITER = 1;
const int Simu::DEFAULT_NUM_ASSEMBLY_THREADS = 0;
const int Simu::DEFAULT_NUM_MATRIX_SOLVER_THREADS = 0;
Simu::Simu():
    PotCons(Off),
    TargetPotCons(1e-3),
    MaxError(DEFAULT_MAX_ERROR),
    CurrentTime(0),
    CurrentChange(0),
    TargetdQ(DEFAULT_TARGET_DQ),
    Maxdt(DEFAULT_MAX_DT),
    EndCriterion(Time),
    LoadQ(DEFAULT_LOAD_Q),
    CurrentDir(""),
    SaveDir(DEFAULT_SAVE_DIR),
    LoadDir(""),
    EndValue(DEFAULT_END_VALUE),
    EndValue_orig(0),
    CurrentIteration(0),

    AssembleMatrix(true),
    MeshModified(true),
    MeshNumber(0),
    OutputEnergy(0),
    OutputFormat(SIMU_OUTPUT_FORMAT_BINARY),
    SaveIter(0),
    SaveFormat(LCview),
    numAsseblyThreads(0),       // 0 MEANS USE ALL AVAILABLE, AND IS DEFAULT
    numMatrixSolverThreads(0),
    MeshName(""),
    dt(DEFAULT_DT) {
    // TODO: make dtLimits a std::vector
    std::copy(DEFAULT_DT_LIMITS.begin(), DEFAULT_DT_LIMITS.end(), dtLimits);
    // TODO: make dtFunction a std::vector;
    std::copy(DEFAULT_DT_FUNCTION.begin(), DEFAULT_DT_FUNCTION.end(), dtFunction);
    // TODO make StretchVector a std::vector
    std::copy(DEFAULT_STRETCH_VECTOR.begin(), DEFAULT_STRETCH_VECTOR.end(), StretchVector);
    // TODO : make RegularGridSize a std::vector
    std::copy(DEFAULT_REGULAR_GRID_SIZE.begin(), DEFAULT_REGULAR_GRID_SIZE.end(), RegularGridSize);
    restrictedTimeStep     = false;
    OutputFormat        = SIMU_OUTPUT_FORMAT_BINARY;
#ifndef DEBUG
    numAsseblyThreads = omp_get_max_threads();
    numMatrixSolverThreads = omp_get_max_threads();
#endif
    // SET MATRIX SOLVER TO AUTO. THIS WILL BE CHANGED LATER
    // IN THE PROGRAM, DEPENDING ON THE PROPERTIES OF THE MATRIX
    // OR AS SPECIFIED IN THE SETTINGS FILE
    QMatrixSolver = Auto;
}
void Simu::PrintSimu() {}

void Simu::setSaveDir(string savedir) {
    SaveDir = savedir;
}
void Simu::setLoadDir(string loaddir) {
    LoadDir = loaddir;
}
void Simu::setMaxError(double me) {
    if (me > 0)
        MaxError = me;
    else {
        printf("error - Simu:setMaxError - negative MaxError - bye!");
        exit(1);
    }
}
void Simu::setdt(const double &td) {
    // clamp time step between min max values if not goinf for steady state (dt=0)
    if (td > 0) {
        dt = td;
        if (td < dtLimits[0])
            dt = dtLimits[0];
        if (td > dtLimits[1])
            dt = dtLimits[1];
    } else {
        dt = 0;
    }
}

void Simu::setdtForced(const double &dt) {
    // FORCES VALUE OF DT
#ifdef DEBUG
    if (dt <= 0) {
        printf("warning in %s, trying to set dt to %e\n", __func__, dt);
        exit(1);
    }
#endif
    this->dt = dt;
}


void Simu::setdtLimits(const vector<double> &vec2) {
    double min = vec2.at(0);
    double max = vec2.at(1);
    if ((min > max) || (min <= 0)) {
        cerr << "error - Simu::setdtLimits - invalid dtLimits - bye! \n" << endl;
        exit(qlc3d_GLOBALS::ERROR_CODE_BAD_SETTINGS_FILE);
    }
    dtLimits[0] = min;
    dtLimits[1] = max;
}

void Simu::setTargetdQ(const double &dq) {
    if (dq <= 0) {
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

void Simu::getdtFunction(double *f) {
    f[0] = dtFunction[0];
    f[1] = dtFunction[1];
    f[2] = dtFunction[2];
    f[3] = dtFunction[3];
}

void Simu::setCurrentIteration(int i) {
    CurrentIteration = i;
}
void Simu::setMaxdt(double maxdt)   {
    Maxdt = maxdt;
}
void Simu::setEndCriterion(string ec) {
    //StringEnum<Simu::EndCriteria> validator("EndCriterion", "Iterations,Time,Change");
    StringEnum<Simu::EndCriteria> validator( SFK_END_CRITERION, Simu::VALID_END_CRITERIA);
    try {
        EndCriterion = validator.getEnumValue(ec);
    } catch (...) {
        validator.printErrorMessage(ec);
        exit(1);
    }
}
void Simu::setEndValue(double ev) {
    // this should prevent end-refinement from overwriting original value of endvalue
    // when resetEndCriterion is called
    if (EndValue == 0)
        EndValue_orig = ev;
    EndValue = ev;
}
void Simu::resetEndCriterion() {
    /*! Resets simu varibales so that simulation is started from beginning, while continuing
    with the current result. This is used after an end-refinement is performed to force simulation to do
    at least one more step.

    For example, is end criterion is MaxChange, Current chancge is
    set to something large so that simlation is effectively restarted
    */
    if (EndCriterion == Change) {
        // Force more iterations by multiplying current change by 1000.
        this->setCurrentChange(this->getEndValue() * 1000.f);
    } else if (EndCriterion == Iterations) {
        this->setEndValue(this->getEndValue() + this->EndValue_orig);   // end counter by original value defined i nsettings file
    } else if (EndCriterion == Time) {
        this->setEndValue(this->EndValue + this->EndValue_orig);
    }
}

void Simu::setCurrentTime(double ct)        {
    CurrentTime = ct;
}
void Simu::setAssembleMatrix(bool yn)       {
    AssembleMatrix = yn;
}
void Simu::setLoadQ(string qbackup)         {
    LoadQ  = qbackup;
}
void Simu::setCurrentChange(double ch)      {
    CurrentChange = ch;
}
void Simu::setMeshName(string meshname)     {
    MeshName = meshname;
}
void Simu::setOutputEnergy(int ope) {
    OutputEnergy = ope;
}
void Simu::setOutputFormat(int opf) {
    if ((opf == SIMU_OUTPUT_FORMAT_BINARY) ||
            (opf == SIMU_OUTPUT_FORMAT_TEXT)) {
        OutputFormat = opf;
    } else {
        printf("error - Simu::setOutputFormat, %i is not a valid output format - bye!\n", opf);
        fflush(stdout);
        exit(1);
    }
}// end setOutputFormat

void Simu::setAsseblyThreadCount(unsigned int numT) {
    numAsseblyThreads = numT;
#ifndef DEBUG
    if (numT == 0)
        numAsseblyThreads = omp_get_max_threads();
#endif
}

void Simu::setMatrixSolverThreadCount(unsigned int numT) {
    numMatrixSolverThreads = numT;
#ifndef DEBUG
    if (numT == 0)
        numMatrixSolverThreads = omp_get_max_threads();
#endif
}

void Simu::setQMatrixSolver(const string &solver) {
    StringEnum<Simu::QMatrixSolvers> validator(SFK_Q_MATRIX_SOLVER, Simu::VALID_Q_MATRIX_SOLVERS);
    try {
        QMatrixSolver = validator.getEnumValue(solver);
    } catch (...) {
        validator.printErrorMessage(solver);
        exit(1);
    }
}

void Simu::setStretchVector(const std::vector<double> vec3) {
    for (int i = 0; i < 3; i++)
        StretchVector[i] = vec3.at(i);
}




string Simu::getMeshName() const {
    // returns mesh name with mesh number appended
    std::string meshname(this->MeshName);
    std::stringstream ss;
    ss << MeshNumber;
    std::string num;
    ss >> num;
    size_t pos = meshname.find_last_of("."); // position of separator point
    meshname.insert(pos, num);
    return meshname;
}

string Simu::getMeshFileNameOnly() {
    // GET FILE NAME OF CURRENT MESH, REMOVE DIRECTORY INFORMATION
    std::string meshname(this->getMeshName());
    size_t pos = meshname.find_last_of("/");
    if (pos != std::string::npos)
        meshname = meshname.substr(pos + 1 , meshname.length() + 1 - pos);
    return meshname;
}

void Simu::IncrementCurrentTime() {
    CurrentTime += dt;
}
void Simu::IncrementCurrentIteration()      {
    CurrentIteration ++;
}
bool Simu::IsRunning()const {
    if (EndCriterion == Iterations) {
        if (CurrentIteration > (int) EndValue) return false;
    } else if (EndCriterion == Time) {
        if (CurrentTime > EndValue) return false;
    } else if (EndCriterion == Change) {
        if (dt > 0) {
            printf("\tdQ / dt = %1.3e , EndValue = %1.3e\n", fabs(getCurrentChange()), EndValue);
            if (fabs(getCurrentChange()) <= EndValue)    // if dQ / dt < maxchange
                return false;
        }
        if (fabs(getCurrentChange()) <= EndValue)
            return false;
    } else {
        printf("error - Simu::IsRunning() - unknowns EndCriterion - bye!\n");
        exit(1);
    }
    return true;
}

void Simu::addSaveFormat(std::string format) {
    /*! Add output file firmat to write*/
    StringEnum<Simu::SaveFormats> validator(SFK_SAVE_FORMAT,
                                            Simu::VALID_SAVE_FORMATS,
                                            {None, LCview, RegularVTK, RegularVecMat,
                                             DirStackZ, LCviewTXT});
    try {
        SaveFormat = SaveFormat | validator.getEnumValue(format);
        std::cout << "added SaveFormat : " << format << std::endl;
    } catch (...) {
        validator.printErrorMessage(format);
        exit(1);
    }
}

void Simu::setSaveFormats(const std::vector<string> saveFormats) {
    for (const auto &s : saveFormats)
        this->addSaveFormat(s);
}

