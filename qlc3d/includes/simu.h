#ifndef SIMU_H
#define SIMU_H
#include <string>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <set>
#include <globals.h>
#define SIMU_END_SIMULATION -20000000 // magic minus bignum
#define SIMU_OUTPUT_FORMAT_BINARY 	0
#define SIMU_OUTPUT_FORMAT_TEXT		1

class Reader; // forward declaration of settings file reader
enum PotentialConsistency {Off, Loop};
enum SimulationMode {TimeStepping, SteadyState};

using namespace std;
class Simu
{
    friend void readSimu(Simu& , Reader&);
public:
    enum SaveFormats {None,
                      LCview,
                      RegularVTK,
                      RegularVecMat, // REGULAR GRID, VECTORS, MATLAB FILE
                      DirStackZ, // REGULAR GRID, CSV-FILE WHERE EA
                      LCviewTXT  // LCview ASCII text file
                     };
    // RegularDirStackZ is written in a CommaSeparatedValue (CSV) text file where each
    // Row is a stack of directors along the z-axis.
    // The director components in each stack are interleaved in order nx,ny,nz, nx,ny,nz, ... as z-increases
    // The stacks are ordered in rows along the x-axis.
    //
    // Possible matrix solver for Q-tensor
    enum QMatrixSolvers { Auto  = 0,
                          PCG   = 1,
                          GMRES = 2};

    // NOTE: the order should match with VALID_END_CRITERIA
    enum EndCriteria {
        Iterations = 0, Time = 1, Change = 2
    };

    //
    // Declare default values for parameters in Simu
    const static vector<string> VALID_END_CRITERIA;
    const static vector<string> VALID_SAVE_FORMATS;
    const static vector<string> VALID_Q_MATRIX_SOLVERS;
    const static string DEFAULT_LOAD_Q;
    const static string DEFAULT_SAVE_DIR;
    const static string DEFAULT_Q_MATRIX_SOLVER;
    const static vector<string> DEFAULT_SAVE_FORMATS;
    const static double DEFAULT_END_VALUE;
    const static double DEFAULT_DT;
    const static double DEFAULT_TARGET_DQ;
    const static double DEFAULT_MAX_DT;
    const static double DEFAULT_MAX_ERROR;
    // int default values can be defined here
    const static int DEFAULT_OUTPUT_ENERGY;
    const static int DEFAULT_OUTPUT_FORMAT;
    const static int DEFAULT_SAVE_ITER;
    const static int DEFAULT_NUM_ASSEMBLY_THREADS;
    const static int DEFAULT_NUM_MATRIX_SOLVER_THREADS;
    // default vectors
    const static vector<double> DEFAULT_STRETCH_VECTOR;
    const static vector<double> DEFAULT_DT_LIMITS;
    const static vector<double> DEFAULT_DT_FUNCTION;
    const static vector<idx>    DEFAULT_REGULAR_GRID_SIZE;
    // default enums
    const static Simu::EndCriteria DEFAULT_END_CRITERION;

private:
    const std::string meshName_;
    //! initial time step is set by user in configuration, the actual time step may change_ during the simulation.
    const double initialTimeStep_;
    const QMatrixSolvers QMatrixSolver_;

    // ADAPTIVE TIME-STEPPING VARIABLES
    const double maxError_;
    const double TargetdQ_; // do newton iterations until this precision
    const double dtLimits_[2];
    const double dtFunction_[4];

    const EndCriteria endCriterion_;
    const std::string loadQ_;
    const std::string saveDir_;        // directory where results are saved
    const double endValue_;

    //
    const double  stretchVector_[3];
    const size_t regularGridSize_[3];  // NUMBER OF NODES IN X,Y AND Z-DIRECTIONS
    const int	outputEnergy_;	// boolean whether or not to calculate energy
    const int	outputFormat_;   // 0/1 -> binary/text (for SaveFormat = LCview) // TODO should be part of list of save formats? Looks like not used anywhere
    const int	saveIter_;       // determines frequency of saving intermediate result files !! CAN THIS BE REMOVED FROM SIMU??
    //const size_t saveFormat_;  // bit field with different save formats TODO has to go. replace e.g. with set/list of enum
    const set<Simu::SaveFormats> saveFormat_;

    const unsigned int numAsseblyThreads_;
    const unsigned int numMatrixSolverThreads_;

    // TODO: mutable state!
    int MeshNumber;     // counts number of modifications. This number is appended to the end of mesh name
    double  EndValue_orig;   // original end values as defined in settings file. this is needed when end refinement is used
    bool MeshModified;
    bool restrictedTimeStep_; // flag to allow/disallow adapting time step size (e.g. just after potential switching)
    //===================================================
    // SOLVER AND ASSEMBLY (METHOD) SPECIFIC VARIABLES
public:
    Simu() = delete; // private, use SimuBuilder to create default valued Simu
    Simu(const std::string &meshName,  double initialTimeStep,
         QMatrixSolvers solver, double maxError, double targetDQ,
         const double dtLimits[2], const double dtFunction[4],
         EndCriteria endCriterion, const std::string &loadQ,
         const std::string &saveDir, double endValue,
         const double stretchVector[3], const size_t regularGridSize[3],
         int outputEnergy, int outputFormat, int saveIter,
         const set<Simu::SaveFormats> saveFormat,
         unsigned int numAsseblyThreads, unsigned int numMatrixSolverThreads

         ): meshName_(meshName), initialTimeStep_(initialTimeStep),QMatrixSolver_(solver),
         maxError_(maxError), TargetdQ_(targetDQ),
         dtLimits_{dtLimits[0], dtLimits[1]},
         dtFunction_{dtFunction[0], dtFunction[1], dtFunction[2], dtFunction[3]},
         endCriterion_(endCriterion), loadQ_(loadQ), saveDir_(saveDir), endValue_(endValue),
         stretchVector_{stretchVector[0], stretchVector[1], stretchVector[2]},
         regularGridSize_{regularGridSize[0], regularGridSize[1], regularGridSize[2]},
         outputEnergy_(outputEnergy), outputFormat_(outputFormat),
         saveIter_(saveIter), saveFormat_(saveFormat), numAsseblyThreads_(numAsseblyThreads),
         numMatrixSolverThreads_(numMatrixSolverThreads)
    {}


    double initialTimeStep() const { return initialTimeStep_; }
    SimulationMode simulationMode() const { return initialTimeStep_ > 0 ? TimeStepping : SteadyState; }
    double getMaxError() const { return maxError_; }
    void getdtFunction(double* f );

    // METHOD VARIABLE ACCESS
    unsigned int getAssemblyThreadCount()const {return numAsseblyThreads_;}
    unsigned int getMatrixSolverThreadCount()const {return numMatrixSolverThreads_;}
    QMatrixSolvers getQMatrixSolver()const {return QMatrixSolver_;}

    string  getLoadQ() const {return loadQ_;}
    string  getSaveDir()const {return saveDir_;}
    EndCriteria  getEndCriterion()const {return endCriterion_;}
    //! returns mesh filename, with MeshNumber appended TODO: delete this
    string  getMeshName()const;
    string  getMeshFileNameOnly(); // returns mesh filename, without directory

    const std::string &meshName() const { return meshName_; }

    //inline double getdt()const {return dt;}
    inline double getTargetdQ()const{return TargetdQ_;}
    inline double getMaxdt()const{return dtLimits_[1];}
    inline double getMindt()const{return dtLimits_[0];}
    double getEndValue()const {return endValue_;}

    double getStretchVectorX()const {return stretchVector_[0];}
    double getStretchVectorY()const {return stretchVector_[1];}
    double getStretchVectorZ()const {return stretchVector_[2];}

    int getSaveIter() const{ return saveIter_;}

    int getOutputEnergy()const{return outputEnergy_;}
    int getOutputFormat()const{return outputFormat_;}
    int  getMeshNumber()const{ return MeshNumber;}
    void setMeshNumber(const int& n) {MeshNumber = n;}
    void IncrementMeshNumber(){ MeshNumber++;}
    const set<Simu::SaveFormats> &getSaveFormat() const { return saveFormat_; }
// REGULAR GRID SIZE
    //void setRegularGridSize(const vector<unsigned int>& vec3);
    size_t getRegularGridXCount(){return regularGridSize_[0];}
    size_t getRegularGridYCount(){return regularGridSize_[1];}
    size_t getRegularGridZCount(){return regularGridSize_[2];}

    // TODO: MUTABLE STATE METHODS. Move to e.g. SimulationState class
    bool restrictedTimeStep() const { return restrictedTimeStep_; }
    void restrictedTimeStep(bool restrictedTimeStep) { restrictedTimeStep_ = restrictedTimeStep; }
    bool IsMeshModified()const{return MeshModified;}
    void setMeshModified(const bool& mm){MeshModified = mm;}
};

class SimuBuilder {
    std::string meshFileName_;
    double initialTimeStep_;
    Simu::QMatrixSolvers qMatrixSolver_;
    double targetDQ_;
    double maxError_;
    double dtLimits_[2];
    double dtFunction_[4];
    Simu::EndCriteria endCriterion_;
    std::string loadQ_;
    std::string saveDir_;
    double endValue_;
    double stretchVector_[3];
    size_t regularGridSize_[3];
    int outputEnergy_;
    int outputFormat_;
    int saveIter_;
    //size_t saveFormat_;
    set<Simu::SaveFormats> saveFormat_;
    unsigned int numAssemblyThreads_;
    unsigned int numMatrixSolverThreads_;

public:
    SimuBuilder():
    meshFileName_(""), initialTimeStep_(Simu::DEFAULT_DT),
    qMatrixSolver_(Simu::QMatrixSolvers::Auto), targetDQ_(Simu::DEFAULT_TARGET_DQ),
    maxError_(Simu::DEFAULT_MAX_ERROR), dtLimits_{Simu::DEFAULT_DT, Simu::DEFAULT_MAX_DT},
    dtFunction_{Simu::DEFAULT_DT_FUNCTION[0], Simu::DEFAULT_DT_FUNCTION[1], Simu::DEFAULT_DT_FUNCTION[2], Simu::DEFAULT_DT_FUNCTION[3]},
    endCriterion_(Simu::DEFAULT_END_CRITERION), loadQ_(""),
    saveDir_(Simu::DEFAULT_SAVE_DIR), endValue_(Simu::DEFAULT_END_VALUE),
    stretchVector_{1., 1., 1.}, regularGridSize_{0, 0, 0},
    outputEnergy_(Simu::DEFAULT_OUTPUT_ENERGY),
    outputFormat_(Simu::DEFAULT_OUTPUT_FORMAT), saveIter_(Simu::DEFAULT_SAVE_ITER),
    saveFormat_{}, numAssemblyThreads_(Simu::DEFAULT_NUM_ASSEMBLY_THREADS),
    numMatrixSolverThreads_(Simu::DEFAULT_NUM_MATRIX_SOLVER_THREADS)
    {}

    SimuBuilder &meshFileName(const std::string &name) { meshFileName_ = name; return *this; };
    SimuBuilder &initialTimeStep(double dt);
    SimuBuilder &qMatrixSolver(const std::string &solverName);
    SimuBuilder &targetDQ(double targetDQ);
    SimuBuilder &maxError(double maxError);
    SimuBuilder &dtLimits(double low, double high);
    SimuBuilder &dtFunction(double v1, double v2, double v3, double v4);
    SimuBuilder &endCriterion(const std::string &name);
    SimuBuilder &loadQ(const std::string &loadQ);
    SimuBuilder &saveDir(const std::string &saveDir);
    SimuBuilder &endValue(double endValue);
    SimuBuilder &stretchVector(double x, double y, double z);
    SimuBuilder &regularGridSize(size_t x, size_t y, size_t z);
    SimuBuilder &outputEnergy(int outputEnergy);
    SimuBuilder &outputFormat(int outputFormat);
    SimuBuilder &saveIter(int saveIter);
    SimuBuilder &saveFormat(const set<std::string> &saveFormats);
    SimuBuilder &numAssemblyThreads(unsigned int n);
    SimuBuilder &numMatrixSolverThreads(unsigned int n);

    Simu* build() const;
};

#endif
