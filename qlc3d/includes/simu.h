#ifndef SIMU_H
#define SIMU_H
#include <string>
#include <sstream>
#include <iostream>
#include <filesystem>
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
class Vec3;
enum PotentialConsistency {Off, Loop};
enum SimulationMode {TimeStepping, SteadyState};

using namespace std;
class Simu {
public:
    // SAVE FORMATS OPTIONS BITFIELDS - MUST BE POWERS OF 2 AS WILL BE USED AS BITFIELDS!!!
    // REMEMBER TO ADD to "validSaveFormatStrings" IN CONSTRUCTOR
    // IF/WHEN ADDING NEW SAVE FORMATS!
    enum SaveFormats {
                      LCview,
                      RegularVTK,
                      RegularVecMat,    // REGULAR GRID, VECTORS, MATLAB FILE
                      DirStackZ,        // REGULAR GRID, CSV-FILE WHERE EA
                      LCviewTXT,        // LCview ASCII text file
                      CsvUnstructured,   // Unstructured csv file, compatible with paraview
                      VTKUnstructuredAsciiGrid
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

    // Declare default values for parameters in Simu
    const static vector<string> VALID_END_CRITERIA;
    const static vector<string> VALID_SAVE_FORMATS;
    const static vector<string> VALID_Q_MATRIX_SOLVERS;
    const static string DEFAULT_LOAD_Q;
    const static string DEFAULT_SAVE_DIR;
    const static string DEFAULT_Q_MATRIX_SOLVER;
    const static double DEFAULT_END_VALUE;
    const static double DEFAULT_DT;
    const static double DEFAULT_TARGET_DQ;
    const static double DEFAULT_MAX_DT;
    const static double DEFAULT_MAX_ERROR;
    // int default values can be defined here
    const static int DEFAULT_OUTPUT_ENERGY;
    const static int DEFAULT_OUTPUT_FORMAT;
    const static int DEFAULT_SAVE_ITER;
    const static double DEFAULT_SAVE_TIME;
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
    /** convergence for Newton iterations in time stepping */
    const double maxError_;
    const double TargetdQ_;             // do newton iterations until this precision
    const double dtLimits_[2];          // min/max time-step values used in adaptive time stepping
    const double dtFunction_[4];        // function coefficients that controll the time step size adjustment
    const EndCriteria endCriterion_;    // how we determine that a simulation has ended
    const std::string loadQ_;
    const std::string saveDir_;         // directory where results are saved
    const std::filesystem::path saveDirAbsolutePath_;
    const double endValue_;
    const double  stretchVector_[3];
    const size_t regularGridSize_[3];   // NUMBER OF NODES IN X,Y AND Z-DIRECTIONS
    const int	outputEnergy_;	        // boolean whether or not to calculate energy
    const int	outputFormat_;          // 0/1 -> binary/text (for SaveFormat = LCview) // TODO should be part of list of save formats? Looks like not used anywhere
    const int	saveIter_;              // determines frequency of saving intermediate result files !! CAN THIS BE REMOVED FROM SIMU??
    const double saveTime_;             // determines frequence of saving intermediate result file, frequency expressed in untis of time
    const set<Simu::SaveFormats> saveFormat_;
    const unsigned int numAsseblyThreads_;
    const unsigned int numMatrixSolverThreads_;
public:
    Simu() = delete; // private, use SimuBuilder to create default valued Simu
    Simu(const std::string &meshName,  double initialTimeStep,
         QMatrixSolvers solver, double maxError, double targetDQ,
         const double dtLimits[2], const double dtFunction[4],
         EndCriteria endCriterion, const std::string &loadQ,
         const std::string &saveDir, double endValue,
         const double stretchVector[3], const size_t regularGridSize[3],
         int outputEnergy, int outputFormat, int saveIter, double saveTime,
         const set<Simu::SaveFormats> saveFormat,
         unsigned int numAsseblyThreads, unsigned int numMatrixSolverThreads,
         const std::filesystem::path &saveDirAbsolutePath

         ): meshName_(meshName), initialTimeStep_(initialTimeStep),QMatrixSolver_(solver),
         maxError_(maxError), TargetdQ_(targetDQ),
         dtLimits_{dtLimits[0], dtLimits[1]},
         dtFunction_{dtFunction[0], dtFunction[1], dtFunction[2], dtFunction[3]},
         endCriterion_(endCriterion), loadQ_(loadQ), saveDir_(saveDir), saveDirAbsolutePath_(saveDirAbsolutePath), endValue_(endValue),
         stretchVector_{stretchVector[0], stretchVector[1], stretchVector[2]},
         regularGridSize_{regularGridSize[0], regularGridSize[1], regularGridSize[2]},
         outputEnergy_(outputEnergy), outputFormat_(outputFormat),
         saveIter_(saveIter), saveTime_(saveTime), saveFormat_(saveFormat), numAsseblyThreads_(numAsseblyThreads),
         numMatrixSolverThreads_(numMatrixSolverThreads)
    {}


    [[nodiscard]] SimulationMode simulationMode() const { return initialTimeStep_ > 0 ? TimeStepping : SteadyState; }

    void getdtFunction(double f[4] ) const;
    [[nodiscard]] std::vector<double> getdtFunction() const;
    // METHOD VARIABLE ACCESS
    [[nodiscard]] unsigned int getAssemblyThreadCount()const {return numAsseblyThreads_;}
    [[nodiscard]] unsigned int getMatrixSolverThreadCount()const {return numMatrixSolverThreads_;}

    [[nodiscard]] const std::string &getLoadQ() const {return loadQ_;}
    /**
     * @brief Returns the directory where results are saved. NOTE: this path is relative to the working directory of
     * the process, not absolute. TODO deprecated, use getSaveDirAbsolutePath() instead
     */
    [[nodiscard]] const std::string &getSaveDir()const {return saveDir_;}
    [[nodiscard]] const std::filesystem::path &getSaveDirAbsolutePath() const { return saveDirAbsolutePath_; }
    [[nodiscard]] const std::string &meshName() const { return meshName_; }

    [[nodiscard]] double initialTimeStep() const { return initialTimeStep_; }
    [[nodiscard]] double getMaxError() const { return maxError_; }
    [[nodiscard]] double getTargetdQ()const{return TargetdQ_;}
    [[nodiscard]] double getMaxdt()const{return dtLimits_[1];}
    [[nodiscard]] double getMindt()const{return dtLimits_[0];}
    [[nodiscard]] double getEndValue()const {return endValue_;}
    [[nodiscard]] Vec3 getStretchVector() const;

    [[nodiscard]] int getSaveIter() const{ return saveIter_;}
    [[nodiscard]] double getSaveTime() const { return saveTime_; }
    [[nodiscard]] int getOutputEnergy()const{return outputEnergy_;}
    [[nodiscard]] int getOutputFormat()const{return outputFormat_;}
    [[nodiscard]] const set<Simu::SaveFormats> &getSaveFormat() const { return saveFormat_; }
    [[nodiscard]] const std::vector<std::string> getSaveFormatStrings() const;
    [[nodiscard]] EndCriteria  getEndCriterion()const {return endCriterion_;}
    [[nodiscard]] QMatrixSolvers getQMatrixSolver()const {return QMatrixSolver_;}
    [[nodiscard]] size_t getRegularGridXCount()const{return regularGridSize_[0];}
    [[nodiscard]] size_t getRegularGridYCount()const{return regularGridSize_[1];}
    [[nodiscard]] size_t getRegularGridZCount()const{return regularGridSize_[2];
    }
// REGULAR GRID SIZE
    size_t getRegularGridXCount(){return regularGridSize_[0];}
    size_t getRegularGridYCount(){return regularGridSize_[1];}
    size_t getRegularGridZCount(){return regularGridSize_[2];}
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
    double saveTime_;
    //size_t saveFormat_;
    set<Simu::SaveFormats> saveFormat_;
    unsigned int numAssemblyThreads_;
    unsigned int numMatrixSolverThreads_;
    std::filesystem::path workingDir_;

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
            outputFormat_(Simu::DEFAULT_OUTPUT_FORMAT), saveIter_(Simu::DEFAULT_SAVE_ITER), saveTime_(Simu::DEFAULT_SAVE_TIME),
            saveFormat_{}, numAssemblyThreads_(Simu::DEFAULT_NUM_ASSEMBLY_THREADS),
            numMatrixSolverThreads_(Simu::DEFAULT_NUM_MATRIX_SOLVER_THREADS),
            workingDir_(std::filesystem::current_path())
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
    SimuBuilder &saveTime(double saveTime);
    SimuBuilder &saveFormat(const set<std::string> &saveFormats);
    SimuBuilder &numAssemblyThreads(unsigned int n);
    SimuBuilder &numMatrixSolverThreads(unsigned int n);

    Simu* build() const;
};

#endif
