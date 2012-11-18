#ifndef SIMU_H
#define SIMU_H
#include <string>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define SIMU_N  0
#define SIMU_END_SIMULATION -20000000 // magig minus bignum

#define SIMU_OUTPUT_FORMAT_BINARY 	0
#define SIMU_OUTPUT_FORMAT_TEXT		1

enum PotentialConsistency {Off, Loop};


using namespace std;
class Simu
{

public:
    // SAVE FORMATS OPTIONS BITFIELDS
    enum SaveFormats {None          = 0,
                      LCview        = 1,
                      RegularVTK    = 2,
                      RegularVecMat = 4     // REGULAR GRID, VECTORS, MATLAB FILE
                     };
private:

    PotentialConsistency PotCons;
    double TargetPotCons; // minimum potential consistency when PotCons == Loop
    double MaxError;
    double CurrentTime;
    double CurrentChange; // change in Q, for determining convergence

    // ADAPTIVE TIME-STEPPING VARIABLES
    double  TargetdQ;
    double  dtLimits[2];
    double  dtFunction[4];
    double  Maxdt;
    string  EndCriterion;
    string  LoadQ;

    string  CurrentDir;      // working directory of qlc3d.exe
    string  SaveDir;        // directory where results are saved
    string  LoadDir;        // directory from where starting results are loaded

    double  EndValue;
    double  EndValue_orig;   // original end values as defined in settings file. this is needed when end refinement is used
    double  StretchVector[3];
    double  EnergyRegion[3]; // x,y,z coordinates for determining regions above/below which to calculate energy
    int CurrentIteration;
    bool AssembleMatrix;

    bool MeshModified;
    int MeshNumber;     // counts number of modifications. This number is appended to the end of mesh name
    int	OutputEnergy;	// boolean whether or not to calculate energy
    int	OutputFormat;   // 0/1 -> binary/text (for SaveFormat = LCview)
    int	SaveIter;       // determines frequency of saving intermediate result files !! CAN THIS BE REMOVED FROM SIMU??
    size_t SaveFormat;  // bit field with different save formats
    size_t RegularGridSize[3];  // NUMBER OF NODES IN X,Y AND Z-DIRECTIONS


public:

    // SAVE FORMAT ID STRINGS
    static const char* SF_LCVIEW;
    static const char* SF_REGULAR_VTK;
    static const char* SF_REGULAR_VECTOR_MATLAB;
    static const char* SF_LCVIEW_TXT;

    string MeshName;
    double dt;
    bool restrictedTimeStep; // flag to allow/disallow adapting time step size (e.g. just after potential switching)

    Simu();
    void PrintSimu();
    void WriteSimu(FILE* fid); // writes simu settings structure to file

    void setMeshName(string meshname);
    void setCurrentDir( const string& curdir){CurrentDir = curdir;}
    void setSaveDir(string savedir);
    void setLoadDir(string loaddir);
    void setMaxError(double me);
    void setEndCriterion( string ec);
    void setEndValue(double ev);
    void resetEndCriterion();       // starts simulation from beginning
    void setdt(double td);  // set dt, but clamps between min-max values
    void setdtForced(const double& dt); // force-sets dt, does not care about min-max values
    void setdtLimits(const double& min, const double& max);
    void setdtFunction(double* f4); // array of length 4
    void getdtFunction(double* f );
    void setCurrentTime(double ct);
    void setCurrentIteration( int i);
    void setAssembleMatrix(bool yn);
    void setCurrentChange(double ch);
    void setLoadQ(string qbackup);
    void setMaxdt(double maxdt);
    void setTargetdQ(const double& dq);


    void setOutputEnergy(int ope);
    void setOutputFormat(int opf);
    void setStretchVectorX(double sx);
    void setStretchVectorY(double sy);
    void setStretchVectorZ(double sz);
    inline void setEnergyRegionX(const double& x) { EnergyRegion[0] = x;}
    inline void setEnergyRegionY(const double& y) { EnergyRegion[1] = y;}
    inline void setEnergyRegionZ(const double& z) { EnergyRegion[2] = z;}
    inline void setSaveIter(const int& saveIter) {this->SaveIter = saveIter;}

    inline void setPotCons(PotentialConsistency pc) {PotCons = pc;}
    inline void setTargetPotCons(const double& tpc) {TargetPotCons = tpc;}

    PotentialConsistency getPotCons() {return PotCons;}
    double getTargetPotCons() {return TargetPotCons;}

    string  getLoadQ() const {return LoadQ;}
    string  getCurrentDir()const {return CurrentDir; }
    string  getLoadDir()const {return LoadDir;}
    string  getSaveDir()const {return SaveDir;}
    string  getEndCriterion()const {return EndCriterion;}
    string  getMeshName()const;  // returns mesh filename, with MeshNumber appended
    string  getMeshFileNameOnly(); // returns mesh filename, without directory

    double getMaxError()const {return MaxError;}
    inline double getdt()const {return dt;}
    inline double getTargetdQ()const{return TargetdQ;}
    inline double getMaxdt()const{return dtLimits[1];}
    inline double getMindt()const{return dtLimits[0];}
    double getCurrentTime() const {return CurrentTime; }
    double getCurrentChange()const {return CurrentChange;}
    double getEndValue()const {return EndValue;}

    double getStretchVectorX()const {return StretchVector[0];}
    double getStretchVectorY()const {return StretchVector[1];}
    double getStretchVectorZ()const {return StretchVector[2];}

    inline double getEnergyRegionX(){return EnergyRegion[0];}
    inline double getEnergyRegionY(){return EnergyRegion[1];}
    inline double getEnergyRegionZ(){return EnergyRegion[2];}
    int getSaveIter() const{ return SaveIter;}

    int getCurrentIteration() const {return CurrentIteration;}
    int getOutputEnergy()const{return OutputEnergy;}
    int getOutputFormat()const{return OutputFormat;}

    void IncrementCurrentTime();
    void IncrementCurrentIteration();
    bool IsRunning() const;
    inline bool IsAssembleMatrix()const {return AssembleMatrix;}

    inline bool IsMeshModified()const{return MeshModified;}
    inline void setMeshModified(const bool& mm){MeshModified = mm;}
    inline int  getMeshNumber()const{ return MeshNumber;}
    inline void setMeshNumber(const int& n) {MeshNumber = n;}
    inline void IncrementMeshNumber(){ MeshNumber++;}

    size_t getSaveFormat()const {return SaveFormat;}
    void clearSaveFormat(){SaveFormat = None; } // CLEARS ALL SAVE TYPES
    void addSaveFormat(std::string format); // ADDS A SAVE FORMAT

// REGULAR GRID SIZE
    void setRegularGridSize(const size_t& nx,
                            const size_t& ny,
                            const size_t& nz){RegularGridSize[0] = nx;
                                              RegularGridSize[1] = ny;
                                              RegularGridSize[2]= nz;}
    size_t getRegularGridXCount(){return RegularGridSize[0];}
    size_t getRegularGridYCount(){return RegularGridSize[1];}
    size_t getRegularGridZCount(){return RegularGridSize[2];}

};
#endif
