
#include <simulation-container.h>
#include <configuration.h>
#include <simu.h>
#include <electrodes.h>
#include <box.h>
#include <alignment.h>
#include <meshrefinement.h>
#include <regulargrid.h>
#include <resultio.h>
#include <qlc3d.h>
#include <inits.h>
#include <eventhandler.h>
#include <util/logging.h>
#include <util/exception.h>
#include <util/stringutil.h>
#include <io/result-output.h>
#include <potential/potential-solver.h>
#include <lc/lc-solver.h>

#include <geom/vec3.h>
#include <spamtrix_ircmatrix.hpp>
#include <util/stringutil.h>

namespace fs = std::filesystem;

SimulationContainer::SimulationContainer(Configuration &config, ResultOutput &resultOut, std::shared_ptr<PotentialSolver> potentialSolver, ILCSolver &lcSolver) :
        configuration(config),
        resultOutput(resultOut),
        potentialSolver(potentialSolver),
        lcSolver(lcSolver),
        electrodes(nullptr),
        boxes(new Boxes()),
        alignment(new Alignment()),
        regGrid(new RegularGrid()),
        eventList(new EventList()) {

    Energy_fid = nullptr;
}

void SimulationContainer::initialise() {
    simulationState_.state(RunningState::INITIALISING);
    simu = configuration.getSimu();
    lc = configuration.getLC();
    electrodes = configuration.getElectrodes();
    // CHANGE CURRENT DIR TO WORKING DIRECTORY
    std::error_code ec;
    fs::current_path(configuration.currentDirectory(), ec);
    if (ec.value()) {
        RUNTIME_ERROR(fmt::format("Could not set working directory to {}", configuration.currentDirectory()));
    }
    Log::info("current working directory is {}", std::filesystem::current_path());

    // Create result output directory if it does not exist
    if (!std::filesystem::exists(simu->getSaveDirAbsolutePath()) && !std::filesystem::create_directories(simu->getSaveDirAbsolutePath())) {
            RUNTIME_ERROR(fmt::format("Could not create directory {}", simu->getSaveDirAbsolutePath().string()));
    }
    Log::info("output and result files will be written into {}", simu->getSaveDirAbsolutePath().string());

    if (simu->getSaveFormat().empty()) {
      Log::warn("No save format specified, no results will be saved. Valid save formats are " + StringUtil::toString(Simu::VALID_SAVE_FORMATS) );
    } else {
      auto saveFormats = simu->getSaveFormatStrings();
      Log::info("results will be saved every {} iterations in {} formats {}", simu->getSaveIter(), saveFormats.size(), StringUtil::toString(saveFormats));
    }

    if (resultOutput.isRegularGridRequired() && (simu->getRegularGridXCount() == 0 || simu->getRegularGridYCount() == 0 || simu->getRegularGridZCount() == 0)) {
      throw std::runtime_error("Regular grid is required by at least one result file format, but it has not been defined in the settings file.");
    }

    eventList->setSaveIter(simu->getSaveIter());
    eventList->setSaveTime(simu->getSaveTime());

    createMeshRefinementEvents(*configuration.refinement(), *eventList);
    createElectrodeSwitchingEvents(*electrodes, *eventList);

    // read missing configuration from file. TODO: all parameters to be provided in Configuration, they should not be read in here
    ReadSettings(configuration.settingsFile(),
                 *boxes,
                 *alignment,
                 *eventList);

    simulationState_.change(0);
    simulationState_.currentTime(0);
    simulationState_.dt(simu->initialTimeStep());
    simulationState_.currentIteration(0);

    // Create a backup settings file in the results directory
    std::filesystem::path settingsBackup = simu->getSaveDirAbsolutePath() / "settings.qfg";
    Log::info("Creating backup of settings file in {}", settingsBackup.string());
    // explicit deletion required due to mingw bug with replace on copy https://sourceforge.net/p/mingw-w64/bugs/852/
    if (fs::exists(settingsBackup) && !fs::remove(settingsBackup)) {
        RUNTIME_ERROR(fmt::format("Unable to delete old settings file backup file {}", settingsBackup));
    }

    if (!std::filesystem::copy_file(configuration.settingsFile(),
                                    simu->getSaveDirAbsolutePath() / "settings.qfg",
                                    std::filesystem::copy_options::overwrite_existing)) {
      RUNTIME_ERROR(fmt::format("Could not back up settings file {} to {}", configuration.settingsFile(), settingsBackup));
    }

    Energy_fid = nullptr; // file for outputting free energy
    // ================================================================
    //	CREATE GEOMETRY
    //	NEED 3 GEOMETRY OBJECTS WHEN USING MESH REFINEMENT
    // ================================================================
    std::filesystem::path meshName = configuration.currentDirectory() / simu->meshName();
    // mesh file is read and geometry is loaded in this function (in inits.cpp)
    prepareGeometry(geom_orig,
                    meshName,
                    *electrodes,
                    *alignment,
                    simu->getStretchVector(),
                    simu->getRegularGridXCount(),
                    simu->getRegularGridYCount(),
                    simu->getRegularGridZCount());
    geom1.setTo(&geom_orig);            // in the beginning working geometry is original

    // SET CONVENIENCE STRUCT OF POINTERS
    //Geometries geometries;
    geometries.geom = &geom1;
    geometries.geom_orig = &geom_orig;
    // ==============================================
    //
    //	POTENTIAL SOLUTION DATA
    //
    //================================================
    Log::info("creating initial electric potential");
    v = SolutionVector((idx) geom1.getnp(), 1);
    v.allocateFixedNodesArrays(geom1);
    v.setPeriodicEquNodes(geom1); // periodic nodes

    // =============================================================
    //
    // 	SET LC INITIAL CONDITIONS
    //
    //==============================================================
    // create vector for 5 * npLC Q-tensor components
    Log::info("Creating initial Q tensor");
    q = SolutionVector(geom1.getnpLC(), 5);    //  Q-tensor for current time step
    qn = SolutionVector(geom1.getnpLC(), 5);   //  Q-tensor from previous time step

    initialiseLcSolutionVector(q, *simu, *lc, *boxes.get(), *alignment, geom1);

    qn = q;  // q-previous = q-current in first iteration

    // SET CONVENIENCE POINTERS STRUCTURE
    solutionVectors.q = &q;
    solutionVectors.qn = &qn;
    solutionVectors.v = &v;

    Log::info("Creating matrix for Q-tensor");
    Kq = createQMatrix(geom1, q, MAT_DOMAIN1);
    //********************************************************************
    //*
    //*		Save Initial configuration and potential
    //*
    //********************************************************************
    Log::info("Saving starting configuration");
    Energy_fid = createOutputEnergyFile(*simu); // done in inits

    handleInitialEvents(simulationState_,
                        *eventList,
                        *electrodes,
                        *alignment,
                        *simu,
                        geometries,
                        solutionVectors,
                        *lc,
                        Kq,
                        resultOutput,
                        *potentialSolver);

    simulationState_.currentIteration(1);
    simulationState_.currentTime(0);
    simulationState_.change(0);
    simulationState_.dt(simu->initialTimeStep());
    startInstant = std::chrono::steady_clock::now();
    maxdq = 0;
}

bool SimulationContainer::hasIteration() const {
    auto end = simu->getEndCriterion();
    double endValue = simu->getEndValue();
    if (end == Simu::Iterations) {
        return simulationState_.currentIteration() <= (int) endValue;
    } else if (end == Simu::Time) {
        return simulationState_.currentTime() <= endValue;
    } else if (end == Simu::Change) {
        // Change is unknown until at leas one iteration has run, so at lest run one iteration.
        if (simulationState_.currentIteration() == 1) {
            return true;
        }

        double currentChange = simulationState_.change();
        if (simu->simulationMode() == TimeStepping) {
            Log::info("|dQ| = {}, EndValue = {}", fabs(currentChange), endValue);
        }
        if (currentChange <= endValue) {
          Log::info("Change {} is smaller than end value {}, ending simulation.", currentChange, endValue);
          return false;
        }
    } else {
        assert(false);
    }
    return true;
}

void SimulationContainer::runIteration() {
    simulationState_.state(RunningState::RUNNING);
    std::chrono::duration<double> elapsedSeconds = std::chrono::steady_clock::now() - startInstant;

    Log::clearIndent();
    Log::info("Iteration {}, Time = {:e}s. Real time = {:.3}s. dt = {}s.",
           simulationState_.currentIteration(),
           simulationState_.currentTime(),
           elapsedSeconds.count(),
           simulationState_.dt());

    Log::incrementIndent();
    // mve this to event handling/result output
    if (simu->getOutputEnergy()) {
        CalculateFreeEnergy(Energy_fid,
                            simulationState_.currentIteration(),
                            simulationState_.currentTime(),
                            *lc,
                            &geom1,
                            &v,
                            &q);
    }

    // CALCULATES Q-TENSOR AND POTENTIAL
    maxdq = updateSolutions();
    Log::info("maxdq = {}", maxdq);
    /// UPDATE CURRENT TIME
    simulationState_.change(maxdq);
    simulationState_.currentTime(simulationState_.currentTime() + simulationState_.dt());

    /// CALCULATE NEW TIME STEP SIZE BASED ON maxdq
    adjustTimeStepSize();
    /// INCREMENT ITERATION COUNTER

    ///EVENTS (ELECTRODES SWITCHING, RESULT OUTPUT ETC.)
    handleEvents(*eventList,
                 *electrodes,
                 *alignment,
                 *simu,
                 simulationState_,
                 geometries,
                 solutionVectors,
                 *lc,
                 Kq,
                 resultOutput,
                 *potentialSolver);

    // TODO: should this be done together with incrementing time?
    simulationState_.currentIteration(simulationState_.currentIteration() + 1);
}

void SimulationContainer::postSimulationTasks() {
    simulationState_.state(RunningState::COMPLETED);
    resultOutput.writeResults(*geometries.geom, v, q, simulationState_);
}

const SimulationState &SimulationContainer::currentState() const {
    return simulationState_;
}

double SimulationContainer::updateSolutions() {
    double maxdq = 0;
    potentialSolver->solvePotential(v, q, geom1);
    int QSolver = configuration.getSolverSettings()->getQ_Solver();

    switch (QSolver) {
        case Q_Solver_PCG: // TODO: cleanup. Exactly same calcQ3d function is called in both cases.
            //maxdq = calcQ3d(&q, &qn, &v, geom1, lc.get(), simu.get(), simulationState_, Kq, configuration.getSolverSettings().get(), alignment.get());
            maxdq = lcSolver.solve(q, v, geom1, simulationState_);
            break;
        case Q_Solver_GMRES:
            //maxdq = calcQ3d(&q, &qn, &v, geom1, lc.get(), simu.get(), simulationState_, Kq, configuration.getSolverSettings().get(), alignment.get());
          maxdq = lcSolver.solve(q, v, geom1, simulationState_);
            break;
        case Q_Solver_Explicit:
            RUNTIME_ERROR("Q_Solver_Explicit is not implemented yet.");
        default:
            RUNTIME_ERROR("Unknown Q_Solver value");
    }

    return maxdq;
}

void SimulationContainer::adjustTimeStepSize() {
    if (simu->simulationMode() == SteadyState) {
        return;
    }
    // if electrode switching etc. has just happened, dont adapt time step this time
    // but set switch to false to allow adjustment starting next step
    if (simulationState_.restrictedTimeStep()) {
        simulationState_.restrictedTimeStep(false);
        return;
    }
    double dt = simulationState_.dt();


    // INCREMENT CURRENT TIME BEFORE CHANGING TIME STEP SIZE
    // dt SHOULD BE CORRECT HERE TO MAKE SURE THAT A TIME STEP
    // COINCIDES WITH AN EVENT, IF ANY.

    /// simu.IncrementCurrentTime();

    // ADAPT TIME STEP ACCORDING TO CONVERGENCE

    double R = simu->getTargetdQ() / fabs(maxdq);

    double Rf[4];
    simu->getdtFunction(&Rf[0]);

    double Rmin = Rf[0]; // S = 0
    double RLmin = Rf[1]; // S = R (min)
    double RLmax = Rf[2]; // S = R (max)
    double Rmax = Rf[3]; // S = 2

    double Smax = 2;
    double S = 1;

    // LINEAR REGION
    if ((R < RLmax) && (R > RLmin)) {
        S = R;
    } else
        // BELOW LINEAR
    if (R <= RLmin) {
        double k = RLmin / (RLmin - Rmin);
        double dR = (RLmin - R);
        S = RLmin - k * dR;
    } else if (R >= RLmax) {// above LINEAR
        double k = (Smax - RLmax) / (Rmax - RLmax);
        double dr = R - RLmax;
        S = RLmax + k * dr;
        if (S > Smax) S = Smax;
    }

    Log::info("Scaling dt by {}.", S);
    double newdt = dt * S;
    if (newdt < simu->getMindt()) {
      Log::info("limiting time step to min dt = {}", simu->getMindt());
      newdt = simu->getMindt();
    }
    else if (newdt > simu->getMaxdt()) {
      Log::info("limiting time step to max dt = {}", simu->getMaxdt());
      newdt = simu->getMaxdt();
    }

    simulationState_.dt(newdt);
}
