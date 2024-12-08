
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
#include <simulation-adaptive-time-step.h>
#include <spamtrix_ircmatrix.hpp>
#include "util/stopwatch.h"
#include "geom/periodicity.h"

namespace fs = std::filesystem;

SimulationContainer::SimulationContainer(Configuration &config,
                                         ResultOutput &resultOut,
                                         std::shared_ptr<PotentialSolver> potentialSolver,
                                         ILCSolver &lcSolver,
                                         EventList &eventList,
                                         SimulationState &simulationState,
                                         SimulationAdaptiveTimeStep &adaptiveTimeStep) :
        configuration(config),
        resultOutput(resultOut),
        potentialSolver(potentialSolver),
        lcSolver(lcSolver),
        electrodes(nullptr),
        boxes(config.getInitialVolumeOrientation()),
        alignment(*config.getAlignment()),
        regGrid(new RegularGrid()),
        eventList(eventList),
        simulationState(simulationState),
        adaptiveTimeStep(adaptiveTimeStep) {

    Energy_fid = nullptr;
}

void SimulationContainer::initialise() {
    simulationState.state(RunningState::INITIALISING);
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

    eventList.setSaveIter(simu->getSaveIter());
    eventList.setSaveTime(simu->getSaveTime());

    createMeshRefinementEvents(*configuration.refinement(), eventList);
    createElectrodeSwitchingEvents(*electrodes, eventList);

    simulationState.change(0);
    simulationState.setCurrentTime(0);
    simulationState.dt(simu->initialTimeStep());
    simulationState.currentIteration(0);

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
                    alignment,
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
    v.initialisePotentialBoundaries(electrodes->getCurrentPotentials(simulationState.currentTime().getTime()), geom1);

    // =============================================================
    //
    // 	SET LC INITIAL CONDITIONS
    //
    //==============================================================
    // create vector for 5 * npLC Q-tensor components
    Log::info("Creating initial Q tensor");
    q = SolutionVector(geom1.getnpLC(), 5);    //  Q-tensor for current time step

    initialiseLcSolutionVector(q, *simu, *lc, *boxes, alignment, geom1);

    geom1.clearPeriodicNodesMapping(); // release resources that are not needed anymore

    // SET CONVENIENCE POINTERS STRUCTURE
    solutionVectors.q = &q;
    solutionVectors.v = &v;

    //********************************************************************
    //*
    //*		Save Initial configuration and potential
    //*
    //********************************************************************
    Log::info("Saving starting configuration");
    Energy_fid = createOutputEnergyFile(*simu); // done in inits

    handleInitialEvents(simulationState,
                        eventList,
                        *electrodes,
                        alignment,
                        *simu,
                        geometries,
                        solutionVectors,
                        *lc,
                        resultOutput,
                        *potentialSolver,
                        adaptiveTimeStep);

    simulationState.currentIteration(1);
    simulationState.setCurrentTime(0);
    simulationState.change(0);
    simulationState.dt(simu->initialTimeStep());
    startInstant = std::chrono::steady_clock::now();
    maxdq = 0;
}

bool SimulationContainer::hasIteration() const {
    auto end = simu->getEndCriterion();
    double endValue = simu->getEndValue();
    if (end == Simu::Iterations) {
        return simulationState.currentIteration() <= (int) endValue;
    } else if (end == Simu::Time) {
        return simulationState.currentTime().equals(endValue) || simulationState.currentTime().lessThan(endValue);
    } else if (end == Simu::Change) {
        // Change is unknown until at leas one iteration has run, so at lest run one iteration.
        if (simulationState.currentIteration() == 1) {
            return true;
        }

        double currentChange = simulationState.change();
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
  Stopwatch stopwatch;
  stopwatch.start();
  simulationState.state(RunningState::RUNNING);
  //adjustTimeStepSize(); // calculate time step size for this iteration.
  adaptiveTimeStep.calculateTimeStep(simulationState);

  std::chrono::duration<double> elapsedSeconds = std::chrono::steady_clock::now() - startInstant;

  Log::clearIndent();
  Log::info("Iteration {}, Time = {:e}s. Real time = {:.3}s. dt = {}s.",
            simulationState.currentIteration(),
            simulationState.currentTime().getTime(),
            elapsedSeconds.count(),
            simulationState.dt());

  Log::incrementIndent();

  // mve this to event handling/result output
  if (simu->getOutputEnergy()) {
    CalculateFreeEnergy(Energy_fid,
                        simulationState.currentIteration(),
                        simulationState.currentTime().getTime(),
                        *lc,
                        &geom1,
                        &v,
                        &q);
  }

  // CALCULATES Q-TENSOR AND POTENTIAL
  auto solverResult = updateSolutions();
  simulationState.change(solverResult.dq);

  simulationState.currentTime().increment(simulationState.dt());

  // Note: currently event must be handled after incrementing time, but before incrementing iteration,
  // otherwise we can miss save event defined by iteration.
  handleEvents(eventList,
               *electrodes,
               alignment,
               *simu,
               simulationState,
               geometries,
               solutionVectors,
               *lc,
               resultOutput,
               *potentialSolver,
               adaptiveTimeStep);

  simulationState.currentIteration(simulationState.currentIteration() + 1);

  Log::info("Total iteration time={:.3}s, LC assembly time = {:.3}s, LC solver time = {:.3}s, |dQ| = {:.3}",
  stopwatch.elapsedSeconds(), solverResult.elapsedTimes.assemblyTimeSeconds, solverResult.elapsedTimes.solveTimeSeconds, solverResult.dq);
}

void SimulationContainer::postSimulationTasks() {
    simulationState.state(RunningState::COMPLETED);
    resultOutput.writeResults(*geometries.geom, v, q, simulationState);
}

const SimulationState &SimulationContainer::currentState() const {
    return simulationState;
}

LCSolverResult SimulationContainer::updateSolutions() {
    potentialSolver->solvePotential(v, q, geom1);
    return lcSolver.solve(q, v, geom1, simulationState);
}
