//
// Created by eero on 03/04/2021.
//
#include <simulation-container.h>
#include <configuration.h>
#include <simu.h>
#include <electrodes.h>
#include <lc.h>

#include <box.h>
#include <alignment.h>
#include <meshrefinement.h>
#include <regulargrid.h>
#include <resultio.h>

#include <filesysfun.h>
#include <qlc3d.h>
#include <inits.h>
#include <calcpot3d.h> // TODO: create PotentialSolver class?
#include <eventhandler.h>
#include <util/logging.h>
#include <util/exception.h>

#include <spamtrix_ircmatrix.hpp>

SimulationContainer::SimulationContainer(Configuration &config) :
        configuration(config),
        electrodes(new Electrodes()),
        boxes(new Boxes()),
        alignment(new Alignment()),
        regGrid(new RegularGrid()),
        eventList(new EventList()),
        settings(new Settings()) {

    Energy_fid = nullptr;
}

void SimulationContainer::initialise() {
    simulationState_.state(RunningState::INITIALISING);
    // CHANGE CURRENT DIR TO WORKING DIRECTORY
    if (!FilesysFun::setCurrentDirectory(configuration.currentDirectory())) {
        RUNTIME_ERROR("Could not set working directory to " + configuration.currentDirectory());
    }
    Log::info("current working directory is {}", FilesysFun::getCurrentDirectory());

    // get parameters provided in configuration
    simu = configuration.simu();
    lc = configuration.lc();

    eventList->setSaveIter(simu->getSaveIter());
    eventList->setSaveTime(simu->getSaveTime());

    createMeshRefinementEvents(*configuration.refinement(), *eventList);

    // read missing configuration from file. TODO: all parameters to be provided in configuration
    ReadSettings(configuration.settingsFileName(),
                 *boxes,
                 *alignment,
                 *electrodes,
                 *eventList,
                 *settings);

    simulationState_.change(0);
    simulationState_.currentTime(0);
    simulationState_.dt(simu->initialTimeStep());
    simulationState_.currentIteration(0);

    // Create a backup settings file in the results directory
    // TODO: this should contain all the data, not jst a simple copy of the file. Full data, including default values,
    // will make results more repeatable in future and debugging easier.
    if (!FilesysFun::copyFile(configuration.settingsFileName(), simu->getSaveDir(), "settings.qfg")) {
        RUNTIME_ERROR("Could not back up settings file.");
    }

    Energy_fid = nullptr; // file for outputting free energy
    // ================================================================
    //	CREATE GEOMETRY
    //	NEED 3 GEOMETRY OBJECTS WHEN USING MESH REFINEMENT
    // ================================================================
    std::string meshName = configuration.currentDirectory() + "/" + simu->meshName(); // TODO
    // mesh file is read and geometry is loaded in this function (in inits.cpp)
    prepareGeometry(geom_orig,
                    meshName,
                    *simu,
                    *alignment,
                    *electrodes);
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
    v.setPeriodicEquNodes(&geom1); // periodic nodes

    // =============================================================
    //
    // 	SET LC INITIAL CONDITIONS
    //
    //==============================================================
    // create vector for 5 * npLC Q-tensor components
    Log::info("Creating initial Q tensor");
    q = SolutionVector(geom1.getnpLC(), 5);    //  Q-tensor for current time step
    qn = SolutionVector(geom1.getnpLC(), 5);   //  Q-tensor from previous time step
    SetVolumeQ(&q, lc->S0(), boxes.get(), geom1.getPtrTop());
    setSurfacesQ(&q, alignment.get(), lc->S0(), &geom1);

    //  LOAD Q FROM RESULT FILE
    if (!simu->getLoadQ().empty()) {
        ResultIO::ReadResult(*simu, q);
        setStrongSurfacesQ(&q, alignment.get(), lc->S0(), &geom1); // over writes surfaces with user specified values
    }
    q.setFixedNodesQ(alignment.get(), geom1.e);  // set fixed surface anchoring
    q.setPeriodicEquNodes(&geom1);          // periodic nodes
    q.EnforceEquNodes(geom1);                // makes sure values at periodic boundaies match
    qn = q;                                   // q-previous = q-current in first iteration

    // SET CONVENIENCE POINTERS STRUCTURE
    //solutionvectors;
    solutionVectors.q = &q;
    solutionVectors.qn = &qn;
    solutionVectors.v = &v;

    // Make matrices for potential and Q-tensor..
    Log::info("Creating sparse matrix for potential");
    Kpot = createPotentialMatrix(geom1, v, 0, *electrodes);

    Log::info("Creating matrix for Q-tensor");
    Kq = createQMatrix(geom1, q, MAT_DOMAIN1);
    //********************************************************************
    //*
    //*		Save Initial configuration and potential
    //*
    //********************************************************************
    Log::info("Saving starting configuration");
    ResultIO::CreateSaveDir(*simu);
    Energy_fid = createOutputEnergyFile(*simu); // done in inits

    handleInitialEvents(simulationState_,
                        *eventList,
                        *electrodes,
                        *alignment,
                        *simu,
                        geometries,
                        solutionVectors,
                        *lc,
                        *settings,
                        Kpot,
                        Kq
    );

    simulationState_.currentIteration(1);
    simulationState_.currentTime(0);
    simulationState_.change(0);
    simulationState_.dt(simu->initialTimeStep());
    //time_t t1, t2;
    time(&t1);
    time(&t2);
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
            return false;
        }
    } else {
        assert(false);
    }
    return true;
}

void SimulationContainer::runIteration() {
    simulationState_.state(RunningState::RUNNING);
    time(&t2);
    Log::info("Iteration {}, Time = {}s. Real time = {}s. dt = {}s.",
           simulationState_.currentIteration(),
           simulationState_.currentTime(),
           (float) t2 - t1,
           simulationState_.dt());

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

    /// UPDATE CURRENT TIME
    simulationState_.change(maxdq);
    simulationState_.currentTime(simulationState_.currentTime() + simulationState_.dt());

    /// CALCULATE NEW TIME STEP SIZE BASED ON maxdq
    adjustTimeStepSize();
    /// INCREMENT ITERATION COUNTER

    ///EVENTS (ELECTRODES SWITCHING ETC.)
    handleEvents(*eventList,
                 *electrodes,
                 *alignment,
                 *simu,
                 simulationState_,
                 geometries,
                 solutionVectors,
                 *lc,
                 *settings,
                 Kpot,
                 Kq);

    // TODO: should this be done together with incrementing time?
    simulationState_.currentIteration(simulationState_.currentIteration() + 1);
}

void SimulationContainer::postSimulationTasks() {
    simulationState_.state(RunningState::COMPLETED);

    handleResultOutput(simulationState_,
                       *simu.get(),
                       lc->S0(),
                       *geometries.geom,
                       v,
                       q);
}

const SimulationState &SimulationContainer::currentState() const {
    return simulationState_;
}

double SimulationContainer::updateSolutions() {
    double maxdq = 0;

    calcpot3d(Kpot, &v, &q, *lc.get(), geom1, settings.get(), electrodes.get());

    int QSolver = settings->getQ_Solver();
    switch (QSolver) {
        case Q_Solver_PCG: // TODO: cleanup. Exactly same calcQ3d function is called in both cases.
            maxdq = calcQ3d(&q, &qn, &v, geom1, lc.get(), simu.get(), simulationState_, Kq, settings.get(), alignment.get());
            break;
        case Q_Solver_GMRES:
            maxdq = calcQ3d(&q, &qn, &v, geom1, lc.get(), simu.get(), simulationState_, Kq, settings.get(), alignment.get());
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
        newdt = simu->getMindt();
    }

    simulationState_.dt(newdt);
}
