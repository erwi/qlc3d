//
// Created by eero on 03/04/2021.
//
#ifndef PROJECT_QLC3D_SIMULATION_CONTAINER_H
#define PROJECT_QLC3D_SIMULATION_CONTAINER_H

#include <simu.h> // TODO: remove all of these
#include <electrodes.h>
#include <lc.h>
#include <simu.h>
#include <electrodes.h>
#include <lc.h>
#include <solver-settings.h>
#include <box.h>
#include <alignment.h>
#include <meshrefinement.h>
#include <regulargrid.h>

#include <ctime>
#include <memory>
#include <spamtrix_ircmatrix.hpp>
#include <eventhandler.h>

#include <simulation-state.h>
#include "lc/lc-solver.h"

class Configuration;
class SimulationState;
class ResultOutput;
class PotentialSolver;
class ILCSolver;

class SimulationContainer {
    Configuration &configuration;
    ResultOutput &resultOutput;
    std::shared_ptr<PotentialSolver> potentialSolver;
    ILCSolver &lcSolver;
    std::shared_ptr<Simu> simu;
    std::shared_ptr<Electrodes> electrodes;
    std::shared_ptr<LC> lc;
    std::unique_ptr<Boxes> boxes;
    std::unique_ptr<Alignment> alignment;
    std::unique_ptr<RegularGrid> regGrid;
    std::unique_ptr<EventList> eventList;

    // state related internal variables. TODO clean them up
    std::chrono::steady_clock::time_point startInstant;
    double maxdq;
    FILE *Energy_fid;

    // geometries and mesh
    Geometry geom1;
    Geometry geom_orig;
    Geometries geometries;

    // solution vectors
    SolutionVectors solutionVectors;
    SolutionVector q;
    SolutionVector qn;
    SolutionVector v;

    // potential consistency arrays // TODO: probably never used
    std::unique_ptr<double[]> v_cons;
    std::unique_ptr<double[]> q_cons;
    std::unique_ptr<double[]> qn_cons;

    // sparse matrices
    SpaMtrix::IRCMatrix Kq;

    SimulationState simulationState_;

    // private methods
    double updateSolutions();
    void adjustTimeStepSize();

public:
    SimulationContainer(Configuration &config, ResultOutput &resultOut, std::shared_ptr<PotentialSolver> potentialSolver, ILCSolver &lcSolver);
    /*!
     * Sets up simulation state. Reads configuration, loads mesh geometry etc.
     */
    void initialise();

    /*!
     * Whether the simulation has additional work to do (= more iterations to process)
     * @return true if there are additional simulation iterations, false if there are no more simulation iterations.
     */
    bool hasIteration() const;

    /*!
     * Runs the next iteration of the current simulation
     */
    void runIteration();

    /*!
     * Performs any additional tasks after all simulation iterations have been completed.
     */
    void postSimulationTasks();

    /*!
     * Returns the current simulation state.
     * @return current simulation state object
     */
    const SimulationState& currentState() const;
};
#endif //PROJECT_QLC3D_SIMULATION_CONTAINER_H
