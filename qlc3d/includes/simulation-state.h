//
// Created by eero on 03/04/2021.
//
#ifndef PROJECT_QLC3D_SIMULATION_STATE_H
#define PROJECT_QLC3D_SIMULATION_STATE_H
#include <cassert>
enum class RunningState {
    INITIALISING, RUNNING, COMPLETED
};

struct Progress {
    //! the current iteration
    int iteration_;
    //! the current simulation time
    double simulationTime_;
    //! the current simulation time step
    double dt_;
    //! the currently elapsed real time
    double realTime_;
    //! the change_ in Q-tensor at current iteration
    double change_;
    //! the current running state
    RunningState runningState_;
};

struct Events {
    /*!
     * Indicates that the current time step has been restricted (set to minimum) because an event that requires it
     * has occurred. When this is false, the time step may not be adapted during this iteration.
     * TODO: rename to mayAdaptTimeStep ?
     */
    bool restrictedTimeStep_;
    /*!
     * Indicates that the current mesh has been modified and a new mesh file must be created.
     * TODO: consider writing the mesh as part of mesh refinement instead?
     */
    bool meshModified_;

    //! counter to keep track of number of mesh modifications
    unsigned int meshNumber_;

    Events(): restrictedTimeStep_{false}, meshModified_{false}, meshNumber_{0}
    { }
};

/*!
 * A class for storing the current simulation state. This is mutable and evolves with each iteration.
 */
class SimulationState {
    Progress progress_;
    Events events_;

public:
    void currentIteration(int iteration) { progress_.iteration_ = iteration; }
    [[nodiscard]] int currentIteration() const { return progress_.iteration_; }

    [[nodiscard]] double currentTime() const { return progress_.simulationTime_; }
    void currentTime(const double &time) { progress_.simulationTime_ = time; }

    [[nodiscard]] double dt() const { return progress_.dt_; }
    void dt(const double& dt) {
        assert(dt >= 0);
        progress_.dt_ = dt;
    }

    [[nodiscard]] double change() const { return progress_.change_; }
    void change(const double &change) { progress_.change_ = change; }

    [[nodiscard]] RunningState state() const { return progress_.runningState_; }
    void state(RunningState runningState) { progress_.runningState_ = runningState; }

    [[nodiscard]] bool restrictedTimeStep() const { return events_.restrictedTimeStep_; }
    void restrictedTimeStep(bool isRestricted) { events_.restrictedTimeStep_ = isRestricted; }

    [[nodiscard]] bool meshModified() const { return events_.meshModified_; }
    void meshModified(bool isModified) { events_.meshModified_ = isModified; }
    [[nodiscard]] unsigned int meshNumber() const { return events_.meshNumber_; }
    void incrementMeshNumber() { events_.meshNumber_++; }
};
#endif //PROJECT_QLC3D_SIMULATION_STATE_H
