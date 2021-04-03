//
// Created by eero on 03/04/2021.
//
#ifndef PROJECT_QLC3D_SIMULATION_STATE_H
#define PROJECT_QLC3D_SIMULATION_STATE_H
#include <cassert>
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
};

/*!
 * A class for storing the current simulation state. This is mutable and evolves with each iteration.
 */
class SimulationState {
    Progress progress_;

public:
    void currentIteration(int iteration) { progress_.iteration_ = iteration; }
    int currentIteration() const { return progress_.iteration_; }

    double currentTime() const { return progress_.simulationTime_; }
    void currentTime(const double &time) { progress_.simulationTime_ = time; }

    double dt() const { return progress_.dt_; }
    void dt(const double& dt) {
        assert(dt >= 0);
        progress_.dt_ = dt;
    }

    double change() const { return progress_.change_; }
    void change(const double &change) { progress_.change_ = change; }
};
#endif //PROJECT_QLC3D_SIMULATION_STATE_H
