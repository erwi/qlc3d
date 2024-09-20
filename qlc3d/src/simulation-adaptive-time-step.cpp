#include <simulation-adaptive-time-step.h>
#include <simu.h>
#include <simulation-state.h>
#include <eventlist.h>
#include <util/logging.h>
#include "util/exception.h"

SimulationAdaptiveTimeStep::SimulationAdaptiveTimeStep(SimulationAdaptiveTimeStepParameters parameters): parameters(parameters) {
  Log::info("Create adaptive time stepping with parameters: {}", parameters);
  assert(parameters.minTimeStep > 0);
  assert(parameters.maxTimeStep > parameters.minTimeStep);
}

void SimulationAdaptiveTimeStep::calculateTimeStep(SimulationState &simulationState) {
  if (parameters.isSteadyState) {
    return;
  }

  if (simulationState.currentIteration() <= 1) {
    return;
  }

  if (simulationState.change() <= 0 || simulationState.dt() <= 0) {
    auto msg = fmt::format("change={} or currentTimeStep={} is 0 or negative.", simulationState.change(), simulationState.dt());
    RUNTIME_ERROR(msg);
  }

  if (simulationState.restrictedTimeStep()) { // time step was restricted in previous iteration, most likely due to an event
    Log::info("Limiting time step to min dt = {}", parameters.minTimeStep);
    simulationState.restrictedTimeStep(false);
    simulationState.dt(parameters.minTimeStep);
  }

  double newdt = interpolateNewTimeStep(simulationState.dt(), simulationState.change());

  // Clamp time step to min and max values
  if (newdt < parameters.minTimeStep) {
    Log::info("limiting time step to min dt = {}", parameters.minTimeStep);
    newdt = parameters.minTimeStep;
  }
  else if (newdt > parameters.maxTimeStep) {
    Log::info("limiting time step to max dt = {}", parameters.maxTimeStep);
    newdt = parameters.maxTimeStep;
  }

  // make sure this time step size does not cause skipping an event
  newdt = dontMissEvent(newdt, simulationState);

  simulationState.dt(newdt);
}

double SimulationAdaptiveTimeStep::dontMissEvent(double newdt, const SimulationState &simulationState) const {
  // make sure this time step size does not cause skipping an event
  SimulationTime currentTime = simulationState.currentTime();
  SimulationTime nextTime = SimulationTime(currentTime).increment(newdt);

  if (nextTime.greaterThan(this->nextEventTime_)) { // would skip over next event, so reduce time step so next iteration coincides with event
    // todo: save the unreduced time step and revert back to it next time round if possible.
    newdt = nextEventTime_.getTime() - currentTime.getTime();
    Log::info("Current time is {}. Limiting time step to {} to accommodate timing of next event at {}", currentTime, newdt, this->nextEventTime_);

    if (newdt < parameters.minTimeStep) {
      Log::warn("dt={} is less than minTimeStep {} to accommodate timing of next event.", newdt, parameters.minTimeStep);
    }
    if (newdt <= 0) {
      Log::error("newdt={} is less than or equal to 0. parameters={}", newdt, parameters);
      RUNTIME_ERROR(fmt::format("newdt={} is less than or equal to 0.", newdt));
    }
  }
  return newdt;
}

double SimulationAdaptiveTimeStep::interpolateNewTimeStep(double prevTimeStep, double maxdq) const {
  // Ratio of target dQ to actual dQ. The aim is to keep this ratio close to 1
  double R = parameters.targetChange / fabs(maxdq);

  // calculate the next time step as dtNext = dt * S(R), where S is a function with three linear regions.
  // We tend to more aggressively reduce the time step than increase it.
  // In the middle region around R = 1, S = R. Increase/decrease time step linearly with R.
  // In the higher region R > RLmax, S > 1, but capped at 2. Increase time step gently
  // In the lower region R < RLmin, S < 1, but capped at 0. Decrease time step aggressively

  //double Rf[4];
  //simu.getdtFunction(&Rf[0]);
  double Rmin = parameters.interpolationFunction[0]; //Rf[0]; // S = 0
  double RLmin = parameters.interpolationFunction[1]; //Rf[1]; // S = R (min)
  double RLmax = parameters.interpolationFunction[2]; //Rf[2]; // S = R (max)
  double Rmax = parameters.interpolationFunction[3]; //Rf[3]; // S = 2

  const double Smax = 2;
  double S = 1;

  // Change is close to target.
  if ((R < RLmax) && (R > RLmin)) {
    S = R;
  }
  // Change is too great compared to target.
  else if (R <= RLmin) {
    double k = RLmin / (RLmin - Rmin); // gradient dS/dR, where S in range [0 -> RLmin]
    double dR = (RLmin - R);
    S = RLmin - k * dR;
  }
  // Change is too small compared to target.
  else if (R >= RLmax) {
    double k = (Smax - RLmax) / (Rmax - RLmax); // gradient dS/dR, where S in range [RLMax -> 2]
    double dr = R - RLmax;
    S = RLmax + k * dr;
    if (S > Smax) S = Smax;
  }

  return prevTimeStep * S;
}

void SimulationAdaptiveTimeStep::setNextEventTime(const SimulationTime &nextEventTime) {
  if (!nextEventTime.equals(nextEventTime_)) {
    Log::info("Next event time is set to {}", nextEventTime);
  }

  nextEventTime_.setTime(nextEventTime);
}