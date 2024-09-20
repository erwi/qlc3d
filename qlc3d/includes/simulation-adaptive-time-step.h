#ifndef PROJECT_QLC3D_SIMULATION_ADAPTIVE_TIME_STEP_H
#define PROJECT_QLC3D_SIMULATION_ADAPTIVE_TIME_STEP_H
#include <fmt/format.h>
#include <vector>
#include <simulation-time.h>

struct SimulationAdaptiveTimeStepParameters {
  /** Is the simulation in steady state? If yes, the time step is not modified as it's ignored anyway. */
  const bool isSteadyState;
  /** 
   * Minimum time step size, seconds. A value of e.g. 1e-9. In some case an even smaller time step 
   * may be produced in order to avoid missing  an event.
   */
  const double minTimeStep;
  /** Maximum time step size, seconds. A value of e.g. 1e-4*/
  const double maxTimeStep;
  
  /** Target correction in first iteration of predictor/corrector implicit time stepping. A value of e.g. 0.01 */
  const double targetChange;
  
  /** Parameter for the poorly documented piecewise linear interpolation function. Values of e.g. {0.5, 0.8, 1.2, 10} seem to work OK*/
  const std::vector<double> interpolationFunction;
};

template <>
class fmt::formatter<SimulationAdaptiveTimeStepParameters> {
public:
  constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
  template <typename Context>
  constexpr auto format (SimulationAdaptiveTimeStepParameters const& p, Context& ctx) const {
    return format_to(ctx.out(), "isSteadyState={}, minTimeStep={}, maxTimeStep={}, targetChange={}, interpolationFunction=[{},{},{},{}]", 
                     p.isSteadyState, p.minTimeStep, p.maxTimeStep, p.targetChange, p.interpolationFunction[0], p.interpolationFunction[1], p.interpolationFunction[2], p.interpolationFunction[3]);
  }
};

class SimulationState;

class SimulationAdaptiveTimeStep {
  const SimulationAdaptiveTimeStepParameters parameters;

  // default is infinity, i.e. no event is expected until they are set using setNextEventTime
  SimulationTime nextEventTime_ = SimulationTime(std::numeric_limits<double>::max());

  double interpolateNewTimeStep(double prevTimeStep, double maxdq) const ;

  double dontMissEvent(double newdt, const SimulationState &simulationState) const;

public:
  SimulationAdaptiveTimeStep(SimulationAdaptiveTimeStepParameters parameters);

  /** Calculate next time step size based on current time step and change in Q-tensor */
  void calculateTimeStep(SimulationState &simulationState);

  /** Set the time of next event. This is used to avoid skipping an event by reducing the time step size. */
  void setNextEventTime(const SimulationTime &nextEventTime);
};

#endif //PROJECT_QLC3D_SIMULATION_ADAPTIVE_TIME_STEP_H
