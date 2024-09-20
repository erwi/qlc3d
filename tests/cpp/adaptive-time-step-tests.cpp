#include <catch.h>
#include <simu.h>
#include <simulation-state.h>
#include <eventlist.h>
#include <simulation-adaptive-time-step.h>

const double MIN_TIME_STEP = 1e-9;
const double MAX_TIME_STEP = 1e-4;
const double TARGET_CHANGE = 0.01;
const double INITIAL_TIME_STEP = 1e-6;

const SimulationAdaptiveTimeStepParameters STEADY_STATE_PARAMETERS = {true, MIN_TIME_STEP, MAX_TIME_STEP, TARGET_CHANGE, {0.5, 0.8, 1.2, 10}};
const SimulationAdaptiveTimeStepParameters TIME_STEPPING_PARAMETERS = {false, MIN_TIME_STEP, MAX_TIME_STEP, TARGET_CHANGE, {0.5, 0.8, 1.2, 10}};

TEST_CASE("Create adaptive time stepping") {
  SimulationAdaptiveTimeStepParameters parameters = {true, MIN_TIME_STEP, MAX_TIME_STEP, TARGET_CHANGE, {0.5, 0.8, 1.2, 10}};
  SimulationAdaptiveTimeStep adaptiveTimeStep(parameters);
  REQUIRE(&adaptiveTimeStep != nullptr);
}

TEST_CASE("Dont change time step in steady state") {
  SimulationAdaptiveTimeStep adaptiveTimeStep(STEADY_STATE_PARAMETERS);

  SimulationState simulationState;
  simulationState.dt(INITIAL_TIME_STEP);
  simulationState.change(0.1 * TARGET_CHANGE);

  adaptiveTimeStep.calculateTimeStep(simulationState);

  REQUIRE(simulationState.dt() == INITIAL_TIME_STEP); // unchanged
}

TEST_CASE("Dont change time step in first iteration") {
  SimulationAdaptiveTimeStep adaptiveTimeStep(TIME_STEPPING_PARAMETERS);

  SimulationState simulationState;
  simulationState.dt(INITIAL_TIME_STEP);
  simulationState.change(0.1 * 0.01);
  simulationState.currentIteration(1);

  // WHEN: adapt time step during first iteration
  adaptiveTimeStep.calculateTimeStep(simulationState);

  // THEN: time step should be unchanged
  REQUIRE(simulationState.dt() == INITIAL_TIME_STEP); // unchanged

  // WHEN: next iteration
  simulationState.currentIteration(2);
  adaptiveTimeStep.calculateTimeStep(simulationState);

  // THEN: time step should be increased
  REQUIRE(simulationState.dt() > INITIAL_TIME_STEP); // changed
}

TEST_CASE("Decrease time step when change is too great") {
  SimulationAdaptiveTimeStep adaptiveTimeStep(TIME_STEPPING_PARAMETERS);
  SimulationState simulationState;
  simulationState.currentIteration(2);
  simulationState.dt(INITIAL_TIME_STEP);
  simulationState.change(1.5 * TARGET_CHANGE); // change is too great

  // WHEN
  adaptiveTimeStep.calculateTimeStep(simulationState);

  // THEN: time step should be decreased
  REQUIRE(simulationState.dt() < INITIAL_TIME_STEP); // time step has decreased
}

TEST_CASE("Dont change time step if change matches target") {
  SimulationAdaptiveTimeStep adaptiveTimeStep(TIME_STEPPING_PARAMETERS);
  SimulationState simulationState;
  simulationState.currentIteration(2);
  simulationState.dt(INITIAL_TIME_STEP);
  simulationState.change(TARGET_CHANGE); // change matches target

  // WHEN
  adaptiveTimeStep.calculateTimeStep(simulationState);

  // THEN: time step should be unchanged
  REQUIRE(simulationState.dt() == INITIAL_TIME_STEP); // time step is unchanged
}

TEST_CASE("Reduce time step to avoid skipping event") {
  SimulationAdaptiveTimeStep adaptiveTimeStep(TIME_STEPPING_PARAMETERS);
  SimulationState simulationState;
  simulationState.currentIteration(2);
  simulationState.dt(INITIAL_TIME_STEP);
  simulationState.change(TARGET_CHANGE); // time step matches target, so should not be changed
  simulationState.currentTime().setTime(1e-3);

  // WHEN
  adaptiveTimeStep.setNextEventTime(SimulationTime(1e-3 + 0.5678 * INITIAL_TIME_STEP)); // next event is in half a time step
  adaptiveTimeStep.calculateTimeStep(simulationState);

  // THEN: time step should be reduced
  REQUIRE(simulationState.dt() == Approx(0.5678 * INITIAL_TIME_STEP).epsilon(1e-9)); // time step has decreased
}

TEST_CASE("Restring time step for one iteration if simulation state requires it") {
  SimulationAdaptiveTimeStep adaptiveTimeStep(TIME_STEPPING_PARAMETERS);
  SimulationState simulationState;
  simulationState.currentIteration(2);
  simulationState.dt(INITIAL_TIME_STEP);
  simulationState.change(TARGET_CHANGE); // time step matches target, so should not be changed
  simulationState.restrictedTimeStep(true);

  // WHEN
  adaptiveTimeStep.calculateTimeStep(simulationState);

  // THEN: time step should be reduced to minimum, but restriction should be lifted
  REQUIRE(simulationState.dt() == MIN_TIME_STEP); // time step has decreased
  REQUIRE(simulationState.restrictedTimeStep() == false); // restriction has been lifted
}