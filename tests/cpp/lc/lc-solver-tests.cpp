#include <catch.h>
#include <memory>
#include <lc/lc-solver.h>
#include <lc.h>
#include <solver-settings.h>
#include <electrodes.h>
#include <geometry.h>
#include <test-util.h>
#include <inits.h>
#include <lc-representation.h>
#include <solutionvector.h>
#include <simulation-state.h>
#include <spamtrix_ircmatrix.hpp>
#include "qlc3d.h"
#include "util/logging.h"
#include "spamtrix_matrixmaker.hpp"
#include "spamtrix_vector.hpp"
#include "io/lcview-result-output.h"
#include <geom/coordinates.h>

//<editor-fold desc="TestUtil">
struct TestData {
  unique_ptr<Geometry> geom;
  unique_ptr<SolutionVector> q;
  unique_ptr<SolutionVector> v;
};

TestData setUp1DGeometry(Alignment &alignmentIn, const LC &lc, double easyTopTiilt, double easyBottomTilt) {
  auto *geom = new Geometry();
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {0, 0});
  prepareGeometry(*geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, alignmentIn, {1, 1, 1});

  auto *v = new SolutionVector(geom->getnp(), 1);
  v->allocateFixedNodesArrays(*geom);
  v->setPeriodicEquNodes(*geom);
  v->setFixedNodesPot(electrodes->getCurrentPotentials(0));

  auto* q = new SolutionVector(geom->getnpLC(), 5);
  // set volume orientation
  //auto topSurface = alignmentIn.getSurface(0);
  //auto bottomSurface = alignmentIn.getSurface(1);


  for (idx i = 0; i < geom->getnpLC(); i++) {
    Vec3 p = geom->getCoordinates().getPoint(i);
    double tiltDegrees = easyBottomTilt + (easyTopTiilt - easyBottomTilt) * p.z();
    auto director = qlc3d::Director::fromDegreeAngles(tiltDegrees, 1, lc.S0());
    q->setValue(i, director);
  }

  setSurfacesQ(*q, alignmentIn, lc.S0(), *geom);
  q->setFixedNodesQ(alignmentIn, geom->getTriangles());
  q->setPeriodicEquNodes(*geom);
  q->EnforceEquNodes(*geom);

  return {std::unique_ptr<Geometry>(geom),
          std::unique_ptr<SolutionVector>(q),
          std::unique_ptr<SolutionVector>(v),
  };
}

void steadyStateSolve(const LC &lc, Alignment &alignment, SolutionVector &q, SolutionVector &v, const Geometry &geom, int maxIter) {
  SimulationState simulationState;
  simulationState.dt(0);

  auto solverSettings = std::make_shared<SolverSettings>();
  solverSettings->setV_GMRES_Toler(1e-9);
  SteadyStateLCSolver solver(lc, *solverSettings, alignment);
  int iter = 0;
  for (iter = 0; iter < maxIter; iter++) {
    LCSolverResult solverResult = solver.solve(q, v, geom, simulationState);

    Log::info("iter={}, dq={}", iter, solverResult.dq);

    if (solverResult.converged && solverResult.dq < 1e-9) {
      Log::info("converged at iter={}", iter);
      return;
    }

    REQUIRE(solverResult.converged);
    REQUIRE(solverResult.solverType == LCSolverType::STEADY_STATE);
    REQUIRE(solverResult.iterations == 1);
  }
  FAIL("Did not converge in given number of iterations " + std::to_string(maxIter));
}

qlc3d::Director findDirectorAtZ(const Geometry &geom, SolutionVector &q, double zLevel) {
  // find director at minimum z value
  double minDeltaZ = std::numeric_limits<double>::max();
  auto coordinates = geom.getCoordinates();
  qlc3d::Director directorOut(1, 0, 0, 0.5);
  for (unsigned int i = 0; i < coordinates.size(); i++) {
    Vec3 p = coordinates.getPoint(i);

    double deltaZ = std::abs(p.z() - zLevel);

    if (deltaZ < minDeltaZ) {
      minDeltaZ = deltaZ;
      directorOut = q.getDirector(i);
    }
  }

  // reverse director if it points in negative x direction
  if (directorOut.vector().dot({1, 0, 0}) < 0) {
    directorOut = {directorOut.vector().x() * -1, directorOut.vector().y() * -1, directorOut.vector().z() * -1, directorOut.S()};
  }

  return directorOut;
}

//</editor-fold>

TEST_CASE("Create Solver") {
  auto lc = std::shared_ptr<LC>(LCBuilder().build());
  auto settings = std::make_shared<SolverSettings>();
  Alignment alignment;
  SteadyStateLCSolver solver(*lc, *settings, alignment);
}

TEST_CASE("[SteadyState] Relax elastic distortions with strong anchoring") {
  // ARRANGE
  // Set up LC with uniform distortion with -45 degrees tilt at bottom and +45 degrees tilt at top
  // Apply no electric field. Anchoring is trong on both top and bottom surfaces.
  auto lc = std::shared_ptr<LC>(LCBuilder().build());

  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {0, 0});

  // Set LC director to uniform vertical direction
  const double bottomTilt = -45;
  const double topTilt = 45;
  const double twistDegrees = 1;
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, topTilt, twistDegrees));
  alignment.addSurface(Surface::ofStrongAnchoring(2, bottomTilt, twistDegrees));

  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, alignment, {1, 1, 1});

  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  SolutionVector q(geom.getnpLC(), 5);
  SolutionVector qn(geom.getnpLC(), 5);

  // set volume orientation
  for (idx i = 0; i < geom.getnpLC(); i++) {
    Vec3 p = geom.getCoordinates().getPoint(i);
    double tiltDegrees = bottomTilt + (topTilt - bottomTilt) * p.z();
    auto director = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, lc->S0());
    q.setValue(i, director);
  }

  setSurfacesQ(q, alignment, lc->S0(), geom);

  q.setFixedNodesQ(alignment, geom.getTriangles());
  q.setPeriodicEquNodes(geom);
  q.EnforceEquNodes(geom);

  SimulationState simulationState;
  simulationState.dt(0);

  auto solverSettings = std::make_shared<SolverSettings>();
  solverSettings->setV_GMRES_Toler(1e-9);
  SteadyStateLCSolver solver(*lc, *solverSettings, alignment);

  // ACT
  // solve to tolerance of 1e-9
  int iter = 0;
  for (iter = 0; iter < 11; iter++) {
    LCSolverResult solverResult = solver.solve(q, v, geom, simulationState);


    if (solverResult.dq < 1e-9) {
      Log::info("converged at iter={}", iter);
      break;
    }

    REQUIRE(solverResult.converged);
    REQUIRE(solverResult.solverType == LCSolverType::STEADY_STATE);
    REQUIRE(solverResult.iterations == 1);
    REQUIRE(solverResult.elapsedTimes.solveTimeSeconds > 0);
    REQUIRE(solverResult.elapsedTimes.assemblyTimeSeconds > 0);
  }

  //writer.write("/home/eero/Desktop/after.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);
  // ASSERT
  // Converge in 6 iterations
  REQUIRE(iter <= 6);

  // LC director orientation should be unchanged for original
  for (idx i = 0; i < geom.getnpLC(); i++) {
    Vec3 p = geom.getCoordinates().getPoint(i);
    double expectedTiltDegrees = bottomTilt + (topTilt - bottomTilt) * p.z();
    auto expectedDirector = qlc3d::Director::fromDegreeAngles(expectedTiltDegrees, twistDegrees, lc->S0());
    double dot = expectedDirector.vector().dot(q.getDirector(i).vector());

    REQUIRE(std::abs(dot) == Approx(1).margin(1e-6));
  }
}

TEST_CASE("[SteadyState] Relax elastic distortions with weak anchoring") {
  // ARRANGE
  // Set up LC with uniform distortion with -45 degrees tilt at bottom and +45 degrees tilt at top
  // Apply no electric field. Anchoring is weak on both top and bottom surfaces.
  auto lc = std::shared_ptr<LC>(LCBuilder()
          .K11(1e-11)
          .K22(1e-11)
          .K33(1e-11)
          .build());

  const double easyTopTilt = 45;
  const double easyBottomTilt = -easyTopTilt;
  const double easyTwistDegrees = 1;
  const double Wexpected = 1e-4;
  Alignment alignment;
  alignment.addSurface(Surface::ofWeakAnchoring(1, easyTopTilt, easyTwistDegrees, Wexpected, 1, 1));
  alignment.addSurface(Surface::ofWeakAnchoring(2, easyBottomTilt, easyTwistDegrees, Wexpected, 1, 1));

  auto data = setUp1DGeometry(alignment, *lc, easyTopTilt, easyBottomTilt);
  SolutionVector &q = *data.q;
  SolutionVector &v = *data.v;
  Geometry &geom = *data.geom;

  steadyStateSolve(*lc, alignment, q, v, geom, 10);

  //writer.write("/home/eero/Desktop/after.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);
  // ASSERT
  auto bottomDir = findDirectorAtZ(geom, q, 0);
  auto topDir = findDirectorAtZ(geom, q, 1);
  Log::info("bottomDir={}, tilt={}, twist={}", bottomDir.vector(), bottomDir.tiltDegrees(), bottomDir.twistDegrees());
  Log::info("topDir={}, tilt={}, twist={}", topDir.vector(), topDir.tiltDegrees(), topDir.twistDegrees());

  double deltaSurfaceTilt = easyTopTilt - topDir.tiltDegrees();
  double deltaTilt = topDir.tiltDegrees() - bottomDir.tiltDegrees();
  double K = lc->K11();
  double d = 1 * 1e-6; // cell thickness in metres
  // convert degrees to radians
  deltaSurfaceTilt *= M_PI / 180.;
  deltaTilt *= M_PI / 180.;

  // The basic torque balance equation in Oseen-Frank vector formulation is given by
  // W = K * deltaTilt / (d * sin(2 * deltaSurfaceTilt))
  // where detaTilt is the total tilt difference between top and bottom surfaces
  // and deltaSurfaceTilt is the difference between the easy tilt angle and the actual tilt angle at the top surface.
  // This assumes that bottom and top surfaces have the same anchoring strengths but opposing easy tilt angles, so that
  // deltaSurfaceTilt is the same at both surfaces.
  //
  // In the Q-tensor formulation, the equation is modified to include a factor of (3 S0 / 2)
  // W = 2 * K * deltaTilt / (3 * S0 * d * sin(2 * deltaSurfaceTilt))

  double Weffective = (2. * K * deltaTilt) / (3 * lc->S0() * d * sin(2 * deltaSurfaceTilt));
  double R = Weffective / Wexpected;
  Log::info("W={}, R={}", Weffective, Weffective / Wexpected);

  REQUIRE(R == Approx(1).margin(1e-3));
}

TEST_CASE("[SteadyState] Relax elastic distortions with weak homeotropic anchoring") {
  auto lc = std::unique_ptr<LC>(LCBuilder()
                                        .K11(1e-11)
                                        .K22(1e-11)
                                        .K33(1e-11)
                                        .build());

  const double easyTopTilt = 90;
  const double easyBottomTilt = 0;
  const double easyTwistDegrees = 1;
  const double Wexpected = 1e-4;
  Alignment alignment;
  alignment.addSurface(Surface::ofWeakHomeotropic(1, Wexpected));
  alignment.addSurface(Surface::ofStrongAnchoring(2, easyBottomTilt, easyTwistDegrees));

  auto data = setUp1DGeometry(alignment, *lc, easyTopTilt, easyBottomTilt);
  SolutionVector &q = *data.q;
  SolutionVector &v = *data.v;
  Geometry &geom = *data.geom;

  vtkIOFun::UnstructuredGridWriter writer;
  writer.write("/home/eero/Desktop/before.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);

  // ACT
  steadyStateSolve(*lc, alignment, q, v, geom, 100);
  writer.write("/home/eero/Desktop/after.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);

  auto topDir = findDirectorAtZ(geom, q, 1);
  Log::info("topDir={}, tilt={}, twist={}", topDir.vector(), topDir.tiltDegrees(), topDir.twistDegrees());
}

TEST_CASE("[SteadyState] Relax elastic distortions with planar degenerate anchoring") {
  auto lc = std::unique_ptr<LC>(LCBuilder()
                                        .K11(1e-11)
                                        .K22(1e-11)
                                        .K33(1e-11)
                                        .build());

  const double easyTopTilt = 0;
  const double easyBottomTilt = 0;
  const double easyTwistDegrees = 1;
  const double Wexpected = 1e-4;
  Alignment alignment;
  alignment.addSurface(Surface::ofPlanarDegenerate(1, Wexpected));
  alignment.addSurface(Surface::ofStrongAnchoring(2, easyBottomTilt, easyTwistDegrees));

  auto data = setUp1DGeometry(alignment, *lc, easyTopTilt, easyBottomTilt);
  SolutionVector &q = *data.q;
  SolutionVector &v = *data.v;
  Geometry &geom = *data.geom;

  //vtkIOFun::UnstructuredGridWriter writer;
  //writer.write("/home/eero/Desktop/before.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);

  // ACT
  steadyStateSolve(*lc, alignment, q, v, geom, 100);
  //writer.write("/home/eero/Desktop/after.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);

  auto topDir = findDirectorAtZ(geom, q, 1);
  Log::info("topDir={}, tilt={}, twist={}", topDir.vector(), topDir.tiltDegrees(), topDir.twistDegrees());
}

TEST_CASE("[SteadyState] Electric switching with applied potential and three elastic constants") {
  // ARRANGE
  // Solve for steady state switching with uniform e-field. The expected mid-plane tilt angle is
  // assumed to be correct, determined at a time when the "examples/steady-state-switching-1d" example
  // is giving good agreement between qlc3d and lc3k results.
  auto lc = std::shared_ptr<LC>(LCBuilder()
          .K11(6.2e-12)
          .K22(3.9e-12)
          .K33(8.2e-12)
          .A(-0.0867e5)
          .B(-2.133e6)
          .C(1.733e6)
          .eps_par(18.5)
          .eps_per(7.0)
          .build());

  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {0, 0});
  // Set LC director to uniform vertical direction
  const double expectedMidTilt = 84.470529;
  const double bottomTilt = 5;
  const double midTilt = expectedMidTilt - bottomTilt;
  const double twistDegrees = 0;

  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, bottomTilt, twistDegrees));
  alignment.addSurface(Surface::ofStrongAnchoring(2, bottomTilt, twistDegrees));

  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, alignment, {1, 1, 1});

  const double topPotential = 2.0;
  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  SolutionVector q(geom.getnpLC(), 5);
  SolutionVector qn(geom.getnpLC(), 5);

  // set volume orientation
  for (idx i = 0; i < geom.getnpLC(); i++) {
    Vec3 p = geom.getCoordinates().getPoint(i);
    double tiltDegrees = bottomTilt + midTilt * p.z() * (1 - p.z()) * 4;
    auto director = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, lc->S0());
    q.setValue(i, director);

    // set potential to a uniform e-field
    double pot = p.z() * topPotential;
    v.setValue(i, 0, pot);
  }

  setSurfacesQ(q, alignment, lc->S0(), geom);

  q.setFixedNodesQ(alignment, geom.getTriangles());
  q.setPeriodicEquNodes(geom);
  q.EnforceEquNodes(geom);

  SimulationState simulationState;
  simulationState.dt(0);

  auto solverSettings = std::make_shared<SolverSettings>();
  solverSettings->setV_GMRES_Toler(1e-9);
  SteadyStateLCSolver solver(*lc, *solverSettings, alignment);

  // ACT
  // solve to tolerance of 1e-9
  int iter = 0;
  for (iter = 0; iter < 11; iter++) {
    auto solverResult = solver.solve(q, v, geom, simulationState);
    Log::info("iter={}, dq={}", iter, solverResult.dq);

    REQUIRE(solverResult.solverType == LCSolverType::STEADY_STATE);
    REQUIRE(solverResult.converged == true);
    REQUIRE(solverResult.iterations == 1);

    if (solverResult.dq < 1e-9) {
      Log::info("converged at iter={}", iter);
      break;
    }
  }

  // ASSERT
  // Mid-plane tilt angle (= maximum tilt angle) should match the expected value
  double maxTilt = 0;
  for (int i = 0; i < geom.getnpLC(); i++) {
    auto director = q.getDirector(i);
    maxTilt = std::max(director.tiltDegrees(), maxTilt);
  }
  Log::info("max tilt: {}", maxTilt);
  REQUIRE(maxTilt == Approx(expectedMidTilt).margin(1e-6));
}

TEST_CASE("[Dynamic] Switching dynamics with applied potential and three elastic constants") {
  // ARRANGE
  // Solve for steady state switching with uniform e-field. The expected mid-plane tilt angle is
  // assumed to be correct, determined at a time when the "examples/steady-state-switching-1d" example
  // is giving good agreement between qlc3d and lc3k results.
  auto lc = std::shared_ptr<LC>(LCBuilder()
                                        .K11(6.2e-12)
                                        .K22(3.9e-12)
                                        .K33(8.2e-12)
                                        .A(-0.0867e5)
                                        .B(-2.133e6)
                                        .C(1.733e6)
                                        .eps_par(18.5)
                                        .eps_per(7.0)
                                        .build());

  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {0, 0});
  // Set LC director to uniform vertical direction
  const double expectedMidTilt = 84.470529;
  const double bottomTilt = 5;
  const double midTilt = expectedMidTilt - bottomTilt;
  const double twistDegrees = 0;

  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, bottomTilt, twistDegrees));
  alignment.addSurface(Surface::ofStrongAnchoring(2, bottomTilt, twistDegrees));

  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, alignment, {1, 1, 1});

  const double topPotential = 2;
  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  SolutionVector q(geom.getnpLC(), 5);
  SolutionVector qn(geom.getnpLC(), 5);

  // set volume orientation
  for (idx i = 0; i < geom.getnpLC(); i++) {
    Vec3 p = geom.getCoordinates().getPoint(i);
    double tiltDegrees = bottomTilt + midTilt * p.z() * (1 - p.z()) * 4;
    auto director = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, lc->S0());
    q.setValue(i, director);

    // set potential to a uniform e-field
    double pot = p.z() * topPotential;
    v.setValue(i, 0, pot);
  }

  setSurfacesQ(q, alignment, lc->S0(), geom);

  q.setFixedNodesQ(alignment, geom.getTriangles());
  q.setPeriodicEquNodes(geom);
  q.EnforceEquNodes(geom);

  SimulationState simulationState;
  simulationState.dt(1e-4);

  auto solverSettings = std::make_shared<SolverSettings>();
  solverSettings->setV_GMRES_Toler(1e-9);

  TimeSteppingLCSolver solver(*lc, *solverSettings, 1e-6, alignment, 10);

  // ACT
  // solve to tolerance of 1e-9
  auto solverResult = solver.solve(q, v, geom, simulationState);

  // ASSERT
  // Solver should have converged to required tolerance in given number of iterations. Required number of iterations
  // is observation and may change with changes to the solver.
  REQUIRE(solverResult.solverType == LCSolverType::TIME_STEPPING);
  REQUIRE(solverResult.converged == true);
  REQUIRE(solverResult.iterations <= 8);
  REQUIRE(solverResult.maxIterationsReached == false);
}

TEST_CASE("[Dynamic] Abort Newton iterations if convergence is not reached") {
  // ARRANGE
  // Solve for steady state switching with uniform e-field. The expected mid-plane tilt angle is
  // assumed to be correct, determined at a time when the "examples/steady-state-switching-1d" example
  // is giving good agreement between qlc3d and lc3k results.
  auto lc = std::shared_ptr<LC>(LCBuilder()
                                        .K11(6.2e-12)
                                        .K22(3.9e-12)
                                        .K33(8.2e-12)
                                        .A(-0.0867e5)
                                        .B(-2.133e6)
                                        .C(1.733e6)
                                        .eps_par(18.5)
                                        .eps_per(7.0)
                                        .build());

  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {0, 0});
  // Set LC director to uniform vertical direction
  const int maxNewtonIterations = 3;
  const double expectedMidTilt = 84.470529;
  const double bottomTilt = 5;
  const double midTilt = expectedMidTilt - bottomTilt;
  const double twistDegrees = 0;

  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, bottomTilt, twistDegrees));
  alignment.addSurface(Surface::ofStrongAnchoring(2, bottomTilt, twistDegrees));

  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, alignment, {1, 1, 1});

  const double topPotential = 2;
  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  SolutionVector q(geom.getnpLC(), 5);
  SolutionVector qn(geom.getnpLC(), 5);

  // set volume orientation
  for (idx i = 0; i < geom.getnpLC(); i++) {
    Vec3 p = geom.getCoordinates().getPoint(i);
    double tiltDegrees = bottomTilt + midTilt * p.z() * (1 - p.z()) * 4;
    auto director = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, lc->S0());
    q.setValue(i, director);

    // set potential to a uniform e-field
    double pot = p.z() * topPotential;
    v.setValue(i, 0, pot);
  }

  setSurfacesQ(q, alignment, lc->S0(), geom);

  q.setFixedNodesQ(alignment, geom.getTriangles());
  q.setPeriodicEquNodes(geom);
  q.EnforceEquNodes(geom);

  SimulationState simulationState;
  simulationState.dt(1e-4);

  auto solverSettings = std::make_shared<SolverSettings>();
  solverSettings->setV_GMRES_Toler(1e-9);

  TimeSteppingLCSolver solver(*lc, *solverSettings, 1e-6, alignment, maxNewtonIterations);

  // ACT
  // solve to tolerance of 1e-9
  auto solverResult = solver.solve(q, v, geom, simulationState);

  // ASSERT
  // Solver should have converged to required tolerance in given number of iterations. Required number of iterations
  // is observation and may change with changes to the solver.
  REQUIRE(solverResult.solverType == LCSolverType::TIME_STEPPING);
  REQUIRE(solverResult.converged == true);
  REQUIRE(solverResult.iterations == maxNewtonIterations);
  REQUIRE(solverResult.maxIterationsReached == true);
  REQUIRE(solverResult.elapsedTimes.solveTimeSeconds > 0);
  REQUIRE(solverResult.elapsedTimes.assemblyTimeSeconds > 0);
}

