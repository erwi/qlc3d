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
#include "spamtrix_diagpreconditioner.hpp"
#include "spamtrix_iterativesolvers.hpp"
#include "lc/time-stepping-lc-solver.h"
#include <qassembly_macros.h>
#include <geom/coordinates.h>

TEST_CASE("Create Solver") {
  auto lc = std::shared_ptr<LC>(LCBuilder().build());
  auto settings = std::make_shared<SolverSettings>();
  LCSolver solver(*lc, *settings);
}

TEST_CASE("Relax elastic distortions") {
  // ARRANGE
  // Set up LC with uniform distortion with -45 degrees tilt at bottom and +45 degrees tilt at top
  // Apply no electric field. Anchoring is trong on both top and bottom surfaces.
  auto lc = std::shared_ptr<LC>(LCBuilder().build());

  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {0, 0});
  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, {1, 1, 1});

  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  // Set LC director to uniform vertical direction
  const double bottomTilt = -45;
  const double topTilt = 45;
  const double twistDegrees = 1;
  Alignment alignment;
  alignment.addSurface(1, "Strong", 1, {topTilt, twistDegrees, 0}, 1, 1, {});
  alignment.addSurface(2, "Strong", 1, {bottomTilt, twistDegrees, 0}, 1, 1, {});

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
  LCSolver solver(*lc, *solverSettings);

  //vtkIOFun::UnstructuredGridWriter writer;
  //writer.write("/home/eero/Desktop/before.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);

  // ACT
  // solve to tolerance of 1e-9
  int iter = 0;
  for (iter = 0; iter < 11; iter++) {
    double dq = solver.solve(q, v, geom, simulationState);


    if (dq < 1e-9) {
      Log::info("converged at iter={}", iter);
      break;
    }
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

TEST_CASE("Steady state switching with applied potential and three elastic constants") {
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
  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, {1, 1, 1});

  const double topPotential = 2.0;
  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  // Set LC director to uniform vertical direction
  const double expectedMidTilt = 84.470529;
  const double bottomTilt = 5;
  const double midTilt = expectedMidTilt - bottomTilt;
  const double twistDegrees = 0;

  Alignment alignment;
  alignment.addSurface(1, "Strong", 1, {bottomTilt, twistDegrees, 0}, 1, 1, {});
  alignment.addSurface(2, "Strong", 1, {bottomTilt, twistDegrees, 0}, 1, 1, {});

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
  LCSolver solver(*lc, *solverSettings);

  // ACT
  // solve to tolerance of 1e-9
  int iter = 0;
  for (iter = 0; iter < 11; iter++) {
    double dq = solver.solve(q, v, geom, simulationState);
    if (dq < 1e-9) {
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

TEST_CASE("Switching dynamics with applied potential and three elastic constants") {
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
  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, {1, 1, 1});

  const double topPotential = 2.0;
  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  // Set LC director to uniform vertical direction
  const double expectedMidTilt = 84.470529;
  const double bottomTilt = 5;
  const double midTilt = expectedMidTilt - bottomTilt;
  const double twistDegrees = 0;

  Alignment alignment;
  alignment.addSurface(1, "Strong", 1, {bottomTilt, twistDegrees, 0}, 1, 1, {});
  alignment.addSurface(2, "Strong", 1, {bottomTilt, twistDegrees, 0}, 1, 1, {});

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
  simulationState.dt(1e-9);

  auto solverSettings = std::make_shared<SolverSettings>();
  solverSettings->setV_GMRES_Toler(1e-9);

  TimeSteppingLCSolver solver(*lc, *solverSettings);

  // ACT
  // solve to tolerance of 1e-9
  int iter = 0;
  for (iter = 0; iter < 11; iter++) {
    double dq = solver.solve(q, v, geom, simulationState);
    if (dq < 1e-9) {
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
