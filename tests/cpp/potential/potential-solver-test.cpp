#include <catch.h>
#include <potential/potential-solver.h>
#include <electrodes.h>
#include <inits.h>
#include <test-util.h>
#include <memory>
#include <geom/coordinates.h>
#include <geometry.h>
#include <lc.h>
#include <solver-settings.h>
#include <solutionvector.h>
#include <spamtrix_ircmatrix.hpp>
#include <lc-representation.h>
#include <qlc3d.h>

TEST_CASE("Create solver") {
  auto electrodes = std::make_shared<Electrodes>();
  auto lc = std::shared_ptr<LC>(LCBuilder().build());
  auto settings = std::make_shared<SolverSettings>();
  PotentialSolver solver(electrodes, lc, settings);
}

TEST_CASE("Solve potential 1D mesh - Expect v = z") {
  // ARRANGE
  // Read thin "1D" mesh with bottom Electrode1 at z=0 and Electrode2 at z=1
  // Set the fixed potential of Electrode1 to 1 and Electrode2 to 0.
  // The LC material is set to uniform director in the vertical direction so that dielectric anisotropy
  // has no effect on the potential solution and the potential should vary linearly from 0 to 1 w.r.t. the mesh z-coordinate
  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, {1, 1, 1});

  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  // Set LC director to uniform vertical direction
  SolutionVector q(geom.getnpLC(), 5);
  auto director = qlc3d::Director(0, 0, 1, 0.6);
  for (idx i = 0; i < geom.getnpLC(); i++) {
    q.setValue(i, director);
  }

  auto lc = std::shared_ptr<LC>(LCBuilder().build());
  auto solverSettings = std::make_shared<SolverSettings>();
  PotentialSolver solver(electrodes, lc, solverSettings);

  // ACT
  solver.solvePotential(v, q, geom);

  // ASSERT
  // expect that the potential solution varies linearly from 0 to 1 w.r.t. the mesh z-coordinate
  // i.e. potential(z) = z for every point in the mesh
  double maxDiff = 0;
  for (idx i = 0; i < v.getnDoF(); i++) {
    double meshZ = geom.getCoordinates().getPoint(i).z();
    double pot = v.getValue(i);
    double diff = std::abs(meshZ - pot);
    maxDiff = std::max(maxDiff, diff);
  }

  REQUIRE(maxDiff < 1e-6);
}

TEST_CASE("Solve pseudo 2D mesh with Neumann boundaries") {
  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  prepareGeometry(geom, TestUtil::RESOURCE_PSEUDO_2D_NEUMANN_GMSH_MESH, *electrodes, {1, 1, 1});

  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  // Set LC director to uniform 45 degree tilt angle
  SolutionVector q(geom.getnpLC(), 5);
  auto director = qlc3d::Director::fromDegreeAngles(45, 0, 0.5);
  const double maxZ = geom.getZmax();
  for (idx i = 0; i < geom.getnpLC(); i++) {
    q.setValue(i, director);
  }

  auto lc = std::shared_ptr<LC>(LCBuilder()
          .eps_par(1)
          .eps_per(5)
          .build());

  auto solverSettings = std::make_shared<SolverSettings>();
  PotentialSolver solver(electrodes, lc, solverSettings);

  // ACT
  solver.solvePotential(v, q, geom);

  // ASSERT
  // TODO: check that potential values on any x-y plane are same. This is not the case currently.

  //vtkIOFun::UnstructuredGridWriter writer;
  //writer.write("/home/eero/Desktop/pseudo2d.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);
}