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
  auto alignment = Alignment();
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, alignment, {1, 1, 1});

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
  solverSettings->setV_GMRES_Toler(1e-9);
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
  auto alignment = Alignment();
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

  prepareGeometry(geom, TestUtil::RESOURCE_PSEUDO_2D_NEUMANN_GMSH_MESH, *electrodes, alignment, {1, 1, 1});

  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  // Set LC director to uniform 45 degree tilt angle
  SolutionVector q(geom.getnpLC(), 5);
  auto director = qlc3d::Director::fromDegreeAngles(45, 0, 0.5);
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

  //vtkIOFun::UnstructuredGridWriter writer;
  //writer.write("/home/eero/Desktop/pseudo2d.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);

  // ASSERT
  // Check that potential values equal the z-coordinate value everywhere
  for (unsigned int i = 0; i < geom.getnp(); i++) {
    double z = geom.getCoordinates().getPoint(i).z();
    double pot = v.getValue(i);
    REQUIRE(pot == Approx(z).margin(1e-6));
  }
}

TEST_CASE("Set uniform Electric field along z-axis") {
  // ARRANGE: minimal set-up required. Presence of electric field in electrodes is sufficient.
  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  electrodes->setElectricField({0, 0, 1});
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));

  prepareGeometry(geom, TestUtil::RESOURCE_THIN_GID_MESH, *electrodes, alignment, {1, 1, 1});
  auto lc = std::shared_ptr<LC>(LCBuilder().build());
  auto solverSettings = std::make_shared<SolverSettings>();

  SolutionVector v(geom.getnp(), 1);
  SolutionVector q(geom.getnpLC(), 5);

  PotentialSolver solver(electrodes, lc, solverSettings);

  // ACT
  solver.solvePotential(v, q, geom);

  // ASSERT
  for (unsigned int i = 0; i < v.getnDoF(); i++) {
    double pot = v.getValue(i);
    Vec3 p = geom.getCoordinates().getPoint(i);

    double expectedPot = p.z() - 0.5; // expect 0 potential at centre of mesh, with z-bounds 0 to 1
    REQUIRE(pot == Approx(expectedPot).margin(1e-12));
  }
}

TEST_CASE("Solve potential - mesh with dielectric layer and Neumann boundaries") {
  // ARRANGE:
  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  electrodes->setDielectricPermittivities({1});
  Alignment alignment;
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
  prepareGeometry(geom, TestUtil::RESOURCE_UNIT_CUBE_DIELECTRIC_NEUMAN_GMSH_MESH, *electrodes, alignment, {1, 1, 1});

  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  // Set LC director to uniform 45 degree tilt angle
  SolutionVector q(geom.getnpLC(), 5);
  auto director = qlc3d::Director::fromDegreeAngles(45, 0, 0.5);
  for (idx i = 0; i < geom.getnpLC(); i++) {
    q.setValue(i, director);
  }

  auto solverSettings = std::make_shared<SolverSettings>();

  SECTION("LC material with permittivity to match the dielectric layer") {
    // LC material that matches permittivity of dielectric layer
    auto lc = std::shared_ptr<LC>(LCBuilder()
                                          .eps_par(1)
                                          .eps_per(1)
                                          .build());

    // ACT
    PotentialSolver solver(electrodes, lc, solverSettings);
    solver.solvePotential(v, q, geom);

    // ASSERT
    // check that potential value is 0.5 * z for every point, since mesh ranges from 0 to 2 along z-axis
    for (unsigned int i = 0; i < geom.getCoordinates().size(); i++) {
      double z = geom.getCoordinates().getPoint(i).z();
      double pot = v.getValue(i);
      REQUIRE(pot == Approx(0.5 * z).margin(3e-4));
    }
  }

  SECTION("LC material with permittivity 2x that of the dielectric layer") {
    // LC material that matches permittivity of dielectric layer
    auto lc = std::shared_ptr<LC>(LCBuilder()
                                          .eps_par(2)
                                          .eps_per(2)
                                          .build());

    // ACT
    PotentialSolver solver(electrodes, lc, solverSettings);
    solver.solvePotential(v, q, geom);

    // ASSERT
    // check that potential value is 0.5 * z for every point, since mesh ranges from 0 to 2 along z-axis
    for (unsigned int i = 0; i < geom.getCoordinates().size(); i++) {
      double z = geom.getCoordinates().getPoint(i).z();
      double pot = v.getValue(i);
      if (z < 1) { // dielectric region
        // about double the gradient so ranging from 0 to 0.333
        double expected = 2 * z / 3.;
        REQUIRE(pot == Approx(expected).margin(3e-4));
      } else if (z > 1) {
        // about half the gradient, so ranging from 0.666 to 1.0;
        double expected = 2./ 3 + (z - 1) / 3;
        REQUIRE(pot == Approx(expected).margin(3e-4));
      }
    }
  }
}

TEST_CASE("Convenience debugging set-up, not a test!") {
  return;

  // set the path to an existing mesh file to calculate potential on it
  auto path = std::filesystem::path("/path/to/mesh/file.msh");
  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2, 3, 4, 5, 6}, {5, 0, 0, 0, 0, 0});
  electrodes->setDielectricPermittivities({1});

  Alignment alignment;
  prepareGeometry(geom, path, *electrodes, alignment, {1, 1, 1});

  SolutionVector v(geom.getnp(), 1);
  v.allocateFixedNodesArrays(geom);
  v.setPeriodicEquNodes(geom);
  v.setFixedNodesPot(electrodes->getCurrentPotentials(0));

  // Set LC director to uniform 45 degree tilt angle
  SolutionVector q(geom.getnpLC(), 5);
  auto director = qlc3d::Director::fromDegreeAngles(3, 45, 0.5);
  for (idx i = 0; i < geom.getnpLC(); i++) {
    q.setValue(i, director);
  }

  auto defaultSolverSettings = std::make_shared<SolverSettings>();

  // LC material that matches permittivity of dielectric layer
  auto lc = std::shared_ptr<LC>(LCBuilder()
                                        .eps_par(23.1)
                                        .eps_per(6.1)
                                        .build());

  // ACT
  PotentialSolver solver(electrodes, lc, defaultSolverSettings);
  solver.solvePotential(v, q, geom);

  // Write result to file for visualisation
  vtkIOFun::UnstructuredGridWriter writer;
  writer.write("/home/eero/Desktop/pseudo2d.vtk", geom.getnpLC(), geom.getCoordinates(), *geom.t, v, q);
}
