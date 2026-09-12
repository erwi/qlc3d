#include <catch.h>
#include <io/orientation-loader.h>
#include <simu.h>
#include <solutionvector.h>
#include <lc-representation.h>
#include <geom/coordinates.h>
#include <test-util.h>
#include <fmt/format.h>

// Unit tests for qlc3d::loadInitialOrientation - the orchestration function that resolves the
// loadQ/loadOrientation settings precedence and dispatches to the reader/assignment abstraction.

namespace {
  // Writes a minimal, valid LCView text result file with the given directors/S values, matching the
  // fixture format used in orientation-reader-tests.cpp.
  TestUtil::TemporaryFile writeLcViewTextFixture(const std::vector<qlc3d::Director> &directors) {
    std::string contents = "** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmesh.txt\n";
    for (size_t i = 0; i < directors.size(); i++) {
      const auto &d = directors[i];
      contents += fmt::format("{} {:f} {:f} {:f} {:f} {:f} {:f}\n",
                               i + 1, d.nx(), d.ny(), d.nz(), 0., d.S(), d.S());
    }
    return TestUtil::TemporaryFile::withContents(contents);
  }
}

TEST_CASE("loadInitialOrientation does nothing when neither loadQ nor loadOrientation is set") {
  Simu *simu = SimuBuilder().meshFileName("mesh.txt").build();
  SolutionVector q(1, 5);
  Coordinates coords;

  // ARRANGE: initialise q with a known value, distinct from anything a loader would write
  qlc3d::Director initial = qlc3d::Director::fromDegreeAngles(0, 0, 0.3);
  q.setValue(0, initial);

  // ACT
  qlc3d::loadInitialOrientation(*simu, 0.6, coords, q);

  // ASSERT: value is unchanged
  auto actual = q.getDirector(0);
  REQUIRE(actual.S() == Approx(initial.S()).margin(1e-6));
  delete simu;
}

TEST_CASE("loadInitialOrientation loads from loadOrientation using the LCView text reader") {
  // GIVEN a fixture file referenced via the new loadOrientation setting. Uses two rows so the file has
  // enough lines (>=5 total, matching the file-content sniffing loop in createLcViewReader) - see the
  // equivalent comment in orientation-reader-tests.cpp.
  std::vector<qlc3d::Director> directors = {
    qlc3d::Director::fromDegreeAngles(30, 60, 0.65),
    qlc3d::Director::fromDegreeAngles(10, 20, 0.55)
  };
  auto file = writeLcViewTextFixture(directors);

  Simu *simu = SimuBuilder().meshFileName("mesh.txt").loadOrientation(file.name().string()).build();
  SolutionVector q(2, 5);
  Coordinates coords;

  // WHEN
  qlc3d::loadInitialOrientation(*simu, 0.5, coords, q);

  // THEN the SolutionVector matches the fixture
  for (size_t i = 0; i < directors.size(); i++) {
    auto actual = q.getDirector((idx) i);
    REQUIRE(actual.S() == Approx(directors[i].S()).margin(1e-6));
    REQUIRE(actual.nx() == Approx(directors[i].nx()).margin(1e-6));
    REQUIRE(actual.ny() == Approx(directors[i].ny()).margin(1e-6));
    REQUIRE(actual.nz() == Approx(directors[i].nz()).margin(1e-6));
  }
  delete simu;
}

TEST_CASE("loadInitialOrientation loads from legacy loadQ (regression, deprecated)") {
  // GIVEN the same fixture, but referenced via the legacy loadQ setting
  std::vector<qlc3d::Director> directors = {
    qlc3d::Director::fromDegreeAngles(30, 60, 0.65),
    qlc3d::Director::fromDegreeAngles(10, 20, 0.55)
  };
  auto file = writeLcViewTextFixture(directors);

  Simu *simu = SimuBuilder().meshFileName("mesh.txt").loadQ(file.name().string()).build();
  SolutionVector q(2, 5);
  Coordinates coords;

  // WHEN - should not throw, and results should be identical to the loadOrientation case
  qlc3d::loadInitialOrientation(*simu, 0.5, coords, q);

  for (size_t i = 0; i < directors.size(); i++) {
    auto actual = q.getDirector((idx) i);
    REQUIRE(actual.S() == Approx(directors[i].S()).margin(1e-6));
    REQUIRE(actual.nx() == Approx(directors[i].nx()).margin(1e-6));
    REQUIRE(actual.ny() == Approx(directors[i].ny()).margin(1e-6));
    REQUIRE(actual.nz() == Approx(directors[i].nz()).margin(1e-6));
  }
  delete simu;
}

TEST_CASE("loadInitialOrientation throws when both loadQ and loadOrientation are set") {
  Simu *simu = SimuBuilder().meshFileName("mesh.txt")
      .loadQ("some-file.abc")
      .loadOrientation("some-file2.abc")
      .build();
  SolutionVector q(1, 5);
  Coordinates coords;

  REQUIRE_THROWS_WITH(qlc3d::loadInitialOrientation(*simu, 0.5, coords, q),
                      Catch::Contains("loadQ") && Catch::Contains("loadOrientation"));
  delete simu;
}

TEST_CASE("loadInitialOrientation current S0 mode preserves file orientation and replaces scalar order") {
  // GIVEN a stored LCView file whose tensor magnitude differs from the active material S0
  std::vector<qlc3d::Director> directors = {
      qlc3d::Director::fromDegreeAngles(30, 60, 0.65),
      qlc3d::Director::fromDegreeAngles(10, 20, 0.55)
  };
  auto file = writeLcViewTextFixture(directors);

  Simu *simu = SimuBuilder().meshFileName("mesh.txt")
      .loadOrientation(file.name().string())
      .loadInitialOrientationS0Mode("current")
      .build();
  SolutionVector q(2, 5);
  Coordinates coords;

  // WHEN the active material S0 is different from the file's stored scalar magnitude
  qlc3d::loadInitialOrientation(*simu, 0.35, coords, q);

  // THEN the loaded director direction is preserved while the scalar order is replaced by the current S0
  for (size_t i = 0; i < directors.size(); i++) {
    auto actual = q.getDirector((idx) i);
    REQUIRE(actual.S() == Approx(0.35).margin(1e-6));
    REQUIRE(actual.nx() == Approx(directors[i].nx()).margin(1e-6));
    REQUIRE(actual.ny() == Approx(directors[i].ny()).margin(1e-6));
    REQUIRE(actual.nz() == Approx(directors[i].nz()).margin(1e-6));
  }
  delete simu;
}

TEST_CASE("loadInitialOrientation current S0 mode overrides CSV row s while keeping director orientation") {
  // GIVEN a director CSV whose row scalar differs from the current material S0
  std::string contents = "x,y,z,nx,ny,nz,s\n0,0,0,0,0,1,0.65\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  Simu *simu = SimuBuilder().meshFileName("mesh.txt")
      .loadOrientation(file.name().string())
      .loadInitialOrientationS0Mode("current")
      .stretchVector(1, 1, 1)
      .build();
  SolutionVector q(1, 5);
  Coordinates coords(std::vector<Vec3>{Vec3(0, 0, 0)});

  // WHEN the current material S0 differs from the file value
  qlc3d::loadInitialOrientation(*simu, 0.45, coords, q);

  // THEN the CSV director direction is preserved but the scalar order is replaced by the current S0
  auto actual = q.getDirector(0);
  REQUIRE(actual.S() == Approx(0.45).margin(1e-6));
  REQUIRE(actual.nx() == Approx(0.).margin(1e-6));
  REQUIRE(actual.ny() == Approx(0.).margin(1e-6));
  REQUIRE(actual.nz() == Approx(1.).margin(1e-6));
  delete simu;
}

TEST_CASE("loadInitialOrientation: single director CSV point applies to every mesh node") {
  // GIVEN a single-row director CSV, and mesh nodes at several different locations
  std::string contents = "x,y,z,nx,ny,nz,s\n0,0,0,1,0,0,0.65\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  Simu *simu = SimuBuilder().meshFileName("mesh.txt")
      .loadOrientation(file.name().string())
      .stretchVector(1, 1, 1)
      .build();
  SolutionVector q(3, 5);
  Coordinates coords(std::vector<Vec3>{Vec3(0, 0, 0), Vec3(5, 5, 5), Vec3(-3, 2, 1)});

  // WHEN
  qlc3d::loadInitialOrientation(*simu, 0.5, coords, q);

  // THEN every node's director/S matches the single CSV row
  for (idx i = 0; i < 3; i++) {
    auto actual = q.getDirector(i);
    REQUIRE(actual.S() == Approx(0.65).margin(1e-6));
    REQUIRE(actual.nx() == Approx(1.).margin(1e-6));
    REQUIRE(actual.ny() == Approx(0.).margin(1e-6));
    REQUIRE(actual.nz() == Approx(0.).margin(1e-6));
  }
  delete simu;
}

TEST_CASE("loadInitialOrientation: multiple director CSV points, nearest wins") {
  // GIVEN two CSV rows at distinct locations, with distinct orientations
  std::string contents = "x,y,z,nx,ny,nz,s\n";
  contents += "0,0,0,1,0,0,0.5\n";   // near origin
  contents += "10,0,0,0,1,0,0.7\n";  // near x=10
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  Simu *simu = SimuBuilder().meshFileName("mesh.txt")
      .loadOrientation(file.name().string())
      .stretchVector(1, 1, 1)
      .build();
  SolutionVector q(2, 5);
  // node 0 is closer to the (0,0,0) sample; node 1 is closer to the (10,0,0) sample
  Coordinates coords(std::vector<Vec3>{Vec3(1, 0, 0), Vec3(9, 0, 0)});

  // WHEN
  qlc3d::loadInitialOrientation(*simu, 0.6, coords, q);

  // THEN
  auto node0 = q.getDirector(0);
  REQUIRE(node0.S() == Approx(0.5).margin(1e-6));
  REQUIRE(node0.nx() == Approx(1.).margin(1e-6));

  auto node1 = q.getDirector(1);
  REQUIRE(node1.S() == Approx(0.7).margin(1e-6));
  REQUIRE(node1.ny() == Approx(1.).margin(1e-6));
  delete simu;
}

TEST_CASE("loadInitialOrientation: StretchVector is applied to director CSV locations before matching") {
  // GIVEN a CSV point at unstretched x=1, and a StretchVector of {2,1,1}, so the CSV point's location in
  // stretched (mesh) space is x=2. A mesh node sits at stretched x=2 (near the scaled CSV point), and another
  // sits far away at stretched x=100. Without applying the stretch scaling to the CSV location, the CSV point
  // (unstretched x=1) would appear closer to a hypothetical node at x=1 than to the true match at x=2 -- this
  // test constructs the mesh node at the *stretched* location that only becomes nearest once the CSV location
  // is correctly scaled by StretchVector.
  std::string contents = "x,y,z,nx,ny,nz,s\n1,0,0,1,0,0,0.55\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  Simu *simu = SimuBuilder().meshFileName("mesh.txt")
      .loadOrientation(file.name().string())
      .stretchVector(2, 1, 1)
      .build();
  SolutionVector q(1, 5);
  Coordinates coords(std::vector<Vec3>{Vec3(2, 0, 0)}); // stretched location of the CSV point

  // WHEN
  qlc3d::loadInitialOrientation(*simu, 0.5, coords, q);

  // THEN the single mesh node picks up the (correctly stretch-scaled) CSV sample
  auto actual = q.getDirector(0);
  REQUIRE(actual.S() == Approx(0.55).margin(1e-6));
  REQUIRE(actual.nx() == Approx(1.).margin(1e-6));
  delete simu;
}
