#include <catch.h>
#include <io/orientation-assignment.h>
#include <solutionvector.h>
#include <lc-representation.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>

// Unit tests for qlc3d::NearestNeighborAssignment - the point-cloud assignment strategy that assigns each
// LC mesh node the tensor value of the nearest OrientationSample by Euclidean distance. Tested in isolation
// against synthetic OrientationSample lists with locations, independent of any file format.

namespace {
  qlc3d::OrientationSample makeSample(double s, double x, double y, double z) {
    return qlc3d::OrientationSample{
      qlc3d::TTensor::fromDirector(qlc3d::Director::fromDegreeAngles(0, 0, s)),
      Vec3(x, y, z)
    };
  }
}

TEST_CASE("NearestNeighborAssignment assigns single sample to all mesh nodes regardless of node count") {
  // GIVEN a single located sample and 5 mesh nodes at different coordinates
  std::vector<qlc3d::OrientationSample> samples = { makeSample(0.6, 1.0, 2.0, 3.0) };

  Coordinates coords(std::vector<Vec3>{
    Vec3(0, 0, 0), Vec3(10, 0, 0), Vec3(0, 10, 0), Vec3(0, 0, 10), Vec3(5, 5, 5)
  });
  SolutionVector q(5, 5);

  // WHEN assigning
  qlc3d::NearestNeighborAssignment assignment;
  assignment.assign(samples, coords, q);

  // THEN every node receives the single sample's tensor value
  auto expected = samples[0].tensor.toDirector();
  for (idx i = 0; i < q.getnDoF(); i++) {
    auto actual = q.getDirector(i);
    REQUIRE(actual.S() == Approx(expected.S()).margin(1e-6));
    REQUIRE(actual.nx() == Approx(expected.nx()).margin(1e-6));
    REQUIRE(actual.ny() == Approx(expected.ny()).margin(1e-6));
    REQUIRE(actual.nz() == Approx(expected.nz()).margin(1e-6));
  }
}

TEST_CASE("NearestNeighborAssignment picks the nearest sample per mesh node") {
  // GIVEN mesh nodes and samples arranged along the x-axis so the expected nearest sample per node is
  // unambiguous:
  //   samples at x = 0 (S=0.5) and x = 10 (S=0.9)
  //   nodes at x = 1 (nearer to sample 0), x = 4 (nearer to sample 0), x = 6 (nearer to sample 1),
  //   x = 9 (nearer to sample 1)
  std::vector<qlc3d::OrientationSample> samples = {
    makeSample(0.5, 0, 0, 0),
    makeSample(0.9, 10, 0, 0)
  };

  Coordinates coords(std::vector<Vec3>{
    Vec3(1, 0, 0), Vec3(4, 0, 0), Vec3(6, 0, 0), Vec3(9, 0, 0)
  });
  SolutionVector q(4, 5);

  qlc3d::NearestNeighborAssignment assignment;
  assignment.assign(samples, coords, q);

  REQUIRE(q.getDirector(0).S() == Approx(0.5).margin(1e-6));
  REQUIRE(q.getDirector(1).S() == Approx(0.5).margin(1e-6));
  REQUIRE(q.getDirector(2).S() == Approx(0.9).margin(1e-6));
  REQUIRE(q.getDirector(3).S() == Approx(0.9).margin(1e-6));
}

TEST_CASE("NearestNeighborAssignment breaks exact distance ties by first-in-list-order") {
  // GIVEN two samples equidistant from a mesh node
  std::vector<qlc3d::OrientationSample> samples = {
    makeSample(0.5, -1, 0, 0), // first in list
    makeSample(0.9, 1, 0, 0)   // same distance from node at origin
  };

  Coordinates coords(std::vector<Vec3>{ Vec3(0, 0, 0) });
  SolutionVector q(1, 5);

  qlc3d::NearestNeighborAssignment assignment;
  assignment.assign(samples, coords, q);

  // THEN the first sample in list order wins the tie
  REQUIRE(q.getDirector(0).S() == Approx(0.5).margin(1e-6));
}

TEST_CASE("NearestNeighborAssignment throws on empty sample list") {
  std::vector<qlc3d::OrientationSample> samples;
  Coordinates coords(std::vector<Vec3>{ Vec3(0, 0, 0) });
  SolutionVector q(1, 5);

  qlc3d::NearestNeighborAssignment assignment;
  REQUIRE_THROWS_AS(assignment.assign(samples, coords, q), std::runtime_error);
}

TEST_CASE("NearestNeighborAssignment throws if any sample lacks a location") {
  std::vector<qlc3d::OrientationSample> samples = {
    makeSample(0.5, 0, 0, 0),
    qlc3d::OrientationSample{qlc3d::TTensor::fromDirector(qlc3d::Director::fromDegreeAngles(0, 0, 0.7)), std::nullopt}
  };
  Coordinates coords(std::vector<Vec3>{ Vec3(0, 0, 0) });
  SolutionVector q(1, 5);

  qlc3d::NearestNeighborAssignment assignment;
  REQUIRE_THROWS_AS(assignment.assign(samples, coords, q), std::runtime_error);
}
