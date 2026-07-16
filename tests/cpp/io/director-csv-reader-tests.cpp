#include <catch.h>
#include <io/director-csv-reader.h>
#include <io/orientation-reader.h>
#include <lc-representation.h>
#include <geom/vec3.h>
#include <test-util.h>
#include <fmt/format.h>

// Unit tests for qlc3d::DirectorCsvReader - the director CSV point-cloud format parser introduced in WP4,
// plus tests for qlc3d::createOrientationReader's extension-based dispatch.

TEST_CASE("DirectorCsvReader parses header columns in shuffled order") {
  // GIVEN a header with columns in a non-canonical order
  std::string contents = "S,nz,nx,z,ny,x,y\n";
  contents += "0.6,0,1,3,0,1,2\n"; // S=0.6, nz=0, nx=1, z=3, ny=0, x=1, y=2
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  qlc3d::DirectorCsvReader reader;
  auto samples = reader.read(file.name().string(), /*s0=*/0.9);

  REQUIRE(samples.size() == 1);
  REQUIRE(samples[0].location.has_value());
  REQUIRE(samples[0].location->x() == Approx(1.));
  REQUIRE(samples[0].location->y() == Approx(2.));
  REQUIRE(samples[0].location->z() == Approx(3.));

  auto dir = samples[0].tensor.toDirector();
  REQUIRE(dir.S() == Approx(0.6));
  REQUIRE(dir.nx() == Approx(1.).margin(1e-6));
  REQUIRE(dir.ny() == Approx(0.).margin(1e-6));
  REQUIRE(dir.nz() == Approx(0.).margin(1e-6));
}

TEST_CASE("DirectorCsvReader uses s0 parameter when S column is missing") {
  std::string contents = "x,y,z,nx,ny,nz\n";
  contents += "0,0,0,1,0,0\n";
  contents += "1,1,1,0,1,0\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  qlc3d::DirectorCsvReader reader;
  auto samples = reader.read(file.name().string(), /*s0=*/0.7);

  REQUIRE(samples.size() == 2);
  for (const auto &sample : samples) {
    REQUIRE(sample.tensor.toDirector().S() == Approx(0.7));
  }
}

TEST_CASE("DirectorCsvReader uses per-row S column values, ignoring s0") {
  std::string contents = "x,y,z,nx,ny,nz,s\n";
  contents += "0,0,0,1,0,0,0.5\n";
  contents += "1,1,1,0,1,0,0.65\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  qlc3d::DirectorCsvReader reader;
  auto samples = reader.read(file.name().string(), /*s0=*/0.1);

  REQUIRE(samples.size() == 2);
  REQUIRE(samples[0].tensor.toDirector().S() == Approx(0.5));
  REQUIRE(samples[1].tensor.toDirector().S() == Approx(0.65));
}

TEST_CASE("DirectorCsvReader normalizes un-normalized director vectors") {
  std::string contents = "x,y,z,nx,ny,nz\n";
  contents += "0,0,0,2,0,0\n"; // un-normalized (2,0,0)
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  qlc3d::DirectorCsvReader reader;
  auto samples = reader.read(file.name().string(), 0.6);

  REQUIRE(samples.size() == 1);
  REQUIRE(samples[0].tensor.toDirector().vector().norm() == Approx(1.));
}

TEST_CASE("DirectorCsvReader throws on zero-length director row") {
  std::string contents = "x,y,z,nx,ny,nz\n";
  contents += "0,0,0,0,0,0\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  qlc3d::DirectorCsvReader reader;
  REQUIRE_THROWS(reader.read(file.name().string(), 0.6));
}

TEST_CASE("DirectorCsvReader throws naming missing required column") {
  // Missing "nz"
  std::string contents = "x,y,z,nx,ny\n";
  contents += "0,0,0,1,0\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  qlc3d::DirectorCsvReader reader;
  REQUIRE_THROWS_WITH(reader.read(file.name().string(), 0.6), Catch::Contains("nz"));
}

TEST_CASE("DirectorCsvReader throws on unrecognized extra column") {
  std::string contents = "x,y,z,nx,ny,nz,extra\n";
  contents += "0,0,0,1,0,0,1\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  qlc3d::DirectorCsvReader reader;
  REQUIRE_THROWS_WITH(reader.read(file.name().string(), 0.6), Catch::Contains("extra"));
}

TEST_CASE("DirectorCsvReader throws on empty file") {
  auto file = TestUtil::TemporaryFile::withContents("", ".csv");

  qlc3d::DirectorCsvReader reader;
  REQUIRE_THROWS(reader.read(file.name().string(), 0.6));
}

TEST_CASE("DirectorCsvReader throws on header-only file (no data rows)") {
  auto file = TestUtil::TemporaryFile::withContents("x,y,z,nx,ny,nz\n", ".csv");

  qlc3d::DirectorCsvReader reader;
  REQUIRE_THROWS(reader.read(file.name().string(), 0.6));
}

TEST_CASE("DirectorCsvReader single data row returns one sample and producesLocations is true") {
  std::string contents = "x,y,z,nx,ny,nz\n0,0,0,1,0,0\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  qlc3d::DirectorCsvReader reader;
  REQUIRE(reader.producesLocations());

  auto samples = reader.read(file.name().string(), 0.6);
  REQUIRE(samples.size() == 1);
}

TEST_CASE("createOrientationReader dispatches to DirectorCsvReader for .csv extension") {
  std::string contents = "x,y,z,nx,ny,nz\n0,0,0,1,0,0\n";
  auto file = TestUtil::TemporaryFile::withContents(contents, ".csv");

  auto reader = qlc3d::createOrientationReader(file.name().string());
  REQUIRE(dynamic_cast<qlc3d::DirectorCsvReader*>(reader.get()) != nullptr);
}

TEST_CASE("createOrientationReader falls back to LCView sniffing for non-csv extension") {
  // Reuse the LCView text fixture format (>= 5 lines total for the marker-sniffing loop).
  std::string contents = "** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmesh.txt\n";
  for (int i = 0; i < 3; i++) {
    contents += fmt::format("{} 1.000000 0.000000 0.000000 0.000000 0.600000 0.600000\n", i + 1);
  }
  auto file = TestUtil::TemporaryFile::withContents(contents);

  auto reader = qlc3d::createOrientationReader(file.name().string());
  REQUIRE(dynamic_cast<qlc3d::LcViewTextReader*>(reader.get()) != nullptr);
}
