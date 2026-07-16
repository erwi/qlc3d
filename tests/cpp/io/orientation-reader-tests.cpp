#include <catch.h>
#include <io/orientation-reader.h>
#include <io/orientation-assignment.h>
#include <solutionvector.h>
#include <lc-representation.h>
#include <geom/coordinates.h>
#include <test-util.h>
#include <fmt/format.h>
#include <cstdio>
#include <cstring>

// Unit tests for the orientation reader/assignment abstraction introduced to replace
// ResultIO::ReadResult/readTextLcViewResultFile/readBinaryLcViewResultFile. These readers must not touch
// SolutionVector/mesh directly - they return plain OrientationSample data - so they are tested here in
// isolation, without needing a real mesh.

namespace {
  // Writes a minimal, valid LCView text result file with the given directors/S values.
  TestUtil::TemporaryFile writeLcViewTextFixture(const std::vector<qlc3d::Director> &directors) {
    std::string contents = "** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmesh.txt\n";
    for (size_t i = 0; i < directors.size(); i++) {
      const auto &d = directors[i];
      contents += fmt::format("{} {:f} {:f} {:f} {:f} {:f} {:f}\n",
                               i + 1, d.nx(), d.ny(), d.nz(), 0., d.S(), d.S());
    }
    return TestUtil::TemporaryFile::withContents(contents);
  }

  // Writes a minimal, valid LCView binary result file containing npLC nodes worth of TTensor values.
  TestUtil::TemporaryFile writeLcViewBinaryFixture(const std::vector<qlc3d::TTensor> &tensors) {
    TestUtil::TemporaryFile file = TestUtil::TemporaryFile::empty();
    FILE *fid = fopen(file.name().string().c_str(), "wb");
    REQUIRE(fid != nullptr);
    // 5 header lines are discarded by the reader, matching the legacy writer's output (which emits an
    // extra blank line after the "Result Time" line).
    fprintf(fid, "%s\n", "** Result Time :   0.000000000");
    fprintf(fid, "\n");
    fprintf(fid, "** z Compression Ratio :  1.00000\n");
    fprintf(fid, "%s\n", "mesh.txt");
    fprintf(fid, "RAW FLOAT TRI - S0, np, nsols\n");
    int np = (int) tensors.size();
    fprintf(fid, "%g %d %d\r\n", 0.5, np, 5);
    for (auto &t : tensors) {
      float q1 = (float) t.t1(), q2 = (float) t.t2(), q3 = (float) t.t3(),
            q4 = (float) t.t4(), q5 = (float) t.t5();
      fwrite(&q1, sizeof(float), 1, fid);
      fwrite(&q2, sizeof(float), 1, fid);
      fwrite(&q3, sizeof(float), 1, fid);
      fwrite(&q5, sizeof(float), 1, fid); // note: file stores q5 before q4, matching legacy writer/reader
      fwrite(&q4, sizeof(float), 1, fid);
    }
    fclose(fid);
    return file;
  }
}

TEST_CASE("LcViewTextReader reads samples in file order without locations") {
  std::vector<qlc3d::Director> directors = {
    qlc3d::Director::fromDegreeAngles(0, 0, 0.6),
    qlc3d::Director::fromDegreeAngles(45, 90, 0.7)
  };
  auto file = writeLcViewTextFixture(directors);

  qlc3d::LcViewTextReader reader;
  REQUIRE_FALSE(reader.producesLocations());

  auto samples = reader.read(file.name().string(), /*s0=*/0.5);
  REQUIRE(samples.size() == 2);
  for (size_t i = 0; i < directors.size(); i++) {
    REQUIRE_FALSE(samples[i].location.has_value());
    auto dir = samples[i].tensor.toDirector();
    REQUIRE(dir.S() == Approx(directors[i].S()).margin(1e-6));
    REQUIRE(dir.nx() == Approx(directors[i].nx()).margin(1e-6));
    REQUIRE(dir.ny() == Approx(directors[i].ny()).margin(1e-6));
    REQUIRE(dir.nz() == Approx(directors[i].nz()).margin(1e-6));
  }
}

TEST_CASE("LcViewTextReader stops reading at all-zero director row") {
  // GIVEN a text file where npLC (2) is less than the total number of rows (3), the writer pads
  // remaining rows with all-zero director vectors to mark "end of LC region".
  std::string contents = "** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmesh.txt\n";
  contents += "1 1.000000 0.000000 0.000000 0.000000 0.600000 0.600000\n";
  contents += "2 0.000000 1.000000 0.000000 0.000000 0.600000 0.600000\n";
  contents += "3 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000\n";
  auto file = TestUtil::TemporaryFile::withContents(contents);

  qlc3d::LcViewTextReader reader;
  auto samples = reader.read(file.name().string(), 0.5);

  // THEN only the two non-zero rows are returned; the padding row is not.
  REQUIRE(samples.size() == 2);
}

TEST_CASE("LcViewBinaryReader reads samples without locations") {
  std::vector<qlc3d::TTensor> tensors = {
    qlc3d::TTensor{0.1, 0.2, 0.3, 0.4, 0.5},
    qlc3d::TTensor{-0.1, 0.15, -0.2, 0.25, -0.3}
  };
  auto file = writeLcViewBinaryFixture(tensors);

  qlc3d::LcViewBinaryReader reader;
  REQUIRE_FALSE(reader.producesLocations());

  auto samples = reader.read(file.name().string(), 0.5);
  REQUIRE(samples.size() == 2);
  for (size_t i = 0; i < tensors.size(); i++) {
    REQUIRE_FALSE(samples[i].location.has_value());
    REQUIRE(samples[i].tensor.t1() == Approx(tensors[i].t1()).margin(1e-6));
    REQUIRE(samples[i].tensor.t2() == Approx(tensors[i].t2()).margin(1e-6));
    REQUIRE(samples[i].tensor.t3() == Approx(tensors[i].t3()).margin(1e-6));
    REQUIRE(samples[i].tensor.t4() == Approx(tensors[i].t4()).margin(1e-6));
    REQUIRE(samples[i].tensor.t5() == Approx(tensors[i].t5()).margin(1e-6));
  }
}

TEST_CASE("createLcViewReader dispatches based on file content") {
  SECTION("Binary marker selects LcViewBinaryReader") {
    auto file = writeLcViewBinaryFixture({qlc3d::TTensor{0.1, 0.2, 0.3, 0.4, 0.5}});
    auto reader = qlc3d::createLcViewReader(file.name().string());
    REQUIRE(dynamic_cast<qlc3d::LcViewBinaryReader*>(reader.get()) != nullptr);
  }

  SECTION("No binary marker selects LcViewTextReader") {
    // Use several rows so the file has enough lines (>=5) for the marker-sniffing loop to safely
    // read without hitting EOF, mirroring realistic result files.
    auto file = writeLcViewTextFixture({
      qlc3d::Director::fromDegreeAngles(0, 0, 0.6),
      qlc3d::Director::fromDegreeAngles(10, 20, 0.6),
      qlc3d::Director::fromDegreeAngles(20, 40, 0.6)
    });
    auto reader = qlc3d::createLcViewReader(file.name().string());
    REQUIRE(dynamic_cast<qlc3d::LcViewTextReader*>(reader.get()) != nullptr);
  }

  SECTION("Non-existent file throws") {
    REQUIRE_THROWS_AS(qlc3d::createLcViewReader("/no/such/file.dat"), std::invalid_argument);
  }
}

TEST_CASE("ExactOrderAssignment assigns samples 1:1 onto SolutionVector in order") {
  std::vector<qlc3d::OrientationSample> samples = {
    qlc3d::OrientationSample{qlc3d::TTensor::fromDirector(qlc3d::Director::fromDegreeAngles(0, 0, 0.6)), std::nullopt},
    qlc3d::OrientationSample{qlc3d::TTensor::fromDirector(qlc3d::Director::fromDegreeAngles(30, 60, 0.7)), std::nullopt}
  };
  SolutionVector q(2, 5);
  Coordinates coords;

  qlc3d::ExactOrderAssignment assignment;
  assignment.assign(samples, coords, q);

  for (size_t i = 0; i < samples.size(); i++) {
    auto expected = samples[i].tensor.toDirector();
    auto actual = q.getDirector((idx) i);
    REQUIRE(actual.S() == Approx(expected.S()).margin(1e-6));
    REQUIRE(actual.nx() == Approx(expected.nx()).margin(1e-6));
    REQUIRE(actual.ny() == Approx(expected.ny()).margin(1e-6));
    REQUIRE(actual.nz() == Approx(expected.nz()).margin(1e-6));
  }
}

TEST_CASE("ExactOrderAssignment throws on sample count mismatch") {
  std::vector<qlc3d::OrientationSample> samples = {
    qlc3d::OrientationSample{qlc3d::TTensor::fromDirector(qlc3d::Director::fromDegreeAngles(0, 0, 0.6)), std::nullopt}
  };
  SolutionVector q(2, 5); // expects 2 samples, only 1 given
  Coordinates coords;

  qlc3d::ExactOrderAssignment assignment;
  REQUIRE_THROWS_WITH(assignment.assign(samples, coords, q),
                      Catch::Contains("The loaded result file size 1 does not match the expected size 2"));
}
