#include <catch.h>
#include <test-util.h>
#include <reader.h>

TEST_CASE("Read key value pair") {
  std::string contents = "key = value";
  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  Reader r;
  r.readSettingsFile(settingsFile.name());

  REQUIRE(r.get<std::string>("key", "not found") == "value");
}

TEST_CASE("Return default value when key not found") {
  std::string contents = "key = value";
  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  Reader r;
  r.readSettingsFile(settingsFile.name());

  REQUIRE(r.get<std::string>("notfound", "default") == "default");
}

TEST_CASE("Check if value type is array") {
  std::string contents = "key = [1, 2, 3]";
  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  Reader r;
  r.readSettingsFile(settingsFile.name());

  REQUIRE(r.isValueArray("key") == true);
}
