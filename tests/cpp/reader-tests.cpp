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

TEST_CASE("Read array of numbers") {
  std::string contents = "key = [1, 2, 3]";
  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  Reader r;
  r.readSettingsFile(settingsFile.name());

  auto array = r.get<std::vector<double>>("key", {});
  REQUIRE(array == std::vector<double>{1, 2, 3});
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

TEST_CASE("Read string value with spaces") {
  std::string contents = "key = \"value with spaces\"";
  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  Reader r;
  r.readSettingsFile(settingsFile.name());

  REQUIRE(r.get<std::string>("key", "") == "value with spaces");
}

TEST_CASE("Key should not contain double quote characters") {
  std::string contents = "\"key\" = woo";
  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  try {
    Reader r;
    r.readSettingsFile(settingsFile.name());
  } catch (ReaderError &re) {
    REQUIRE(re.errorMessage == "Key contains invalid character(s)");
    return;
  }
  FAIL("Expected ReaderError to be thrown");
}

TEST_CASE("Detect invalid string value with one double-quote") {
  std::string contents = "key = \"value";
  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  try {
    Reader r;
    r.readSettingsFile(settingsFile.name());
    std::string val = r.get<std::string>("key", "");
  } catch (ReaderError &re) {
    REQUIRE(re.lineNumber == 1);
    return;
  }
  FAIL("Expected ReaderError to be thrown");
}

TEST_CASE("Check if value type is array of strings") {
  std::string contents = "key = [\"one\",    \"two\",\t   \"three\"]";
  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  Reader r;
  r.readSettingsFile(settingsFile.name());

  REQUIRE(r.isValueArrayOfStrings("key") == true);
}
