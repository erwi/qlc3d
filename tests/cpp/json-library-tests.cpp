#include <catch.h>
#include <../../extern/json/json.hpp>

TEST_CASE("nlohmann json header can be included and used") {
  nlohmann::json payload = {
      {"name", "qlc3d"},
      {"version", 1},
      {"enabled", true}
  };

  REQUIRE(payload.at("name").get<std::string>() == "qlc3d");
  REQUIRE(payload.at("version").get<int>() == 1);
  REQUIRE(payload.at("enabled").get<bool>());
}
