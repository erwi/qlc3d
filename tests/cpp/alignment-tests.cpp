#include <catch.h>
#include <alignment.h>

TEST_CASE("Strong anchoring properties") {
  auto s = Surface::ofStrongAnchoring(1, 45, 23);

  REQUIRE(s.getAnchoringType() == AnchoringType::Strong);
  REQUIRE(s.getStrength() == std::numeric_limits<double>::infinity());
  REQUIRE(std::isnan(s.getK1()) == true);
  REQUIRE(std::isnan(s.getK2()) == true);
  REQUIRE(s.getEasyTilt() == 45);
  REQUIRE(s.getEasyTwist() == 23);
  REQUIRE(s.usesSurfaceNormal() == false);
  REQUIRE(s.isStrong() == true);
  REQUIRE(s.getOverrideVolume() == true);
  REQUIRE(s.getFixLCNumber() == 1);
  REQUIRE(s.toString().find("Anchoring:Strong") != std::string::npos);
}

TEST_CASE("Strong anchoring properties with overrideVolume set to false") {
  auto s = Surface::ofStrongAnchoring(1, 45, 23, false);

  REQUIRE(s.getAnchoringType() == AnchoringType::Strong);
  REQUIRE(s.getStrength() == std::numeric_limits<double>::infinity());
  REQUIRE(std::isnan(s.getK1()) == true);
  REQUIRE(std::isnan(s.getK2()) == true);
  REQUIRE(s.getEasyTilt() == 45);
  REQUIRE(s.getEasyTwist() == 23);
  REQUIRE(s.usesSurfaceNormal() == false);
  REQUIRE(s.isStrong() == true);
  REQUIRE(s.getOverrideVolume() == false);
  REQUIRE(s.getFixLCNumber() == 1);
  REQUIRE(s.toString().find("Anchoring:Strong") != std::string::npos);
}


TEST_CASE("Weak anchoring properties") {
  auto s = Surface::ofWeakAnchoring(1, 45, 23, 0.5, 0.1, 0.2);

  REQUIRE(s.getAnchoringType() == AnchoringType::Weak);
  REQUIRE(s.getStrength() == 0.5);
  REQUIRE(s.getK1() == 0.1);
  REQUIRE(s.getK2() == 0.2);
  REQUIRE(s.getEasyTilt() == 45);
  REQUIRE(s.getEasyTwist() == 23);
  REQUIRE(s.usesSurfaceNormal() == false);
  REQUIRE(s.isStrong() == false);
  REQUIRE(s.getOverrideVolume() == true);
  REQUIRE(s.getFixLCNumber() == 1);
  REQUIRE(s.toString().find("Anchoring:Weak") != std::string::npos);
}

TEST_CASE("Homeotropic anchoring properties") {
  auto s = Surface::ofStrongHomeotropic(1);

  REQUIRE(s.getAnchoringType() == AnchoringType::Homeotropic);
  REQUIRE(s.getStrength() == std::numeric_limits<double>::infinity());
  REQUIRE(s.getK1() == 0);
  REQUIRE(s.getK2() == 1);
  REQUIRE(s.getEasyTilt() == 90);  // changes at every node depending on surface normal
  REQUIRE(std::isnan(s.getEasyTwist()));
  REQUIRE(s.usesSurfaceNormal() == true);
  REQUIRE(s.isStrong() == true);
  REQUIRE(s.getOverrideVolume() == true);
  REQUIRE(s.getFixLCNumber() == 1);
  REQUIRE(s.toString().find("Anchoring:Homeotropic") != std::string::npos);
}

TEST_CASE("Weak homeotropic anchoring properties") {
  auto s = Surface::ofWeakHomeotropic(1, 0.5);

  REQUIRE(s.getAnchoringType() == AnchoringType::WeakHomeotropic);
  REQUIRE(s.getStrength() == 0.5);
  REQUIRE(s.getK1() == 0);
  REQUIRE(s.getK2() == 1);
  REQUIRE(s.getEasyTilt() == 90);  // changes at every node depending on surface normal
  REQUIRE(std::isnan(s.getEasyTwist()));
  REQUIRE(s.usesSurfaceNormal() == true);
  REQUIRE(s.isStrong() == false);
  REQUIRE(s.getOverrideVolume() == true);
  REQUIRE(s.getFixLCNumber() == 1);
  REQUIRE(s.toString().find("Anchoring:WeakHomeotropic") != std::string::npos);
}

TEST_CASE("Degenerate anchoring properties") {
  auto s = Surface::ofPlanarDegenerate(1, 0.5);

  REQUIRE(s.getAnchoringType() == AnchoringType::Degenerate);
  REQUIRE(s.getStrength() == 0.5);
  REQUIRE(s.getK1() == 0);
  REQUIRE(s.getK2() == 1);
  REQUIRE(s.getEasyTilt() == 0);
  REQUIRE(std::isnan(s.getEasyTwist()));
  REQUIRE(s.usesSurfaceNormal() == true);
  REQUIRE(s.isStrong() == false);
  REQUIRE(s.getOverrideVolume() == true);
  REQUIRE(s.getFixLCNumber() == 1);
  REQUIRE(s.toString().find("Anchoring:Degenerate") != std::string::npos);
}

TEST_CASE("Freeze anchoring properties") {
  auto s = Surface::ofFreeze(1);

  REQUIRE(s.getAnchoringType() == AnchoringType::Freeze);
  REQUIRE(s.getStrength() == std::numeric_limits<double>::infinity());
  REQUIRE(std::isnan(s.getK1()) == true);
  REQUIRE(std::isnan(s.getK2()) == true);
  REQUIRE(std::isnan(s.getEasyTilt()) == true);
  REQUIRE(std::isnan(s.getEasyTwist()) == true);
  REQUIRE(s.usesSurfaceNormal() == false);
  REQUIRE(s.isStrong() == true);
  REQUIRE(s.getOverrideVolume() == false);
  REQUIRE(s.getFixLCNumber() == 1);
  REQUIRE(s.toString().find("Anchoring:Freeze") != std::string::npos);
}