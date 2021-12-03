//
// Created by eero on 10/04/2021.
//
#include <catch.h>
#include <lc-representation.h>
const double MARGIN = 1e-12;

TEST_CASE("Detect non-unit-length Director creation") {
    REQUIRE_THROWS_AS(qlc3d::Director(2, 0, 0, 1), std::invalid_argument);
}

TEST_CASE("Create Directors by tilt, twist angles") {
    using namespace qlc3d;
    SECTION("tilt = 0, twist = 0") {
        auto n = Director::fromDegreeAngles(0, 0, 1);
        REQUIRE(n.nx() == 1);
        REQUIRE(n.ny() == 0);
        REQUIRE(n.nz() == 0);
        REQUIRE(n.S() == 1);
    }

    SECTION("tilt = 90 degrees, twist > 0") {
        auto n = Director::fromDegreeAngles(90, 1, 1);
        REQUIRE(n.nx() == Approx(0).margin(MARGIN));
        REQUIRE(n.ny() == Approx(0).margin(MARGIN));
        REQUIRE(n.nz() == Approx(1).margin(MARGIN));
    }

    SECTION("tilt = 45 degees, twist = -45 degrees") {
        auto n = Director::fromDegreeAngles(45, -45, 1);
        REQUIRE(n.nx() == Approx(0.5).margin(MARGIN));
        REQUIRE(n.ny() == Approx(-0.5).margin(MARGIN));
        REQUIRE(n.nz() == Approx(sqrt(2) / 2).margin(MARGIN));
    }
}

TEST_CASE("Convert between director and tensor representations") {
    using namespace qlc3d;

    SECTION("back and forth director/Q-tensor conversion for n=[1, 0, 0]") {
        Director dir{1, 0, 0, 0.5};
        QTensor q = QTensor::fromDirector(dir);
        Director dir2 = q.toDirector();

        REQUIRE(dir.nx() == Approx(dir2.nx()).margin(MARGIN));
        REQUIRE(dir.ny() == Approx(dir2.ny()).margin(MARGIN));
        REQUIRE(dir.nz() == Approx(dir2.nz()).margin(MARGIN));
        REQUIRE(dir.S() == Approx(dir2.S()).margin(MARGIN));
    }

    SECTION("back and forth director/Q-tensor conversion for tilt = 45degrees, twist = -45 degrees") {
        auto dir = Director::fromDegreeAngles(45, -45, 0.5);
        auto q = QTensor::fromDirector(dir);
        auto dir2 = q.toDirector();

        REQUIRE(dir.nx() == Approx(dir2.nx()).margin(MARGIN));
        REQUIRE(dir.ny() == Approx(dir2.ny()).margin(MARGIN));
        REQUIRE(dir.nz() == Approx(dir2.nz()).margin(MARGIN));
        REQUIRE(dir.S() == Approx(dir2.S()).margin(MARGIN));
    }

    SECTION("bac and forth director/T-tensor conversion for tilt = 45 degrees, twist = -45 degrees") {
        auto d = Director::fromDegreeAngles(45, -45, 0.5);
        auto t = TTensor::fromDirector(d);
        auto d2 = t.toDirector();

        REQUIRE(d.nx() == Approx(d2.nx()).margin(MARGIN));
        REQUIRE(d.ny() == Approx(d2.ny()).margin(MARGIN));
        REQUIRE(d.nz() == Approx(d2.nz()).margin(MARGIN));
        REQUIRE(d.S() == Approx(d2.S()).margin(MARGIN));
    }
}
