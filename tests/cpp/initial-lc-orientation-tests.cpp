#include <catch.h>
#include <qlc3d.h>
#include "lc-representation.h"

const double MARGIN = 1e-12;

TEST_CASE("Set initial LC orientation") {

    double S0 = 0.55;

    Boxes boxes;

    SECTION("Set director to uniform 90 degrees inside Normal box") {
        SolutionVector q(2, 5); // 2 nodes
        double points[6] =  {
                0, 0, 0,    // inside box
                1, 0, 0};   // outside box

        // box with 90-degree tilt
        boxes.addBox(1, "Normal", {}, {0, 0.5}, {0, 1}, {0, 1}, {90, 0}, {0, 0});

        SetVolumeQ(&q, S0, &boxes, points);

        qlc3d::Director d1 = q.getDirector(0);
        qlc3d::Director d2 = q.getDirector(1);

        // d1 should be {0, 0, 1}
        REQUIRE(d1.nx() == Approx(0).margin(MARGIN));
        REQUIRE(d1.ny() == Approx(0).margin(MARGIN));
        REQUIRE(d1.nz() == Approx(1).margin(MARGIN));
        REQUIRE(d1.S() == Approx(S0).margin(MARGIN));

        // d2 should be {the default 1, 0, 0}
        REQUIRE(d2.nx() == Approx(1).margin(MARGIN));
        REQUIRE(d2.ny() == Approx(0).margin(MARGIN));
        REQUIRE(d2.nz() == Approx(0).margin(MARGIN));
        REQUIRE(d2.S() == Approx(S0).margin(MARGIN));
    }

    SECTION("Set director to linearly increasing 0-90 degree tilt, twist inside Normal box") {
        SolutionVector q(3, 5); // 3 nodes: bottom, mid, top
        double points[9] = {
                0, 0, 0,    // bottom
                0, 0, 0.5,  // mid
                0, 0, 1};   // top

        boxes.addBox(1, "Normal", {},
                     {0, 1},
                     {0, 1},
                     {0, 1},
                     {0, 90}, // tilt increasing from 0 to 90
                     {0, 90});  // twist increasing from 0 to 90
        SetVolumeQ(&q, S0, &boxes, points);

        qlc3d::Director bottom = q.getDirector(0);
        qlc3d::Director mid = q.getDirector(1);
        qlc3d::Director top = q.getDirector(2);

        // bottom should be {1, 0, 0}
        REQUIRE(bottom.S() == Approx(S0).margin(MARGIN));
        REQUIRE(bottom.nx() == Approx(1).margin(MARGIN));
        REQUIRE(bottom.ny() == Approx(0).margin(MARGIN));
        REQUIRE(bottom.nz() == Approx(0).margin(MARGIN));

        // mid should have 45-degree tilt,twist
        REQUIRE(mid.S() == Approx(S0).margin(MARGIN));
        REQUIRE(mid.tiltDegrees() == Approx(45).margin(MARGIN));
        REQUIRE(mid.twistDegrees() == Approx(45).margin(MARGIN));

        // top should have 90-degree tilt,twist
        REQUIRE(top.S() == Approx(S0).margin(MARGIN));
        REQUIRE(top.nz() == Approx(1).margin(MARGIN));
        REQUIRE(top.tiltDegrees() == Approx(90).margin(MARGIN));

    }
}
