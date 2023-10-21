#include <catch.h>
#include <fe/shapefunction3d.h>
#include <fe/fe-util.h>
#include <geom/vec3.h>
#include <geom/coordinates.h>
#include <geometry.h>
#include <electrodes.h>
#include <test-util.h>
#include <inits.h>
#include <util/logging.h>
#include "fe/shapefunction2d.h"

TEST_CASE("Linear tet 3D shape function") {
  ShapeFunction s;


  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  prepareGeometry(geom, TestUtil::RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH, *electrodes, {1, 1, 1});
  auto tets = geom.getTetrahedra();
  auto coords = geom.getCoordinates();
  idx elemNodes[4] = {0, 0, 0, 0};
  Vec3 elemCoords[4] = {Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)};

  SECTION("Check gaussian integration parameters") {
    REQUIRE(11 == s.numGaussPoints());

    double sumWeight = 0;
    for(;s.hasNextPoint(); s.nextPoint()) {
      sumWeight += s.getWeight();
    }

    REQUIRE(sumWeight == Approx(1.0 / 6.).margin(1e-12)); // sum of weights should be 1/6 because Tetrahedron volume is 6x determinant
  }

  SECTION("Integrate total volume of unit cube") {
    double totalVolume = 0;

    for (int iTet = 0; iTet < tets.getnElements(); iTet++) {
      tets.loadNodes(iTet, elemNodes);
      coords.loadCoordinates(&elemNodes[0], &elemNodes[4], elemCoords);

      double determinant = tets.getDeterminant(iTet);
      s.initialiseElement(elemCoords, determinant);

      for (; s.hasNextPoint(); s.nextPoint()) {

        double mul = s.getWeight() * determinant;
        for (int i = 0; i < 4; i++) {
          totalVolume += mul * s.N(i);
        }
      }
    }

    totalVolume *= 1e18;
    Log::info("Total volume: {}", totalVolume);

    REQUIRE(totalVolume == Approx(1.0).margin(1e-9));
  }

  SECTION("Volume integral of f(x,y,z) = x^2 * y * z in a unit cube") {
    double totalIntegral = 0;
    for (int iTet = 0; iTet < tets.getnElements(); iTet++) {
      tets.loadNodes(iTet, elemNodes);
      coords.loadCoordinates(&elemNodes[0], &elemNodes[4], elemCoords);

      double determinant = tets.getDeterminant(iTet) * 1e18;
      s.initialiseElement(elemCoords, determinant);

      double x[] = {elemCoords[0].x(), elemCoords[1].x(), elemCoords[2].x(), elemCoords[3].x()};
      double y[] = {elemCoords[0].y(), elemCoords[1].y(), elemCoords[2].y(), elemCoords[3].y()};
      double z[] = {elemCoords[0].z(), elemCoords[1].z(), elemCoords[2].z(), elemCoords[3].z()};

      for (;s.hasNextPoint(); s.nextPoint()) {

        double xlocal = s.calculate(x);
        double ylocal = s.calculate(y);
        double zlocal = s.calculate(z);

        double mul = s.getWeight() * determinant;
        totalIntegral += mul * xlocal * xlocal * ylocal * zlocal;
      }
    }
    Log::info("xIntegral: {}", totalIntegral);
    REQUIRE(totalIntegral == Approx(0.0833333333).margin(1e-9));
  }

  SECTION("Evaluate gradients") {
    for (int iTet = 0; iTet < tets.getnElements(); iTet++) {
      tets.loadNodes(iTet, elemNodes);
      coords.loadCoordinates(&elemNodes[0], &elemNodes[4], elemCoords);

      double determinant = tets.getDeterminant(iTet) * 1e18;
      s.initialiseElement(elemCoords, determinant);

      double x[] = {elemCoords[0].x(), elemCoords[1].x(), elemCoords[2].x(), elemCoords[3].x()};
      double y[] = {elemCoords[0].y(), elemCoords[1].y(), elemCoords[2].y(), elemCoords[3].y()};
      double z[] = {elemCoords[0].z(), elemCoords[1].z(), elemCoords[2].z(), elemCoords[3].z()};

      for (; s.hasNextPoint(); s.nextPoint()) {
        double xx = s.gradientX(x);
        double xy = s.gradientY(x);
        double xz = s.gradientZ(x);

        double yx = s.gradientX(y);
        double yy = s.gradientY(y);
        double yz = s.gradientZ(y);

        double zx = s.gradientX(z);
        double zy = s.gradientY(z);
        double zz = s.gradientZ(z);

        REQUIRE(xx == Approx(1).margin(1e-12));
        REQUIRE(xy == Approx(0).margin(1e-12));
        REQUIRE(xz == Approx(0).margin(1e-12));

        REQUIRE(yx == Approx(0).margin(1e-12));
        REQUIRE(yy == Approx(1).margin(1e-12));
        REQUIRE(yz == Approx(0).margin(1e-12));

        REQUIRE(zx == Approx(0).margin(1e-12));
        REQUIRE(zy == Approx(0).margin(1e-12));
        REQUIRE(zz == Approx(1).margin(1e-12));
      }
    }
  }
}

TEST_CASE("3D boundary integral in a unit cube") {
  ShapeSurf4thOrder s;

  SECTION("Check Gaussian intergration parameters") {
    REQUIRE(7 == s.numGaussPoints());

    double weightSum = 0;
    for (; s.hasNextPoint(); s.nextPoint()) {
      weightSum += s.getWeight();
    }
    REQUIRE(weightSum == Approx(0.5).margin(1e-12));
  }



  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  prepareGeometry(geom, TestUtil::RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH, *electrodes, {1, 1, 1});
  auto tets = geom.getTetrahedra();
  auto tris = geom.getTriangles();
  auto coords = geom.getCoordinates();
  idx tetNodes[4] = {0, 0, 0, 0};
  idx triNodes[3] = {0, 0, 0};

  Vec3 tetCoords[4] = {Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)};
  const Vec3 *triCoord = &tetCoords[0];

  SECTION("Integrate total surface area of Neumann boundaries in a unit cube") {
    double totalArea = 0;
    for (int iTri = 0; iTri < tris.getnElements(); iTri++) {
      unsigned int material = tris.getMaterialNumber(iTri);
      if (MAT_NEUMANN != material) {
        continue;
      }

      unsigned int iTet = tris.getConnectedVolume(iTri);

      double tetDet = tets.getDeterminant(iTet);
      double triDet = tris.getDeterminant(iTri);

      tris.loadNodes(iTri, triNodes);
      tets.loadNodes(iTet, tetNodes);

      // reorder tet nodes so that the triangle is the first three nodes
      reorderTetNodes(tetNodes, triNodes);

      coords.loadCoordinates(&tetNodes[0], &tetNodes[4], tetCoords);
      tetCoords[0] *= 1e-6;
      tetCoords[1] *= 1e-6;
      tetCoords[2] *= 1e-6;
      tetCoords[3] *= 1e-6;

      // calculate triangle area of the Numan boundary
      Vec3 v1 = triCoord[1] - triCoord[0];
      Vec3 v2 = triCoord[2] - triCoord[0];
      double area = 0.5 * v1.cross(v2).norm();

      // calculate triangle area using Gauss integration over the element
      double triArea = 0;
      s.initialiseElement(tetCoords, tetDet);
      for (; s.hasNextPoint(); s.nextPoint()) {
        double mul = s.getWeight() * triDet;
        for (int i = 0; i < 4; i++) {
          triArea += mul * s.N(i);
        }
      }

      REQUIRE(triArea == Approx(area).margin(1e-12));
      totalArea += triArea;
    }

    totalArea *= 1e12;
    Log::info("Total area: {}", totalArea);

    REQUIRE(totalArea == Approx(4.0).margin(1e-9));
  }

  SECTION("Integrate 3D gradients over Neumann boundaries in a unit cube") {

    for (int iTri = 0; iTri < tris.getnElements(); iTri++) {
      unsigned int material = tris.getMaterialNumber(iTri);
      if (MAT_NEUMANN != material) {
        continue;
      }

      Vec3 surfaceNorma = tris.getSurfaceNormal(iTri);

      unsigned int iTet = tris.getConnectedVolume(iTri);
      double tetDet = tets.getDeterminant(iTet);
      double triDet = tris.getDeterminant(iTri);

      tris.loadNodes(iTri, triNodes);
      tets.loadNodes(iTet, tetNodes);

      // reorder tet nodes so that the triangle is the first three nodes
      reorderTetNodes(tetNodes, triNodes);

      coords.loadCoordinates(&tetNodes[0], &tetNodes[4], tetCoords);
      tetCoords[0] *= 1e-6;
      tetCoords[1] *= 1e-6;
      tetCoords[2] *= 1e-6;
      tetCoords[3] *= 1e-6;

      double x[4] = {tetCoords[0].x(), tetCoords[1].x(), tetCoords[2].x(), tetCoords[3].x()};
      double y[4] = {tetCoords[0].y(), tetCoords[1].y(), tetCoords[2].y(), tetCoords[3].y()};
      double z[4] = {tetCoords[0].z(), tetCoords[1].z(), tetCoords[2].z(), tetCoords[3].z()};


      // calculate triangle area using Gauss integration over the element
      double triArea = 0;
      s.initialiseElement(tetCoords, tetDet);

      double xx = s.gradientX(x);
      double xy = s.gradientY(x);
      double xz = s.gradientZ(x);

      double totalXx = 0;

      for (; s.hasNextPoint(); s.nextPoint()) {
        double mul = s.getWeight() * triDet;
        for (int i = 0; i < 4; i++) {
          //totalXx += mul *
        }
      }
    }
  }
}