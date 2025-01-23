#include <catch.h>
#include <fe/fe-util.h>
#include <geom/vec3.h>
#include <geom/coordinates.h>
#include <geometry.h>
#include <electrodes.h>
#include <test-util.h>
#include <inits.h>
#include <util/logging.h>
#include <fe/gaussian-quadrature.h>
#include <fe/keast.h>

TEST_CASE("Linear tet 3D shape function") {
  GaussianQuadratureTet<11> g = gaussQuadratureTet4thOrder();


  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  auto alignment = Alignment();
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
  prepareGeometry(geom, TestUtil::RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH, electrodes, alignment, {1, 1, 1});
  auto tets = geom.getTetrahedra();
  auto coords = geom.getCoordinates();
  idx elemNodes[4] = {0, 0, 0, 0};
  Vec3 elemCoords[4] = {Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)};

  SECTION("Check gaussian integration parameters") {
    REQUIRE(11 == g.numGaussPoints());

    double sumWeight = 0;
    for(;g.hasNextPoint(); g.nextPoint()) {
      sumWeight += g.weight();
    }

    REQUIRE(sumWeight == Approx(1.0 / 6.).margin(1e-12)); // sum of weights should be 1/6 because Tetrahedron volume is 6x determinant
  }

  SECTION("Integrate total volume of unit cube") {
    double totalVolume = 0;

    for (int iTet = 0; iTet < tets.getnElements(); iTet++) {
      tets.loadNodes(iTet, elemNodes);
      coords.loadCoordinates(&elemNodes[0], &elemNodes[4], elemCoords);

      double determinant = tets.getDeterminant(iTet);
      g.initialiseElement(elemCoords, determinant);

      for (; g.hasNextPoint(); g.nextPoint()) {

        double mul = g.weight() * determinant;
        for (int i = 0; i < 4; i++) {
          totalVolume += mul * g.N(i);
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
      g.initialiseElement(elemCoords, determinant);

      double x[] = {elemCoords[0].x(), elemCoords[1].x(), elemCoords[2].x(), elemCoords[3].x()};
      double y[] = {elemCoords[0].y(), elemCoords[1].y(), elemCoords[2].y(), elemCoords[3].y()};
      double z[] = {elemCoords[0].z(), elemCoords[1].z(), elemCoords[2].z(), elemCoords[3].z()};

      for (;g.hasNextPoint(); g.nextPoint()) {

        double xlocal = g.sample(x);
        double ylocal = g.sample(y);
        double zlocal = g.sample(z);

        double mul = g.weight() * determinant;
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
      g.initialiseElement(elemCoords, determinant);

      double x[] = {elemCoords[0].x(), elemCoords[1].x(), elemCoords[2].x(), elemCoords[3].x()};
      double y[] = {elemCoords[0].y(), elemCoords[1].y(), elemCoords[2].y(), elemCoords[3].y()};
      double z[] = {elemCoords[0].z(), elemCoords[1].z(), elemCoords[2].z(), elemCoords[3].z()};

      for (; g.hasNextPoint(); g.nextPoint()) {
        double xx = g.sampleX(x);
        double xy = g.sampleY(x);
        double xz = g.sampleZ(x);

        double yx = g.sampleX(y);
        double yy = g.sampleY(y);
        double yz = g.sampleZ(y);

        double zx = g.sampleX(z);
        double zy = g.sampleY(z);
        double zz = g.sampleZ(z);

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

  SECTION("Sampling 6 nodal values (permittivity tensor)") {
    double nodalValues[4][6] = {
            {1, 2, 3, 4, 5, 6}, // values at node 1
            {1, 2, 3, 4, 5, 6}, // values at node 2
            {1, 2, 3, 4, 5, 6}, // values at node 3
            {1, 2, 3, 4, 5, 6}  // values at node 4
            };

    double v1, v2, v3, v4, v5, v6;
    g.sampleAll(nodalValues, v1, v2, v3, v4, v5, v6);

    REQUIRE(v1 == Approx(1).margin(1e-12));
    REQUIRE(v2 == Approx(2).margin(1e-12));
    REQUIRE(v3 == Approx(3).margin(1e-12));
    REQUIRE(v4 == Approx(4).margin(1e-12));
    REQUIRE(v5 == Approx(5).margin(1e-12));
    REQUIRE(v6 == Approx(6).margin(1e-12));
  }

  SECTION("Sampling 5 nodal values (Q-tensor)") {
    double nodalValues[4][5] = {
            {1, 2, 3, 4, 5}, // values at node 1
            {1, 2, 3, 4, 5}, // values at node 2
            {1, 2, 3, 4, 5}, // values at node 3
            {1, 2, 3, 4, 5}  // values at node 4
    };

    double v1, v2, v3, v4, v5;
    g.sampleAll(nodalValues, v1, v2, v3, v4, v5);

    REQUIRE(v1 == Approx(1).margin(1e-12));
    REQUIRE(v2 == Approx(2).margin(1e-12));
    REQUIRE(v3 == Approx(3).margin(1e-12));
    REQUIRE(v4 == Approx(4).margin(1e-12));
    REQUIRE(v5 == Approx(5).margin(1e-12));
  }
}

TEST_CASE("New tet 3D shape function - linear tet") { // TODO: repeat this with quadratic element
  TetShapeFunction shape(1); // = createLinearTetShapeFunction();

  shape.initialise(Keast4);
  //shape.initialise(Keast0);


  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  auto alignment = Alignment();
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
  prepareGeometry(geom, TestUtil::RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH, electrodes, alignment, {1, 1, 1});
  auto tets = geom.getTetrahedra();
  auto coords = geom.getCoordinates();
  idx elemNodes[4] = {0, 0, 0, 0};
  Vec3 elemCoords[4] = {Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)};



  SECTION("Check gaussian integration parameters") {
    //REQUIRE(11 == shape.getNumGaussPoints());

    double sumWeight = 0;
    for (; shape.hasNextPoint(); shape.nextPoint()) {
      sumWeight += shape.getWeight();
    }

    REQUIRE(sumWeight == Approx(1.0 / 6.).margin(1e-12)); // sum of weights should be 1/6 because Tetrahedron volume is 6x determinant
  }

  SECTION("Integrate total volume of unit cube") {
    double totalVolume = 0;

    for (int iTet = 0; iTet < tets.getnElements(); iTet++) {
      tets.loadNodes(iTet, elemNodes);
      coords.loadCoordinates(&elemNodes[0], &elemNodes[4], elemCoords);

      double determinant = tets.getDeterminant(iTet);


      for (; shape.hasNextPoint(); shape.nextPoint()) {
        shape.initialiseElement(elemCoords, determinant);
        double mul = shape.getWeight() * determinant;
        for (int i = 0; i < 4; i++) {
          totalVolume += mul * shape.N(i);
        }
      }
    }

    totalVolume *= 1e18;
    Log::info("Total volume: {}", totalVolume);

    REQUIRE(totalVolume == Approx(1.0).margin(1e-9));
  }

  SECTION("Evaluate gradients") {
    for (int iTet = 0; iTet < tets.getnElements(); iTet++) {
      tets.loadNodes(iTet, elemNodes);
      coords.loadCoordinates(&elemNodes[0], &elemNodes[4], elemCoords);

      double determinant = tets.getDeterminant(iTet) * 1e18;
      double x[] = {elemCoords[0].x(), elemCoords[1].x(), elemCoords[2].x(), elemCoords[3].x()};
      double y[] = {elemCoords[0].y(), elemCoords[1].y(), elemCoords[2].y(), elemCoords[3].y()};
      double z[] = {elemCoords[0].z(), elemCoords[1].z(), elemCoords[2].z(), elemCoords[3].z()};

      for (; shape.hasNextPoint(); shape.nextPoint()) {
        shape.initialiseElement(elemCoords, determinant);
        double xx = shape.sampleX(x);
        double xy = shape.sampleY(x);
        double xz = shape.sampleZ(x);

        double yx = shape.sampleX(y);
        double yy = shape.sampleY(y);
        double yz = shape.sampleZ(y);

        double zx = shape.sampleX(z);
        double zy = shape.sampleY(z);
        double zz = shape.sampleZ(z);

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

TEST_CASE("New tet 3D shape function - quadratic tet") {
  TetShapeFunction shape(2); // = createLinearTetShapeFunction();
  shape.initialise(Keast4);


  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  auto alignment = Alignment();
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
  prepareGeometry(geom, TestUtil::RESOURCE_SMALL_CUBE_QUADRATIC_GMSH_MESH, electrodes, alignment, {1, 1, 1}, 0, 0, 0);
  auto tets = geom.getTetrahedra();
  auto coords = geom.getCoordinates();
  idx elemNodes[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Vec3 elemCoords[10] = {Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0),
                         Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)};

  SECTION("Check gaussian integration parameters") {
    //REQUIRE(11 == shape.getNumGaussPoints());

    double sumWeight = 0;
    for (; shape.hasNextPoint(); shape.nextPoint()) {
      sumWeight += shape.getWeight();
    }

    REQUIRE(sumWeight == Approx(1.0 / 6.).margin(
            1e-12)); // sum of weights should be 1/6 because Tetrahedron volume is 6x determinant
  }

  SECTION("Integrate total volume of unit cube") {
    double totalVolume = 0;

    for (int iTet = 0; iTet < tets.getnElements(); iTet++) {
      tets.loadNodes(iTet, elemNodes);
      coords.loadCoordinates(&elemNodes[0], &elemNodes[tets.getnNodes()], elemCoords);

      double determinant = tets.getDeterminant(iTet);
      shape.initialiseElement(elemCoords, determinant);

      for (; shape.hasNextPoint(); shape.nextPoint()) {

        double mul = shape.getWeight() * determinant;
        for (int i = 0; i < tets.getnNodes(); i++) {
          totalVolume += mul * shape.N(i);
        }
      }
    }

    totalVolume *= 1e18;
    Log::info("Total volume: {}", totalVolume);

    REQUIRE(totalVolume == Approx(1.0).margin(1e-9));
  }

  SECTION("Evaluate gradients") {
    for (int iTet = 0; iTet < tets.getnElements(); iTet++) {
      tets.loadNodes(iTet, elemNodes);

      // swap last nodes because gmsh order is different
      // TODO: this should be done during reading/initialisation for each element, if detected?
      std::swap(elemNodes[8], elemNodes[9]);

      coords.loadCoordinates(&elemNodes[0], &elemNodes[10], elemCoords);
      double x[] = {elemCoords[0].x(), elemCoords[1].x(), elemCoords[2].x(), elemCoords[3].x(), elemCoords[4].x(), elemCoords[5].x(), elemCoords[6].x(), elemCoords[7].x(), elemCoords[8].x(), elemCoords[9].x()};
      double y[] = {elemCoords[0].y(), elemCoords[1].y(), elemCoords[2].y(), elemCoords[3].y(), elemCoords[4].y(), elemCoords[5].y(), elemCoords[6].y(), elemCoords[7].y(), elemCoords[8].y(), elemCoords[9].y()};
      double z[] = {elemCoords[0].z(), elemCoords[1].z(), elemCoords[2].z(), elemCoords[3].z(), elemCoords[4].z(), elemCoords[5].z(), elemCoords[6].z(), elemCoords[7].z(), elemCoords[8].z(), elemCoords[9].z()};

      double determinant = tets.getDeterminant(iTet) * 1e18;

      for (; shape.hasNextPoint(); shape.nextPoint()) {
        shape.initialiseElement(elemCoords, determinant); // TODO: call this inside nextPoint()?
        double xx = shape.sampleX(x);
        double xy = shape.sampleY(x);
        double xz = shape.sampleZ(x);

        double yx = shape.sampleX(y);
        double yy = shape.sampleY(y);
        double yz = shape.sampleZ(y);

        double zx = shape.sampleX(z);
        double zy = shape.sampleY(z);
        double zz = shape.sampleZ(z);

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
  GaussianQuadratureTet<7> g = gaussQuadratureTetBoundaryIntegral4thOrder();
  SECTION("Check Gaussian intergration parameters") {
    REQUIRE(7 == g.numGaussPoints());

    double weightSum = 0;
    for (; g.hasNextPoint(); g.nextPoint()) {
      weightSum += g.weight();
    }
    REQUIRE(weightSum == Approx(0.5).margin(1e-12));
  }

  Geometry geom;
  auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
  auto alignment = Alignment();
  alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
  alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
  prepareGeometry(geom, TestUtil::RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH, electrodes, alignment, {1, 1, 1});
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
      reorderBoundaryTetNodes(tetNodes, triNodes);

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
      g.initialiseElement(tetCoords, tetDet);
      for (; g.hasNextPoint(); g.nextPoint()) {
        double mul = g.weight() * triDet;
        for (int i = 0; i < 4; i++) {
          triArea += mul * g.N(i);
        }
      }

      REQUIRE(triArea == Approx(area).margin(1e-12));
      totalArea += triArea;
    }

    totalArea *= 1e12;
    Log::info("Total area: {}", totalArea);

    REQUIRE(totalArea == Approx(4.0).margin(1e-9));
  }
}

TEST_CASE("Linear triangle 2D share function") {
  GaussianQuadratureTri<7> g = gaussianQuadratureTri4thOrder();

  SECTION("Sum of weights should equal 0.5") {
    double sumWeight = 0;
    for (; g.hasNextPoint(); g.nextPoint()) {
      sumWeight += g.weight();
    }
    REQUIRE(sumWeight == Approx(0.5).margin(1e-12)); // because determinant is 2x triangle area
  }

  SECTION("Integrate total area of a triangle") {

    // triangle with coordinates (0, 0, 0), (1, 0, 0), (0, 1, 0)
    Vec3 elemCoordinates[3] = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0)};
    const double area = 0.5;
    double determinant = 2 * area; // determinant is half triangle are;

    g.initialiseElement(elemCoordinates, determinant);

    double totalArea = 0;
    for (; g.hasNextPoint(); g.nextPoint()) {
      double mul = g.weight() * determinant;
      for (int i = 0; i < 3; i++) {
       totalArea += mul * g.N(i);
      }
    }

    REQUIRE(totalArea == Approx(area).margin(1e-12));
  }

  SECTION("Integrate total surface area of a unit cube mesh") {
    Geometry geom;
    auto electrodes = Electrodes::withInitialPotentials({1, 2}, {1, 0});
    auto alignment = Alignment();
    alignment.addSurface(Surface::ofStrongAnchoring(1, 0, 0));
    alignment.addSurface(Surface::ofStrongAnchoring(2, 0, 0));
    prepareGeometry(geom, TestUtil::RESOURCE_UNIT_CUBE_NEUMANN_GMSH_MESH, electrodes, alignment, {1, 1, 1});
    auto tris = geom.getTriangles();
    auto coords = geom.getCoordinates();
    idx elemNodes[3] = {0, 0, 0};
    Vec3 elemCoords[3] = {Vec3(0, 0, 0), Vec3(0, 0, 0), Vec3(0, 0, 0)};

    double totalArea = 0;
    for (idx it = 0; it < tris.getnElements(); it++) {
      tris.loadNodes(it, elemNodes);
      coords.loadCoordinates(&elemNodes[0], &elemNodes[3], elemCoords);

      double determinant = tris.getDeterminant(it);
      g.initialiseElement(elemCoords, determinant);

      double area = 0;
      for (; g.hasNextPoint(); g.nextPoint()) {
        double mul = g.weight() * determinant;
        for (int i = 0; i < 3; i++) {
          area += mul * g.N(i);
        }
      }
      totalArea += area;
    }

    // convert coordinates from microns to metres
    totalArea *= 1e12;
    REQUIRE(totalArea == Approx(6.0).margin(1e-9));
  }
}