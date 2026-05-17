// Combined LC energy integrator - implements both volume and surface integration.
// Volume integration uses tetrahedral Gauss quadrature (Keast8).
// Surface integration uses triangular Gauss quadrature (Tri4thOrder).
// Both integrals are also individually callable for independent testing.

#include <energy/lc-energy-integrator.h>
#include <energy/lc-energy-density.h>
#include <fe/gaussian-quadrature.h>
#include <fe/keast.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <solutionvector.h>
#include <mesh/mesh.h>
#include <lc.h>
#include <globals.h>
#include <material_numbers.h>
#include <lc-representation.h>
#include <alignment.h>
#include <util/logging.h>
#include <cmath>

// ============================================================================
// Volume integration
// ============================================================================

EnergyResult integrateVolumeEnergy(const LC &lc, const Geometry &geom,
                                   const SolutionVector &v, const SolutionVector &q) {
    const EnergyMaterialParams mat = lcToEnergyMaterialParams(lc);
    const Mesh &tets = geom.getTetrahedra();
    const Coordinates &coords = geom.getCoordinates();

    const ElementType elementType = tets.getElementType();
    const unsigned int elementOrder = getElementOrder(elementType);
    const unsigned int npe = getNodesPerElement(elementType);

    TetShapeFunction shapes(elementOrder);
    shapes.setIntegrationPoints(Keast8);

    EnergyResult result{};

    std::vector<unsigned int> tetNodes(npe, 0);
    std::vector<Vec3> tetCoords(npe, Vec3());
    std::vector<qlc3d::TTensor> qNodal(npe, qlc3d::TTensor());
    std::vector<double> potential(npe, 0.0);

    const unsigned int elementCount = tets.getnElements();

    for (unsigned int indTet = 0; indTet < elementCount; indTet++) {
        if (tets.getMaterialNumber(indTet) > MAT_DOMAIN7) {
            continue;
        }

        // Load node indices, coordinates (converting µm → m), Q-tensor and potential values
        tets.loadNodes(indTet, tetNodes.data());
        coords.loadCoordinates(tetNodes.data(), tetNodes.data() + npe, tetCoords.data());
        for (auto &c : tetCoords) {
            c *= qlc3d::units::MICROMETER_TO_METER;
        }
        q.loadQtensorValues(tetNodes.data(), tetNodes.data() + npe, qNodal.data());
        v.loadValues(tetNodes.data(), tetNodes.data() + npe, potential.data());

        const double det = tets.getDeterminant(indTet);

        for (; shapes.hasNextPoint(); shapes.nextPoint()) {
            shapes.initialiseElement(tetCoords.data(), det);

            // Interpolate Q-tensor components and their spatial derivatives at this Gauss point
            double q1, q2, q3, q4, q5;
            double q1x, q2x, q3x, q4x, q5x;
            double q1y, q2y, q3y, q4y, q5y;
            double q1z, q2z, q3z, q4z, q5z;
            shapes.sampleQ(qNodal, q1, q2, q3, q4, q5);
            shapes.sampleQX(qNodal, q1x, q2x, q3x, q4x, q5x);
            shapes.sampleQY(qNodal, q1y, q2y, q3y, q4y, q5y);
            shapes.sampleQZ(qNodal, q1z, q2z, q3z, q4z, q5z);

            // Electric field E = -∇φ (potential gradient is ∂φ/∂x_i → E_i = -∂φ/∂x_i)
            double Vx = shapes.sampleX(potential.data());
            double Vy = shapes.sampleY(potential.data());
            double Vz = shapes.sampleZ(potential.data());

            GaussPointData gp;
            gp.q1 = q1; gp.q2 = q2; gp.q3 = q3; gp.q4 = q4; gp.q5 = q5;
            gp.q1x = q1x; gp.q1y = q1y; gp.q1z = q1z;
            gp.q2x = q2x; gp.q2y = q2y; gp.q2z = q2z;
            gp.q3x = q3x; gp.q3y = q3y; gp.q3z = q3z;
            gp.q4x = q4x; gp.q4y = q4y; gp.q4z = q4z;
            gp.q5x = q5x; gp.q5y = q5y; gp.q5z = q5z;
            // Store potential gradient as Ex=-Vx etc. (E = -∇φ)
            gp.Ex = -Vx; gp.Ey = -Vy; gp.Ez = -Vz;

            const double weight = shapes.getWeight() * det;
            result.elastic      += weight * elasticEnergyDensity(gp, mat);
            result.thermotropic += weight * thermotropicEnergyDensity(gp, mat);
            result.electric     += weight * electricEnergyDensity(gp, mat);
        }
    }

    return result;
}

// ============================================================================
// Surface integration
// ============================================================================


double integrateSurfaceEnergy(const Alignment &alignment, const Geometry &geom,
                              const SolutionVector &q, double S0) {
    const Mesh &tris = geom.getTriangles();
    const unsigned int elementCount = tris.getnElements();

    const ElementType elementType = tris.getElementType();
    const unsigned int elementOrder = getElementOrder(elementType);
    const unsigned int npe = tris.getnNodes();

    // Get weak surface map: FIXLC number → Surface
    std::unordered_map<unsigned int, Surface> weakSurfaces = alignment.getWeakSurfacesByFixLcNumber();
    if (weakSurfaces.empty()) {
        return 0.0;
    }

    TriShapeFunction shapes(elementOrder);
    shapes.setIntegrationPoints(Tri4thOrder);

    std::vector<unsigned int> triNodes(npe, 0);
    std::vector<qlc3d::TTensor> qNodal(npe, qlc3d::TTensor());
    std::vector<Vec3> triCoords(npe, Vec3());
    std::vector<Vec3> vec1(npe, Vec3());
    std::vector<Vec3> vec2(npe, Vec3());

    double totalSurfaceEnergy = 0.0;

    for (unsigned int indTri = 0; indTri < elementCount; indTri++) {
        // Identify which anchoring condition applies via FIXLC number
        unsigned int fixLcNumber = tris.getFixLCNumber(indTri);
        auto it = weakSurfaces.find(fixLcNumber);
        if (it == weakSurfaces.end()) {
            continue; // not a weak-anchoring surface
        }
        const Surface &surface = it->second;

        tris.loadNodes(indTri, triNodes.data());
        geom.getCoordinates().loadCoordinates(triNodes.data(), triNodes.data() + npe, triCoords.data());
        for (auto &c : triCoords) {
            c *= qlc3d::units::MICROMETER_TO_METER;
        }
        q.loadQtensorValues(triNodes.data(), triNodes.data() + npe, qNodal.data());

        // Determine anchoring principal axes at each node
        if (surface.usesSurfaceNormal()) {
            // Homeotropic-style: v2 is the surface normal, v1 is zero (K1=0 enforced)
            for (unsigned int i = 0; i < npe; i++) {
                vec1[i] = Vec3(0, 0, 0);
                vec2[i] = geom.getNodeNormal(triNodes[i]);
            }
        } else {
            for (unsigned int i = 0; i < npe; i++) {
                vec1[i] = surface.getV1();
                vec2[i] = surface.getV2();
            }
        }

        // W is negative for WeakHomeotropic (repulsive) to match solver convention
        const double W = surface.getAnchoringType() == WeakHomeotropic
                         ? -surface.getStrength()
                         :  surface.getStrength();
        const double K1 = surface.getK1();
        const double K2 = surface.getK2();

        const double triDet = tris.getDeterminant(indTri);

        for (; shapes.hasNextPoint(); shapes.nextPoint()) {
            double q1, q2, q3, q4, q5;
            shapes.sampleQ(qNodal.data(), q1, q2, q3, q4, q5);

            Vec3 v1, v2;
            shapes.sample(vec1, v1);
            shapes.sample(vec2, v2);

            double density = surfaceEnergyDensity(q1, q2, q3, q4, q5, v1, v2, W, K1, K2, S0);
            totalSurfaceEnergy += shapes.getWeight() * triDet * density;
        }
    }
    return totalSurfaceEnergy;
}

// ============================================================================
// Combined integrator
// ============================================================================

EnergyResult integrateEnergy(const LC &lc, const Geometry &geom,
                             const SolutionVector &v, const SolutionVector &q,
                             const Alignment &alignment) {
    EnergyResult result = integrateVolumeEnergy(lc, geom, v, q);
    result.surface = integrateSurfaceEnergy(alignment, geom, q, lc.S0());
    return result;
}

