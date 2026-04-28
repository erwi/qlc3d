#include <refinement/q-tensor-interpolator.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <solutionvector.h>
#include <material_numbers.h>
#include <util/exception.h>
#include <mesh/mesh.h>
#include <algorithm>
#include <cassert>

// MAT_DOMAIN1 is defined in lc.h
#include <lc.h>

double interpolateScalar(const double loc[4], const double S[4]) {
    // Linearly interpolate a scalar field at a point given by barycentric coordinates.
    // Result = sum_i(loc[i] * S[i]) over the four tet corners.
    double val = 0;
    for (int i = 0; i < 4; i++) {
        val += loc[i] * S[i];
    }
    return val;
}

void evaluateTet10ShapeFunctions(const double loc[4], double N[10]) {
    // Evaluate the 10 quadratic TET10 shape functions at a point given by
    // barycentric (volume) coordinates loc[0..3].
    //
    // Volume coordinate mapping (identical to gaussianquadrature.h initialiseQuadraticTet):
    //   N1 = loc[0] = 1 - r - s - t
    //   N2 = loc[1] = r
    //   N3 = loc[2] = s
    //   N4 = loc[3] = t
    //
    // Node ordering (Gmsh TET10):
    //   [0] corner A  [1] corner B  [2] corner C  [3] corner D
    //   [4] mid AB    [5] mid BC    [6] mid AC    [7] mid AD
    //   [8] mid CD    [9] mid BD
    const double N1 = loc[0];
    const double N2 = loc[1];
    const double N3 = loc[2];
    const double N4 = loc[3];

    // Corner nodes: Ni*(2*Ni - 1)
    N[0] = N1 * (2.0 * N1 - 1.0);
    N[1] = N2 * (2.0 * N2 - 1.0);
    N[2] = N3 * (2.0 * N3 - 1.0);
    N[3] = N4 * (2.0 * N4 - 1.0);

    // Mid-edge nodes: 4*Ni*Nj
    N[4] = 4.0 * N1 * N2; // AB
    N[5] = 4.0 * N2 * N3; // BC
    N[6] = 4.0 * N3 * N1; // AC
    N[7] = 4.0 * N1 * N4; // AD
    N[8] = 4.0 * N3 * N4; // CD
    N[9] = 4.0 * N2 * N4; // BD
}

SolutionVector interpolateQTensor(Geometry &geomNew,
                                  Geometry &geomOld,
                                  const SolutionVector &qOld) {
    // Map each new LC node to the old tet that contains it using the regular
    // grid for fast spatial lookup.
    std::vector<idx> pInTet;
    geomOld.genIndToTetsByCoords(pInTet,
                                 geomNew.getCoordinates(),
                                 true,   // error if coordinate is not found
                                 true    // only consider LC elements
                                 );

    SolutionVector qNew;
    qNew.Resize(geomNew.getnpLC(), 5);

    auto &oldTets = geomOld.getTetrahedra();
    const ElementType srcType = oldTets.getElementType();
    const size_t npLC = static_cast<size_t>(geomNew.getnpLC());

    if (srcType != ElementType::LINEAR_TETRAHEDRON &&
        srcType != ElementType::QUADRATIC_TETRAHEDRON) {
        RUNTIME_ERROR("interpolateQTensor: unsupported source element type " + toString(srcType));
    }

    for (size_t ind = 0; ind < npLC; ind++) {
        // Only LC volumes are interpolated.
        assert(geomOld.getTetrahedra().getMaterialNumber(pInTet[ind]) == MAT_DOMAIN1);

        // Compute barycentric coordinates of the new node inside the old tet.
        double loc[4];
        Vec3 targetPoint = geomNew.getCoordinates().getPoint(ind);
        oldTets.calcLocCoords(pInTet[ind], geomOld.getCoordinates(), targetPoint, loc);

        // Read the first 4 (corner) node indices — used by both paths and the
        // fast-path copy below.
        size_t n[4];
        for (int k = 0; k < 4; k++) {
            n[k] = oldTets.getNode(pInTet[ind], k);
        }

        // Corner-node fast path: if the largest barycentric coordinate is
        // close to 1.0 the new node coincides with an existing corner node,
        // so copy its Q values directly without interpolation.  This applies
        // to both TET4 and TET10 sources because the barycentric coordinates
        // always refer to the four corner nodes.
        const size_t ind_max = std::max_element(loc, loc + 4) - loc;
        if (loc[ind_max] >= 0.99999) {
            for (int i = 0; i < 5; i++) {
                qNew.setValue(ind, i, qOld.getValue(n[ind_max], i));
            }
        } else if (srcType == ElementType::QUADRATIC_TETRAHEDRON) {
            // TET10 quadratic path: evaluate all 10 shape functions and
            // accumulate contributions from all 10 nodes.  This reproduces
            // polynomial fields up to degree 2 exactly.
            double N[10];
            evaluateTet10ShapeFunctions(loc, N);

            // Mid-edge fast path: if any mid-edge shape function is ~1 the new
            // node coincides with that mid-edge node, so copy directly.
            // Mid-edge nodes are indices 4..9.
            int midEdgeFastNode = -1;
            for (int k = 4; k < 10; k++) {
                if (N[k] >= 0.99999) {
                    midEdgeFastNode = static_cast<int>(k);
                    break;
                }
            }

            if (midEdgeFastNode >= 0) {
                const size_t srcNode = oldTets.getNode(pInTet[ind], midEdgeFastNode);
                for (int i = 0; i < 5; i++) {
                    qNew.setValue(ind, i, qOld.getValue(srcNode, i));
                }
            } else {
                for (int i = 0; i < 5; i++) {
                    double val = 0.0;
                    for (int k = 0; k < 10; k++) {
                        val += N[k] * qOld.getValue(oldTets.getNode(pInTet[ind], k), i);
                    }
                    qNew.setValue(ind, i, val);
                }
            }
        } else {
            // LINEAR_TETRAHEDRON path: linearly interpolate over the 4 corner nodes.
            for (int i = 0; i < 5; i++) {
                double qo[4] = {
                    qOld.getValue(n[0], i),
                    qOld.getValue(n[1], i),
                    qOld.getValue(n[2], i),
                    qOld.getValue(n[3], i)
                };
                qNew.setValue(ind, i, interpolateScalar(loc, qo));
            }
        }
    }

    return qNew;
}

