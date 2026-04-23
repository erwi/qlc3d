#include <refinement/q-tensor-interpolator.h>
#include <geometry.h>
#include <geom/coordinates.h>
#include <solutionvector.h>
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
    const size_t npLC = static_cast<size_t>(geomNew.getnpLC());

    for (size_t ind = 0; ind < npLC; ind++) {
        // Only LC volumes are interpolated.
        assert(geomOld.getTetrahedra().getMaterialNumber(pInTet[ind]) == MAT_DOMAIN1);

        // Compute barycentric coordinates of the new node inside the old tet.
        double loc[4];
        Vec3 targetPoint = geomNew.getCoordinates().getPoint(ind);
        oldTets.calcLocCoords(pInTet[ind], geomOld.getCoordinates(), targetPoint, loc);

        size_t n[4];
        for (int k = 0; k < 4; k++) {
            n[k] = oldTets.getNode(pInTet[ind], k);
        }

        // Corner-node fast path: if the largest barycentric coordinate is
        // close to 1.0 the new node coincides with an existing corner node,
        // so copy its Q values directly without interpolation.
        size_t ind_max = std::max_element(loc, loc + 4) - loc;
        if (loc[ind_max] >= 0.99999) {
            for (int i = 0; i < 5; i++) {
                qNew.setValue(ind, i, qOld.getValue(n[ind_max], i));
            }
        } else {
            // General path: interpolate each of the 5 Q-tensor components
            // independently as a scalar field using barycentric weights.
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

