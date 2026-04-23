#include <refinement/element-selector.h>
#include <refinement/refinement-spec.h>
#include <geometry.h>
#include <solutionvector.h>
#include <globals.h>
#include <refinement.h>  // RED_TET, MAT_DIELECTRIC1, MAT_DOMAIN7
#include <stdexcept>
#include <algorithm>
#include <set>

#include "geom/coordinates.h"

namespace qlc3d::refinement {

// ── helpers shared between selectors ─────────────────────────────────────────

namespace {

/** Mark all LC tets that contain any node from @p indp as RED_TET. */
void markTetsByNodeIndex(const std::vector<idx>& indp,
                         const Geometry& geom,
                         std::vector<idx>& i_tet)
{
    auto& tets = geom.getTetrahedra();
    std::vector<std::set<idx>> p_to_t;
    tets.gen_p_to_elem(p_to_t);

    for (idx node : indp) {
        for (idx tet : p_to_t[node]) {
            if (tets.getMaterialNumber(tet) <= MAT_DOMAIN7) {
                i_tet[tet] = RED_TET;
            }
        }
    }
}

/** Return the max absolute change across all 5 Q-tensor components within element @p elem. */
double maxDeltaQ(idx elem, const Geometry& geom, const SolutionVector& q)
{
    auto& t = geom.getTetrahedra();
    const idx numNodes = t.getnNodes();
    double maxdq = 0.0;
    for (idx dim = 0; dim < 5; ++dim) {
        double qe[4] = {0, 0, 0, 0};
        for (idx j = 0; j < numNodes; ++j) {
            qe[j] = q.getValue(t.getNode(elem, j), dim);
        }
        double mxq = *std::max_element(qe, qe + 4);
        double mnq = *std::min_element(qe, qe + 4);
        if ((mxq - mnq) > maxdq) maxdq = mxq - mnq;
    }
    return maxdq;
}

} // anonymous namespace

// ── ChangeSelector ────────────────────────────────────────────────────────────

class ChangeSelector : public ElementSelector {
    const RefinementSpec& spec_;
public:
    explicit ChangeSelector(const RefinementSpec& spec) : spec_(spec) {}

    [[nodiscard]] bool needsQTensor() const override { return true; }
    [[nodiscard]] unsigned int getRefIter() const override { return spec_.getRefIter(); }

    void selectTets(const Geometry& geom,
                    const SolutionVector* q,
                    int refIter,
                    std::vector<idx>& i_tet) const override
    {
        if (!q) throw std::invalid_argument("ChangeSelector requires a non-null SolutionVector");
        double threshold = spec_.getValue(static_cast<size_t>(refIter));
        auto& tets = geom.getTetrahedra();
        const idx numTets = tets.getnElements();
        for (idx i = 0; i < numTets; ++i) {
            if (tets.getMaterialNumber(i) >= MAT_DIELECTRIC1) continue;
            if (maxDeltaQ(i, geom, *q) >= threshold) {
                i_tet[i] = RED_TET;
            }
        }
    }
};

// ── SphereSelector ────────────────────────────────────────────────────────────

class SphereSelector : public ElementSelector {
    const RefinementSpec& spec_;
public:
    explicit SphereSelector(const RefinementSpec& spec) : spec_(spec) {}

    [[nodiscard]] unsigned int getRefIter() const override { return spec_.getRefIter(); }

    void selectTets(const Geometry& geom,
                    const SolutionVector* /*q*/,
                    int refIter,
                    std::vector<idx>& i_tet) const override
    {
        double rad = spec_.getValue(static_cast<size_t>(refIter));
        double radSq = rad * rad;
        const Vec3 centre(spec_.getX()[0], spec_.getY()[0], spec_.getZ()[0]);

        std::vector<idx> nodesInside;
        auto& coords = geom.getCoordinates();
        for (idx i = 0; i < geom.getnpLC(); ++i) {
            if (coords.getPoint(i).distanceSquared(centre) <= radSq) {
                nodesInside.push_back(i);
            }
        }
        if (!nodesInside.empty()) {
            markTetsByNodeIndex(nodesInside, geom, i_tet);
        }
    }
};

// ── BoxSelector ───────────────────────────────────────────────────────────────

class BoxSelector : public ElementSelector {
    const RefinementSpec& spec_;
public:
    explicit BoxSelector(const RefinementSpec& spec) : spec_(spec) {}

    [[nodiscard]] unsigned int getRefIter() const override { return spec_.getRefIter(); }

    void selectTets(const Geometry& geom,
                    const SolutionVector* /*q*/,
                    int refIter,
                    std::vector<idx>& i_tet) const override
    {
        auto& x = spec_.getX();
        auto& y = spec_.getY();
        auto& z = spec_.getZ();

        double xMin = x[refIter * 2],     xMax = x[refIter * 2 + 1];
        double yMin = y[refIter * 2],     yMax = y[refIter * 2 + 1];
        double zMin = z[refIter * 2],     zMax = z[refIter * 2 + 1];

        std::vector<idx> nodesInside;
        auto& coords = geom.getCoordinates();
        for (idx i = 0; i < geom.getnpLC(); ++i) {
            auto& p = coords.getPoint(i);
            if (p.x() >= xMin && p.x() <= xMax &&
                p.y() >= yMin && p.y() <= yMax &&
                p.z() >= zMin && p.z() <= zMax)
            {
                nodesInside.push_back(i);
            }
        }
        if (!nodesInside.empty()) {
            markTetsByNodeIndex(nodesInside, geom, i_tet);
        }
    }
};

// ── Factory methods ───────────────────────────────────────────────────────────

std::unique_ptr<ElementSelector> ElementSelector::makeChange(const RefinementSpec& spec)
{
    return std::make_unique<ChangeSelector>(spec);
}

std::unique_ptr<ElementSelector> ElementSelector::makeSphere(const RefinementSpec& spec)
{
    return std::make_unique<SphereSelector>(spec);
}

std::unique_ptr<ElementSelector> ElementSelector::makeBox(const RefinementSpec& spec)
{
    return std::make_unique<BoxSelector>(spec);
}

std::unique_ptr<ElementSelector> ElementSelector::fromSpec(const RefinementSpec& spec)
{
    switch (spec.getType()) {
        case RefinementSpec::Type::Change: return makeChange(spec);
        case RefinementSpec::Type::Sphere: return makeSphere(spec);
        case RefinementSpec::Type::Box:    return makeBox(spec);
        default:
            throw std::invalid_argument("ElementSelector::fromSpec: unrecognised RefinementSpec type");
    }
}

} // namespace qlc3d::refinement

