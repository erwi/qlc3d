#ifndef QLC3D_REFINEMENT_MIDPOINT_NODE_LOOKUP_H
#define QLC3D_REFINEMENT_MIDPOINT_NODE_LOOKUP_H

#include <globals.h>
#include <unordered_map>
#include <utility>

namespace qlc3d::refinement {

using MidpointEdge = std::pair<idx, idx>;

namespace detail {

/** Returns the canonical ordering for an undirected edge. */
inline MidpointEdge canonicalEdge(idx n1, idx n2) {
    return n1 < n2 ? MidpointEdge{n1, n2} : MidpointEdge{n2, n1};
}

/** Hash functor for undirected midpoint-edge keys. */
struct MidpointEdgeHash {
    std::size_t operator()(const MidpointEdge &edge) const noexcept;
};

} // namespace detail

/**
 * @brief Abstract edge-to-midpoint lookup.
 *
 * Implementations must treat edges as undirected, so lookup(1, 2) and
 * lookup(2, 1) return the same midpoint node index.
 */
class MidpointNodeLookup {
public:
    virtual ~MidpointNodeLookup() = default;

    /**
     * @brief Look up the midpoint node index for an edge.
     *
     * @param n1 First endpoint of the edge.
     * @param n2 Second endpoint of the edge.
     * @return The node index of the inserted midpoint.
     */
    virtual idx lookup(idx n1, idx n2) const = 0;
};


/**
 * @brief Midpoint lookup backed by an unordered map.
 *
 * This implementation is convenient for unit tests and small synthetic cases.
 */
class MapMidpointNodeLookup final : public MidpointNodeLookup {
public:
    /** Register the midpoint node for an undirected edge. */
    void registerEdge(idx n1, idx n2, idx midpointNode);

    /**
     * @brief Look up the midpoint node index for an edge.
     *
     * @param n1 First endpoint of the edge.
     * @param n2 Second endpoint of the edge.
     * @return The node index of the inserted midpoint.
     */
    idx lookup(idx n1, idx n2) const override;

private:
    std::unordered_map<MidpointEdge, idx, detail::MidpointEdgeHash> midpoints_;
};

} // namespace qlc3d::refinement

#endif // QLC3D_REFINEMENT_MIDPOINT_NODE_LOOKUP_H

