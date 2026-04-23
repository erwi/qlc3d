#ifndef PERIODIC_EDGE_EXPANDER_H
#define PERIODIC_EDGE_EXPANDER_H

#include <line.h>
#include <geometry.h>
#include <vector>
#include <memory>

/**
 * @brief Catalogues the periodic boundary edges of a mesh and can expand
 *        a set of bisected edges with their periodic mirror counterparts.
 *
 * Constructed once from a Geometry.  Call expand() one or more times to
 * add mirror edges to a bisected-edge list.
 */
class PeriodicEdgeExpander {
public:
    /**
     * @brief Constructs the expander and catalogues all boundary edges.
     * @param geom  The geometry whose periodic boundary triangles are scanned.
     */
    explicit PeriodicEdgeExpander(Geometry& geom);

    ~PeriodicEdgeExpander();

    /**
     * @brief Appends mirror edges for every edge in @p lines that lies on a
     *        periodic boundary face or corner.
     *        Mirror edges are appended; duplicates are removed before returning.
     * @param lines  In/out vector of bisected edges.
     */
    void expand(std::vector<Line>& lines) const;

    /** @return Edges lying on the front (min-Y) periodic face. */
    const std::vector<Line>& frontEdges()  const;
    /** @return Edges lying on the back (max-Y) periodic face. */
    const std::vector<Line>& backEdges()   const;
    /** @return Edges lying on the left (min-X) periodic face. */
    const std::vector<Line>& leftEdges()   const;
    /** @return Edges lying on the right (max-X) periodic face. */
    const std::vector<Line>& rightEdges()  const;
    /** @return Edges lying on the top (max-Z) periodic face. */
    const std::vector<Line>& topEdges()    const;
    /** @return Edges lying on the bottom (min-Z) periodic face. */
    const std::vector<Line>& bottomEdges() const;

private:
    struct Catalogue;
    std::unique_ptr<Catalogue> cat_;
    Geometry& geom_;

    /**
     * @brief Finds the translational mirror of @p edge in @p candidates.
     * @param edge        The edge to match.
     * @param candidates  Pool of candidate mirror edges.
     * @param geom        Geometry used for coordinate lookups.
     * @param dir         3-element mask: 1 = coordinate must match, 0 = may differ.
     * @return Pointer to the matching edge, or nullptr if not found.
     */
    static const Line* findTranslation(const Line& edge,
                                       const std::vector<Line>& candidates,
                                       Geometry& geom,
                                       const double dir[3]);
};

#endif // PERIODIC_EDGE_EXPANDER_H

