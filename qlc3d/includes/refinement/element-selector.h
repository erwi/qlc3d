#ifndef ELEMENT_SELECTOR_H
#define ELEMENT_SELECTOR_H

#include <globals.h>
#include <memory>
#include <vector>

class Geometry;
class SolutionVector;
class RefinementSpec;

namespace qlc3d::refinement {

/**
 * Abstract interface for selecting tetrahedra that require refinement.
 *
 * Implementations mark entries in @p i_tet as RED_TET (value 6) for
 * tetrahedra that meet the selector's criteria.
 *
 * Use the static factory methods (makeChange, makeSphere, makeBox, fromSpec)
 * to create instances rather than constructing concrete subclasses directly.
 */
class ElementSelector {
public:
    virtual ~ElementSelector() = default;

    /**
     * Mark tetrahedra in @p i_tet that satisfy this selector's criterion.
     *
     * @param geom     Current working geometry.
     * @param q        Q-tensor solution vector. May be nullptr for selectors
     *                 that do not require it (Sphere, Box). ChangeSelector
     *                 requires a non-null pointer.
     * @param refIter  Zero-based refinement sub-iteration index.
     * @param i_tet    In/out vector of per-tet refinement type codes. Selected
     *                 tets are set to RED_TET; others are left unchanged.
     */
    virtual void selectTets(const Geometry& geom,
                            const SolutionVector* q,
                            int refIter,
                            std::vector<idx>& i_tet) const = 0;

    /**
     * @return true if this selector needs a valid Q-tensor solution vector.
     */
    [[nodiscard]] virtual bool needsQTensor() const { return false; }

    /**
     * @return The number of refinement sub-iterations this selector covers.
     */
    [[nodiscard]] virtual unsigned int getRefIter() const = 0;

    /**
     * Create a ChangeSelector that marks tets where the max change across any
     * Q-tensor component exceeds a threshold.
     *
     * @param spec  RefinementSpec of type Change.
     */
    static std::unique_ptr<ElementSelector> makeChange(const RefinementSpec& spec);

    /**
     * Create a SphereSelector that marks tets containing any LC node inside a
     * sphere defined by @p spec.
     *
     * @param spec  RefinementSpec of type Sphere.
     */
    static std::unique_ptr<ElementSelector> makeSphere(const RefinementSpec& spec);

    /**
     * Create a BoxSelector that marks tets containing any LC node inside an
     * axis-aligned bounding box defined by @p spec.
     *
     * @param spec  RefinementSpec of type Box.
     */
    static std::unique_ptr<ElementSelector> makeBox(const RefinementSpec& spec);

    /**
     * Create an ElementSelector of the appropriate concrete type based on
     * @p spec's type.
     *
     * @param spec  RefinementSpec whose type determines the selector.
     * @throws std::invalid_argument if the type is unrecognised.
     */
    static std::unique_ptr<ElementSelector> fromSpec(const RefinementSpec& spec);
};

} // namespace qlc3d::refinement

#endif // ELEMENT_SELECTOR_H

