#ifndef REFINEMENT_SPEC_H
#define REFINEMENT_SPEC_H

#include <memory>
#include <string>
#include <vector>

/**
 * RefinementSpec is the unified descriptor for a mesh refinement region, replacing both
 * RefinementConfig (settings-file representation) and RefInfo (runtime representation).
 *
 * Instances are created via the static factory methods makeExplicit() or makePeriodic().
 */
class RefinementSpec {
public:
    enum class Type { Change, Sphere, Box };

    /**
     * Constructs a spec for a refinement event that fires at a specific iteration or time.
     * Exactly one of @p iteration or @p time must be positive; the other must be ≤ 0.
     *
     * @param type     Refinement type string: "Change", "Sphere", or "Box" (case-insensitive).
     * @param iteration Fire at this iteration number (> 0), or ≤ 0 if time-based.
     * @param time      Fire at this simulation time (> 0), or ≤ 0 if iteration-based.
     * @param values    Refinement threshold values.
     * @param x         X-coordinates of refinement region.
     * @param y         Y-coordinates of refinement region.
     * @param z         Z-coordinates of refinement region.
     * @return A unique_ptr owning the new RefinementSpec.
     * @throws std::invalid_argument if arguments are inconsistent or validation fails.
     */
    static std::unique_ptr<RefinementSpec> makeExplicit(
        const std::string &type,
        long iteration,
        double time,
        std::vector<double> values,
        std::vector<double> x,
        std::vector<double> y,
        std::vector<double> z);

    /**
     * Constructs a spec for a periodically re-occurring refinement.
     *
     * @param type   Refinement type string: "Change", "Sphere", or "Box" (case-insensitive).
     * @param values Refinement threshold values.
     * @param x      X-coordinates of refinement region.
     * @param y      Y-coordinates of refinement region.
     * @param z      Z-coordinates of refinement region.
     * @return A unique_ptr owning the new RefinementSpec.
     * @throws std::invalid_argument if validation fails.
     */
    static std::unique_ptr<RefinementSpec> makePeriodic(
        const std::string &type,
        std::vector<double> values,
        std::vector<double> x,
        std::vector<double> y,
        std::vector<double> z);

    /** @return The type of refinement region. */
    [[nodiscard]] Type getType() const { return type_; }

    /** @return true if this spec fires periodically (not tied to a specific iteration or time). */
    [[nodiscard]] bool isPeriodic() const { return periodic_; }

    /**
     * @return The iteration at which this event fires, or 0 if periodic or time-based.
     */
    [[nodiscard]] long getIteration() const { return iteration_; }

    /**
     * @return The simulation time at which this event fires, or 0.0 if periodic or iteration-based.
     */
    [[nodiscard]] double getTime() const { return time_; }

    /**
     * @return Number of refinement sub-iterations this spec describes.
     *         For Change/Sphere: values.size(). For Box: x.size() / 2.
     */
    [[nodiscard]] unsigned int getRefIter() const { return refIter_; }

    /**
     * @param i Index into the values array.
     * @return Value at index i.
     */
    [[nodiscard]] double getValue(size_t i) const;

    /** @return X-coordinates of refinement region. */
    [[nodiscard]] const std::vector<double>& getX() const { return x_; }

    /** @return Y-coordinates of refinement region. */
    [[nodiscard]] const std::vector<double>& getY() const { return y_; }

    /** @return Z-coordinates of refinement region. */
    [[nodiscard]] const std::vector<double>& getZ() const { return z_; }

    /**
     * Create an independent copy of this spec.
     * @return A new unique_ptr owning a deep copy.
     */
    [[nodiscard]] std::unique_ptr<RefinementSpec> clone() const;

private:
    RefinementSpec() = default;

    static Type parseType(const std::string &type);
    static unsigned int calcRefIter(Type type, const std::vector<double> &values, const std::vector<double> &x);
    static void validateFields(Type type,
                                const std::vector<double> &values,
                                const std::vector<double> &x,
                                const std::vector<double> &y,
                                const std::vector<double> &z);

    Type type_{};
    bool periodic_{false};
    long iteration_{0};
    double time_{0.0};
    unsigned int refIter_{0};
    std::vector<double> values_;
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> z_;
};

#endif // REFINEMENT_SPEC_H


