#ifndef MESHREFINEMENT_H
#define MESHREFINEMENT_H
#include <list>
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <material_numbers.h>
#include <refinement/refinement-spec.h>

/**
 * RefinementConfig represents the data in a settings file for a single REFINEMENT object. These
 * are later translated into refinement events.
 */
class RefinementConfig {
    /** throws exception if an invalid RefinementConfig is detected */
    void validate();
public:
    const std::string type_;
    const std::vector<int> iterations_;
    const std::vector<double> times_;
    const std::vector<double> values_;
    const std::vector<double> x_;
    const std::vector<double> y_;
    const std::vector<double> z_;

    RefinementConfig(std::string type, std::vector<int> iterations, std::vector<double> times,
                     std::vector<double> values,
                     std::vector<double> x, std::vector<double> y, std::vector<double> z) :
            type_ { type },
            iterations_ { iterations },
            times_ { times },
            values_ { values },
            x_ { x }, y_ { y }, z_ { z } {
        validate();
    }

    /**
     * if no iteration or time is defined, then this is assumed to occur periodically, with a period defined as
     * repRefIter and/or repRefTime
     */
    bool occursPeriodically() const { return iterations_.empty() && times_.empty(); }

    /**
     * Convert this RefinementConfig to a RefinementSpec, expanding each explicit
     * iteration or time into a separate spec.  If periodic, creates a single periodic spec.
     * @return vector of unique_ptr<RefinementSpec>
     */
    [[nodiscard]] std::vector<std::unique_ptr<RefinementSpec>> toSpecs() const;
};

class MeshRefinement {
    unsigned int repRefIter_ = 0;    // iteration period for repeating mesh refinement
    std::vector<RefinementConfig> refinementConfigs_;

public:
    void setRepRefIter(unsigned int iter) { repRefIter_ = iter; }
    [[nodiscard]] unsigned int getRepRefIter() const { return repRefIter_; }


    void setRefinementConfig(std::vector<RefinementConfig> &&ref);
    [[nodiscard]] const std::vector<RefinementConfig>& getRefinementConfig() const { return refinementConfigs_; }
};
#endif
