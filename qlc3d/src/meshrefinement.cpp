#include <meshrefinement.h>
#include <refinement/refinement-spec.h>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <util/logging.h>

//<editor-fold RefinementConfig>
void RefinementConfig::validate() {
    // Normalize type to lowercase for comparison
    std::string typeLower = type_;
    std::transform(typeLower.begin(), typeLower.end(), typeLower.begin(), ::tolower);

    if (typeLower == "change") {
        if (values_.empty()) {
            throw std::invalid_argument("Refinement type 'Change' requires non-empty values.");
        }
    } else if (typeLower == "sphere") {
        if (values_.empty() || x_.empty() || y_.empty() || z_.empty()) {
            throw std::invalid_argument("Refinement type 'Sphere' requires non-empty values and x, y, z coordinates.");
        }
    } else if (typeLower == "box") {
        size_t nx = x_.size(), ny = y_.size(), nz = z_.size();
        if (nx == 0 || ny == 0 || nz == 0) {
            throw std::invalid_argument("Refinement type 'Box' requires non-empty x, y, z coordinates.");
        }
        if (nx != ny || nx != nz) {
            throw std::invalid_argument("Refinement type 'Box' requires equal numbers of x, y, z coordinates.");
        }
        if (nx % 2 != 0) {
            throw std::invalid_argument("Refinement type 'Box' requires an even number of x, y, z coordinates (pairs per box).");
        }
    } else {
        throw std::invalid_argument(fmt::format("Unknown refinement type '{}'.", type_));
    }
}

void MeshRefinement::setRefinementConfig(std::vector<RefinementConfig> &&ref) {
    refinementConfigs_ = std::move(ref);
}

std::vector<std::unique_ptr<RefinementSpec>> RefinementConfig::toSpecs() const {
    std::vector<std::unique_ptr<RefinementSpec>> specs;
    if (occursPeriodically()) {
        specs.push_back(RefinementSpec::makePeriodic(type_, values_, x_, y_, z_));
    } else {
        for (int iter : iterations_) {
            specs.push_back(RefinementSpec::makeExplicit(type_, static_cast<long>(iter), -1.0, values_, x_, y_, z_));
        }
        for (double time : times_) {
            specs.push_back(RefinementSpec::makeExplicit(type_, 0, time, values_, x_, y_, z_));
        }
    }
    return specs;
}
//</editor-fold RefinementConfig>
