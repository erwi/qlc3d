#include <meshrefinement.h>
#include <sstream>

void RefinementConfig::validate() {
    // TODO:
    // check that Type is known
    // check that correct number of values are given
    // etc.
}

MeshRefinement::MeshRefinement() {
    refiter = 0;
    needs_new_mesh = true;
}

void MeshRefinement::setRefinementConfig(std::vector<RefinementConfig> &&ref) {
    refinementConfigs_ = std::move(ref);
}
