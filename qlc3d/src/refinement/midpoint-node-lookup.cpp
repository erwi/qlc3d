#include <refinement/midpoint-node-lookup.h>
#include <stdexcept>

namespace qlc3d::refinement {

std::size_t detail::MidpointEdgeHash::operator()(const MidpointEdge &edge) const noexcept {
	const auto h1 = std::hash<idx>{}(edge.first);
	const auto h2 = std::hash<idx>{}(edge.second);
	return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6U) + (h1 >> 2U));
}


void MapMidpointNodeLookup::registerEdge(idx n1, idx n2, idx midpointNode) {
	midpoints_[detail::canonicalEdge(n1, n2)] = midpointNode;
}

idx MapMidpointNodeLookup::lookup(idx n1, idx n2) const {
	const auto edge = detail::canonicalEdge(n1, n2);
	const auto it = midpoints_.find(edge);
	if (it == midpoints_.end()) {
		throw std::out_of_range("No midpoint node registered for edge.");
	}
	return it->second;
}

} // namespace qlc3d::refinement

