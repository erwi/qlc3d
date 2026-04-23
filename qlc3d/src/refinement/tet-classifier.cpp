#include <refinement/tet-classifier.h>
#include <refinement/periodic-edge-expander.h>
#include <refinement.h>
#include <geometry.h>
#include <geom/periodicity.h>
#include <util/exception.h>
#include <fmt/format.h>
#include <algorithm>
#include <set>

namespace qlc3d::refinement {

std::vector<Line> collectRedTetEdges(const Mesh& tets, const std::vector<idx>& i_tet) {
    std::vector<Line> lines;
    for (idx i = 0; i < (idx) i_tet.size(); i++) {
        if (i_tet[i] == RED_TET) {
            idx n[4] = { tets.getNode(i, 0), tets.getNode(i, 1),
                         tets.getNode(i, 2), tets.getNode(i, 3) };
            lines.push_back(Line(n[0], n[1]));
            lines.push_back(Line(n[0], n[2]));
            lines.push_back(Line(n[0], n[3]));
            lines.push_back(Line(n[1], n[2]));
            lines.push_back(Line(n[1], n[3]));
            lines.push_back(Line(n[2], n[3]));
        }
    }
    std::sort(lines.begin(), lines.end());
    auto u = std::unique(lines.begin(), lines.end());
    lines.resize(u - lines.begin());
    return lines;
}

void countEdgesPerElement(const Mesh& elements,
                          const std::vector<std::set<idx>>& p_to_elem,
                          const std::vector<Line>& lines,
                          std::vector<idx>& i_elem,
                          std::vector<std::set<idx>>& m_to_l) {
    const idx nElem = elements.getnElements();
    i_elem.assign(nElem, 0);
    m_to_l.assign(nElem, std::set<idx>{});

    const int numPoints = (int) p_to_elem.size();
    for (idx l = 0; l < (idx) lines.size(); l++) {
        int n1 = lines[l].L[0];
        int n2 = lines[l].L[1];
        if (n1 < numPoints && n2 < numPoints) {
            std::vector<idx> shared;
            std::set_intersection(p_to_elem[n1].begin(), p_to_elem[n1].end(),
                                  p_to_elem[n2].begin(), p_to_elem[n2].end(),
                                  std::back_inserter(shared));
            for (idx e : shared) {
                i_elem[e]++;
                m_to_l[e].insert(l);
            }
        }
    }
}

Num_Ref_Tet assignRefinementTypes(std::vector<idx>& i_tet) {
    Num_Ref_Tet nrt;
    for (idx& val : i_tet) {
        if (val == 0) continue;
        else if (val == GREEN1_TET) nrt.green1++;
        else if (val == GREEN2_TET) nrt.green2++;
        else if (val == GREEN3_TET) nrt.green3++;
        else if (val > GREEN3_TET) {
            val = RED_TET;
            nrt.red++;
        }
    }
    return nrt;
}

void resolveGreen3RedAmbiguity(std::vector<idx>& i_tet,
                               const std::vector<Line>& lines,
                               const std::vector<std::set<idx>>& t_to_l) {
    for (idx i = 0; i < (idx) i_tet.size(); i++) {
        if (i_tet[i] == GREEN3_TET) {
            std::set<idx> nodes;
            for (idx li : t_to_l[i]) {
                nodes.insert(lines[li].L[0]);
                nodes.insert(lines[li].L[1]);
            }
            if (nodes.size() == 4) {
                i_tet[i] = RED_TET;
            } else if (nodes.size() != 3) {
                RUNTIME_ERROR(fmt::format("Expected 3 nodes, got {}.", nodes.size()));
            }
        }
    }
}

ClassificationResult classifyRefinement(Geometry& geom, std::vector<idx> i_tet_in) {
    ClassificationResult result;
    result.i_tet = std::move(i_tet_in);

    auto& tets = geom.getTetrahedra();
    std::vector<std::set<idx>> p_to_t;
    tets.gen_p_to_elem(p_to_t);

    PeriodicEdgeExpander expander(geom);

    Num_Ref_Tet nrt = assignRefinementTypes(result.i_tet);

    if (nrt.red == 0) {
        // Nothing to do – set up empty output
        result.i_tri.assign(geom.getTriangles().getnElements(), 0);
        return result;
    }

    unsigned int prevRed = 0;
    do {
        prevRed = nrt.red;
        result.lines = collectRedTetEdges(tets, result.i_tet);
        expander.expand(result.lines);
        countEdgesPerElement(tets, p_to_t, result.lines, result.i_tet, result.t_to_l);
        resolveGreen3RedAmbiguity(result.i_tet, result.lines, result.t_to_l);
        nrt = assignRefinementTypes(result.i_tet);
    } while (nrt.red > prevRed);

    // Triangle classification
    std::vector<std::set<idx>> p_to_e;
    geom.getTriangles().gen_p_to_elem(p_to_e);
    result.i_tri.assign(geom.getTriangles().getnElements(), 0);
    countEdgesPerElement(geom.getTriangles(), p_to_e,
                         result.lines, result.i_tri, result.e_to_l);
    return result;
}

} // namespace qlc3d::refinement


