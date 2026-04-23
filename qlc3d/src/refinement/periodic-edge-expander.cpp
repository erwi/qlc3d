/**
 * Implementation of PeriodicEdgeExpander.
 *
 * Strategy:
 *  1. Catalogue  – The constructor scans all MAT_PERIODIC triangles and
 *     classifies their edges into per-face and per-corner bucket vectors
 *     (front/back, left/right, top/bottom, and 12 corner lines).
 *
 *  2. Match  – findTranslation() is a static predicate that iterates a
 *     candidate list and returns the first edge that is a translational
 *     copy of the query edge (i.e. the same shape, shifted in the
 *     coordinate direction(s) indicated by the dir[] mask).
 *
 *  3. Expand – expand() walks every bisected edge in the input list,
 *     detects which periodic face or corner it belongs to, calls
 *     findTranslation() against the opposite bucket, and appends the
 *     mirror edge.  Duplicates are removed before returning.
 *
 * Corner-line special cases
 * -------------------------
 * For left/right periodicity the four vertical corner lines (c0..c3,
 * along Z) each need 3 mirrors – one in each of the other 3 corners.
 * For fully periodic (top/bottom also active) the 8 horizontal corner
 * lines (ca..cd along X, cA..cD along Y) are handled analogously, but
 * top-bottom periodicity currently throws RUNTIME_ERROR because it has
 * not been fully implemented in the refinement pipeline.
 */

#include <refinement/periodic-edge-expander.h>
#include <material_numbers.h>
#include <geom/periodicity.h>
#include <util/exception.h>
#include <line.h>
#include <geometry.h>

// ---- Catalogue ----------------------------------------------------------------

struct PeriodicEdgeExpander::Catalogue {
    // face edges
    std::vector<Line> front, back, left, right, top, bottom;
    // vertical corners (along Z): c0=xmin/ymin, c1=xmax/ymin, c2=xmax/ymax, c3=xmin/ymax
    std::vector<Line> c0, c1, c2, c3;
    // horizontal corners along X: ca=ymin/zmin, cb=ymax/zmin, cc=ymax/zmax, cd=ymin/zmax
    std::vector<Line> ca, cb, cc, cd;
    // horizontal corners along Y: cA=xmin/zmin, cB=xmax/zmin, cC=xmax/zmax, cD=xmin/zmax
    std::vector<Line> cA, cB, cC, cD;
};

// ---- Constructor (Catalogue step) --------------------------------------------

PeriodicEdgeExpander::PeriodicEdgeExpander(Geometry& geom)
    : cat_(std::make_unique<Catalogue>()), geom_(geom)
{
    PeriodicityType periodicType(geom.getTriangles());
    bool peritype[3] = {
        periodicType.isFrontBackPeriodic(),
        periodicType.isLeftRightPeriodic(),
        periodicType.isTopBottomPeriodic()
    };

    for (idx i = 0; i < geom.getTriangles().getnElements(); i++) {
        if (geom.getTriangles().getMaterialNumber(i) != MAT_PERIODIC) continue;

        idx n[3];
        geom.getTriangles().loadNodes(i, n);
        Line lines[3] = {Line(n[0], n[1]), Line(n[0], n[2]), Line(n[1], n[2])};

        for (idx l = 0; l < 3; l++) {
            // CASE 1: front/back only
            if (peritype[0] && !peritype[1] && !peritype[2]) {
                if      (lines[l].isOnFrontSurface(&geom)) cat_->front.push_back(lines[l]);
                else if (lines[l].isOnBackSurface(&geom))  cat_->back.push_back(lines[l]);
            }
            // CASE 2: front/back + left/right
            else if (peritype[0] && peritype[1] && !peritype[2]) {
                if      (lines[l].isOnFrontSurface(&geom)) cat_->front.push_back(lines[l]);
                else if (lines[l].isOnBackSurface(&geom))  cat_->back.push_back(lines[l]);

                if      (lines[l].isOnLeftSurface(&geom))  cat_->left.push_back(lines[l]);
                else if (lines[l].isOnRightSurface(&geom)) cat_->right.push_back(lines[l]);

                // vertical corner classification
                if      (lines[l].isCorn0(&geom)) cat_->c0.push_back(lines[l]);
                else if (lines[l].isCorn1(&geom)) cat_->c1.push_back(lines[l]);
                else if (lines[l].isCorn2(&geom)) cat_->c2.push_back(lines[l]);
                else if (lines[l].isCorn3(&geom)) cat_->c3.push_back(lines[l]);
            }
            // CASE 3: fully periodic (front/back + left/right + top/bottom)
            else if (peritype[0] && peritype[1] && peritype[2]) {
                if      (lines[l].isOnFrontSurface(&geom))  cat_->front.push_back(lines[l]);
                else if (lines[l].isOnBackSurface(&geom))   cat_->back.push_back(lines[l]);

                if      (lines[l].isOnLeftSurface(&geom))   cat_->left.push_back(lines[l]);
                else if (lines[l].isOnRightSurface(&geom))  cat_->right.push_back(lines[l]);

                if      (lines[l].isOnTopSurface(&geom))    cat_->top.push_back(lines[l]);
                else if (lines[l].isOnBottomSurface(&geom)) cat_->bottom.push_back(lines[l]);

                if      (lines[l].isCorn0(&geom)) cat_->c0.push_back(lines[l]);
                else if (lines[l].isCorn1(&geom)) cat_->c1.push_back(lines[l]);
                else if (lines[l].isCorn2(&geom)) cat_->c2.push_back(lines[l]);
                else if (lines[l].isCorn3(&geom)) cat_->c3.push_back(lines[l]);

                if      (lines[l].isCorna(&geom)) cat_->ca.push_back(lines[l]);
                else if (lines[l].isCornb(&geom)) cat_->cb.push_back(lines[l]);
                else if (lines[l].isCornc(&geom)) cat_->cc.push_back(lines[l]);
                else if (lines[l].isCornd(&geom)) cat_->cd.push_back(lines[l]);

                if      (lines[l].isCornA(&geom)) cat_->cA.push_back(lines[l]);
                else if (lines[l].isCornB(&geom)) cat_->cB.push_back(lines[l]);
                else if (lines[l].isCornC(&geom)) cat_->cC.push_back(lines[l]);
                else if (lines[l].isCornD(&geom)) cat_->cD.push_back(lines[l]);
            }
        }
    }

    // deduplicate all buckets
    uniquefy_line_vector(cat_->front);
    uniquefy_line_vector(cat_->back);
    uniquefy_line_vector(cat_->left);
    uniquefy_line_vector(cat_->right);
    uniquefy_line_vector(cat_->top);
    uniquefy_line_vector(cat_->bottom);
    uniquefy_line_vector(cat_->c0);
    uniquefy_line_vector(cat_->c1);
    uniquefy_line_vector(cat_->c2);
    uniquefy_line_vector(cat_->c3);
    uniquefy_line_vector(cat_->ca);
    uniquefy_line_vector(cat_->cb);
    uniquefy_line_vector(cat_->cc);
    uniquefy_line_vector(cat_->cd);
    uniquefy_line_vector(cat_->cA);
    uniquefy_line_vector(cat_->cB);
    uniquefy_line_vector(cat_->cC);
    uniquefy_line_vector(cat_->cD);
}

PeriodicEdgeExpander::~PeriodicEdgeExpander() = default;

// ---- Accessors ---------------------------------------------------------------

const std::vector<Line>& PeriodicEdgeExpander::frontEdges()  const { return cat_->front; }
const std::vector<Line>& PeriodicEdgeExpander::backEdges()   const { return cat_->back; }
const std::vector<Line>& PeriodicEdgeExpander::leftEdges()   const { return cat_->left; }
const std::vector<Line>& PeriodicEdgeExpander::rightEdges()  const { return cat_->right; }
const std::vector<Line>& PeriodicEdgeExpander::topEdges()    const { return cat_->top; }
const std::vector<Line>& PeriodicEdgeExpander::bottomEdges() const { return cat_->bottom; }

// ---- findTranslation ---------------------------------------------------------

const Line* PeriodicEdgeExpander::findTranslation(const Line& edge,
                                                   const std::vector<Line>& candidates,
                                                   Geometry& geom,
                                                   const double dir[3])
{
    // We need a non-const copy to call the non-const Line methods.
    Line mutableEdge = edge;
    double mutableDir[3] = {dir[0], dir[1], dir[2]};
    for (const Line& c : candidates) {
        Line mutableC = c;
        if (mutableEdge.isTranslationOf(mutableC, &geom, mutableDir)) {
            return &c;
        }
    }
    return nullptr;
}

// ---- expand ------------------------------------------------------------------

void PeriodicEdgeExpander::expand(std::vector<Line>& lines) const {
    PeriodicityType periodicityType(geom_.getTriangles());
    if (!periodicityType.isAnyPeriodic()) return;

    std::vector<Line> newlines;

    for (Line& edge : lines) {
        // --- Front / Back ---
        if (periodicityType.isFrontBackPeriodic()) {
            bool found = true;
            double dir[3] = {1, 0, 1}; // X and Z must match; Y (depth) may shift
            if (edge.isOnFrontSurface(&geom_)) {
                found = false;
                const Line* mirror = findTranslation(edge, cat_->back, geom_, dir);
                if (mirror) { newlines.push_back(*mirror); found = true; }
            } else if (edge.isOnBackSurface(&geom_)) {
                found = false;
                const Line* mirror = findTranslation(edge, cat_->front, geom_, dir);
                if (mirror) { newlines.push_back(*mirror); found = true; }
            }
            if (!found) {
                RUNTIME_ERROR("Could not find periodic front/back face line.")
            }
        }

        // --- Left / Right ---
        if (periodicityType.isLeftRightPeriodic()) {
            bool found = true;
            double dir[3] = {0, 1, 1}; // Y and Z must match; X may shift
            if (edge.isOnLeftSurface(&geom_)) {
                found = false;
                const Line* mirror = findTranslation(edge, cat_->right, geom_, dir);
                if (mirror) { newlines.push_back(*mirror); found = true; }
            } else if (edge.isOnRightSurface(&geom_)) {
                found = false;
                const Line* mirror = findTranslation(edge, cat_->left, geom_, dir);
                if (mirror) { newlines.push_back(*mirror); found = true; }
            }
            if (!found) {
                RUNTIME_ERROR("Could not find periodic left/right face line.")
            }

            // Vertical corner lines need 3 mirrors (one in each other corner)
            if (edge.isTopBottomCornerLine(&geom_)) {
                int cc = 0;
                double cdir[3] = {0, 0, 1}; // Z must match; X,Y may shift
                if (const Line* m = findTranslation(edge, cat_->c0, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->c1, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->c2, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->c3, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (cc != 3) {
                    RUNTIME_ERROR("Problem finding periodic top-bottom corner lines.")
                }
            }
        }

        // --- Top / Bottom ---
        if (periodicityType.isTopBottomPeriodic()) {
            bool found = true;
            double dir[3] = {1, 1, 0}; // X and Y must match; Z may shift
            if (edge.isOnTopSurface(&geom_)) {
                found = false;
                const Line* mirror = findTranslation(edge, cat_->bottom, geom_, dir);
                if (mirror) { newlines.push_back(*mirror); found = true; }
            } else if (edge.isOnBottomSurface(&geom_)) {
                found = false;
                const Line* mirror = findTranslation(edge, cat_->top, geom_, dir);
                if (mirror) { newlines.push_back(*mirror); found = true; }
            }
            if (!found) {
                RUNTIME_ERROR("Could not find periodic top/bottom line.")
            }

            // Front/back corner lines (horizontal, along Y)
            if (edge.isFrontBackCornerLine(&geom_)) {
                int cc = 0;
                double cdir[3] = {0, 1, 0}; // Y must match; X,Z may shift
                if (const Line* m = findTranslation(edge, cat_->cA, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->cB, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->cC, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->cD, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (cc != 3) {
                    RUNTIME_ERROR("Problem finding periodic front-back corner lines.")
                }
            } else if (edge.isLeftRightCornerLine(&geom_)) {
                int cc = 0;
                double cdir[3] = {1, 0, 0}; // X must match; Y,Z may shift
                if (const Line* m = findTranslation(edge, cat_->ca, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->cb, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->cc, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (const Line* m = findTranslation(edge, cat_->cd, geom_, cdir)) { newlines.push_back(*m); cc++; }
                if (cc != 3) {
                    RUNTIME_ERROR("Problem finding periodic left-right corner lines.")
                }
            }
            RUNTIME_ERROR("Top-bottom periodicity not implemented in mesh refinement yet.");
        }
    }

    lines.insert(lines.end(), newlines.begin(), newlines.end());
    uniquefy_line_vector(lines);
}

