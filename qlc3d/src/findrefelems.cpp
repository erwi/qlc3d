#include <refinement.h>
#include <algorithm>
#include <globals.h>
#include "geom/coordinates.h"

double get_elem_maxdQ(const idx elem,        // index to element
                      const Geometry& geom,     // currrent geometry
                      const SolutionVector& q)  // Q-tensor for gurrent geometry
{
// RETURNS ABSOLUTE MAXIMUM CHANGE IN ANY OF THE 5 Q-TENSOR COMPONENTS
// WITHIN ELEMENT elem


    double qe[4] = {0,0,0,0};
    double maxdq = 0;
    auto &t = geom.getTetrahedra();
    const idx numNodes = t.getnNodes();
    for (idx dim = 0 ; dim < 5 ; dim ++) {// for q1 -> q5
        for (idx j = 0; j < numNodes; j++){ // for each node in this tet
            idx nn = t.getNode(elem, j); // node number
            qe[j] = q.getValue( nn , dim ); // get g
        } // end for each node

        double mxq = *max_element(qe , qe+4);
        double mnq = *min_element(qe , qe+4);

        if ( (mxq-mnq) > maxdq ) {maxdq = mxq-mnq;}

    }// end for each dimension

    return maxdq;
}




void findTets_Change( const RefInfo& refinfo,
                      vector <idx>& i_tet,
                      const int refiter,
                      const Geometry& geom,
                      const SolutionVector& q
                      )
{
    // SELECTS ELEMENTS WHERE Q-TENSOR CHANGE IS ABOVE ALLOWED THRESHOLD,
    // AS SPECIFIED IN REFINFO OBJECT

    double thq = refinfo.getValue( refiter ); // GET THRESHOLD
    auto &t = geom.getTetrahedra();
    const idx numTets = t.getnElements();
    for (idx i = 0; i <  numTets; i++) {
        // IF NON-LC TETRAHEDRON: DON'T DO ANYTHING
        if (t.getMaterialNumber(i) >= MAT_DIELECTRIC1 )
            continue;

        double maxdq = get_elem_maxdQ( i , geom, q);
        if ( maxdq >= thq ) { // MARK AS RED IF....
            i_tet[i] = RED_TET;
        }
    }// end for
} // END findTets_Change


void selectTetsByCoordIndex(const std::vector<idx>& indp,
                            const Geometry& geom,
                            vector <idx>& i_tet
                            )
{
// SETS INDEX TO RED TETS TO "RED_TET" FOR THOSE TETRAHEDRAL
// ELEMENTS THAT CONTAIN NODES THAT ARE INDEXED IN IN indp
  auto &tets = geom.getTetrahedra();
    // MAKE P TO TET INDEX
    std::vector< std::set<idx> > p_to_t;
    tets.gen_p_to_elem(p_to_t);
    std::set<idx> :: iterator itr;
    for (idx i = 0; i < (idx) indp.size(); i++) {
        idx node = indp[i];
        for (itr = p_to_t[node].begin(); itr != p_to_t[node].end(); ++itr) {
            idx mat = tets.getMaterialNumber( *itr );
            if (mat <= MAT_DOMAIN7 )
                i_tet[(*itr)] = RED_TET;
        }
    }
}

void findTets_Sphere(const RefInfo& refinfo,
                     vector <idx>& i_tet,
                     const int refiter,
                     const Geometry& geom
                     )
{

    double rad = refinfo.getValue( refiter ); // SPHERE RADIUS
    rad*=rad;                                 // ARDIUS SQUARED
    const Vec3 centre(refinfo.getX()[0], refinfo.getY()[0], refinfo.getZ()[0]);

    // MAKE INDEX OF ALL POINTS THAT ARE SUFFICIENTLY CLOSE
    std::vector<idx> p_close;
    auto &coords = geom.getCoordinates();
    for (idx i = 0 ; i < geom.getnpLC() ; i++) {
      auto p = coords.getPoint(i);
      double distsqr = p.distanceSquared(centre);

      if (distsqr <= rad) {
        p_close.push_back(i);
      }
    }

    if(p_close.empty()) { // NO NODES SELECTED. CAN RETURN
      return;
    }

    selectTetsByCoordIndex(p_close, geom, i_tet);
}

void findTets_Box(const RefInfo& refinfo,
                  vector<idx>& i_tet,
                  const int refiter,
                  const Geometry& geom)
{
  // FINDS INDEX TO ALL TETS THAT ARE WITHIN A BOX DEFINED FOR THIS REFITER

  // GET CORRECT COORDINATE VALUES FOR THIS REFITER
  auto x = refinfo.getX();
  auto y = refinfo.getY();
  auto z = refinfo.getZ();

  // BOX LIMITS [min,max]
  double xLim[2] = {x[refiter*2], x[refiter*2 +1]};
  double yLim[2] = {y[refiter*2], y[refiter*2 +1]};
  double zLim[2] = {z[refiter*2], z[refiter*2 +1]};

  // FIND INDEX TO ALL NODES THAT ARE WITHIN BOX
  std::vector<idx> indp;
  auto &coords = geom.getCoordinates();
  for (idx i = 0; i < geom.getnpLC(); i++) {
    auto &p = coords.getPoint(i);

    if (((xLim[0] <= p.x()) && (xLim[1] >= p.x())) &&
        ((yLim[0] <= p.y()) && (yLim[1] >= p.y())) &&
        ((zLim[0] <= p.z()) && (zLim[1] >= p.z()))) {
      indp.push_back(i);
    }
  }

  if ( indp.empty() )
    return;

  selectTetsByCoordIndex(indp, geom, i_tet); // MARK TETS THAT CONTAIN FOUND NODES

}
