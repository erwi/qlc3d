#include <refinfo.h>
#include <refinement.h>
#include <algorithm>
#include <globals.h>

double get_elem_maxdQ(const idx elem,        // index to element
                      const Geometry& geom,     // currrent geometry
                      const SolutionVector& q)  // Q-tensor for gurrent geometry
{
// RETURNS ABSOLUTE MAXIMUM CHANGE IN ANY OF THE 5 Q-TENSOR COMPONENTS
// WITHIN ELEMENT elem


    double qe[4] = {0,0,0,0};
    double maxdq = 0;
    for (idx dim = 0 ; dim < 5 ; dim ++) // for q1 -> q5
    {
        for ( idx j = 0 ; j < geom.t->getnNodes() ; j++) // for each node in this tet
        {
            idx nn = geom.t->getNode(elem, j); // node number
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

    for ( idx i = 0 ; i <  geom.t->getnElements() ; i++)// for each tetrahedron
    {
        // IF NON-LC TETRAHEDRON: DON'T DO ANYTHING
        if (geom.t->getMaterialNumber(i) >= MAT_DIELECTRIC1 )
            continue;

        double maxdq = get_elem_maxdQ( i , geom, q);
        if ( maxdq >= thq ) // MARK AS RED IF....
        {
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

    // MAKE P TO TET INDEX
    std::vector< std::set<idx> > p_to_t;
    geom.t->gen_p_to_elem(p_to_t);
    std::set<idx> :: iterator itr;
    for (idx i = 0 ; i < (idx) indp.size() ; i++)
    {
        idx node = indp[i];
        for (itr = p_to_t[node].begin() ; itr != p_to_t[node].end() ; ++itr)
        {
            idx mat = geom.t->getMaterialNumber( *itr );
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
    double centre[3] = {0,0,0};
    refinfo.getCoord( centre[0], centre[1], centre[2] );

    // MAKE INDEX OF ALL POINTS THAT ARE SUFFICIENTLY CLOSE
    std::vector<idx> p_close;
    for (idx i = 0 ; i < geom.getnpLC() ; i++)
    {
        double distsqr = geom.getAbsDistSqr( i , centre );

        if ( distsqr <= rad )
            p_close.push_back( i );
    }

    if(p_close.empty() ) // NO NODES SELECTED. CAN RETURN
        return;

    selectTetsByCoordIndex(p_close, geom, i_tet);
}

void findTets_Box(const RefInfo& refinfo,
                  vector<idx>& i_tet,
                  const int refiter,
                  const Geometry& geom)
{
    // FINDS INDEX TO ALL TETS THAT ARE WITHIN A BOX DEFINED FOR THIS REFITER

    // GET CORRECT COORDINATE VALUES FOR THIS REFITER
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    refinfo.getCoords(x,y,z);

    // BOX LIMITS [min,max]
    double xLim[2] = {x[refiter*2], x[refiter*2 +1]};
    double yLim[2] = {y[refiter*2], y[refiter*2 +1]};
    double zLim[2] = {z[refiter*2], z[refiter*2 +1]};

    // FIND INDEX TO ALL NODES THAT ARE WITHIN BOX
    std::vector<idx> indp;
    for (idx i = 0; i < geom.getnpLC() ; i++)
    {
        double px = geom.getpX(i);
        double py = geom.getpY(i);
        double pz = geom.getpZ(i);

        if ( ((xLim[0]<=px) && (xLim[1]>=px)) &&
             ((yLim[0]<=py) && (yLim[1]>=py)) &&
             ((zLim[0]<=pz) && (zLim[1]>=pz)) )
        {
            indp.push_back(i);
        }
    }

    if ( indp.empty() )
        return;

    selectTetsByCoordIndex(indp, geom, i_tet); // MARK TETS THAT CONTAIN FOUND NODES

}
