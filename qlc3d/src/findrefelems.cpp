#include <refinfo.h>
#include <refinement.h>
#include <algorithm>
#include <globals.h>

double get_elem_maxdQ(const size_t elem,
                      const Geometry& geom,
                      const SolutionVector& q)
{
// RETURNS ABSOLUTE MAXIMUM CHANGE IN ANY OF THE 5 Q-TENSOR COMPONENTS
// WITHIN ELEMENT elem


    double qe[4] = {0,0,0,0};
    double maxdq = 0;
    for (unsigned int dim = 0 ; dim < 5 ; dim ++) // for q1 -> q5
    {
        for ( int j = 0 ; j < geom.t->getnNodes() ; j++) // for each node in this tet
        {
            int nn = geom.t->getNode(elem, j); // node number
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

    for ( idx i = 0 ; i < (idx) geom.t->getnElements() ; i++)
    {
        double maxdq = get_elem_maxdQ( i , geom, q);

        if ( maxdq >= thq ) // MARK AS RED IF....
        {
            i_tet[i] = RED_TET;
        }
    }
} // END findTets_Change

void findTets_Sphere(const RefInfo& refinfo,
                     vector <idx>& i_tet,
                     const int refiter,
                     const Geometry& geom
                     )
{

    double rad = refinfo.getValue( refiter ); // SPHERE RADIUS
    rad*=rad;                                 // ARDIUS SQUARED
    double centre[3] = { 0,0,0};
    refinfo.getCoord( centre[0], centre[1], centre[2] );

    // MAKE INDEX OF ALL POINTS THAT ARE SUFFICIENTLY CLOSE
    std::vector<size_t> p_close;
    for (size_t i = 0 ; i < geom.getnp() ; i++)
    {
        double distsqr = geom.getAbsDistSqr( i , centre );

        if ( distsqr <= rad )
            p_close.push_back( i );
    }

    // MAKE p_close INDEX UNIQUE
    std::sort( p_close.begin() , p_close.end() );
    std::vector<size_t>::iterator uitr = unique( p_close.begin(), p_close.end() );
    p_close.erase( uitr , p_close.end() );

    // LOOP OVER EACH ELEMENT
    for (size_t i = 0 ; i < geom.t->getnElements() ; i++)
    {
        int* tt = geom.t->getPtrToElement( i );         // GET SHORTCUT TO ELEMENT NODE INDEX

        for ( size_t j = 0 ; j < p_close.size() ; j++)    // LOOP OVER ALL SELECTED NODES
        {
            if ( ( tt[0] == p_close[j] ) || // IF ELEMENT CONTAINS NODE j
                 ( tt[1] == p_close[j] ) ||
                 ( tt[2] == p_close[j] ) ||
                 ( tt[3] == p_close[j] ) )
            {
                i_tet[i] = RED_TET;         // MARK ELEMENT AS RED
                break;
            }
        }
    }

}
