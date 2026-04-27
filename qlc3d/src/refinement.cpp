
#include <refinement.h>
#include <refinement/tet-classifier.h>
#include <material_numbers.h>
#include <meshrefinement.h>
#include <line.h>
#include <set>
#include <globals.h>
#include <geom/coordinates.h>
#include <util/exception.h>
#include <util/logging.h>
#include "geom/periodicity.h"

using std::set;


void modify_geometry(Geometry& geom,
                     vector <idx>& i_tet, vector <idx>& i_tri,
                     vector <Vec3>& new_p,
                     vector <idx>& new_t, vector <idx>& new_e,
                     vector <idx> new_mat_t ,
                     vector <idx> new_mat_e
                     ){

    if ( i_tet.empty() ) return;

    //geom.appendCoordinates( new_p ); // this is done asap after new coordinate calculation so that new tet node ordering can be calculated to ensure positive jdet
    // MAKE LIST OF TETS AND TRIS TO BE REMOVED. THESE ARE REPLACED BY THE NEWLY
    // CREATED ONES
    set <idx> ind_remove_tets;
    for (idx i = 0 ; i < (idx) i_tet.size() ; i++)
    {
        if ( i_tet[i] >0 )
        {
            ind_remove_tets.insert(i);
        }
    }
    set <idx> ind_remove_tris;
    for( idx i = 0 ; i < (idx) i_tri.size(); i++ )
    {
        if ( i_tri[i] )
        {
            ind_remove_tris.insert(i);
        }
    }

    auto &tets = geom.getTetrahedra();
    auto &tris = geom.getTriangles();
    tets.removeElements( ind_remove_tets );
    tets.appendElements(new_t, new_mat_t);
    tris.removeElements( ind_remove_tris );
    tris.appendElements(new_e, new_mat_e);

    geom.ReorderDielectricNodes(); // Dielectric nodes are moved last
    //geom.t->setMaxNodeNumber( (unsigned int) geom.getnp() );
    tris.setConnectedVolume(&tets);
    tets.calculateDeterminants3D(geom.getCoordinates());
    tets.ScaleDeterminants(qlc3d::units::CUBIC_MICROMETER_TO_CUBIC_METER);
    tris.calculateSurfaceNormals(geom.getCoordinates(), &tets);
    tris.ScaleDeterminants(qlc3d::units::SQUARE_MICROMETER_TO_SQUARE_METER);
 }

void Refine(Geometry& geom,                 // SOURCE (OLD) GEOMETRY
            vector <idx> & i_tet)        // REFINEMENT TYPES VECTOR
{
    using namespace qlc3d::refinement;

    Num_Ref_Tet nrt = assignRefinementTypes(i_tet);
    if (nrt.red == 0) {
        RUNTIME_ERROR("No red tetrahedra selected.")
    }
    Log::info("{} Red tetrahedra selected.", nrt.red);

    auto classResult = classifyRefinement(geom, i_tet);

    vector<Vec3> new_p;
    vector <idx>  new_t;
    vector <idx>  new_e;
    vector <idx> new_mat_t;
    vector <idx> new_mat_e;

    create_new_elements( geom,
                         classResult.i_tet, classResult.i_tri,
                         classResult.lines,
                         classResult.t_to_l, classResult.e_to_l,
                         new_p,
                         new_t, new_mat_t,
                         new_e, new_mat_e);

    modify_geometry(geom,
                    classResult.i_tet, classResult.i_tri,
                    new_p, new_t, new_e,
                    new_mat_t, new_mat_e );
}// end void Refine
