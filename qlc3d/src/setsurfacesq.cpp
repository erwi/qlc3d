#include "../includes/qlc3d.h"
#include <math.h>
#include <globals.h>
#include <util/logging.h>
#include <util/exception.h>

void setGlobalAngles(SolutionVector *q,
                     double S0,
                     double tilt,
                     double twist,
                     vector<idx>* ind_nodes)
{
    // SETS GLOBAL TILT AND TIST ANGLES FOR Q-TENSOR AT NODES SPCIFIED IN ind_nodes

    tilt = tilt * PI / 180.0;
    twist= twist* PI / 180.0;

    double nx = cos(tilt)*cos(twist);
    double ny = cos(tilt)*sin(twist);
    double nz = sin(tilt);

    double a1 = S0 * (3 * nx * nx - 1) / 2.0;
    double a2 = S0 * (3 * ny * ny - 1) / 2.0;
    double a3 = S0 * (3 * nx * ny) / 2.0;
    double a4 = S0 * (3 * ny * nz) / 2.0;
    double a5 = S0 * (3 * nx * nz) / 2.0;

    // convert tensor basis
    vector <idx>::iterator itr;

    for (itr = ind_nodes->begin() ; itr != ind_nodes->end() ; ++itr) {
        q->setValue(*itr,0, (a1+a2)*(sqrt(6)/-2) );
        q->setValue(*itr,1, (a1+(a1+a2)/-2)*sqrt(2) );
        q->setValue(*itr,2, a3*sqrt(2) );
        q->setValue(*itr,3, a4*sqrt(2) );
        q->setValue(*itr,4, a5*sqrt(2) );
    }
}


void setHomeotropic(SolutionVector* q,
                    double S0,
                    vector<idx>* ind_nodes,
                    Geometry* geom){
    vector <idx>::iterator i;
    for (i = ind_nodes->begin(); i != ind_nodes->end() ; ++i) {
        double nx = geom->getNodeNormalsX(*i);
        double ny = geom->getNodeNormalsY(*i);
        double nz = geom->getNodeNormalsZ(*i);

        double a1 = S0 * (3 * nx * nx - 1) / 2.0;
        double a2 = S0 * (3 * ny * ny - 1) / 2.0;
        double a3 = S0 * (3 * nx * ny) / 2.0;
        double a4 = S0 * (3 * ny * nz) / 2.0;
        double a5 = S0 * (3 * nx * nz) / 2.0;

        q->setValue(*i, 0, (a1 + a2) * (sqrt(6) / -2));
        q->setValue(*i, 1, (a1 + (a1 + a2) / -2) * sqrt(2));
        q->setValue(*i, 2, a3 * sqrt(2));
        q->setValue(*i, 3, a4 * sqrt(2));
        q->setValue(*i, 4, a5 * sqrt(2));
    }
}


void setFrozenSurfaces(SolutionVector* q, vector<idx>* ind_nodes){
    // Freezes Q-tensor to current values at this surface
    // frozen surfaces don't need anything done
    q = q; // SUPPRESS WARNINGS
    ind_nodes = ind_nodes;
}

void setPolymeriseSurfaces(SolutionVector*q, double Sth){
    // FIXES Q-TENSOR ON ALL NODES WHERE ORDER PARAMETER IS BELOW DEFINED VALUE Sth
    // NOTHING NEEDS TO BE DONE HERE. SEE SolutionVector::setFixedNodesQ WHERE NODES
    // ARE SELECTED
    q=q;         // SUPPRESS WARNINGS
    Sth = Sth;

}

void setStrongSurfacesQ(SolutionVector *q,
                        Alignment* alignment,
                        double S0,
                        Geometry* geom)
{
    // sets only strong anchoring surfaces
    vector<idx> ind_nodes;
    for (int i = 0 ; i < alignment->getnSurfaces() ; i++ ){
        ind_nodes.clear();
        // creates index of all nodes of this alignment surface type
        // 08/02/12 geom->e->FindIndexToMaterialNodes((i+1)*MAT_FIXLC1, &ind_nodes);
        //geom->e->listNodesOfMaterial( ind_nodes, (i+1)*MAT_FIXLC1 );
        geom->e->listFixLCSurfaces(ind_nodes, i+1);

        if ( !ind_nodes.empty() ) { // if nodes found
            // get type of current surface
            AnchoringType aType = alignment->surface[i]->getAnchoringType();
            // depending on type, do different things...
            if (aType == Strong) {
                double tilt = alignment->surface[i]->getEasyTilt();
                double twist= alignment->surface[i]->getEasyTwist();
                setGlobalAngles(q, S0, tilt, twist, &ind_nodes);
            }
            else if (aType == Homeotropic) {
                setHomeotropic(q, S0, &ind_nodes, geom);
            }
            else if (aType == Freeze) {
                setFrozenSurfaces(q, &ind_nodes);
            }
        }// end if alignment nodes found
    }// end for loop over alignment surfaces
}
//end void setStrongSurfaces


void setManualNodesAnchoring(SolutionVector *q, double S0, Surface& surf){
    /*!
  fixes maually defined nodes to tilt/twist values
  */

    double tilt = surf.getEasyTilt();
    double twist = surf.getEasyTwist();

    // CONVERT SURFACE PARAMS TO VALID NODE INDEX VECTOR
    vector<idx> nodes_idx;
    for (size_t i = 0 ; i < surf.Params.size() ; i++){
        if ( (surf.Params[i] < 0 ) ){
            RUNTIME_ERROR("Negative node index when setting ManualNodesAnchoring.");
        }
        idx nodeIdx = static_cast<idx>(surf.Params[i]);
        nodes_idx.push_back( nodeIdx);
    }
    setGlobalAngles(q, S0, tilt, twist, &nodes_idx);
}


void setSurfacesQ(SolutionVector *q, Alignment* alignment, double S0,  Geometry* geom){
    /*! sets q-tensor values for all surfaces. */
    vector<idx> ind_nodes;
    
    // loop over all surfaces loaded from settings file
    for (int i = 0 ; i < alignment->getnSurfaces() ; i++ ){
        const Surface &surf = alignment->getSurface(i);

        ind_nodes.clear();
        geom->e->listFixLCSurfaces(ind_nodes, i+1);
        // MANUAL NODES HAVE NO SURFACE TRIANGLES IN MESH AND MUST BE HANDLED SEPARATELY
        if (surf.getAnchoringType() == ManualNodes) {
            Log::info("FIXLC{} is manual nodes anchoring.", i + 1);
            // MANUAL NODES SHOULD NOT DE DEFINED FOR A FIXLC# THAT IS PRESENT IN THE MESH
            if (!ind_nodes.empty()) {
                RUNTIME_ERROR(fmt::format("Can not set ManualNodes for FIXLC{}. Manual nodes should not be defined for a "
                              "FIXLC number that is present in the mesh. First available surface for this type is FIXLC{}",
                              i + 1, alignment->getnSurfaces() + 1));
            }
            setManualNodesAnchoring(q, S0, *(alignment->surface[i]));
            continue;
        }
        //
        // creates index of all nodes of this alignment surface type
        if ( !ind_nodes.empty() ) { // if nodes found

            AnchoringType aType = surf.getAnchoringType();
            double strength = surf.getStrength();

            if ((aType == Strong) ||
                (aType == Weak) ||
                (aType == Degenerate)) {
                double tilt = alignment->surface[i]->getEasyTilt();
                double twist= alignment->surface[i]->getEasyTwist();
                setGlobalAngles(q, S0, tilt, twist, &ind_nodes);
            }
            // if homeotropic OR degenerate with negative strength
            else if ((aType == Homeotropic) ||
                     ((aType == Degenerate) && (strength < 0))) {
                setHomeotropic(q, S0, &ind_nodes, geom);
            }
            else if (aType == Freeze) {
                setFrozenSurfaces(q, &ind_nodes);
            }
            else if ( aType == Polymerise) { // FREEZES ALL NODES WHOSE ORDER IS BELOW VALUE DEFINED IN STRENGTH
                double Sth = alignment->surface[i]->getStrength();
                setPolymeriseSurfaces(q,Sth);
            }
            else {
                RUNTIME_ERROR("Unhandled anchoring type " + surf.getAnchoringTypeName() + ".")
            }
        } else {
            RUNTIME_ERROR(fmt::format("FIXLC{} has ben defined in settings file, but no such surface found in "
                                      "the mesh.", i + 1));
    }
}// end for loop over alignment surfaces


}
//end void setSurfacesQ




