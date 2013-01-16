#include "../includes/qlc3d.h"
#include <math.h>
#include <globals.h>

void setGlobalAngles(SolutionVector *q,
                     LC* lc,
                     double tilt,
                     double twist,
                     vector<idx>* ind_nodes)
{
    // SETS GLOBAL TILT AND TIST ANGLES FOR Q-TENSOR AT NODES SPCIFIED IN ind_nodes

    tilt = tilt * PI / 180.0;
    twist= twist* PI / 180.0;

    //printf(" tilt = %f, twist = %i\n",tilt,twist);
    double nx = cos(tilt)*cos(twist);
    double ny = cos(tilt)*sin(twist);
    double nz = sin(tilt);

    double S = lc->getS0();

    double a1 = S*(3* nx * nx - 1)/2.0;
    double a2 = S*(3* ny * ny - 1)/2.0;
    double a3 = S*(3* nx * ny)/2.0;
    double a4 = S*(3* ny * nz)/2.0;
    double a5 = S*(3* nx * nz)/2.0;

    // convert tensor basis
    vector <idx>::iterator itr;

    //printf("setting angles for %i nodes\n", (int) ind_nodes->size() );
    for (itr = ind_nodes->begin() ; itr != ind_nodes->end() ; ++itr)
    {
        q->setValue(*itr,0, (a1+a2)*(sqrt(6)/-2) );
        q->setValue(*itr,1, (a1+(a1+a2)/-2)*sqrt(2) );
        q->setValue(*itr,2, a3*sqrt(2) );
        q->setValue(*itr,3, a4*sqrt(2) );
        q->setValue(*itr,4, a5*sqrt(2) );
    }
}


void setHomeotropic(SolutionVector* q,
                    LC* lc,
                    vector<idx>* ind_nodes,
                    Geometry* geom)
{
    vector <idx>::iterator i;
    double S = lc->getS0();
    for (i = ind_nodes->begin(); i != ind_nodes->end() ; ++i)
    {

        double nx = geom->getNodeNormalsX(*i);
        double ny = geom->getNodeNormalsY(*i);
        double nz = geom->getNodeNormalsZ(*i);

        double a1 = S*(3* nx * nx - 1)/2.0;
        double a2 = S*(3* ny * ny - 1)/2.0;
        double a3 = S*(3* nx * ny)/2.0;
        double a4 = S*(3* ny * nz)/2.0;
        double a5 = S*(3* nx * nz)/2.0;

        q->setValue(*i,0, (a1+a2)*(sqrt(6)/-2) );
        q->setValue(*i,1, (a1+(a1+a2)/-2)*sqrt(2) );
        q->setValue(*i,2, a3*sqrt(2) );
        q->setValue(*i,3, a4*sqrt(2) );
        q->setValue(*i,4, a5*sqrt(2) );
    }
}


void setFrozenSurfaces(SolutionVector* q, vector<idx>* ind_nodes)
{
    // Freezes Q-tensor to current values at this surface
    // frozen surfaces don't need anything done
    q = q; // SUPPRESS WARNINGS
    ind_nodes = ind_nodes;
    printf("frozen surface\n");

}

void setPolymeriseSurfaces(SolutionVector*q, double Sth)
{
// FIXES Q-TENSOR ON ALL NODES WHERE ORDER PARAMETER IS BELOW DEFINED VALUE Sth
// NOTHING NEEDS TO BE DONE HERE. SEE SolutionVector::setFixedNodesQ WHERE NODES
// ARE SELECTED
   q=q;         // SUPPRESS WARNINGS
   Sth = Sth;

}

void setStrongSurfacesQ(SolutionVector *q,
                        Alignment* alignment,
                        LC* lc,
                        Geometry* geom)
{
    // sets only strong anchoring surfaces
    vector<idx> ind_nodes;
    for (int i = 0 ; i < alignment->getnSurfaces() ; i++ ){
        ind_nodes.clear();
        // creates index of all nodes of this alignment surface type
        // 08/02/12 geom->e->FindIndexToMaterialNodes((i+1)*MAT_FIXLC1, &ind_nodes);
        geom->e->listNodesOfMaterial( ind_nodes, (i+1)*MAT_FIXLC1 );


        if ( !ind_nodes.empty() ){ // if nodes found
            string AnchoringType = alignment->surface[i]->getAnchoringType();

            if ( (AnchoringType.compare("Strong") == 0) ){
                double tilt = alignment->surface[i]->getEasyTilt();
                double twist= alignment->surface[i]->getEasyTwist();
                setGlobalAngles(q,lc,tilt,twist,&ind_nodes);
            }

            // homeotropic anchoring counts as strong for now...
            else if ( (AnchoringType.compare("Homeotropic")==0)  ){
                setHomeotropic(q,lc,&ind_nodes,geom);
            }
            // frozen surfaces are also strong
            else if ( AnchoringType.compare("Freeze") == 0 )
            {
                setFrozenSurfaces(q, &ind_nodes);
            }
        }// end if alignment nodes found
    }// end for loop over alignment surfaces
}
//end void setStrongSurfaces




void setSurfacesQ(SolutionVector *q, Alignment* alignment, LC* lc,  Geometry* geom){
    /*! sets q-tensor values for all surfaces. */
    vector<idx> ind_nodes;
    
    // loop over all surfaces loaded from settings file
    for (int i = 0 ; i < alignment->getnSurfaces() ; i++ )
    {
        ind_nodes.clear();
        // creates index of all nodes of this alignment surface type
        // 08/02/12 geom->e->FindIndexToMaterialNodes((i+1)*MAT_FIXLC1, &ind_nodes);
        geom->e->listNodesOfMaterial( ind_nodes , (i+1)*MAT_FIXLC1 );
        if ( !ind_nodes.empty() )
        { // if nodes found

            string AnchoringType = alignment->surface[i]->getAnchoringType();
            //unsigned int AnchoringNum = alignment->surface[i]->getAnchoringNum(); // <--- SHOULD USE THIS
            if ((AnchoringType.compare("Strong") == 0) ||
                    (AnchoringType.compare("Weak") == 0))
            {
                double tilt = alignment->surface[i]->getEasyTilt();
                double twist= alignment->surface[i]->getEasyTwist();
                setGlobalAngles(q,lc,tilt,twist,&ind_nodes);
            }

            // if homeotropic OR degenerate with negative strength
            else if ( (AnchoringType.compare("Homeotropic")==0)  ||
                      (AnchoringType.compare("Degenerate") ==0 &&
                       (alignment->surface[i]->getStrength() < 0) ) )
            {
                setHomeotropic(q,lc,&ind_nodes,geom);
            }
            else if ( AnchoringType.compare("Freeze") == 0 )
            {
                setFrozenSurfaces(q, &ind_nodes);
            }
            else if ( AnchoringType.compare("Polymerise") ) // FREEZES ALL NODES WHOSE ORDER IS BELOW VALUE DEFINED IN STRENGTH
            {
                double Sth = alignment->surface[i]->getStrength();
                setPolymeriseSurfaces(q,Sth);
            }
        }// end if alignment nodes found
    }// end for loop over alignment surfaces


}
//end void setSurfacesQ




