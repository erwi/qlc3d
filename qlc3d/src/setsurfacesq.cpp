#include "../includes/qlc3d.h"
#include <math.h>


void setGlobalAngles(SolutionVector *q, LC* lc,double tilt, double twist, vector<int>* ind_nodes){
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
    vector <int>::iterator itr;
	
	//printf("setting angles for %i nodes\n", (int) ind_nodes->size() );
    for (itr = ind_nodes->begin() ; itr != ind_nodes->end() ; itr++)
    {
        q->setValue(*itr,0, (a1+a2)*(sqrt(6)/-2) );
        q->setValue(*itr,1, (a1+(a1+a2)/-2)*sqrt(2) );
        q->setValue(*itr,2, a3*sqrt(2) );
        q->setValue(*itr,3, a4*sqrt(2) );
        q->setValue(*itr,4, a5*sqrt(2) );
    }
}


void setHomeotropic(SolutionVector* q,LC* lc,vector<int>* ind_nodes,Geometry* geom){
    vector <int>::iterator i;
    double S = lc->getS0();
    for (i = ind_nodes->begin(); i != ind_nodes->end() ; i++){


        //double px = geom->getpX(*i);
        //double py = geom->getpY(*i);
        //double pz = geom->getpZ(*i);

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

// Freezes Q-tensor to current values at this surface
void setFrozenSurfaces(SolutionVector* q, vector<int>* ind_nodes)
{
	// frozen surfaces don't need anything done
	printf("frozen surface\n");
	
}



void setStrongSurfacesQ(SolutionVector *q, Alignment* alignment, LC* lc,  Geometry* geom){
    // sets only strong anchoring surfaces
    vector<int> ind_nodes;
    for (int i = 0 ; i < alignment->getnSurfaces() ; i++ ){
        ind_nodes.clear();
        // creates index of all nodes of this alignment surface type
        geom->e->FindIndexToMaterialNodes((i+1)*MAT_FIXLC1, &ind_nodes);
		
        if (ind_nodes.size()>0){ // if nodes found
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
    vector<int> ind_nodes;
    
    // loop over all surfaces loaded from settings file
    for (int i = 0 ; i < alignment->getnSurfaces() ; i++ )
    {
        ind_nodes.clear();
        // creates index of all nodes of this alignment surface type
        geom->e->FindIndexToMaterialNodes((i+1)*MAT_FIXLC1, &ind_nodes);

        if (ind_nodes.size()>0)
        { // if nodes found
    
            string AnchoringType = alignment->surface[i]->getAnchoringType();
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
        }// end if alignment nodes found
    }// end for loop over alignment surfaces

	
}
//end void setSurfacesQ




