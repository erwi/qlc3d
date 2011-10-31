#include <meshrefinement.h>
#include <refinement.h>
#include <vector>
#include <stdio.h>
using std::vector;


double getMaxChange(Geometry& geom, SolutionVector& v){
/*! Loops through mesh and checks for largest change in variable v*/
    double dv = 0;
    int e = 0;
    int npt = geom.t->getnNodes();
    //double* vals = new double[ npt ];
    double vals[4] = {0,0,0,0};
    for (unsigned int elem = 0 ; elem < (unsigned int) geom.t->getnElements() ; elem++){
        for (int dim = 0 ; dim <  v.getnDimensions() ; dim++){
            for (int node = 0 ; node < npt ; node++){
                int nn = geom.t->getNode(elem, node);
                vals[node] = v.getValue(nn,dim);
            }
                double max = *max_element( vals , vals+ npt ); // find max and min values in this element
                double min = *min_element( vals , vals+ npt );
                double dvl = max - min; // local dv

                if ( dvl > dv ){
                    dv = dvl;
                    e = elem;
                    printf(" elem %i, dim %i\n", elem, dim);

                }
                    }// end for dim
    }// end for elem
    printf(" maddq in element %i\n", e);

    //delete vals;
    return dv;
}

/*
bool needsRefinement(Geometry &geom, SolutionVector &q, MeshRefinement &meshrefinement){
// checks whether additional end refinement is needed. returns true if needed.
//

    double minchange = 1e50;
    vector  < AutoRef > :: iterator itr = meshrefinement.AutoRefinement.begin();

    // minimum allowed dq is assumed to be the last specified value in MaxValue vector in AUTOREF settings structure
    for (itr; itr!= meshrefinement.AutoRefinement.end() ; itr++){
        int nitr = itr->getNumIterations();
        minchange = itr->getMaxValue( nitr-1 );
    }

    if (meshrefinement.AutoRefinement.size() > 0 ){
        printf("min allowed dq = %f\n", minchange);
        printf("checking for maximum change in Q per element...");
        double dq = getMaxChange( geom , q);
        printf("max dq = %f\n", dq);
        return ( minchange < dq );
    }

    // return false by default if AUTOREF is not defined
        return false;


}
*/
