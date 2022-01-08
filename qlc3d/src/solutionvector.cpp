#include <solutionvector.h>
#include <material_numbers.h>
#include <algorithm>
#include <stdlib.h>
#include <omp.h>
#include <lc-representation.h>
const double SolutionVector::BIGNUM = 1e99;


// HACK DECLARATION OF TENSORTOVECTOR, NEEDED FOR FIXING POLYMERISED NODES
//
double *tensortovector(double *a, int npLC);


SolutionVector::~SolutionVector() {
    ClearFixed();
    if (Values) {
        free(Values);
    }
    Values = NULL;
    if (Elim) {
        free(Elim);
    }
    Elim = NULL;
    if (EquNodes != NULL) {
        free(EquNodes);
    }
    EquNodes = NULL;
}
void SolutionVector::ClearFixed() {
    if (FixedNodes && FixedValues) {
        free(FixedNodes);
        free(FixedValues);
        FixedNodes = NULL;
        FixedValues = NULL;
        setnFixed(0);
    }
    if (IsFixed) {
        free(IsFixed);
        IsFixed = NULL;
    }
    if (FixedNodeMaterial) {
        free(FixedNodeMaterial);
        FixedNodeMaterial = NULL;
    }
}

SolutionVector &SolutionVector::operator=(const SolutionVector &r) {
    if (this == &r)
        return *this;
    // CLEAR ALL DATA
    ClearAll();
    // REALLOCATE AND COPY
    nDoF    = r.nDoF;
    nFixed  = r.nFixed;
    nDimensions = r.nDimensions;
    nFreeNodes  = r.nFreeNodes;
    // VALUES
    if (r.Values) {
        Values = (double *) malloc(nDoF * nDimensions * sizeof(double));
        memcpy(Values, r.Values , nDoF * nDimensions * sizeof(double));
    }
    // FIXED NODED
    //ClearFixed();
    if (r.FixedValues) {
        FixedValues = (double *) malloc(nFixed * nDimensions * sizeof(double));
        memcpy(FixedValues, r.FixedValues, nFixed * nDimensions * sizeof(double));
    }
    if (r.FixedNodes) {
        FixedNodes  = (idx *)    malloc(nFixed * nDimensions * sizeof(idx));
        memcpy(FixedNodes, r.FixedNodes, nFixed * nDimensions * sizeof(idx));
    }
    if (r.FixedNodeMaterial) {
        FixedNodeMaterial = (idx *) malloc(nFixed * sizeof(idx));
        memcpy(FixedNodeMaterial, r.FixedNodeMaterial, nFixed * sizeof(idx));
    }
    // PERIODIC EQU NODES AND ELIM
    if (r.Elim) {
        Elim    = (idx *)    malloc(nDoF * nDimensions *  sizeof(idx));
        memcpy(Elim, r.Elim , nDoF * nDimensions * sizeof(idx));
    }
    if (r.EquNodes) {
        EquNodes    = (idx *)    malloc(nDoF * nDimensions *  sizeof(idx));
        memcpy(EquNodes, r.EquNodes, nDoF * nDimensions * sizeof(idx));
    }
    if (r.IsFixed) {
        IsFixed = (bool *)   malloc(nDoF * nDimensions *  sizeof(bool));
        memcpy(IsFixed, r.IsFixed, nDoF * nDimensions * sizeof(bool));
    }
    return *this;
}

SolutionVector::SolutionVector():
    nDoF(0),
    nFixed(0),
    nDimensions(0),
    FixedNodeMaterial(NULL),
    IsFixed(NULL),
    Elim(NULL),
    EquNodes(NULL),
    nFreeNodes(0),
    FixedNodes(NULL),
    FixedValues(NULL),
    Values(NULL) {
}

SolutionVector::SolutionVector(idx np, idx dim):
    nDoF(np),
    nFixed(0),
    nDimensions(dim),
    FixedNodeMaterial(NULL),
    IsFixed(NULL),
    Elim(NULL),
    EquNodes(NULL),
    nFreeNodes(np),
    FixedNodes(NULL),
    FixedValues(NULL),
    Values(NULL) {
    Allocate((unsigned int) np, (unsigned int) dim);
}

void SolutionVector::Allocate(const idx np, const idx ndim) {
    Resize(np , ndim); // OOPS. TWO NAMES, ONE FUNCTION, SOUNDS FAMILIAR...
}

void SolutionVector::setnFixed(idx n) {
    nFixed = n;
}

void SolutionVector::setValuesTo(const double &value) {
    // ALL VALUES ARE SET TO THAT OF THE INPUT VALUE
#ifndef DEBUG
    #pragma omp parallel for
#endif
    for (size_t i = 0 ; i < (size_t) nDoF ; i ++)
        Values[i] = value;
}// end void setValuesTo

void SolutionVector::setValuesTo(const double *values) {
    // ALL VALUES ARE SET TO THOSE OF THE INPUT VECTOR
    // VECTOR LENGHTS MUCH MATCH, NO CHECKING IS PERFORMED HERE!!
    idx n = nDoF * nDimensions;
    memcpy(Values, values, n * sizeof(double));
    //for (idx i = 0 ; i < n ; i++)
    //    this->Values[i] = values[i];
}
void SolutionVector::setValuesTo(const SolutionVector &other) {
    this->setValuesTo(other.Values);
}
void SolutionVector::Resize(const unsigned int &n, const unsigned int &dim) {
    ClearAll();
    nDoF = n;
    nFreeNodes = n;     // all nodes are free until set fixed
    nDimensions = dim;
    Values = (double *) malloc(nDoF * nDimensions * sizeof(double));
    memset(Values, 0, nDoF * nDimensions * sizeof(double));
}

void SolutionVector::ClearAll() {
    // CLEAR ALL DATA
    nDoF = 0;
    nFixed = 0;
    nDimensions = 0;
    nFreeNodes = 0;
    if (IsFixed)    free(IsFixed);      IsFixed     = NULL;
    if (Elim)       free(Elim);         Elim        = NULL;
    if (EquNodes)   free(EquNodes);     EquNodes    = NULL;
    if (FixedNodes) free(FixedNodes);   FixedNodes  = NULL;
    if (FixedValues)free(FixedValues);  FixedValues = NULL;
    if (Values)     free(Values);       Values      = NULL;
    if (FixedNodeMaterial) free(FixedNodeMaterial); FixedNodeMaterial = NULL;
}

void SolutionVector::setFixedNodesQ(Alignment *alignment, Mesh *e) {
    /*! sets fixed nodes index list and value lists, assuming values are currently correct
    FixedNodes is index list to all nodes that are fixied
    FixedValues are the corresponding individual values that are fixed
    */
    // 1. First get index to all strong anchoring nodes
    vector <idx> ind_to_nodes;
    for (int i = 0 ; i < alignment->getnSurfaces() ; i ++) {
        auto surf = alignment->getSurface(i); // ref to i'th surface

        if (surf.isStrong()) {
            printf("FIXLC%i is strong\n", i + 1);
            vector <idx> temp_index;
            // IN CASE OF MANUAL NODES ANCHORING, THE NODE INDEXES ARE
            // ALREADY KNOWN AND LISTED IN THE Surface.Params VECTROR
            if (surf.getAnchoringType() == ManualNodes) {
                // Cast params vector <double> to indices
                for (size_t i = 0 ; i < surf.Params.size() ; i++)
                    temp_index.push_back(static_cast<idx>(surf.Params[i]));
            }
            // HACK TO FIX POLYMERISED NODES
            /*
            if ( alignment->getAnchoringNum(i) == ANCHORING_POLYMERISE ) {
                temp_index.clear();

                // Sth IS THRESHOLD VALUE FOR CHOOSING POLYMERISED NODES
                double Sth = alignment->surface[i]->getStrength();

                // CONVERT QURRENT Q-TENSOR TO DIRECTOR AND ORDER PARAMETERS
                idx npLC = getnDoF();
                double *n = tensortovector(Values, npLC ); // get vector data

                // FIND INDEX TO ALL NODES WHERE ORDER IS LESS OR EQUAL TO Sth
                for (idx i = 3*npLC ; i < 4*npLC; i++){ // FOR EACH ORDER PARAMETER
                    if ( n[i] <= Sth )
                        temp_index.push_back(i - 3*npLC);
                }
                double minS = *std::min_element(n+3*npLC, n+4*npLC);
                double maxS = *std::max_element(n+3*npLC, n+4*npLC);
                printf("S range = [%f,%f], Sth = %f, Polymerised %u nodes", minS, maxS,Sth,(idx)temp_index.size() );

                delete [] n;

                // CHECK THAT NOT ALL NODES WERE CHOSEN
                if (temp_index.size() >= (size_t) npLC) {
                    printf("error in %s, %s, setting Polymerised fixed nodes",__func__,__FILE__);
                    printf("too many nodes Polymerised. Try setting the threshold to a larger value - bye!\n");
                    exit(1);
                }
            }
            */
            //else{    // IF NOT POLYMERISED SURFACE, SELECT BY SURFACE MATERIAL NUMBER
            //e->listNodesOfMaterial( temp_index , curFixLC_MAT_NUMBER );
            else {
                e->listFixLCSurfaces(temp_index, i + 1);
            }
            ind_to_nodes.insert(ind_to_nodes.end(), temp_index.begin(), temp_index.end());
        } else {
            printf("FIXLC%i is not strong, it is %s\n", i + 1, alignment->surface[i]->getAnchoringTypeName().c_str());
        }
    }// end for i
    // INDEX TO FIXED NODES MAY CONTAIN DUPLICATED ENTRIES, IF TWO FIXED
    // LC SURFACES ARE NEXT TO EACHOTHER. REMOVE THESE
    sort(ind_to_nodes.begin() , ind_to_nodes.end());
    vector <idx>::iterator itr;
    itr = unique(ind_to_nodes.begin() , ind_to_nodes.end());
    ind_to_nodes.erase(itr, ind_to_nodes.end());
    //2. allocate memory for index and fixed values
    if (FixedNodes != NULL) free(FixedNodes);
    if (FixedValues != NULL) free(FixedValues);
    nFixed = ind_to_nodes.size();
    // IF NO FIXED NODES EXIST, LEAVE
    if (!nFixed) {
        setBooleanFixedNodeList();  // SET ALL TO NON-FIXED
        return;
    }
    FixedNodes = (idx *) malloc(getnDimensions() * nFixed * sizeof(idx));
    if (FixedNodes == NULL) {
        printf("error - SolutionVector::setFixedNodesQ() - could not allocate memory for fixed nodes, bye!");
        exit(1);
    }
    FixedValues = (double *) malloc(getnDimensions() * nFixed * sizeof(double));
    if (FixedValues == NULL) {
        printf("error - SolutionVector::setFixedNodesQ() - could not allocate memory for fixed values, bye!");
        exit(1);
    }
    // 3. copy index and values to arrays. This assumes that current Q-tensor values are correct
    // and will be fixed at these values (frozen to current)
    idx i = 0;
    //  nFreeNodes = nDoF - nFixed;
    for (itr = ind_to_nodes.begin() ; itr != ind_to_nodes.end() ; ++itr, ++i) {
        idx fixedNodeNumber = *itr;
        FixedNodes[i] = fixedNodeNumber;
        FixedNodes[i + 1 * nFixed] = fixedNodeNumber + 1 * getnDoF();
        FixedNodes[i + 2 * nFixed] = fixedNodeNumber + 2 * getnDoF();
        FixedNodes[i + 3 * nFixed] = fixedNodeNumber + 3 * getnDoF();
        FixedNodes[i + 4 * nFixed] = fixedNodeNumber + 4 * getnDoF();
        FixedValues[i]            = getValue(*itr, 0); //q1
        FixedValues[i + 1 * nFixed] = getValue(*itr, 1); //q2
        FixedValues[i + 2 * nFixed] = getValue(*itr, 2); //q3;
        FixedValues[i + 3 * nFixed] = getValue(*itr, 3); //q4;
        FixedValues[i + 4 * nFixed] = getValue(*itr, 4); //q5;
    }
    setBooleanFixedNodeList();
}

/*
void SolutionVector::setFixedNodes(vector<int> *Material, vector<double> *val ,
                                   idx *Elem,
                                   idx *Mat,
                                   idx nElem,
                                   idx nNodes)
{// sets fixed node lists with node numbers and corresponding fixed values
    // *Material = fixed material numbers defined in settings file
    // *val = values corresponding to *Material
    // *Elem = array of element node numbers
    // *Mat  = array of material numbers corresponding to elem
    // nElem and nNodes = dimensions of element node arrays

    if (Material->size()!=val->size())
    {
        printf("\nerror-solutionvector::SetFixedNodes\nnumber of fixed values must equal number of materials - bye!");
        exit(1);
    }

    nFixed = 0; // reset number of fixed nodes
    if (FixedNodes != NULL) free(FixedNodes);   // these may have been set earlier, and should be reset
    if (FixedValues!= NULL) free(FixedValues);

    vector<int> vec_FixedNodes;
    vector<double> vec_FixedValues;


    //int n_fix_mat = Material->size();
    //printf("%i fixed materials\n",n_fix_mat);

    vector<int>::iterator i;
    vector<int>::iterator j;
    vector<double>::iterator Val; // fixed value

    for ( i = Material->begin(), Val =val->begin() ; i != Material->end() ; i++,Val++ ){

        // find index to all elements of material i
        vector<int> fix_elem;
        vector<int> fix_nodes;
        for (idx x = 0 ; x < nElem ; x ++)
        {
            if ( (Mat[x] & 31*MAT_ELECTRODE1) == *i )  // bitwise testing of material number
            {
                fix_elem.push_back(x);
                //printf("found\n");
            }
        }

        // ...then get all nodes from the fixed elements
        for (j = fix_elem.begin(); j != fix_elem.end() ; j ++)
        {
            for (int x = 0 ; x < nNodes ; x ++)
            {
                fix_nodes.push_back(Elem[(*j)*nNodes+x]);
            }
        }
        // remove repeated nodes from list
        sort(fix_nodes.begin(), fix_nodes.end());
        vector<int>::iterator sorted;
        sorted = unique(fix_nodes.begin(), fix_nodes.end());

        //add nodes and values to final vector
        for (j = fix_nodes.begin() ; j != sorted ; j++)
        {
            //printf("%i ",*j);
            vec_FixedNodes.push_back(*j);
            vec_FixedValues.push_back(*Val);
        }
    }

    // Allocate memory for fixed node lists
    nFixed = vec_FixedNodes.size();
    FixedNodes = (int*)malloc(nFixed * sizeof(int));
    FixedValues= (double*)malloc(nFixed * sizeof(double));

    // nFreeNodes = nDoF - nFixed; // update free nodes count
    // Copy values form vectors to memory
    int c= 0;
    for (i = vec_FixedNodes.begin(), Val = vec_FixedValues.begin(), c = 0;
         i != vec_FixedNodes.end() ; i ++, Val++, c++)
    {
        FixedNodes[c] = *i;
        FixedValues[c] = *Val;
    }



}
*/
void SolutionVector::allocateFixedNodesArrays(Geometry &geom) {
    // ALLOCATES ARRAYS FOR MANAGING FIXED NODES.
    //
    // CALCULATES NUMBER OF FIXED NODES nFixed
    // SETS FIXED NODES INDEX ARRAY FixedNodes VALUES
    // ALLOCATES MEMORY FOR FIXED VALUES, BUT DOES NOT SET VALUES
    // SETS FIXED NODES BOOLEAN FLAG ARRAY (CALLS setBoleanFixedNodesList)
    // THIS IS AN INITIALIZATION ONLY FUNCTION.
    // RESET ALL FIRST, IF RESIZING SOLUTION VECTOR
    // (E.G. AFTER MESH REFIENEMENT)
    if (nFixed > 0) {
        printf("error in %s, fixed arrays are already initialised - bye!\n", __func__);
        exit(1);
    }
    // SEPARATE INDEX TO ALL FIXED NODES
    // (SOME WILL BE REPEATED)
    Mesh &e = *geom.e;   // PTR TO SURFACE MESH
    std::vector<SolutionVectorNameSpace::node> fixed_nodes;
    for (idx i = 0 ; i < e.getnElements() ; i++) {
        int mat = e.getMaterialNumber(i);
        size_t indE = MATNUM_TO_ELECTRODE_NUMBER((size_t) mat);
        if (indE) { // IF ELECTRODE ELEMENT
            for (idx j = 0 ; j < e.getnNodes() ; j++)
                fixed_nodes.push_back(
                    SolutionVectorNameSpace::
                    node(e.getNode(i, j), mat));
        }
    }
    // REMOVE REPEATED NODE INDEXES
    sort(fixed_nodes.begin(), fixed_nodes.end());
    std::vector<SolutionVectorNameSpace::node> :: iterator itr;
    itr = unique(fixed_nodes.begin(), fixed_nodes.end());
    fixed_nodes.erase(itr, fixed_nodes.end());
    nFixed = fixed_nodes.size();
    // CREATE FIXED NODES INDEXES AND VALUES ARRAYS
    FixedNodes = (idx *) malloc(nFixed * sizeof(idx) * nDimensions);
    FixedNodeMaterial = (idx *) malloc(nFixed * sizeof(idx));
    FixedValues = (double *) malloc(nFixed * nDimensions * sizeof(double));
    if ((!FixedNodes) || (!FixedValues) || (!FixedNodeMaterial)) {
        printf("error in %s, malloc returned NULL - bye!\n", __func__);
        exit(1);
    }
    // SET FIXED NODES INDEX ARRAY VALUES
    for (idx i = 0 ; i < nFixed ; i++) {
        for (idx j = 0 ; j < nDimensions ; j++) {
            FixedNodes[i + j * nFixed] = fixed_nodes[i].nodenum;
        }
        FixedNodeMaterial[i] = fixed_nodes[i].mat;
    }
    // ALL FIXED VALUES TO 0
    memset(FixedValues, 0 , nFixed * nDimensions * sizeof(double));
    setBooleanFixedNodeList();  // SET BOOLEAN FLAGS
}

void SolutionVector::setFixedNodesPot(Electrodes *electrodes) {
    // SETS VALUES IN FixedValues
    if (electrodes->isEField())  // DON'T DO ANYTHING IF UNIFORM E-FIELD DEFINED
        return;
    for (idx i = 0 ; i < nFixed ; i++) {
        int mat = FixedNodeMaterial[i]; // GET MATERIAL NUMBER FOR ith FIXED NODE
        size_t indE = MATNUM_TO_ELECTRODE_NUMBER((size_t) mat);
        if (indE) {
            double pot = electrodes->getCurrentElectrodePotential(indE - 1);
            FixedValues[i] = pot;
        }
    }
}

void SolutionVector::setFixedNodesPot(Electrodes &electrodes,
                                      Mesh *surface_mesh)

{
    // sets fixed nodes for potentials,
    // this method is horrible and need to be rewritten
    if (nFixed > 0) { //IF RE-SETTING, CLEAR OLD
        if (FixedNodes) {
            free(FixedNodes);
            FixedNodes = NULL;
        }
        if (FixedValues) {
            free(FixedValues);
            FixedValues = NULL;
        }
        if (IsFixed) {
            free(IsFixed);
            IsFixed = NULL;
        }
        nFixed = 0;
    }
    for (size_t i = 0 ; i < electrodes.getnElectrodes() ; i++) { // for each electrode
        double pot = electrodes.getCurrentElectrodePotential(i);
        // set all fixed nodes for electrode i to its current potential
        AddFixed((i + 1)*MAT_ELECTRODE1,
                 pot ,
                 surface_mesh);
    }// end for each electrode
    setBooleanFixedNodeList();
}
// end void setFixedNodesPot


void SolutionVector::setBooleanFixedNodeList() {
    // A BOOLEAN FLAG FOR EACH NODE, WHETHER IT IS FIXED OR NOT.
    // MAYBE A PARAMETERS VARIABLE WIT BIT MASKS WOULD BE MORE EFFICIENT
    // IF MANY "FLAG" ARRAYS ARE NEEDED
    if (IsFixed != NULL) {
        free(IsFixed);
    }
    // ALLOCATE MEMORY FOR ARRAY
    size_t size = nDimensions * nDoF * sizeof(bool);
    IsFixed = (bool *) malloc(size);
    memset(IsFixed , false , size);  // set all to false
    // SET VALUE TO TRUE/FALSE FOR EACH NODE
    idx numFixedDoFs = getnFixed() * nDimensions;
    for (idx i = 0 ; i < numFixedDoFs ; i ++) { // then set only fixed nodes to true
        idx indToFixed = FixedNodes[i];
        IsFixed[indToFixed] = true;
    }
}


void SolutionVector::PrintFixedNodes() {
    //debugging
    printf("number of fixed nodes is : %i\n", nFixed);
    for (idx i = 0 ; i < nFixed ; i++) {
        for (idx j = 0 ; j < nDimensions ; j++) {
            printf("fixed node %u, node number %u, value %1.3f, ", i, FixedNodes[i], FixedValues[i]);
        }
        if (FixedNodeMaterial) {
            printf("material %u", FixedNodeMaterial[i]);
        }
        printf("\n");
    }
}
void SolutionVector::PrintValues() {
    for (idx i = 0 ; i < nDoF * nDimensions; i++) {
        printf("Value[%i] = %e\n", i, Values[i]);
    }
}
void SolutionVector::PrintElim() {
    if (!Elim) {
        printf("Elim array is NULL\n");
        return;
    }
    printf("\nPrinting for dimension 0 only!:\n");
    for (idx i = 0 ; i < nDoF ; i++) {
        printf("Elim[%i] = %i\n", i, Elim[i]);
    }
}
void SolutionVector::PrintEquNodes() {
    if (!EquNodes) {
        printf("EquNodes array is NULL\n");
        return;
    }
    printf("\nEquNodes:\n");
    for (idx i = 0 ; i < nDoF ; i ++) {
        printf("EquNodes[%i] = %i \n", i, EquNodes[i]);
    }
}
void SolutionVector::PrintIsFixed() {
    printf("%i fixed, %i free nodes:\n", this->getnFixed(), this->getnFreeNodes());
    for (idx i = 0 ; i < this->getnDoF()*this->getnDimensions() ; i++) {
        printf("node[%i] = ", i);
        if (getIsFixed(i)) {
            printf("TRUE\n");
        } else {
            printf("FALSE\n");
        }
    }
}

void SolutionVector::setToFixedValues() {
    // SETS ALL VALUES TO CORRECT FIXED VALUES
    // MAKE SURE ARRAYS HAVE BEEN INITIALISED
    if (((FixedNodes == NULL) || (FixedValues == NULL)) &&
            (nFixed > 0)) {
        printf("error - SolutionVector::setToFixedValues, NULL pointer - bye!\n");
        exit(1);
    }
    for (idx i = 0 ; i < nFixed * getnDimensions(); i ++) {
        int ind = FixedNodes[i];
        double val = FixedValues[i];
        Values[ ind ] = val;
    }
}




void SolutionVector::setPeriodicEquNodes(Geometry *geom) {
    /*!
    SET VALUES IN THE ELIM ARRAY. THE ELIM ARRAY CONTAINS
    MAPPINGS FORM A NODE NUMBER TO ITS ACTUAL
    DEGREE OF FREEDOM (ITS ROW/COL POSITION IN THE GLOBAL MATRIX)
    */
    // IF NO PERIODIC NODES PRESENT, DON'T GENERATE EQUIVALENT NODES INDEXES
    if (!geom->getleft_right_is_periodic() &&
            !geom->gettop_bottom_is_periodic() &&
            !geom->getfront_back_is_periodic() &&
            (this->nFixed == 0)
       ) {
        return; // no periodic boundaries, can return
    }
    if (Elim != NULL) free(Elim);   // allocate memory for equivalent nodes
    Elim  = (idx *) malloc(nDoF * nDimensions * sizeof(idx));
    // NODAL EQUIVALENCIES HAVE BEEN SET.
    // REPLACE DEPENDENT NODES WITH THEIR
    // INDEPENDENT EQUIVALENT NODES
    std::vector <idx> elim(nDoF, 0);    // convenience working copy of Elim
    for (idx i = 0 ; i < (idx) this->getnDoF() ; i++) {
        elim[i] =  geom->getPeriodicEquNode(i) ;
    }
    // MARK FIXED NODES. THSE WILL BE LATER ON REMOVED FROM
    // FREE DEGREES OF FREEDOM
    for (idx i = 0 ; i < nDoF ; i++) {
        if (this->getIsFixed(i))
            elim[i] = NOT_AN_INDEX;
    }
    nFreeNodes = 0;
    std::vector <idx> elima(nDoF, 0);   // Elim altered
    for (idx i = 0 ; i < nDoF ; i++)    // SET TO 1,2,3...
        elima[i] = i;
    elima.resize(nDoF);
    // LOOP OVER EACH NODE. DECREASE INDEX TO ALL INDEPENDENT DOFs
    // THAT COME AFTER A DEPENDENT NODE (EQUIVALENT TO SHIFTING LEFT
    // ROWS/COLUMNS OF A MATRIX AFTER A COLUMN IS REMOVED)
    for (idx i = 0 ; i < nDoF ; i++) {
        if (elim[i] != i) {  // IF i'th NODE IS DEPENDENT
            for (idx j = i ; j < nDoF ; j++) { // SHIFT DOWN ALL DOF INDEXES AFTER IT
                elima[j] --;
            }
        }
    }
    // SET DEPENDENT VARAIBLE INDEXES TO POINT TO CORRECT
    // INDEPENDENT DOF
    for (idx i = 0 ; i < nDoF ; i++) { // SET CORRECT VALUES
        if ((elim[i] != i) && (elim[i] != NOT_AN_INDEX)) { // IF i'th NODE IS DEPENDENT ( AND NOT FIXED)
            elima[i] = elima[ elim[i] ]; // IT WILL DEPEND ON THE CORRECTED DOF INDEX
        } else if (elim[i] == NOT_AN_INDEX) { // KEEP FIXED NODE FLAGS
            elima[i] = NOT_AN_INDEX;
        }
    }
    // TOTAL NUMBER OF FREE DOFs THAT NEED TO BE SOLVED (PER DIMENSION)
    // nFreeNodes = *max_element(elima.begin(), elima.end() ) + 1;
    for (idx i = 0 ; i < (idx) elima.size() ; i++) {
        if (elima[i] < NOT_AN_INDEX)
            nFreeNodes = std::max(nFreeNodes , elima[i] + 1);
    }
    // COPY BACK VALUES TO Elim ARRAY
    for (idx i = 0 ; i < nDoF ; i++)
        Elim[i] = elima[i];
    // EXPAND Elim IF MORE THAN ONE DIMENSIONS PER NODE
    if (nDimensions > 1) {
        for (idx j = 1; j < nDimensions ; j ++) {
            for (idx i = 0 ; i < nDoF ; i ++) {
                if (Elim[i] == NOT_AN_INDEX)
                    Elim[j * nDoF + i] = NOT_AN_INDEX;
                else
                    Elim[j * nDoF + i] = Elim[i] + j * nFreeNodes;
            }
        }
    }
}// end setPeriondicEquNodes

void SolutionVector::setCornerElim(
    list <int> &corn0, // sets periodic equivalent nodes for 4 corners
    list <int> &corn1, // corn1[i] = corn0[i]
    list <int> &corn2, // corn2[i] = corn0[i]
    list <int> &corn3, // corn3[i] = corn0[i]
    int *Elim,
    const int &dim, // direction of corner 0,1,2 -> x,y,z
    double *p // coordinates
)

{
    /*! This sets Elim vector for periodic corner nodes.
    * Nodes are arranged so that corn1[i] = corn2[i] = corn2[i] = corn0[i]
    * AbsDist is a function pointer to a function that calculates absolute distance along
    * X, Y or Z axes (use geom->getAbsNDist, where N is X, Y, or Z).
    * getCoord is a function pointer to getting X,Y or Z coordinate (replaces geom->getpN)
    * */
    list<int> :: iterator c0;
    list<int> :: iterator c1;
    list<int> :: iterator c2;
    list<int> :: iterator c3;
    double eps = 1e-5;
    for (c0 = corn0.begin() ; c0 != corn0.end() ; ++c0) {
        // outer corner loop - corn0
        bool found = false;
        double dist = 0;
        double p0 = p[(*c0) * 3 + dim ]; // dim coordinate of node c0
        int iNearest = -1;
        double minDist = BIGNUM;
        // COMPARISON WITH C1
        for (c1 = corn1.begin() ; c1 != corn1.end() ; ++c1) {
            // inner corner loop - corn1
            double p1 = p[(*c1) * 3 + dim];
            dist = fabs(p0 - p1);
            if (dist <= minDist) {
                minDist = dist;
                iNearest = *c1;
            }
            if (dist <= eps) { // if same coordinate -> match!
                Elim[*c1] = *c0;
                found = true;
                break;
            }
        }
        if (!found) {
            printf("error - corner node 1 not found - bye\n");
            printf("corner node 0 = %i,\tat [%e,%e,%e]\n", *c0, p[*c0 * 3 + 0], p[*c0 * 3 + 1], p[*c0 * 3 + 2]);
            printf("nearest node  = %i,\tat [%e,%e,%e]\n", iNearest, p[iNearest * 3 + 0], p[iNearest * 3 + 1], p[iNearest * 3 + 2]);
            printf("distance = %e\n", minDist);
            exit(1);
        }
        // COMPARISON WITH C2
        found = false;
        iNearest = -1;
        minDist = BIGNUM;
        for (c2 = corn2.begin() ; c2 != corn2.end() ; ++c2) {
            // inner corner loop - corn2
            double p2 = p[(*c2) * 3 + dim ];
            dist = fabs(p0 - p2);
            if (dist <= minDist) {
                minDist = dist;
                iNearest = *c2;
            }
            if (dist <= eps) { // if same cordinate -> match!
                Elim[*c2] = *c0;
                found = true;
                break;
            }
        }
        if (!found) {
            printf("error - corner node 2 not found - bye\n");
            printf("corner node 0 = %i,\tat [%e,%e,%e]\n", *c0, p[*c0 * 3 + 0], p[*c0 * 3 + 1], p[*c0 * 3 + 2]);
            printf("nearest node  = %i,\tat [%e,%e,%e]\n", iNearest, p[iNearest * 3 + 0], p[iNearest * 3 + 1], p[iNearest * 3 + 2]);
            printf("distance = %e\n", minDist);
            exit(1);
        }
        // COMPARISON WITH C3
        found = false;
        iNearest = -1;
        minDist = BIGNUM;
        for (c3 = corn3.begin() ; c3 != corn3.end() ; ++c3) {
            // inner corner loop - corn2
            double p3 = p[(*c3) * 3 + dim ];
            dist = fabs(p0 - p3);
            if (dist <= minDist) {
                minDist = dist;
                iNearest = *c1;
            }
            if (dist <= eps) {  // if same cordinate -> match!
                Elim[*c3] = *c0;
                found = true;
                break;
            }
        }
        if (!found) {
            printf("error - corner node 3 not found - bye\n");
            printf("corner node 0 = %i,\tat [%e,%e,%e]\n", *c0, p[*c0 * 3 + 0], p[*c0 * 3 + 1], p[*c0 * 3 + 2]);
            printf("nearest node  = %i,\tat [%e,%e,%e]\n", iNearest, p[iNearest * 3 + 0], p[iNearest * 3 + 1], p[iNearest * 3 + 2]);
            printf("distance = %e\n", minDist);
            exit(1);
        }
    }// end for loop over corn0 nodes c0
}// end void setCornerElim



void SolutionVector::setFaceElim(list <int> &face0,  // face1[i] = face0[i]
                                 list <int> &face1,
                                 int *Elim,
                                 const int &norm, // face normal, 0,1,2 -> x,y,z
                                 double *p) { // pointer to node coordinates
    // norm is face normal, need coordinates to vectors parallel to it
    int ind1 = (norm + 1) % 3;   // ind is a pre-calculated offeset to coordinate comparison in p.
    int ind2 = (norm + 2) % 3;   // e.g. norm = 0 -> ind1 = 1, ind2 = 2
    double eps = 1e-5; // accuracy of coordinate comparison
    // SEARCH FOR NODE EQUIVALENCIES BY COMPARING COORDINATES
    // THAT ARE PERPENDICULAR TO FACE NORMAL
    list <int>:: iterator F0;    // FACE 0 NODES
    list <int>:: iterator F1;    // FACE 1 NODES
    int fc, bc; // debug counters
    for (F0 = face0.begin(), fc = 0; F0 != face0.end() ; ++F0, ++fc) { // LOOP OVER FACE 0
        bool found = false;
        double f1 = p[3 * (*F0) + ind1 ]; // coordinates of node F2 in face0
        double f2 = p[3 * (*F0) + ind2 ];
        for (F1 = face1.begin(), bc = 0; F1 != face1.end(); ++F1, ++bc) {
            double fa = p[3 * (*F1) + ind1]; // coordinates of node F1 in face 1
            double fb = p[3 * (*F1) + ind2];
            //compare coordinates
            double dist1 = fabs(f1 - fa); // distances in plane
            double dist2 = fabs(f2 - fb);
            double tdist = dist1 * dist1 + dist2 * dist2;

            if (tdist < eps * eps) { // compare squared distances
                Elim[*F1] = *F0;
                found = true;
                break;
            }
        }// end for B
        if (!found) {
            printf("error - matching periodic faces with normal %i - bye!\n", norm);
            exit(1);
        }
    }//end for F
}

void SolutionVector::setValue(const idx n,
                              const idx dim,
                              const double val) {
    assert(n < getnDoF());
    assert(dim < getnDimensions());
    Values[n + dim * nDoF ] = val;
}

void SolutionVector::setValue(const idx n, const qlc3d::TTensor &t) {
    assert(getnDimensions() == 5); // this should be SolutionVector for Q-tensor, not potential solution
    setValue(n, 0, t.t1());
    setValue(n, 1, t.t2());
    setValue(n, 2, t.t3());
    setValue(n, 3, t.t4());
    setValue(n, 4, t.t5());
}

void SolutionVector::AddFixed(int mat, double val, Mesh *mesh) {
    vector <idx> ind_p; // index to all nodes of material mat
    // 08/02/12 mesh->FindIndexToMaterialNodes(mat,&ind_p);
    mesh->listNodesOfMaterial(ind_p, (idx) mat);
    // Allocate memory for old + new fixed nodes and values
    idx NewFixedSize = (idx) ind_p.size() + nFixed; // number of old + new fixed nodes
    idx *NewFixedNodes = (idx *) malloc(NewFixedSize * sizeof(idx)); // allocate enough memory for old + new
    if (NewFixedNodes == NULL) {
        printf("error - SolutionVector::AddFixed(int,double,Mesh*) - could not allocate memory for fixed nodes, bye!");
        exit(1);
    }
    double *NewFixedValues = (double *) malloc(NewFixedSize * sizeof(double));
    if (NewFixedValues == NULL) {
        printf("error - SolutionVector::AddFixed(int,double,Mesh*) - could not allocate memry for fixed values, bye!");
        exit(1);
    }
    // Copy old fixed nodes and values and add new ones
    vector <idx>::iterator itr;
    itr = ind_p.begin();
    for (idx i = 0 ; i < NewFixedSize ; i++) {
        if (i < nFixed) { // copy old ones
            NewFixedNodes[i] = FixedNodes[i];
            NewFixedValues[i] = FixedValues[i];
        } else { // add new ones
            NewFixedNodes[i] = *itr; //  = iterator to index to material nodes
            ++itr;
            NewFixedValues[i] = val;
        }
    }// end for i
    // free old arrays and set pointer to New arrays
    if (FixedValues != NULL) free(FixedValues);
    if (FixedNodes != NULL) free(FixedNodes);
    FixedValues = NewFixedValues;
    FixedNodes  = NewFixedNodes;
    NewFixedValues = NULL;
    NewFixedNodes  = NULL; // necessary to keep only single pointer to memory location?
    setnFixed(NewFixedSize); // update counter
}
void SolutionVector::EnforceEquNodes(const Geometry &geom) {
    // MAKES SURE THAT VALUES ON PERIODIC SURFACES ARE OK.
    // THIS MAY BE NEEDED e.g. AT THE START OF A SIMULATION
    // OR TO AVOID ACCUMULATION OF NUMERICAL NOISE(?)
    // CALCULATES PERIODIC EQUIVALEN NODE FROM ELIM IN CASES WHERE
    // WHERE MORE THAN 1 DEGREE OF FREEDOM EXISTS (i.e. Q-TENSOR)
    for (size_t i = 0 ; i < (size_t) nDimensions ; i ++) {
        for (size_t j = 0 ; j < (size_t) nDoF; j++) {
            size_t equDof = geom.getPeriodicEquNode(j);
            size_t dep = i * nDoF + j;        // DEPENDENT NODE
            size_t indep = i * nDoF + equDof; // EQUIVALENT INDEPENDENT NODE
            Values[ dep ] = Values[ indep ];
        }
    }
}









