#include <solutionvector.h>
#include <material_numbers.h>
#include <algorithm>
#include <stdlib.h>
#include <omp.h>
#include <lc-representation.h>
#include <util/exception.h>
#include <util/logging.h>

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
            Log::info("FIXLC{} is strong", i + 1);
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
            Log::info("FIXLC{} is not strong, it is {}", i + 1, alignment->surface[i]->getAnchoringTypeName());
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
    if (FixedNodes == nullptr) {
        RUNTIME_ERROR("Could not allocate memory for fixed nodes.");
    }
    FixedValues = (double *) malloc(getnDimensions() * nFixed * sizeof(double));
    if (FixedValues == nullptr) {
        RUNTIME_ERROR("Could not allocate memory for fixed values.");
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
        RUNTIME_ERROR("Fixed arrays are already initialised.");
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
        RUNTIME_ERROR("Allocation failed when nFixed is " + to_string(nFixed));
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

void SolutionVector::setToFixedValues() {
    // SETS ALL VALUES TO CORRECT FIXED VALUES
    // MAKE SURE ARRAYS HAVE BEEN INITIALISED
    if (((FixedNodes == nullptr) || (FixedValues == nullptr)) &&
            (nFixed > 0)) {
        RUNTIME_ERROR("Null pointer for fixed nodes/values.");
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
    if (NewFixedNodes == nullptr) {
        RUNTIME_ERROR(fmt::format("Could not allocate memory for {} fixed nodes.", nFixed));
    }
    double *NewFixedValues = (double *) malloc(NewFixedSize * sizeof(double));
    if (NewFixedValues == nullptr) {
        RUNTIME_ERROR(fmt::format("Could not allocate memory for {} fixed values.", nFixed));
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









