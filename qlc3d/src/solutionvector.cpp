#include <solutionvector.h>
#include <material_numbers.h>
#include <algorithm>
#include <stdlib.h>
#include <omp.h>
#include <lc-representation.h>
#include <util/exception.h>
#include <util/logging.h>
#include <util/hash.h>
const double SolutionVector::BIGNUM = 1e99;

// HACK DECLARATION OF TENSORTOVECTOR, NEEDED FOR FIXING POLYMERISED NODES
//
double *tensortovector(double *a, int npLC);

SolutionVector::~SolutionVector() {
    ClearFixed();
}
void SolutionVector::ClearFixed() {
    setnFixed(0);
    fixedNodeMaterial.clear();
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
    values = r.values;
    isFixed = r.isFixed;
    fixedNodeMaterial = r.fixedNodeMaterial;
    elim = r.elim;
    equNodes = r.equNodes;
    fixedNodes = r.fixedNodes;
    fixedValues = r.fixedValues;
    return *this;
}

SolutionVector::SolutionVector():
    nDoF(0),
    nFixed(0),
    nDimensions(0),
    nFreeNodes(0) { }

SolutionVector::SolutionVector(idx np, idx dim):
    nDoF(np),
    nFixed(0),
    nDimensions(dim),
    nFreeNodes(np) {
    Resize(np , dim);
}

void SolutionVector::setnFixed(idx n) {
    nFixed = n;
}

void SolutionVector::setValuesTo(const double &value) {
    for (size_t i = 0; i < (size_t) values.size(); i++) {
      values[i] = value;
    }
}

void SolutionVector::setValuesTo(const SolutionVector &other) {
    this->values = other.values;
}

void SolutionVector::Resize(const unsigned int &n, const unsigned int &dim) {
    ClearAll();
    nDoF = n;
    nFreeNodes = n;     // all nodes are free until set fixed
    nDimensions = dim;

    values.clear();
    values.resize(nDoF * nDimensions, 0.0);
}

void SolutionVector::ClearAll() {
    // CLEAR ALL DATA
    nDoF = 0;
    nFixed = 0;
    nDimensions = 0;
    nFreeNodes = 0;

    values.clear();
    isFixed.clear();
    fixedNodeMaterial.clear();
    elim.clear();
    equNodes.clear();

    fixedNodes.clear();
    fixedValues.clear();
}

void SolutionVector::setFixedNodesQ(const Alignment &alignment, const Mesh &e) {
    /*! sets fixed nodes index list and value lists, assuming values are currently correct
    FixedNodes is index list to all nodes that are fixed, i.e. strong anchoring
    FixedValues are the corresponding individual values that are fixed
    */
    // 1. First get index to all strong anchoring nodes
    vector <idx> ind_to_nodes;
    for (int i = 0 ; i < alignment.getnSurfaces() ; i ++) {
        auto surf = alignment.getSurface(i); // ref to i'th surface

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
                e.listFixLCSurfaces(temp_index, i + 1);
            }
            ind_to_nodes.insert(ind_to_nodes.end(), temp_index.begin(), temp_index.end());
        } else {
            Log::info("FIXLC{} is not strong, it is {}", i + 1, alignment.surface[i]->getAnchoringTypeName());
        }
    }// end for i
    // INDEX TO FIXED NODES MAY CONTAIN DUPLICATED ENTRIES, IF TWO FIXED
    // LC SURFACES ARE NEXT TO EACHOTHER. REMOVE THESE
    sort(ind_to_nodes.begin() , ind_to_nodes.end());
    vector <idx>::iterator itr;
    itr = unique(ind_to_nodes.begin() , ind_to_nodes.end());
    ind_to_nodes.erase(itr, ind_to_nodes.end());
    //2. allocate memory for index and fixed values
    fixedNodes.clear();
    fixedValues.clear();
    nFixed = ind_to_nodes.size();
    // IF NO FIXED NODES EXIST, LEAVE
    if (!nFixed) {
        setBooleanFixedNodeList();  // SET ALL TO NON-FIXED
        return;
    }
    fixedNodes.resize(getnDimensions() * nFixed, NOT_AN_INDEX);
    fixedValues.resize(getnDimensions() * nFixed, 0.0);
    // 3. copy index and values to arrays. This assumes that current Q-tensor values are correct
    // and will be fixed at these values (frozen to current)
    idx i = 0;
    //  nFreeNodes = nDoF - nFixed;
    for (itr = ind_to_nodes.begin() ; itr != ind_to_nodes.end() ; ++itr, ++i) {
        idx fixedNodeNumber = *itr;
        fixedNodes[i] = fixedNodeNumber;
        fixedNodes[i + 1 * nFixed] = fixedNodeNumber + 1 * getnDoF();
        fixedNodes[i + 2 * nFixed] = fixedNodeNumber + 2 * getnDoF();
        fixedNodes[i + 3 * nFixed] = fixedNodeNumber + 3 * getnDoF();
        fixedNodes[i + 4 * nFixed] = fixedNodeNumber + 4 * getnDoF();
        fixedValues[i]            = getValue(*itr, 0); //q1
        fixedValues[i + 1 * nFixed] = getValue(*itr, 1); //q2
        fixedValues[i + 2 * nFixed] = getValue(*itr, 2); //q3;
        fixedValues[i + 3 * nFixed] = getValue(*itr, 3); //q4;
        fixedValues[i + 4 * nFixed] = getValue(*itr, 4); //q5;
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
    fixedNodes.resize(nFixed * nDimensions, NOT_AN_INDEX);

    fixedNodeMaterial.clear();
    fixedNodeMaterial.resize(nFixed);
    fixedValues.resize(nFixed * nDimensions, 0);
    // SET FIXED NODES INDEX ARRAY VALUES
    for (idx i = 0 ; i < nFixed ; i++) {
        for (idx j = 0 ; j < nDimensions ; j++) {
            fixedNodes[i + j * nFixed] = fixed_nodes[i].nodenum;
        }
        fixedNodeMaterial.at(i) = fixed_nodes[i].mat;
    }
    setBooleanFixedNodeList();  // SET BOOLEAN FLAGS
}

void SolutionVector::setFixedNodesPot(const std::unordered_map<unsigned int, double> &potentialsByElectrode) {
  if (fixedNodes.empty()) {
      RUNTIME_ERROR("no fixed nodes defined");
  }
    for (const auto& [electrodeNumber, potential] : potentialsByElectrode) {
    Log::info("Setting Electrode{} nodes to potential {}", electrodeNumber, potential);
  }

  if (fixedNodes.size() != fixedValues.size()) {
    RUNTIME_ERROR(fmt::format("fixedNodes.size() ({}) != fixedValues.size() ({})", fixedNodes.size(), fixedValues.size()));
  }

  for (idx i = 0; i < nFixed; i++) {
    int mat = fixedNodeMaterial.at(i); // Get material number for ith fixed node
    size_t electrodeNumber = MATNUM_TO_ELECTRODE_NUMBER((size_t) mat); // greater or equal to 1 if electrode, i.e. 0 if not electrode
    if (electrodeNumber > 0) {
      fixedValues.at(i) = potentialsByElectrode.at(electrodeNumber);
      values[fixedNodes[i]] = fixedValues[i];
    }
  }
}

void SolutionVector::setBooleanFixedNodeList() {
    // A BOOLEAN FLAG FOR EACH NODE, WHETHER IT IS FIXED OR NOT.
    // MAYBE A PARAMETERS VARIABLE WIT BIT MASKS WOULD BE MORE EFFICIENT
    // IF MANY "FLAG" ARRAYS ARE NEEDED
    isFixed.clear();
    isFixed.resize(nDimensions * nDoF, false);

    // SET VALUE TO TRUE/FALSE FOR EACH NODE
    idx numFixedDoFs = getnFixed() * nDimensions;
    for (idx i = 0 ; i < numFixedDoFs ; i ++) { // then set only fixed nodes to true
        idx indToFixed = fixedNodes.at(i);
        isFixed.at(indToFixed) = true;
    }
}

void SolutionVector::setToFixedValues() {
    // SETS ALL VALUES TO CORRECT FIXED VALUES
    // MAKE SURE ARRAYS HAVE BEEN INITIALISED
    for (idx i = 0 ; i < nFixed * getnDimensions(); i ++) {
        idx ind = fixedNodes[i];
        double val = fixedValues[i];
        values[ind] = val;
    }
}

void SolutionVector::setPeriodicEquNodes(const Geometry &geom) {
    /*!
    SET VALUES IN THE ELIM ARRAY. THE ELIM ARRAY CONTAINS
    MAPPINGS FORM A NODE NUMBER TO ITS ACTUAL
    DEGREE OF FREEDOM (ITS ROW/COL POSITION IN THE GLOBAL MATRIX)
    */
    // IF NO PERIODIC NODES PRESENT, DON'T GENERATE EQUIVALENT NODES INDEXES
    if (!geom.getleft_right_is_periodic() &&
            !geom.gettop_bottom_is_periodic() &&
            !geom.getfront_back_is_periodic() &&
            (this->nFixed == 0)
       ) {
        return; // no periodic boundaries, can return
    }
    elim.clear();
    elim.resize(nDoF * nDimensions, 0);
    // NODAL EQUIVALENCIES HAVE BEEN SET.
    // REPLACE DEPENDENT NODES WITH THEIR
    // INDEPENDENT EQUIVALENT NODES
    std::vector <idx> elimt(nDoF, 0);    // convenience working copy of Elim
    for (idx i = 0 ; i < (idx) this->getnDoF() ; i++) {
        elimt.at(i) =  geom.getPeriodicEquNode(i) ;
    }
    // MARK FIXED NODES. THSE WILL BE LATER ON REMOVED FROM
    // FREE DEGREES OF FREEDOM
    for (idx i = 0 ; i < nDoF ; i++) {
        if (this->getIsFixed(i)) {
          elimt.at(i) = NOT_AN_INDEX;
        }
    }
    nFreeNodes = 0;
    std::vector <idx> elima(nDoF, 0);   // Elim altered
    for (idx i = 0 ; i < nDoF ; i++) {   // SET TO 1,2,3...
      elima.at(i) = i;
    }
    elima.resize(nDoF);
    // LOOP OVER EACH NODE. DECREASE INDEX TO ALL INDEPENDENT DOFs
    // THAT COME AFTER A DEPENDENT NODE (EQUIVALENT TO SHIFTING LEFT
    // ROWS/COLUMNS OF A MATRIX AFTER A COLUMN IS REMOVED)
    for (idx i = 0 ; i < nDoF ; i++) {
        if (elimt.at(i) != i) {  // IF i'th NODE IS DEPENDENT
            for (idx j = i ; j < nDoF ; j++) { // SHIFT DOWN ALL DOF INDEXES AFTER IT
                elima.at(j) --;
            }
        }
    }
    // SET DEPENDENT VARAIBLE INDEXES TO POINT TO CORRECT
    // INDEPENDENT DOF
    for (idx i = 0 ; i < nDoF ; i++) { // SET CORRECT VALUES
        if ((elimt.at(i) != i) && (elimt.at(i) != NOT_AN_INDEX)) { // IF i'th NODE IS DEPENDENT ( AND NOT FIXED)
            elima.at(i) = elima.at( elimt.at(i)); // IT WILL DEPEND ON THE CORRECTED DOF INDEX
        } else if (elimt.at(i) == NOT_AN_INDEX) { // KEEP FIXED NODE FLAGS
            elima.at(i) = NOT_AN_INDEX;
        }
    }
    // TOTAL NUMBER OF FREE DOFs THAT NEED TO BE SOLVED (PER DIMENSION)
    // nFreeNodes = *max_element(elima.begin(), elima.end() ) + 1;
    for (idx i = 0 ; i < (idx) elima.size() ; i++) {
        if (elima.at(i) < NOT_AN_INDEX)
            nFreeNodes = std::max(nFreeNodes , elima.at(i) + 1);
    }
    // COPY BACK VALUES TO elim ARRAY
    for (idx i = 0 ; i < nDoF ; i++) {
      elim.at(i) = elima.at(i);
    }
    // EXPAND Elim IF MORE THAN ONE DIMENSIONS PER NODE
    if (nDimensions > 1) {
        for (idx j = 1; j < nDimensions ; j ++) {
            for (idx i = 0 ; i < nDoF ; i ++) {
                if (elim.at(i) == NOT_AN_INDEX) {
                  elim.at(j * nDoF + i) = NOT_AN_INDEX;
                } else {
                  elim.at(j * nDoF + i) = elim.at(i) + j * nFreeNodes;
                }
            }
        }
    }
}// end setPeriondicEquNodes

void SolutionVector::setValue(const idx n,
                              const idx dim,
                              const double val) {
    assert(n < getnDoF());
    assert(dim < getnDimensions());
    values[n + dim * nDoF ] = val;
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
            values[dep] = values[indep];
        }
    }
}

void SolutionVector::setValue(const idx n, const qlc3d::TTensor &t) {
    assert(getnDimensions() == 5); // this should be SolutionVector for Q-tensor, not potential solution
    setValue(n, 0, t.t1());
    setValue(n, 1, t.t2());
    setValue(n, 2, t.t3());
    setValue(n, 3, t.t4());
    setValue(n, 4, t.t5());
}

void SolutionVector::setValue(idx i, const qlc3d::Director &d) {
    assert(getnDimensions() == 5); // this should be SolutionVector for director, not potential solution
    setValue(i, qlc3d::TTensor::fromDirector(d));
}

qlc3d::Director SolutionVector::getDirector(idx i) const {
    assert(getnDimensions() == 5); // this should be SolutionVector for director, not potential solution
    return qlc3d::TTensor(
            getValue(i, 0),
            getValue(i, 1),
            getValue(i, 2),
            getValue(i, 3),
            getValue(i, 4)).toDirector();
}

std::vector<qlc3d::Director> SolutionVector::getDirector() const {
    assert(getnDimensions() == 5); // this should be SolutionVector for director, not potential solution
    std::vector<qlc3d::Director> directors;
    directors.reserve(getnDoF());
    for (idx i = 0 ; i < getnDoF() ; i++) {
        directors.push_back(getDirector(i));
    }
    return directors;
}

int64_t SolutionVector::hashCode() {
    return hashCode64(&values[0], &values[values.size()]);
}



