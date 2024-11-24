#include <solutionvector.h>
#include <material_numbers.h>
#include <algorithm>
#include <lc-representation.h>
#include <util/exception.h>
#include <util/logging.h>
#include <util/hash.h>
#include <dofmap.h>

// HACK DECLARATION OF TENSORTOVECTOR, NEEDED FOR FIXING POLYMERISED NODES
//
double *tensortovector(double *a, int npLC);

SolutionVector::~SolutionVector() {
    fixedNodes.clear();
}
//void SolutionVector::ClearFixed() {
    //setnFixed(0);
    //fixedNodes.clear();
//}

SolutionVector &SolutionVector::operator=(const SolutionVector &r) {
    if (this == &r)
        return *this;
    // CLEAR ALL DATA
    ClearAll();
    // REALLOCATE AND COPY
    nDoF    = r.nDoF;
    nDimensions = r.nDimensions;
    values = r.values;

    DofMap *otherMap = r.dofMap.get();
    if (otherMap != nullptr) {
      dofMap = std::make_unique<DofMap>(*otherMap);
    }

    return *this;
}

SolutionVector::SolutionVector():
    nDoF(0),
    nDimensions(0)
    { }

SolutionVector::SolutionVector(idx np, idx dim):
    nDoF(np),
    nDimensions(dim),
    values (nDoF * nDimensions, 0.) {
    //Resize(np , dim);
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
    nDimensions = dim;
    values.clear();
    values.resize(nDoF * nDimensions, 0.0);

}

void SolutionVector::ClearAll() {
    // CLEAR ALL DATA
    nDoF = 0;
    nDimensions = 0;
    values.clear();
    fixedNodes.clear();
    dofMap.reset();
}

void SolutionVector::setFixedLcNodes(const Alignment &alignment, const Mesh &e) {
  std::unordered_map<unsigned int, double> fixedQByNodeNumber;
  for (auto &a : alignment.surface) {
    if (!a.isStrong()) {
      continue;
    }

    auto fixLCNumber = a.getFixLCNumber();
    auto surfaceNodes = e.listFixLCSurfaceNodes(fixLCNumber);


    for (auto &i : surfaceNodes) {
      fixedQByNodeNumber[i] = getValue(i, 0); // q1
      fixedQByNodeNumber[i + getnDoF()] = getValue(i, 1); // q2
      fixedQByNodeNumber[i + 2 * getnDoF()] = getValue(i, 2); // q3
      fixedQByNodeNumber[i + 3 * getnDoF()] = getValue(i, 3); // q4
      fixedQByNodeNumber[i + 4 * getnDoF()] = getValue(i, 4); // q5
    }
  }

  fixedNodes.setFixedNodesAndValues(fixedQByNodeNumber);
}

void SolutionVector::setFixedPotentials(const Mesh &triangles,
                        const std::unordered_map<unsigned int, double> &potentialByElectrode) {
  fixedNodes.setFixedNodesPot(triangles, potentialByElectrode);

  auto &fn = fixedNodes.getFixedValueByNodeIndex();

  for (auto & [nodeIndex, potential] : fn) {
    values[nodeIndex] = potential;
  }
}

void SolutionVector::setPeriodicEquNodes(const Geometry &geom) {
    /*!
    SET VALUES IN THE ELIM ARRAY. THE ELIM ARRAY CONTAINS
    MAPPINGS FORM A NODE NUMBER TO ITS ACTUAL
    DEGREE OF FREEDOM (ITS ROW/COL POSITION IN THE GLOBAL MATRIX)
    */
    dofMap = std::make_unique<DofMap>(nDoF, nDimensions);
    dofMap->calculateMapping(geom, fixedNodes);
    /*
    // IF NO PERIODIC NODES PRESENT, DON'T GENERATE EQUIVALENT NODES INDEXES
    if (!geom.getleft_right_is_periodic() &&
            !geom.gettop_bottom_is_periodic() &&
            !geom.getfront_back_is_periodic() &&
            (this->nFixed == 0) &&
            fixedNodules.getnFixedNodes() == 0
       ) {
        return; // no periodic boundaries, can return
    }
    //elim.clear();
    std::vector<unsigned int> elim;
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
      if (fixedNodules.isFixedNode(i)) {
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

    dofMap = std::make_unique<DofMap>(elim);
     */
}

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

void SolutionVector::loadValues(const idx *start, const idx *end, double *valuesOut) const {
  for (idx* i = const_cast<idx *>(start); i != end; ++i) {
#ifndef NDEBUG
    if (*i >= getnDoF()) {
      RUNTIME_ERROR("index out of bounds " + std::to_string(*i));
    }
#endif
    valuesOut[i - start] = getValue(*i);
  }
}

void SolutionVector::loadQtensorValues(const idx *start, const idx *end, qlc3d::TTensor* tensorOut) const {
  for (idx* i = const_cast<idx *>(start); i != end; ++i) {
#ifndef NDEBUG
    if (*i >= getnDoF()) {
      RUNTIME_ERROR("index out of bounds " + std::to_string(*i));
    }
#endif
    tensorOut[i - start] = qlc3d::TTensor(
            getValue(*i, 0),
            getValue(*i, 1),
            getValue(*i, 2),
            getValue(*i, 3),
            getValue(*i, 4));
  }
}

void SolutionVector::loadEquNodes(const idx* start, const idx* end, idx* equNodesOut) const {
  for (idx* i = const_cast<idx *>(start); i != end; ++i) {
    equNodesOut[i - start] = getEquNode(*i);
  }
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



