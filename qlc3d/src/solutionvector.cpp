#include <solutionvector.h>
#include <material_numbers.h>
#include <algorithm>
#include <lc-representation.h>
#include <util/exception.h>
#include <util/hash.h>
#include <dofmap.h>
#include <unordered_set>

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
    dofMap.reset();
}

void SolutionVector::initialiseLcBoundaries(const Geometry &geom, const Alignment &alignment) {
  auto &triangles = geom.getTriangles();
  std::unordered_set<unsigned int> allFixedNodes;
  for (auto &a : alignment.surface) {
    if (!a.isStrong()) {
      continue;
    }

    auto surfaceNodes = triangles.listFixLCSurfaceNodes(a.getFixLCNumber());
    allFixedNodes.insert(surfaceNodes.begin(), surfaceNodes.end());
  }

  dofMap = std::make_unique<DofMap>(nDoF, nDimensions);
  dofMap->calculateMapping(geom, allFixedNodes);
  numFixedNodes = allFixedNodes.size();
}

void SolutionVector::initialisePotentialBoundaries(const Geometry &geom,
                                                   const std::unordered_map<unsigned int, double> &potentialByElectrode) {

  if (potentialByElectrode.empty()) { //
    numFixedNodes = 0;
    dofMap = std::make_unique<DofMap>(nDoF, nDimensions);
    dofMap->calculateMapping(geom, {});
    return;
  }

  std::unordered_set<unsigned int> allFixedNodes;
  for (auto& [electrodeNumber, potential] : potentialByElectrode) {
    set<unsigned int> nodeIndices = geom.getTriangles().findElectrodeSurfaceNodes(electrodeNumber);

    for (auto &n : nodeIndices) {
      setValue(n, 0, potential);
    }

    allFixedNodes.insert(nodeIndices.begin(), nodeIndices.end());
  }

  dofMap = std::make_unique<DofMap>(nDoF, nDimensions);
  dofMap->calculateMapping(geom, allFixedNodes);
  numFixedNodes = allFixedNodes.size();
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

const DofMap& SolutionVector::getDofMap() const {
  if (dofMap == nullptr) {
    RUNTIME_ERROR("DofMap not initialised");
  }
  return *dofMap;
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



