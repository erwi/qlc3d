#include <dofmap.h>
#include <geometry.h>
#include <fixednodes.h>
#include <geom/periodicity.h>
#include <util/logging.h>

DofMap::DofMap(unsigned int nDof, unsigned int nDimensions): nDof(nDof), nDimensions(nDimensions) {
  dofs.resize(nDof * nDimensions, NOT_DOF);
}

void DofMap::calculateMapping(const std::unordered_set<unsigned int> &fixedNodes,
                              const std::vector<unsigned int> &periodicNodesMapping) {
  Log::info("Calculating mapping of {} nodes, {} fixed nodes, and {} periodic nodes",
            nDof, fixedNodes.size(), periodicNodesMapping.size());
  // TODO: rewrite this
  dofs.clear();
  dofs.resize(nDof * nDimensions, 0);
  // NODAL EQUIVALENCIES HAVE BEEN SET.
  // REPLACE DEPENDENT NODES WITH THEIR
  // INDEPENDENT EQUIVALENT NODES
  std::vector<idx> elimt(nDof, 0);    // convenience working copy of Elim
  for (idx i = 0; i < nDof; i++) {
    elimt.at(i) = periodicNodesMapping[i];
  }

  // MARK FIXED NODES. THESE WILL BE LATER ON REMOVED FROM
  // FREE DEGREES OF FREEDOM
  for (idx i = 0 ; i < nDof ; i++) {
    if (fixedNodes.find(i) != fixedNodes.end()) {
      elimt.at(i) = NOT_DOF;
    }
  }
  nFreeNodes = 0;
  std::vector <idx> elima(nDof, 0);   // Elim altered
  for (idx i = 0; i < nDof; i++) {   // SET TO 1,2,3...
    elima.at(i) = i;
  }
  elima.resize(nDof);
  // LOOP OVER EACH NODE. DECREASE INDEX TO ALL INDEPENDENT DOFs
  // THAT COME AFTER A DEPENDENT NODE (EQUIVALENT TO SHIFTING LEFT
  // ROWS/COLUMNS OF A MATRIX AFTER A COLUMN IS REMOVED)
  idx reduceCounter = 0;
  std::vector<idx> reduced(nDof, 0);


  for (idx i = 0; i < nDof; i++) { // if i'th node is dependent reduce dof for it and all nodes after it in elima
    if (elimt.at(i) != i) {
      reduceCounter ++;
    }
    reduced.at(i) = reduceCounter;
  }
  for (idx i = 0 ; i < nDof ; i++) {
    elima.at(i) -= reduced.at(i);
  }

  // SET DEPENDENT VARAIBLE INDEXES TO POINT TO CORRECT
  // INDEPENDENT DOF
  for (idx i = 0; i < nDof; i++) { // SET CORRECT VALUES
    if ((elimt.at(i) != i) && (elimt.at(i) != NOT_DOF)) { // IF i'th NODE IS DEPENDENT ( AND NOT FIXED)
      elima.at(i) = elima.at( elimt.at(i)); // IT WILL DEPEND ON THE CORRECTED DOF INDEX
    } else if (elimt.at(i) == NOT_DOF) { // KEEP FIXED NODE FLAGS
      elima.at(i) = NOT_DOF;
    }
  }
  // TOTAL NUMBER OF FREE DOFs THAT NEED TO BE SOLVED (PER DIMENSION)
  // nFreeNodes = *max_element(elima.begin(), elima.end() ) + 1;
  for (unsigned int i : elima) {
    if (i < NOT_DOF) {
      nFreeNodes = std::max(nFreeNodes, i + 1);
    }
  }
  // COPY BACK VALUES TO elim ARRAY
  for (idx i = 0 ; i < nDof ; i++) {
    dofs.at(i) = elima.at(i);
  }
  // EXPAND Elim IF MORE THAN ONE DIMENSIONS PER NODE
  if (nDimensions > 1) {
    for (idx j = 1; j < nDimensions ; j ++) {
      for (idx i = 0 ; i < nDof; i ++) {
        if (dofs.at(i) == NOT_DOF) {
          dofs.at(j * nDof + i) = NOT_DOF;
        } else {
          dofs.at(j * nDof + i) = dofs.at(i) + j * nFreeNodes;
        }
      }
    }
  }
}