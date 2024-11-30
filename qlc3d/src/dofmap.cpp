#include <dofmap.h>
#include <geometry.h>
#include <fixednodes.h>

DofMap::DofMap(unsigned int nDof, unsigned int nDimensions): nDof(nDof), nDimensions(nDimensions) {}

DofMap::DofMap(const std::vector<unsigned int> dofs) : dofs(dofs){

}

DofMap::DofMap(const DofMap &other) : dofs(other.dofs) {

}

void DofMap::calculateMapping(const Geometry &geom, const std::unordered_set<unsigned int> &fixedNodes) {
  // TODO: rewrite this
// IF NO PERIODIC NODES PRESENT, DON'T GENERATE EQUIVALENT NODES INDEXES
  if (!geom.getleft_right_is_periodic() &&
      !geom.gettop_bottom_is_periodic() &&
      !geom.getfront_back_is_periodic() &&
      fixedNodes.empty()) {
    return; // no periodic boundaries, can return
  }

  dofs.clear();
  dofs.resize(nDof * nDimensions, 0);
  // NODAL EQUIVALENCIES HAVE BEEN SET.
  // REPLACE DEPENDENT NODES WITH THEIR
  // INDEPENDENT EQUIVALENT NODES
  std::vector <idx> elimt(nDof, 0);    // convenience working copy of Elim
  for (idx i = 0 ; i < nDof ; i++) {
    elimt.at(i) =  geom.getPeriodicEquNode(i) ;
  }
  // MARK FIXED NODES. THSE WILL BE LATER ON REMOVED FROM
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
  for (idx i = 0; i < nDof; i++) {
    if (elimt.at(i) != i) {  // IF i'th NODE IS DEPENDENT
      for (idx j = i; j < nDof; j++) { // SHIFT DOWN ALL DOF INDEXES AFTER IT
        elima.at(j)--;
      }
    }
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