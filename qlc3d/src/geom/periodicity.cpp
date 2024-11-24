#include <geom/periodicity.h>
#include <geom/vec3.h>
#include <geom/coordinates.h>
#include <mesh.h>
#include <material_numbers.h>
#include <util/logging.h>
#include "util/exception.h"

PeriodicityType::PeriodicityType(const Mesh &triangles) {
  left_right_is_periodic = false;
  front_back_is_periodic = false;
  top_bottom_is_periodic = false;

  const double EPS = 1e-7;
  for (unsigned int i = 0; i < triangles.getnElements(); i++) {
    unsigned int material = triangles.getMaterialNumber(i);
    if (material != MAT_PERIODIC) {
      continue;
    }

    auto normal = triangles.getSurfaceNormal(i);

    if (fabs(fabs(normal.x()) - 1.0) < EPS) {
      left_right_is_periodic = true;
    }
    else if (fabs(fabs(normal.y()) - 1.0) < EPS) {
      front_back_is_periodic = true;
    }
    else if (fabs(fabs(normal.z()) - 1.0) < EPS) {
      top_bottom_is_periodic = true;
    }
    else {
      RUNTIME_ERROR(fmt::format("Periodic surface element {} has invalid normal ({})", i, normal));
    }
    // IF ALL SURFACES HAVE ALREADY BEEN IDENTIFIED AS PERIODIC
    // NO NEED TO CHECK FURTHER TRIANGLES
    if (left_right_is_periodic && front_back_is_periodic && top_bottom_is_periodic) {
      break;
    }
  }
}


set<unsigned int> findPeriodicNodes(const Mesh &tris) {

  set<unsigned int> leftFace, rightFace, frontFace, backFace, topFace, bottomFace;


  auto equals = [](double a, double b) { return fabs(a - b) < 1e-7; };
  auto isLeftFace = [equals](const Vec3 &normal) { return equals(normal.x(), 1.); };
  auto isRightFace = [equals](const Vec3 &normal) { return equals(normal.x(), -1.); };
  auto isFrontFace = [equals](const Vec3 &normal) { return equals(normal.y(), 1.); };
  auto isBackFace = [equals](const Vec3 &normal) { return equals(normal.y(), -1.); };
  auto isTopFace = [equals](const Vec3 &normal) { return equals(normal.z(), 1.); };
  auto isBottomFace = [equals](const Vec3 &normal) { return equals(normal.z(), -1.); };

  const Vec3 rightNormal(-1, 0, 0);
  const Vec3 frontNormal(0, 1, 0);
  const Vec3 backNormal(0, -1, 0);
  const Vec3 topNormal(0, 0, 1);
  const Vec3 bottomNormal(0, 0, -1);

  for (unsigned int i = 0; i < tris.getnElements(); i++) {
    unsigned int material = tris.getMaterialNumber(i);
    if (material != MAT_PERIODIC) {
      continue;
    }

    auto normal = tris.getSurfaceNormal(i);

    if (isLeftFace(normal)) {
      leftFace.insert(i);
    } else if (isRightFace(normal)) {
      rightFace.insert(i);
    } else if (isFrontFace(normal)) {
      frontFace.insert(i);
    } else if (isBackFace(normal)) {
      backFace.insert(i);
    } else if (isTopFace(normal)) {
      topFace.insert(i);
    } else if (isBottomFace(normal)) {
      bottomFace.insert(i);
    }
  }
  return leftFace;
}


PeriodicNodesMapping::PeriodicNodesMapping(const Mesh &tris,
                                           const Coordinates &coords,
                                           const PeriodicityType &periodicityType) {

  if (!periodicityType.isAnyPeriodic()) {
    Log::warn("No periodic surfaces detected. Periodic nodes mapping will not be created.");
    return;
  }

  std::set<unsigned int> periodicNodes = findPeriodicNodes(tris);

  if (periodicityType.isFrontBackPeriodic() && !periodicityType.isLeftRightPeriodic() && !periodicityType.isTopBottomPeriodic()) {

  } else if (periodicityType.isFrontBackPeriodic() && periodicityType.isLeftRightPeriodic() && !periodicityType.isTopBottomPeriodic()) {

  } else if (periodicityType.isFrontBackPeriodic() && periodicityType.isLeftRightPeriodic() && periodicityType.isTopBottomPeriodic()) {

  } else {
    RUNTIME_ERROR("Periodicity type not supported.");
  }
}


void PeriodicNodesMapping::matchNodesFrontBack(const set<unsigned int> periNodes, const Coordinates &coords) {
 std::set<unsigned int> front, back;

 for (auto &i : periNodes) {
  // coordinates.get
 }

}


