#include <geom/coordinates.h>
#include <geom/vec3.h>

void Coordinates::setPoints(std::vector<Vec3> &&points_) {
  points = points_;
}

std::shared_ptr<Coordinates> Coordinates::from(const double *p, unsigned int numPoints) {
  auto coords = std::make_shared<Coordinates>();

  return coords;
}
