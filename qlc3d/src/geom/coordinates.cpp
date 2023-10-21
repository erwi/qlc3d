#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <util/exception.h>
#include <util/logging.h>

Coordinates::Coordinates(std::vector<Vec3> &&points) : points(std::move(points)) {}

void Coordinates::setPoints(std::vector<Vec3> &&points_) {
  points = points_;
}

std::shared_ptr<Coordinates> Coordinates::from(const std::vector<double> &p) {
  if (p.size() % 3 != 0) {
    RUNTIME_ERROR("Coordinates::from: p.size() % 3 != 0");
  }
  std::vector<Vec3> points;
  points.reserve(p.size() / 3);
  for (unsigned int i = 0; i < p.size() / 3; ++i) {
    points.emplace_back(p[3 * i], p[3 * i + 1], p[3 * i + 2]);
  }
  return std::make_shared<Coordinates>(std::move(points));
}

unsigned int Coordinates::size() const {
  return points.size();
}

const Vec3 &Coordinates::getPoint(unsigned int i) const {
  return points[i];
}

std::shared_ptr<Coordinates> Coordinates::clone() const {
  std::vector<Vec3> points = this->points; // deep copy
  auto coords = std::make_shared<Coordinates>();
  coords->setPoints(std::move(points));

  return coords;
}

std::vector<Vec3> Coordinates::findBoundingBox() const {
  double xMin = points[0].x();
  double xMax = points[0].x();
  double yMin = points[0].y();
  double yMax = points[0].y();
  double zMin = points[0].z();
  double zMax = points[0].z();

  for (const auto& point : points) {
    xMin = std::min(xMin, point.x());
    xMax = std::max(xMax, point.x());
    yMin = std::min(yMin, point.y());
    yMax = std::max(yMax, point.y());
    zMin = std::min(zMin, point.z());
    zMax = std::max(zMax, point.z());
  }
  return {{xMin, yMin, zMin}, {xMax, yMax, zMax}};
}

void Coordinates::append(const std::vector<double> &p) {
  if (p.size() % 3 != 0) {
    RUNTIME_ERROR("Coordinates::append: p.size() % 3 != 0");
  }
  for (unsigned int i = 0; i < p.size() / 3; ++i) {
    points.emplace_back(p[3 * i], p[3 * i + 1], p[3 * i + 2]);
  }
}

void Coordinates::scale(const Vec3 &scale) {
  for (auto &point : points) {
    point *= scale;
  }
}

void Coordinates::loadCoordinates(const idx *start, const idx *end, Vec3 *coordinatesOut) const {
  for (auto i = start; i != end; ++i) {
    coordinatesOut[i - start] = points[*i];
  }
}

void Coordinates::clear() {
  points.clear();
}