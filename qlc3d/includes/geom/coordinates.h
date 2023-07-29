#ifndef PROJECT_QLC3D_COORDINATES_H
#define PROJECT_QLC3D_COORDINATES_H
#include <vector>
#include <memory>

class Vec3;

class Coordinates {
  std::vector<Vec3> points;

public:
  //explicit Coordinates(std::vector<Vec3> &points) : points(std::move(points)) {}
  Coordinates() = default;


  static std::shared_ptr<Coordinates> from(const double* p, unsigned int numPoints);

  //const Vec3& getPoint(unsigned int i) const;

private:
  void setPoints(std::vector<Vec3> &&points_);
};

#endif //PROJECT_QLC3D_COORDINATES_H
