#ifndef PROJECT_QLC3D_COORDINATES_H
#define PROJECT_QLC3D_COORDINATES_H
#include <vector>
#include <memory>
#include "globals.h"

class Vec3;

class Coordinates {
  std::vector<Vec3> points;

public:
  Coordinates() = default;
  Coordinates(std::vector<Vec3> &&points);

  static std::shared_ptr<Coordinates> from(const std::vector<double> &p);

  /** get reference to ith point */
  [[nodiscard]] const Vec3& getPoint(unsigned int i) const;

  /** return number of points */
  [[nodiscard]] unsigned int size() const;

  [[nodiscard]] std::vector<Vec3> findBoundingBox() const;

  void clear();
  void append(const std::vector<double> &p);
  void scale(const Vec3& scale);
  /** create a deep copy of this object */
  [[nodiscard]] std::shared_ptr<Coordinates> clone() const;
private:
  void setPoints(std::vector<Vec3> &&points_);
};

#endif //PROJECT_QLC3D_COORDINATES_H
