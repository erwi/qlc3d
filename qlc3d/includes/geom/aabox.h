
#ifndef PROJECT_QLC3D_AABOX_H
#define PROJECT_QLC3D_AABOX_H
#include <geom/vec3.h>


/** Axis-aligned bounding box class */
class AABox {
private:
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;

public:
  AABox() : x_min(0), x_max(0), y_min(0), y_max(0), z_min(0), z_max(0) {}
  AABox(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
          : x_min(xMin), x_max(xMax), y_min(yMin), y_max(yMax), z_min(zMin), z_max(zMax) {
    if (xMin > xMax || yMin > yMax || zMin > zMax) {
      throw std::invalid_argument("Invalid bounds: min bound is greater than max bound");
    }
  }

  [[nodiscard]] bool contains(double x, double y, double z) const {
    return (x >= x_min && x <= x_max &&
            y >= y_min && y <= y_max &&
            z >= z_min && z <= z_max);
  }

  [[nodiscard]] double volume() const {
    return (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
  }

  [[nodiscard]] Vec3 center() const {
    return {(x_min + x_max) / 2.0, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0};
  }

  void setBoundsX(double xMin, double xMax) {
    if (xMin > xMax) {
      throw std::invalid_argument("Invalid bounds: min bound is greater than max bound");
    }
    x_min = xMin;
    x_max = xMax;
  }

  void setBoundsY(double yMin, double yMax) {
    if (yMin > yMax) {
      throw std::invalid_argument("Invalid bounds: min bound is greater than max bound");
    }
    y_min = yMin;
    y_max = yMax;
  }

  void setBoundsZ(double zMin, double zMax) {
    if (zMin > zMax) {
      throw std::invalid_argument("Invalid bounds: min bound is greater than max bound");
    }
    z_min = zMin;
    z_max = zMax;
  }

  // Getters for the box boundaries
  [[nodiscard]] double getXMin() const { return x_min; }
  [[nodiscard]] double getXMax() const { return x_max; }
  [[nodiscard]] double getYMin() const { return y_min; }
  [[nodiscard]] double getYMax() const { return y_max; }
  [[nodiscard]] double getZMin() const { return z_min; }
  [[nodiscard]] double getZMax() const { return z_max; }
};

template <>
class fmt::formatter<AABox> {
public:
  constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
  template <typename Context>
  constexpr auto format (AABox const& b, Context& ctx) const {
    return format_to(ctx.out(), "{{xMin={}, xMax={}, yMin={}, yMax={}, zMin={}, zMax={}}}", b.getXMin(), b.getXMax(), b.getYMin(), b.getYMax(), b.getZMin(), b.getZMax());
  }
};


#endif //PROJECT_QLC3D_AABOX_H
