#ifndef PROJECT_QLC3D_VEC3_H
#define PROJECT_QLC3D_VEC3_H

#include <sstream>
#include <cmath>
#include <fmt/format.h>

/** Simple 3D vector class */
class Vec3 {
  double x_ = 0, y_ = 0, z_ = 0;
public:
  Vec3() = default;

  Vec3(double x, double y, double z) : x_(x), y_(y), z_(z) {}

  Vec3(const Vec3 &) = default;

  Vec3 &operator=(const Vec3 &) = default;

  Vec3(Vec3 &&) = default;

  Vec3 &operator=(Vec3 &&) = default;

  ~Vec3() = default;

  [[nodiscard]] const double &x() const {
    return x_;
  }

  [[nodiscard]] const double &y() const {
    return y_;
  }

  [[nodiscard]] const double &z() const {
    return z_;
  }

  /** get dimension dim where dim=0 is x, dim=1 is y, dim=2 = z else throws runtime error */
  [[nodiscard]] double getDimension(unsigned int dim) const {
    switch (dim) {
      case 0:
        return x_;
      case 1:
        return y_;
      case 2:
        return z_;
      default:
        throw std::runtime_error("invalid dimension " + std::to_string(dim));
    }
  }

  [[nodiscard]] double dot(const Vec3 &rhs) const {
    return x_ * rhs.x_ + y_ * rhs.y_ + z_ * rhs.z_;
  }

  /** check whether distance between this and that is less than eps */
  [[nodiscard]] bool equals(const Vec3 &that, double eps) {
    return distanceSquared(that) < (eps * eps);
  }

  [[nodiscard]] double distanceSquared(const Vec3 &that) const {
    return std::pow(x_ - that.x_, 2) + std::pow(y_ - that.y_, 2) + std::pow(z_ - that.z_, 2);
  }

  Vec3 &operator+=(const Vec3 &rhs) {
    x_ += rhs.x_;
    y_ += rhs.y_;
    z_ += rhs.z_;
    return *this;
  }

  Vec3 &operator-=(const Vec3 &rhs) {
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    z_ -= rhs.z_;
    return *this;
  }

  Vec3 &operator*=(double rhs) {
    x_ *= rhs;
    y_ *= rhs;
    z_ *= rhs;
    return *this;
  }

  Vec3 &operator*=(const Vec3 &rhs) {
    x_ *= rhs.x_;
    y_ *= rhs.y_;
    z_ *= rhs.z_;
    return *this;
  }

  Vec3 &operator/=(double rhs) {
    x_ /= rhs;
    y_ /= rhs;
    z_ /= rhs;
    return *this;
  }

  Vec3 operator-() const {
    return Vec3(-x_, -y_, -z_);
  }

  Vec3 operator*(double rhs) const {
    return Vec3(x_ * rhs, y_ * rhs, z_ * rhs);
  }

  Vec3 operator/(double rhs) const {
    return Vec3(x_ / rhs, y_ / rhs, z_ / rhs);
  }

  [[nodiscard]]

  std::string toString() const {
    return std::to_string(x_) + ", " + std::to_string(y_) + ", " + std::to_string(z_);
  }

  double norm() const {
    return sqrt(x_ * x_ + y_ * y_ + z_ * z_);
  }

  /** squared norm */
  [[nodiscard]] double norm2() const {
    return x_ * x_ + y_ * y_ + z_ * z_;
  }

  Vec3 &normalize() {
    double n = norm();
    x_ /= n;
    y_ /= n;
    z_ /= n;
    return *this;
  }

  [[nodiscard]] Vec3 normalized() const {
    Vec3 v(*this);
    v.normalize();
    return v;
  }

  [[nodiscard]] Vec3 cross(const Vec3 &rhs) const {
    return Vec3(y_ * rhs.z_ - z_ * rhs.y_,
                z_ * rhs.x_ - x_ * rhs.z_,
                x_ * rhs.y_ - y_ * rhs.x_);
  }

  friend Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs) {
    return Vec3(lhs.x_ + rhs.x_, lhs.y_ + rhs.y_, lhs.z_ + rhs.z_);
  }

  friend Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs) {
    return Vec3(lhs.x_ - rhs.x_, lhs.y_ - rhs.y_, lhs.z_ - rhs.z_);
  }
};

template <>
class fmt::formatter<Vec3> {
public:
  constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
  template <typename Context>
  constexpr auto format (Vec3 const& v, Context& ctx) const {
    return format_to(ctx.out(), "{}, {}, {}", v.x(), v.y(), v.z());
  }
};

#endif //PROJECT_QLC3D_VEC3_H
