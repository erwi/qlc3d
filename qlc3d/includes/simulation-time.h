#ifndef PROJECT_QLC3D_SIMULATION_TIME_H
#define PROJECT_QLC3D_SIMULATION_TIME_H

#include <cmath>
#include <fmt/format.h>

class SimulationTime {
  static constexpr double RESOLUTION = 1e-15;
  double time;

public:
  explicit SimulationTime(double time) : time(time) {}
  SimulationTime(const SimulationTime &t) = default;

  SimulationTime& increment(double dt) { time += dt; return *this; }
  SimulationTime& setTime(double t) { time = t; return *this; }
  SimulationTime& setTime(const SimulationTime &t) { time = t.time; return *this; }
  /** TODO: we shouldn't need this? */
  [[nodiscard]] double getTime() const { return time; }

  [[nodiscard]] bool equals(double t) const { return std::abs(time - t) < RESOLUTION; }
  [[nodiscard]] bool greaterThan(double t) const { return time - t > RESOLUTION; }
  [[nodiscard]] bool lessThan(double t) const { return t - time > RESOLUTION; }

  [[nodiscard]] bool equals(const SimulationTime &t) const { return std::abs(time - t.time) < RESOLUTION; }
  [[nodiscard]] bool greaterThan(const SimulationTime &t) const { return time - t.time > RESOLUTION; }
  [[nodiscard]] bool lessThan(const SimulationTime &t) const { return t.time - time > RESOLUTION; }
};

template <>
class fmt::formatter<SimulationTime> {
public:
  constexpr auto parse (format_parse_context& ctx) { return ctx.begin(); }
  template <typename Context>
  constexpr auto format (SimulationTime const& t, Context& ctx) const {
    return format_to(ctx.out(), "time={}", t.getTime());
  }
};

#endif //PROJECT_QLC3D_SIMULATION_TIME_H
