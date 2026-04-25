#ifndef GLOBALS_H
#define GLOBALS_H
#include <limits>
// GLOBAL DEFINITIONS ARE MADE HERE
#ifndef PI
#define PI 3.14159265358979323846264338327950288419716939937510
#endif

typedef unsigned int idx; // ALL INDEX VARIABLES ARE 32 BIT UNSIGNED INTS
static const idx NOT_AN_INDEX = std::numeric_limits<idx>::max();

namespace qlc3d_GLOBALS {
    static const double GLOBAL_COORD_EPS = 1e-10;
    static const int ERROR_CODE_BAD_SETTINGS_FILE = 1;
}

/**
 * @brief Shared conversion factors between micron-based mesh data and SI units.
 */
namespace qlc3d::units {
  inline constexpr double MICROMETER_TO_METER = 1e-6;
  inline constexpr double SQUARE_MICROMETER_TO_SQUARE_METER = MICROMETER_TO_METER * MICROMETER_TO_METER;
  inline constexpr double CUBIC_MICROMETER_TO_CUBIC_METER = SQUARE_MICROMETER_TO_SQUARE_METER * MICROMETER_TO_METER;

  inline constexpr double METER_TO_MICROMETER = 1.0 / MICROMETER_TO_METER;
  inline constexpr double SQUARE_METER_TO_SQUARE_MICROMETER = 1.0 / SQUARE_MICROMETER_TO_SQUARE_METER;
  inline constexpr double CUBIC_METER_TO_CUBIC_MICROMETER = 1.0 / CUBIC_MICROMETER_TO_CUBIC_METER;
}
#endif // GLOBALS_H
