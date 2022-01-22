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
#endif // GLOBALS_H
