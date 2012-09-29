#ifndef GLOBALS_H
#define GLOBALS_H
#include <limits>
// GLOBAL DEFINITIONS ARE MADE HERE

typedef unsigned int idx; // ALL INDEX VARIABLES ARE 32 BIT UNSIGNED INTS

static const idx NOT_AN_INDEX = std::numeric_limits<idx>::max();


namespace qlc3d_GLOBALS
{
    extern double GLOBAL_COORD_EPS;
}

#endif // GLOBALS_H

