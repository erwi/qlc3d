// Define GiD material numbers
#include <unordered_map>
#include <string>
#ifndef MATERIAL_NUMBERS_H
#define MATERIAL_NUMBERS_H

#define MAT_INVALID 0
#define MAT_DOMAIN1 4       // bit 3
#define MAT_DOMAIN2 8       // 4
#define MAT_DOMAIN3 12
#define MAT_DOMAIN4 16
#define MAT_DOMAIN5 20
#define MAT_DOMAIN6 24
#define MAT_DOMAIN7 28

#define MAT_DIELECTRIC1 32  // 4
#define MAT_DIELECTRIC2 36  // 5
#define MAT_DIELECTRIC3 40
#define MAT_DIELECTRIC4 44
#define MAT_DIELECTRIC5 48
#define MAT_DIELECTRIC6 52
#define MAT_DIELECTRIC7 56

#define MAT_ELECTRODE1  64  // 6
#define MAT_ELECTRODE2  128 // 7
#define MAT_ELECTRODE3  192 // 6 and 7
#define MAT_ELECTRODE4  256 // 8
#define MAT_ELECTRODE5  320 // 6 and 8
#define MAT_ELECTRODE6  384 // 7 and 8
#define MAT_ELECTRODE7  448 // 6 and 7 and 8
#define MAT_ELECTRODE8  512 // 9
#define MAT_ELECTRODE9  576 // 6 and 9
static const unsigned int MAT_MAX_ELECTRODES_COUNT = 9;
// UP TO MAT_ELECTRODE16 ?
#define MAT_FIXLC1 2048 // 11
#define MAT_FIXLC2 4096 // 12
#define MAT_FIXLC3 6144 // 11 and 12
#define MAT_FIXLC4 8192 // 13
#define MAT_FIXLC5 10240
#define MAT_FIXLC6 12288
#define MAT_FIXLC7 14336
#define MAT_FIXLC8 16384 // 14
#define MAT_FIXLC9 18432 // 14 and 11

#define MAT_PERIODIC 3
#define MAT_NEUMANN 2   // 0x
namespace qlc3d {
    static const std::unordered_map<std::string, int> MATERIAL_NUMBER_BY_NAME = {
            {"domain1",     MAT_DOMAIN1},
            {"periodic",    MAT_PERIODIC},
            {"neumann",     MAT_NEUMANN},
            // dielectric volumes 1 - 7
            {"dielectric1", MAT_DIELECTRIC1},
            {"dielectric2", MAT_DIELECTRIC2},
            {"dielectric3", MAT_DIELECTRIC3},
            {"dielectric4", MAT_DIELECTRIC4},
            {"dielectric5", MAT_DIELECTRIC5},
            {"dielectric6", MAT_DIELECTRIC6},
            {"dielectric7", MAT_DIELECTRIC7},
            // electrodes 1 - 9
            {"electrode1",  MAT_ELECTRODE1},
            {"electrode2",  MAT_ELECTRODE2},
            {"electrode3",  MAT_ELECTRODE3},
            {"electrode4",  MAT_ELECTRODE4},
            {"electrode5",  MAT_ELECTRODE5},
            {"electrode6",  MAT_ELECTRODE6},
            {"electrode7",  MAT_ELECTRODE7},
            {"electrode8",  MAT_ELECTRODE8},
            {"electrode9",  MAT_ELECTRODE9},
            // fixlc alignment layers 1 - 9
            {"fixlc1",      MAT_FIXLC1},
            {"fixlc2",      MAT_FIXLC2},
            {"fixlc3",      MAT_FIXLC3},
            {"fixlc4",      MAT_FIXLC4},
            {"fixlc5",      MAT_FIXLC5},
            {"fixlc6",      MAT_FIXLC6},
            {"fixlc7",      MAT_FIXLC7},
            {"fixlc8",      MAT_FIXLC8},
            {"fixlc9",      MAT_FIXLC9},
    };
}
inline
int FIXLCN_TO_MATNUM(const int &n) {
    /*! converts fixlc number to material number. e.g. 1 -> 2048 etc.*/
    return n * MAT_FIXLC1;
}

inline
size_t MATNUM_TO_ELECTRODE_NUMBER(const size_t &mat) {
/*!
 * Converts a material number to Electrode number, or 0 if
 * input material number does not corrspond any electrode.
 *
 * e.g.:
 *      32 -> 1
 *      36 -> 2
 *      2080 -> 1  (FIXLC1 + ELECTRODE1)
 *      8768 -> 9  (FIXLC4 + ELECTRODE9)
 */
    // CREATE MASK WITH BITS 6,7,8,9 SET
    size_t mask = 64 | 128 | 256 | 512;
    // ONLY KEEP ELECTRODE BITS OF INPUT VARIABLE
    size_t eleBits = mat & mask;
    // GET ELECTRODE INDEX NUMBER
    size_t eleNum = eleBits / MAT_ELECTRODE1;
    return eleNum;
}

inline
size_t MATNUM_TO_FIXLC_NUMBER(const size_t &mat) {
    // RETURNS INDEX NUMBER OF A FIXLC SURFACE
    // CRAETE MASK WITH BITS 11 -> 14 SET
    size_t mask = MAT_FIXLC1 | MAT_FIXLC2 | MAT_FIXLC4 | MAT_FIXLC8;
    size_t fixlcBits = mat & mask;
    size_t fixlcNum = fixlcBits / MAT_FIXLC1;
    return fixlcNum;
}
#endif
