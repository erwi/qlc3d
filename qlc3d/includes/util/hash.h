#ifndef PROJECT_QLC3D_HASH_H
#define PROJECT_QLC3D_HASH_H
#include <cstdint>
#include <geom/vec3.h>
// This contains some hashing functions. Mainly for use in testing/debugging to detect when an array of values
// has changed.
/**
 * Hash of a 32 bit value
 */
template<typename T32>
inline int64_t hashCode32(T32 v) {
  return 31 * (*(int32_t*)&v);
}

/**
 * Hash of a 64 bit value
 */
template<typename T64>
inline int64_t hashCode64(T64 v) {
  return 31 * (*(int64_t*)&v);
}

/**
 * Hash of a range of 64 bit values
 */
template<typename T64>
inline int64_t hashCode64(const T64* begin, const T64* end) {
  int64_t hash = 0;
  int64_t count = 1;
  for(auto it = begin; it != end; ++it) {
    hash += count * hashCode64(*it);
    count ++;
  }
  hash += 31 * count; // all-0 arrays of different length should have different hash codes
  return hash;
}

/**
 * Hash of a range of 32 bit values
 */
template<typename T32>
inline int64_t hashCode32(const T32* begin, const T32* end) {
  int64_t hash = 0;
  int64_t count = 1;
  for(auto it = begin; it != end; ++it) {
    hash += count * hashCode32(*it);
    count ++;
  }
  hash += 31 * count; // all-0 arrays of different length should have different hash codes
  return hash;
}

#endif
