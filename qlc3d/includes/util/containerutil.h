#ifndef PROJECT_QLC3D_CONTAINERUTIL_H
#define PROJECT_QLC3D_CONTAINERUTIL_H
#include <set>
#include <unordered_set>

template <typename T>
bool contains(std::set<T> s, T value) {
  return s.find(value) != s.end();
}

template <typename T>
bool contains(std::unordered_set<T> s, T value) {
  return s.find(value) != s.end();
}

#endif //PROJECT_QLC3D_CONTAINERUTIL_H
