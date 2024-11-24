#ifndef PROJECT_QLC3D_PERIODICITY_H
#define PROJECT_QLC3D_PERIODICITY_H
#include <unordered_map>
#include <set>

class Mesh;
class Coordinates;

class PeriodicityType {
  bool left_right_is_periodic;
  bool front_back_is_periodic;
  bool top_bottom_is_periodic;
public:

  PeriodicityType(const Mesh &triangles);

  [[nodiscard]] bool isAnyPeriodic() const { return isFrontBackPeriodic() || isLeftRightPeriodic() || isTopBottomPeriodic(); };
  [[nodiscard]] bool isFrontBackPeriodic() const { return front_back_is_periodic; };
  [[nodiscard]] bool isLeftRightPeriodic() const { return left_right_is_periodic; };
  [[nodiscard]] bool isTopBottomPeriodic() const { return top_bottom_is_periodic; };
};

class PeriodicNodesMapping {
  std::unordered_map<unsigned int, unsigned int> mapping;

  void matchNodesFrontBack(const std::set<unsigned int> periNodes, const Coordinates &coords);
  //void mathcNodesFrontBackLeftRight(const set<unsigned int> periNodes, const Coordinates &coords);
  //void matchNodesFrontBackLeftRightTopBottom(const set<unsigned int> periNodes, const Coordinates &coords);

public:
  PeriodicNodesMapping(const Mesh &tris, const Coordinates &coords, const PeriodicityType &periodicity);
  [[nodiscard]] unsigned int getPeriodicNode(unsigned int node) const;


};


#endif //PROJECT_QLC3D_PERIODICITY_H
