#ifndef PROJECT_QLC3D_PERIODICITY_H
#define PROJECT_QLC3D_PERIODICITY_H
#include <unordered_map>
#include <set>
#include <vector>
#include <list>

class Mesh;
class Coordinates;

class PeriodicityType {
  bool left_right_is_periodic;
  bool front_back_is_periodic;
  bool top_bottom_is_periodic;
public:

  explicit PeriodicityType(const Mesh &triangles);

  [[nodiscard]] bool isAnyPeriodic() const { return isFrontBackPeriodic() || isLeftRightPeriodic() || isTopBottomPeriodic(); };
  [[nodiscard]] bool isFrontBackPeriodic() const { return front_back_is_periodic; };
  [[nodiscard]] bool isLeftRightPeriodic() const { return left_right_is_periodic; };
  [[nodiscard]] bool isTopBottomPeriodic() const { return top_bottom_is_periodic; };
};

class PeriodicNodesMapping {
  PeriodicityType periodicityType;
  std::vector<unsigned int> periNodes_;

  void makePeriEquNodes(const PeriodicityType &periodicityType,
                        const Mesh &e,
                        const Coordinates &coordinates);

  void setFacePeriNodes(std::list<unsigned int> &face1,
                        std::list<unsigned int> &face2,
                        const int &norm,
                        const Coordinates &coordinates);

  void setEdgePeriNodes(std::list<unsigned int> &edge1,
                        std::list<unsigned int> &edge2,
                        const int &dim,
                        const Coordinates &coordinates);
public:
  PeriodicNodesMapping(const Mesh &tris, const Coordinates &coords);
  [[nodiscard]] unsigned int getPeriodicNode(unsigned int node) const; // { return periNodes_[node]; };

  void initialisePeriodicNodes(const Mesh &e, const Coordinates &coordinates);

  const std::vector<unsigned int> &getPeriNodes() const { return periNodes_; };
  const PeriodicityType &getPeriodicityType() const { return periodicityType; };
};


#endif //PROJECT_QLC3D_PERIODICITY_H
