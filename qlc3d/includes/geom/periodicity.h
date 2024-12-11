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

  void makePeriEquNodes(const Mesh &e,
                        const Coordinates &coordinates);

  void setFacePeriNodes(const std::vector<unsigned int> &face1,
                        const std::vector<unsigned int> &face2,
                        int norm,
                        const Coordinates &coordinates);

  void setEdgePeriNodes(const std::vector<unsigned int> &edge1,
                        const std::vector<unsigned int> &edge2,
                        int dim,
                        const Coordinates &coordinates);
public:
  PeriodicNodesMapping(const Mesh &tris, const Coordinates &coords);
  [[nodiscard]] const PeriodicityType &getPeriodicityType() const { return periodicityType; };
  [[nodiscard]] const std::vector<unsigned int> &getMapping() const { return periNodes_; }
};


#endif //PROJECT_QLC3D_PERIODICITY_H
