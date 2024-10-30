#ifndef PROJECT_QLC3D_RANDOMGENERATOR_H
#define PROJECT_QLC3D_RANDOMGENERATOR_H
#include <random>
#include <chrono>

class RandomGenerator {
private:
  std::mt19937 gen;
  std::uniform_real_distribution<> dis;

public:
  static unsigned long millisSinceEpoch() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()).count();
  }

  // Constructor that takes a seed value
  RandomGenerator(unsigned int seed)
          : gen(seed), dis(0.0, 1.0) {}

  // Constructor that uses the current time as the seed
  RandomGenerator()
          : gen(RandomGenerator::millisSinceEpoch()), dis(0.0, 1.0) {}

  /* get a random double in the range [0, 1.0] */
  double getRandomDouble() {
    return dis(gen);
  }
};
#endif //PROJECT_QLC3D_RANDOMGENERATOR_H
