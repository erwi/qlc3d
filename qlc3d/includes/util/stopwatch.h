#ifndef PROJECT_QLC3D_STOPWATCH_H
#define PROJECT_QLC3D_STOPWATCH_H
#include <chrono>

class Stopwatch {
  bool running;
  long long elapsed_time; // in milliseconds
  std::chrono::steady_clock::time_point start_time;

public:
  Stopwatch() : running(false), elapsed_time(0) {}

  void start() {
    if (!running) {
      start_time = std::chrono::steady_clock::now();
      running = true;
    }
  }

  void stop() {
    if (running) {
      auto end_time = std::chrono::steady_clock::now();
      elapsed_time += std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
      running = false;
    }
  }

  void reset() {
    running = false;
    elapsed_time = 0;
  }

  [[nodiscard]] double elapsedSeconds() const {
    if (running) {
      auto current_time = std::chrono::steady_clock::now();
      return (elapsed_time + std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count()) / 1000.0;
    }
    return elapsed_time / 1000.0;
  }
};
#endif //PROJECT_QLC3D_STOPWATCH_H
