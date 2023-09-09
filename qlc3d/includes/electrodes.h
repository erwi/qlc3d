#ifndef ELECTRODES_H
#define ELECTRODES_H
#include <stdlib.h>
#include <utility>
#include <vector>
#include <unordered_map>
#include <list>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <eventlist.h>
#include <limits>
#include <memory>

class Vec3;
class Event;

class SwitchingInstance {
    /*!A SwitchingInstance contains information needed to describe a single
    electrode switching event*/
public:
    SwitchingInstance(const double t, const double pot, const size_t en):
        time(t), potential(pot), electrodeNumber(en) {
    }
    static const size_t UNIFORM_E_FIELD;
    double time;
    double potential;
    size_t electrodeNumber;
};

class Electrode {
  unsigned int electrodeNumber;
  std::vector<double> times;
  std::vector<double> potentials;
public:
  Electrode(unsigned int electrodeNumber, std::vector<double> times, std::vector<double> potentials);
  //:
  //  electrodeNumber(electrodeNumber),
  //  times(std::move(times)),
  //  potentials(std::move(potentials)) { };

  [[nodiscard]] unsigned int getElectrodeNumber() const { return electrodeNumber; }
  [[nodiscard]] double getPotentialAtTime(double queryTime) const;
  [[nodiscard]] std::vector<Event*> createSwitchingEvents() const;
};


/*!A class with parameters related to potential calculations and electrodes*/
class Electrodes {
private:
  std::vector<double> currentElectrodePotentials;    // keeps current potential values for each electrode

  std::unordered_map<unsigned int, std::shared_ptr<Electrode>> electrodeMap;
  std::shared_ptr<Vec3> electricField; // optional uniform electric field
  std::vector <double> eps_dielectric; // relative dielectric permittivity of dielectric regions
public:

    Electrodes();
    Electrodes(std::vector<std::shared_ptr<Electrode>> electrodes, std::shared_ptr<Vec3> electricField = nullptr) {
      this->electricField = electricField;
      for (auto e : electrodes) {
        electrodeMap[e->getElectrodeNumber()] = e;
      }
    }

    static std::shared_ptr<Electrodes> withInitialPotentials(std::vector<unsigned int> electrodeNumber, std::vector<double> potential);

    /**
     * return electrode with given number. This is Not 0 based, i.e. electrode 1 is the first one, as de
     * fined in material numbers
     */
    [[nodiscard]] std::shared_ptr<Electrode> getElectrode(unsigned int electrode) const;

    double getDielectricPermittivity(int i) const;  // gets relative dielectric permittivity of dielectric#i
    bool getCalcPot() const;
    bool hasElectricField() const;              // returns true if uniform E-field has been defined
    [[nodiscard]] Vec3 getElectricField() const;
    [[nodiscard]] std::unordered_map<unsigned int, double> getCurrentPotentials(double currentTime) const;
    /** Create list of all electrode switching events */
    [[nodiscard]] std::vector<Event*> createSwitchingEvents() const;
    [[nodiscard]] size_t getnElectrodes() const { return electrodeMap.size(); }
};
#endif
