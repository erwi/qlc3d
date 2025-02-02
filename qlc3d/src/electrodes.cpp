#include <electrodes.h>
#include <util/exception.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <geom/vec3.h>
#include <util/logging.h>
#include <eventlist.h>
#include <algorithm>

// SET MAGIC CONSTANT FOR INDICATING UNIFORM ELECTRIC FIELD
const size_t SwitchingInstance::UNIFORM_E_FIELD = numeric_limits<size_t>::max();

Electrode::Electrode(unsigned int electrodeNumber, const std::vector<double> &times, const std::vector<double> &potentials) {
  this->electrodeNumber = electrodeNumber;

  if (times.size() != potentials.size()) {
    throw std::runtime_error(fmt::format("Number of times ({}) and potentials ({}) do not match for electrode {}.", times.size(), potentials.size(), electrodeNumber));
  }

  if (times.empty()) {
    this->times_ = {0};
    this->potentials_ = {0};
  } else if (times[0] > 0) { // make sure we have a default potential and time for start time 0
    Log::warn("No potential defined for electrode {} at time 0. Adding default potential 0.", electrodeNumber);
    this->times_ = {0};
    this->potentials_ = {0};
  }

  if (!std::is_sorted(times.begin(), times.end())) {
    throw std::runtime_error(fmt::format("Times for electrode {} are not sorted.", electrodeNumber));
  }

  for (auto t : times) { this->times_.push_back(t); }
  for (auto p : potentials) { this->potentials_.push_back(p); }
}

double Electrode::getPotentialAtTime(double queryTime) const {
  if (queryTime < 0) {
    RUNTIME_ERROR(fmt::format("Trying to access electrode potential at time {}.", queryTime));
  }

  // find index to first switching time equal or larger than provided time
  auto itr = times_.rbegin();
  int ind = times_.size() - 1;
  while (itr != times_.rend()) {
    if (*itr <= queryTime) {
      return potentials_[ind];
    }
    ind--;
    itr++;
  }
  RUNTIME_ERROR(fmt::format("Trying to access electrode potential at time {} when only potentials for times {} exist.", queryTime, times_));
}

std::vector<Event*> Electrode::createSwitchingEvents() const {
  std::vector<Event*> events;
  for (double time : times_) {
    SwitchingInstance* si = new SwitchingInstance(time, getPotentialAtTime(time), electrodeNumber);

    events.emplace_back(new Event(EVENT_SWITCHING, time, si));
  }
  return events;
}

//===========================================================
//
//		ELECTRODES
//
//

Electrodes::Electrodes() {
    eps_dielectric.push_back(1.0);
}

Electrodes::Electrodes(std::vector<std::shared_ptr<Electrode>> electrodes) {
  this->electricField = nullptr;
  for (auto e : electrodes) {
    electrodeMap[e->getElectrodeNumber()] = e;
  }
}

Electrodes::Electrodes(const Vec3 &electricField) {
  this->electricField = std::make_unique<Vec3>(electricField);
}

Electrodes Electrodes::withConstantElectricField(const Vec3 &electricField) {
  return Electrodes(electricField);
}

Electrodes Electrodes::withElectrodePotentials(std::vector<std::shared_ptr<Electrode>> electrodes) {
  return Electrodes(electrodes);
}

Electrodes Electrodes::withInitialPotentials(const std::vector<unsigned int> &electrodeNumber, const std::vector<double> &potential) {
  std::vector<std::shared_ptr<Electrode>> electrodes;
  for (unsigned int i = 0; i < electrodeNumber.size(); i++) {
    electrodes.emplace_back(std::shared_ptr<Electrode>(new Electrode(electrodeNumber[i], {0}, {potential[i]})));
  }
  return withElectrodePotentials(electrodes);
}

Electrodes Electrodes::withInitialPotentials(const std::vector<unsigned int> &electrodeNumber, double potential) {
  std::vector<double> potentials;
  potentials.resize(electrodeNumber.size(), potential);
  return withInitialPotentials(electrodeNumber, potentials);
}

double Electrodes:: getDielectricPermittivity(int i) const {
#ifndef NDEBUG
  if (i < 0 || i >= (int) eps_dielectric.size()) {
    RUNTIME_ERROR(fmt::format("Trying to access dielectic {} when only {} dielected materials exist.", i,
                              eps_dielectric.size()));
  }
#endif
	return eps_dielectric[i];
}

void Electrodes::setDielectricPermittivities(std::vector<double> eps) {
  eps_dielectric = std::move(eps);
}

bool Electrodes::isPotentialCalculationRequired() const {
/*! Returns whether potential calculations are necessary */
    if (hasElectricField()) { // if have fixed E-field
        return false;
    }
    if (0 == getnElectrodes()) { // if have no electrodes
        return false;
    }
    // Otherwise
    return true;
}

bool Electrodes::hasElectricField() const {
  return electricField != nullptr;
}

Vec3 Electrodes::getElectricField() const {
    if (!hasElectricField()) {
        RUNTIME_ERROR("No uniform E-field has been defined.");
    }
    return {electricField->x(), electricField->y(), electricField->z()};
}

std::unordered_map<unsigned int, double> Electrodes::getCurrentPotentials(double currentTime) const {
  std::unordered_map<unsigned int, double> currentPotentials;

  if (hasElectricField()) { // return empty electrode potentials if a field is defined
    return currentPotentials;
  }

  for (auto &e : electrodeMap) {
    currentPotentials[e.first] = e.second->getPotentialAtTime(currentTime);
  }

  return currentPotentials;
}

std::vector<Event*> Electrodes::createSwitchingEvents() const {
  std::vector<Event*> events;

  for (auto &e : electrodeMap) {
    auto electrodeEvents = e.second->createSwitchingEvents();
    for (auto &ev : electrodeEvents) {
      events.push_back(ev);
    }
  }

  return events;
}