// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "EnergyMeter.h"

#include "Common/Json.h"
#include "Monitoring/Counters.h"
#include "Monitoring/Energy/Discovery.h"
#include "Monitoring/Energy/Sensor.h"
#include "Monitoring/Energy/Sysfs.h"
#include "Monitoring/Unit.h"
#include "Parallel/MPI.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <mpi.h>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <utils/logger.h>

namespace seissol::monitoring {

namespace {

constexpr std::array<energy::Part, 5> Parts{energy::Part::Node,
                                            energy::Part::Cpu,
                                            energy::Part::Memory,
                                            energy::Part::Accelerators,
                                            energy::Part::Devices};
constexpr std::array<const char*, 5> CounterNames{
    "nodeEnergy", "cpuEnergy", "memoryEnergy", "acceleratorEnergy", "deviceEnergy"};
constexpr std::array<const char*, 5> RateNames{
    "nodePower", "cpuPower", "memoryPower", "acceleratorPower", "devicePower"};
constexpr std::array<std::string_view, 7> Sources{
    "cray-pm", "powercap", "hwmon", "nvml", "amdsmi", "rocmsmi", "levelzero"};

//! two times of sensors closer than this, in seconds, are the same
constexpr double SameTime = 1e-6;

/**
 * @brief The number, minimum, maximum, mean and sum of the values of the ranks that have one.
 * Collective; the same on all ranks.
 */
json::Value overHolders(std::optional<double> value) {
  MPI_Comm comm = seissol::Mpi::mpi.comm();
  std::array<double, 2> sums{value.has_value() ? 1.0 : 0.0, value.value_or(0.0)};
  // the minimum as the maximum of the negated values
  constexpr double Lowest = std::numeric_limits<double>::lowest();
  std::array<double, 2> maxima{value.has_value() ? -value.value() : Lowest,
                               value.has_value() ? value.value() : Lowest};
  MPI_Allreduce(MPI_IN_PLACE, sums.data(), 2, MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, maxima.data(), 2, MPI_DOUBLE, MPI_MAX, comm);
  json::Value result;
  result["count"] = sums[0];
  if (sums[0] > 0) {
    result["min"] = -maxima[0];
    result["max"] = maxima[1];
    result["mean"] = sums[1] / sums[0];
    result["sum"] = sums[1];
  }
  return result;
}

} // namespace

const char* EnergyMeter::counterName(energy::Part part) {
  for (std::size_t i = 0; i < Parts.size(); ++i) {
    if (Parts[i] == part) {
      return CounterNames[i];
    }
  }
  return "";
}

void EnergyMeter::setup(Counters& counters, bool enabled, const std::string& root) {
  enabled_ = enabled;
  const bool leader = seissol::Mpi::mpi.sharedMemMpiRank() == 0;
  if (enabled && leader) {
    for (auto& sensor : energy::nodeSensors(root)) {
      const auto first = sensor->read();
      sensors_.push_back(Sensor{.sensor = std::move(sensor), .last = first});
    }
    nodePowerCap_ = energy::crayPmPowerCap(root);
  }

  // every rank registers the counters of the parts that any node measures
  std::array<int, Parts.size()> measured{};
  for (std::size_t i = 0; i < Parts.size(); ++i) {
    for (const auto& sensor : sensors_) {
      if (sensor.sensor->part() == Parts[i]) {
        measured[i] = 1;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,
                measured.data(),
                static_cast<int>(measured.size()),
                MPI_INT,
                MPI_MAX,
                seissol::Mpi::mpi.comm());
  std::string names;
  for (std::size_t i = 0; i < Parts.size(); ++i) {
    if (measured[i] == 0) {
      continue;
    }
    const auto part = Parts[i];
    counters.add(CounterInfo{.name = CounterNames[i],
                             .rateName = RateNames[i],
                             .unit = &UnitEnergy,
                             .rateUnit = &UnitPower,
                             .scope = CounterScope::Node,
                             .perStage = true},
                 [this, part] { return read(part); });
    names += std::string(names.empty() ? "" : ", ") + CounterNames[i];
  }

  if (!enabled) {
    logInfo() << "The energy of the run is not measured (SEISSOL_ENERGY=0).";
  } else if (names.empty()) {
    logInfo() << "No node of the run has energy sensors that can be read; the energy of the run is "
                 "not measured.";
  } else {
    logInfo() << "The energy of the run is measured:" << names;
    if (seissol::Mpi::mpi.rank() == 0) {
      std::string own;
      for (const auto& sensor : sensors_) {
        own += std::string(own.empty() ? "" : ", ") + sensor.sensor->name();
      }
      logInfo() << "Energy sensors of the node of rank 0:" << (own.empty() ? "none" : own);
    }
  }
}

std::optional<CounterValue> EnergyMeter::read(energy::Part part) {
  std::optional<CounterValue> value;
  // the clock of the sensors, if they all have the same time, as the files of an update of the HPE
  // Cray counters have
  std::optional<double> seconds;
  bool shared = true;
  for (auto& entry : sensors_) {
    if (entry.sensor->part() != part) {
      continue;
    }
    // a read that fails takes the last value, which the next read catches up on
    if (const auto reading = entry.sensor->read(); reading.has_value()) {
      entry.last = reading;
    }
    if (!entry.last.has_value()) {
      continue;
    }
    if (!value.has_value()) {
      value = CounterValue{};
    }
    value.value().amount += entry.last.value().joules;
    const auto time = entry.last.value().seconds;
    if (!time.has_value() ||
        (seconds.has_value() && std::abs(seconds.value() - time.value()) > SameTime)) {
      shared = false;
    }
    seconds = time;
  }
  if (value.has_value() && shared) {
    value.value().time = seconds;
  }
  return value;
}

json::Value EnergyMeter::describe() const {
  json::Value result;
  result["enabled"] = enabled_;
  auto& parts = result["parts"] = json::Value::object();
  for (std::size_t i = 0; i < Parts.size(); ++i) {
    unsigned sources = 0;
    std::size_t count = 0;
    std::optional<double> limit;
    for (const auto& sensor : sensors_) {
      if (sensor.sensor->part() != Parts[i]) {
        continue;
      }
      ++count;
      for (std::size_t source = 0; source < Sources.size(); ++source) {
        if (sensor.sensor->source() == Sources[source]) {
          sources |= 1U << source;
        }
      }
      if (const auto watts = sensor.sensor->powerLimit(); watts.has_value()) {
        limit = limit.value_or(0.0) + watts.value();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &sources, 1, MPI_UNSIGNED, MPI_BOR, seissol::Mpi::mpi.comm());
    auto perNode =
        overHolders(count > 0 ? std::optional<double>(static_cast<double>(count)) : std::nullopt);
    auto limits = overHolders(limit);
    if (sources == 0) {
      continue;
    }
    auto& entry = parts[energy::partName(Parts[i])];
    entry["counter"] = CounterNames[i];
    auto& names = entry["sources"] = json::Value::array();
    for (std::size_t source = 0; source < Sources.size(); ++source) {
      if ((sources & (1U << source)) != 0) {
        names.push(Sources[source]);
      }
    }
    entry["nodes"] = perNode["count"];
    entry["sensorsPerNode"] = std::move(perNode);
    if (limits.find("sum") != nullptr) {
      entry["powerLimitPerNode"] = std::move(limits);
    }
  }
  auto cap = overHolders(nodePowerCap_);
  if (cap.find("sum") != nullptr) {
    result["nodePowerCap"] = std::move(cap);
  }
  return result;
}

} // namespace seissol::monitoring
