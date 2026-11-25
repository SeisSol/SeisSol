// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_ENERGYMETER_H_
#define SEISSOL_SRC_MONITORING_ENERGYMETER_H_

#include "Common/Json.h"
#include "Monitoring/Counters.h"
#include "Monitoring/Energy/Sensor.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace seissol::monitoring {

/**
 * @brief The energy that the nodes of the run take, from the sensors they have.
 *
 * The first rank of every node opens the sensors of the node (see energy::nodeSensors); for every
 * part of a node that any node measures, a counter of the run sums its sensors: nodeEnergy,
 * cpuEnergy, memoryEnergy, acceleratorEnergy and deviceEnergy, in joules, with the rates nodePower,
 * cpuPower, memoryPower, acceleratorPower and devicePower, in watts. The accelerators are the GPUs
 * as the node counts them, the devices the GPUs as their management library counts them: two
 * measurements of the same, which are not to be added. The counters are given for every stage of
 * the run.
 */
class EnergyMeter {
  public:
  /**
   * @brief Opens the sensors of the nodes and registers the counters with @p counters. Collective.
   *
   * With @p enabled false, no sensor is opened and no counter registered; @p root is the directory
   * that /sys is below, which is / but in the tests.
   */
  void setup(Counters& counters, bool enabled, const std::string& root = "/");

  /**
   * @brief For the summary: for each part measured, its counter, the sources that measure it and
   * on how many nodes; the power limits of the nodes and of their GPUs. Collective; complete on
   * rank 0.
   */
  [[nodiscard]] json::Value describe() const;

  //! @brief The name of the counter of the energy of @p part, e.g. cpuEnergy.
  static const char* counterName(energy::Part part);

  private:
  struct Sensor {
    std::unique_ptr<energy::Sensor> sensor;
    //! what it read last, for a read that fails
    std::optional<energy::EnergyReading> last;
  };

  //! @brief The sum of the sensors of @p part on this rank, or none if it has none.
  std::optional<CounterValue> read(energy::Part part);

  std::vector<Sensor> sensors_;
  bool enabled_{false};
  std::optional<double> nodePowerCap_;
};

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_ENERGYMETER_H_
