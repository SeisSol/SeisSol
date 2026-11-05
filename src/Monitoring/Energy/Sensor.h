// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_ENERGY_SENSOR_H_
#define SEISSOL_SRC_MONITORING_ENERGY_SENSOR_H_

#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace seissol::monitoring::energy {

//! @brief What of a node a sensor measures.
enum class Part {
  //! the whole node
  Node,
  //! its CPUs
  Cpu,
  //! its memory
  Memory,
  //! its accelerators, as the node counts them
  Accelerators,
  //! its GPUs, as their management library counts them
  Devices
};

//! @brief The name of @p part in the summary: node, cpu, memory, accelerators or devices.
std::string_view partName(Part part);

//! @brief A value of a sensor: the energy it counted since it was opened, in joules, and the time
//! of the value on the clock of the sensor, in seconds, if it has one.
struct EnergyReading {
  double joules{0};
  std::optional<double> seconds;
};

/**
 * @brief The raw values of a hardware energy counter, which may start over at zero, as the energy
 * since the first value.
 *
 * A counter that runs over starts over at the value given as its range; one without a known range
 * that goes back is taken as reset to zero, as by a reload of its driver.
 */
class Accumulator {
  public:
  explicit Accumulator(std::optional<double> range = {}) : range_(range) {}

  //! @brief Takes the next raw value; returns the energy since the first one.
  double update(double raw);

  private:
  std::optional<double> range_;
  std::optional<double> first_;
  double last_{0};
  double offset_{0};
};

/**
 * @brief A source of energy: one file of a node or one GPU.
 */
class Sensor {
  public:
  Sensor(std::string source, std::string name, Part part)
      : source_(std::move(source)), name_(std::move(name)), part_(part) {}
  virtual ~Sensor() = default;

  Sensor(const Sensor&) = delete;
  Sensor(Sensor&&) = delete;
  auto operator=(const Sensor&) = delete;
  auto operator=(Sensor&&) = delete;

  //! @brief What it reads from: cray-pm, powercap, hwmon, nvml, amdsmi, rocmsmi or levelzero.
  [[nodiscard]] const std::string& source() const { return source_; }

  //! @brief What it reads, e.g. cray-pm/accel0_energy or nvml/0000:3b:00.0.
  [[nodiscard]] const std::string& name() const { return name_; }

  [[nodiscard]] Part part() const { return part_; }

  //! @brief The energy since the sensor was opened; none if it cannot be read now.
  virtual std::optional<EnergyReading> read() = 0;

  //! @brief The power the hardware is limited to, in watts, if it says so.
  [[nodiscard]] virtual std::optional<double> powerLimit() const { return {}; }

  private:
  std::string source_;
  std::string name_;
  Part part_;
};

} // namespace seissol::monitoring::energy

#endif // SEISSOL_SRC_MONITORING_ENERGY_SENSOR_H_
