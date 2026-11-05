// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Sysfs.h"

#include "Monitoring/Energy/Sensor.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace seissol::monitoring::energy {

namespace {

std::filesystem::path below(const std::string& root, std::string_view path) {
  return std::filesystem::path(root) / std::filesystem::path(path).relative_path();
}

//! @brief The entries of @p directory, sorted, or none if it cannot be listed.
std::vector<std::filesystem::path> entries(const std::filesystem::path& directory) {
  std::vector<std::filesystem::path> result;
  std::error_code error;
  for (const auto& entry : std::filesystem::directory_iterator(directory, error)) {
    result.push_back(entry.path());
  }
  std::sort(result.begin(), result.end());
  return result;
}

std::optional<double> parseNumber(const std::string& text) {
  std::istringstream stream(text);
  double value = 0;
  if (stream >> value) {
    return value;
  }
  return {};
}

class CrayPmSensor : public Sensor {
  public:
  CrayPmSensor(std::filesystem::path path, Part part)
      : Sensor("cray-pm", "cray-pm/" + path.filename().string(), part), path_(std::move(path)) {}

  std::optional<EnergyReading> read() override {
    const auto text = readSysfsFile(path_.string());
    if (!text.has_value()) {
      return {};
    }
    const auto value = parseCrayPm(text.value());
    if (!value.has_value() || value->unit != "J") {
      return {};
    }
    return EnergyReading{accumulator_.update(value->value), value->seconds};
  }

  private:
  std::filesystem::path path_;
  Accumulator accumulator_;
};

//! @brief A file that holds a counter of microjoules.
class MicrojouleSensor : public Sensor {
  public:
  MicrojouleSensor(std::string source,
                   std::string name,
                   Part part,
                   std::filesystem::path path,
                   std::optional<double> rangeMicrojoules)
      : Sensor(std::move(source), std::move(name), part), path_(std::move(path)),
        accumulator_(rangeMicrojoules.has_value()
                         ? std::optional<double>(rangeMicrojoules.value() * 1e-6)
                         : std::nullopt) {}

  std::optional<EnergyReading> read() override {
    const auto text = readSysfsFile(path_.string());
    if (!text.has_value()) {
      return {};
    }
    const auto microjoules = parseNumber(text.value());
    if (!microjoules.has_value()) {
      return {};
    }
    return EnergyReading{accumulator_.update(microjoules.value() * 1e-6), {}};
  }

  private:
  std::filesystem::path path_;
  Accumulator accumulator_;
};

} // namespace

std::optional<std::string> readSysfsFile(const std::string& path) {
  const std::ifstream stream(path);
  if (!stream) {
    return {};
  }
  std::ostringstream content;
  content << stream.rdbuf();
  if (stream.bad()) {
    return {};
  }
  auto text = content.str();
  const auto begin = text.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) {
    return std::string{};
  }
  const auto end = text.find_last_not_of(" \t\r\n");
  return text.substr(begin, end - begin + 1);
}

std::optional<CrayPmValue> parseCrayPm(std::string_view text) {
  std::istringstream stream{std::string(text)};
  CrayPmValue result;
  if (!(stream >> result.value)) {
    return {};
  }
  stream >> result.unit;
  double timestamp = 0;
  std::string timestampUnit;
  if ((stream >> timestamp >> timestampUnit) && timestampUnit == "us") {
    result.seconds = timestamp * 1e-6;
  }
  return result;
}

std::vector<std::unique_ptr<Sensor>> crayPmSensors(const std::string& root) {
  std::vector<std::unique_ptr<Sensor>> sensors;
  for (const auto& path : entries(below(root, "/sys/cray/pm_counters"))) {
    const auto file = path.filename().string();
    if (!file.ends_with("energy")) {
      continue;
    }
    std::optional<Part> part;
    if (file == "energy") {
      part = Part::Node;
    } else if (file == "cpu_energy") {
      part = Part::Cpu;
    } else if (file == "memory_energy") {
      part = Part::Memory;
    } else if (file.starts_with("accel")) {
      part = Part::Accelerators;
    }
    if (!part.has_value()) {
      continue;
    }
    auto sensor = std::make_unique<CrayPmSensor>(path, part.value());
    if (sensor->read().has_value()) {
      sensors.push_back(std::move(sensor));
    }
  }
  return sensors;
}

std::optional<double> crayPmPowerCap(const std::string& root) {
  const auto text = readSysfsFile(below(root, "/sys/cray/pm_counters/power_cap").string());
  if (!text.has_value()) {
    return {};
  }
  const auto value = parseCrayPm(text.value());
  // a cap of zero is none
  if (!value.has_value() || value->value <= 0) {
    return {};
  }
  return value->value;
}

std::vector<std::unique_ptr<Sensor>> powercapSensors(const std::string& root) {
  std::vector<std::unique_ptr<Sensor>> sensors;
  // the zones of intel-rapl-mmio repeat those of intel-rapl; the zones of the parts of a package
  // (intel-rapl:0:0, ...) are listed next to it, and their name says which part they are
  std::vector<std::filesystem::path> zones;
  for (const auto& path : entries(below(root, "/sys/class/powercap"))) {
    const auto zone = path.filename().string();
    if (zone.starts_with("intel-rapl:")) {
      zones.push_back(path);
    }
  }
  for (const auto& zone : zones) {
    const auto name = readSysfsFile((zone / "name").string());
    if (!name.has_value()) {
      continue;
    }
    std::optional<Part> part;
    if (name->starts_with("package")) {
      part = Part::Cpu;
    } else if (name.value() == "dram") {
      part = Part::Memory;
    }
    if (!part.has_value()) {
      continue;
    }
    std::optional<double> range;
    if (const auto text = readSysfsFile((zone / "max_energy_range_uj").string())) {
      range = parseNumber(text.value());
    }
    auto sensor =
        std::make_unique<MicrojouleSensor>("powercap",
                                           "powercap/" + zone.filename().string() + "/" + *name,
                                           part.value(),
                                           zone / "energy_uj",
                                           range);
    // only root may read the zones on most kernels
    if (sensor->read().has_value()) {
      sensors.push_back(std::move(sensor));
    }
  }
  return sensors;
}

std::string_view vendorName(Vendor vendor) {
  switch (vendor) {
  case Vendor::Nvidia:
    return "nvidia";
  case Vendor::Amd:
    return "amd";
  case Vendor::Intel:
    return "intel";
  }
  return "unknown";
}

std::set<Vendor> gpuVendors(const std::string& root) {
  // the PCI classes of GPUs: VGA, 3D and display controllers, and processing accelerators
  constexpr std::array<std::string_view, 4> GpuClasses{"0x0300", "0x0302", "0x0380", "0x1200"};
  std::set<Vendor> vendors;
  for (const auto& device : entries(below(root, "/sys/bus/pci/devices"))) {
    const auto cls = readSysfsFile((device / "class").string());
    const auto vendor = readSysfsFile((device / "vendor").string());
    if (!cls.has_value() || !vendor.has_value()) {
      continue;
    }
    const bool gpu = std::any_of(GpuClasses.begin(),
                                 GpuClasses.end(),
                                 [&](std::string_view prefix) { return cls->starts_with(prefix); });
    if (!gpu) {
      continue;
    }
    if (vendor.value() == "0x10de") {
      vendors.insert(Vendor::Nvidia);
    } else if (vendor.value() == "0x1002") {
      vendors.insert(Vendor::Amd);
    } else if (vendor.value() == "0x8086") {
      vendors.insert(Vendor::Intel);
    }
  }
  return vendors;
}

std::vector<std::unique_ptr<Sensor>> hwmonGpuSensors(const std::string& root) {
  std::vector<std::unique_ptr<Sensor>> sensors;
  for (const auto& hwmon : entries(below(root, "/sys/class/hwmon"))) {
    const auto driver = readSysfsFile((hwmon / "name").string());
    if (!driver.has_value() || (driver.value() != "i915" && driver.value() != "xe")) {
      continue;
    }
    auto sensor = std::make_unique<MicrojouleSensor>("hwmon",
                                                     "hwmon/" + hwmon.filename().string() + "/" +
                                                         driver.value(),
                                                     Part::Devices,
                                                     hwmon / "energy1_input",
                                                     std::nullopt);
    if (sensor->read().has_value()) {
      sensors.push_back(std::move(sensor));
    }
  }
  return sensors;
}

} // namespace seissol::monitoring::energy
