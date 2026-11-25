// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Vendors.h"

#include "Monitoring/Energy/Library.h"
#include "Monitoring/Energy/Sensor.h"

#include <cstdint>
#include <iomanip>
#include <ios>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// NOLINTBEGIN(readability-identifier-naming,modernize-avoid-c-arrays)
// the declarations of the libraries, with C linkage, as their headers have them; this file alone
// uses them
namespace seissol::monitoring::energy {
namespace {
extern "C" {
namespace nvml {
using Return = int;
using Device = void*;
struct PciInfo {
  char busIdLegacy[16];
  unsigned int domain;
  unsigned int bus;
  unsigned int device;
  unsigned int pciDeviceId;
  unsigned int pciSubSystemId;
  char busId[32];
};
using Init = Return();
using DeviceGetCount = Return(unsigned int*);
using DeviceGetHandleByIndex = Return(unsigned int, Device*);
using DeviceGetTotalEnergyConsumption = Return(Device, unsigned long long*);
using DeviceGetEnforcedPowerLimit = Return(Device, unsigned int*);
using DeviceGetPciInfo = Return(Device, PciInfo*);
} // namespace nvml

namespace amdsmi {
using Status = int;
using Handle = void*;
struct PowerCapInfo {
  std::uint64_t powerCap;
  std::uint64_t defaultPowerCap;
  std::uint64_t dpmCap;
  std::uint64_t minPowerCap;
  std::uint64_t maxPowerCap;
  std::uint64_t reserved[3];
};
using Init = Status(std::uint64_t);
using GetSocketHandles = Status(std::uint32_t*, Handle*);
using GetProcessorHandles = Status(Handle, std::uint32_t*, Handle*);
// the address is a union of bit fields in 64 bits (amdsmi_bdf_t)
using GetGpuDeviceBdf = Status(Handle, std::uint64_t*);
using GetEnergyCount = Status(Handle, std::uint64_t*, float*, std::uint64_t*);
using GetPowerCapInfo = Status(Handle, std::uint32_t, PowerCapInfo*);
} // namespace amdsmi

namespace rsmi {
using Status = int;
using Init = Status(std::uint64_t);
using NumMonitorDevices = Status(std::uint32_t*);
using DevPciIdGet = Status(std::uint32_t, std::uint64_t*);
using DevEnergyCountGet = Status(std::uint32_t, std::uint64_t*, float*, std::uint64_t*);
using DevPowerCapGet = Status(std::uint32_t, std::uint32_t, std::uint64_t*);
} // namespace rsmi

namespace zes {
using Result = int;
using Handle = void*;
struct PciAddress {
  std::uint32_t domain;
  std::uint32_t bus;
  std::uint32_t device;
  std::uint32_t function;
};
struct PciSpeed {
  std::int32_t gen;
  std::int32_t width;
  std::int64_t maxBandwidth;
};
struct PciProperties {
  int stype;
  void* pNext;
  PciAddress address;
  PciSpeed maxSpeed;
  std::uint8_t haveBandwidthCounters;
  std::uint8_t havePacketCounters;
  std::uint8_t haveReplayCounters;
};
struct PowerProperties {
  int stype;
  void* pNext;
  std::uint8_t onSubdevice;
  std::uint32_t subdeviceId;
  std::uint8_t canControl;
  std::uint8_t isEnergyThresholdSupported;
  std::int32_t defaultLimit;
  std::int32_t minLimit;
  std::int32_t maxLimit;
};
struct EnergyCounter {
  std::uint64_t energy;
  std::uint64_t timestamp;
};
using Init = Result(std::uint32_t);
using DriverGet = Result(std::uint32_t*, Handle*);
using DeviceGet = Result(Handle, std::uint32_t*, Handle*);
using DevicePciGetProperties = Result(Handle, PciProperties*);
using DeviceGetCardPowerDomain = Result(Handle, Handle*);
using DeviceEnumPowerDomains = Result(Handle, std::uint32_t*, Handle*);
using PowerGetProperties = Result(Handle, PowerProperties*);
using PowerGetEnergyCounter = Result(Handle, EnergyCounter*);
} // namespace zes
}
} // namespace
} // namespace seissol::monitoring::energy
// NOLINTEND(readability-identifier-naming,modernize-avoid-c-arrays)

namespace seissol::monitoring::energy {

namespace {

// the values of success, and the constants used, of the libraries
constexpr int Success = 0;
constexpr std::uint64_t AmdSmiInitAmdGpus = 1U << 1U;
constexpr int ZesStructurePciProperties = 0x2;
constexpr int ZesStructurePowerProperties = 0xd;

std::string pciAddress(std::uint64_t domain,
                       std::uint64_t bus,
                       std::uint64_t device,
                       std::uint64_t function) {
  std::ostringstream stream;
  stream << std::hex << std::setfill('0') << std::setw(4) << domain << ":" << std::setw(2) << bus
         << ":" << std::setw(2) << device << "." << function;
  return stream.str();
}

class NvmlSensor : public Sensor {
  public:
  NvmlSensor(std::string name,
             nvml::Device device,
             nvml::DeviceGetTotalEnergyConsumption* energy,
             std::optional<double> limit)
      : Sensor("nvml", std::move(name), Part::Devices), device_(device), energy_(energy),
        limit_(limit) {}

  std::optional<EnergyReading> read() override {
    unsigned long long millijoules = 0;
    if (energy_(device_, &millijoules) != Success) {
      return {};
    }
    return EnergyReading{.joules = accumulator_.update(static_cast<double>(millijoules) * 1e-3),
                         .seconds = {}};
  }

  [[nodiscard]] std::optional<double> powerLimit() const override { return limit_; }

  private:
  nvml::Device device_;
  nvml::DeviceGetTotalEnergyConsumption* energy_;
  std::optional<double> limit_;
  Accumulator accumulator_;
};

class AmdSmiSensor : public Sensor {
  public:
  AmdSmiSensor(std::string name,
               amdsmi::Handle handle,
               amdsmi::GetEnergyCount* energy,
               std::optional<double> limit)
      : Sensor("amdsmi", std::move(name), Part::Devices), handle_(handle), energy_(energy),
        limit_(limit) {}

  std::optional<EnergyReading> read() override {
    std::uint64_t accumulator = 0;
    float resolution = 0;
    std::uint64_t timestamp = 0;
    if (energy_(handle_, &accumulator, &resolution, &timestamp) != Success) {
      return {};
    }
    const double joules = static_cast<double>(accumulator) * static_cast<double>(resolution) * 1e-6;
    return EnergyReading{.joules = accumulator_.update(joules),
                         .seconds = static_cast<double>(timestamp) * 1e-9};
  }

  [[nodiscard]] std::optional<double> powerLimit() const override { return limit_; }

  private:
  amdsmi::Handle handle_;
  amdsmi::GetEnergyCount* energy_;
  std::optional<double> limit_;
  Accumulator accumulator_;
};

class RocmSmiSensor : public Sensor {
  public:
  RocmSmiSensor(std::string name,
                std::uint32_t index,
                rsmi::DevEnergyCountGet* energy,
                std::optional<double> limit)
      : Sensor("rocmsmi", std::move(name), Part::Devices), index_(index), energy_(energy),
        limit_(limit) {}

  std::optional<EnergyReading> read() override {
    std::uint64_t accumulator = 0;
    float resolution = 0;
    std::uint64_t timestamp = 0;
    if (energy_(index_, &accumulator, &resolution, &timestamp) != Success) {
      return {};
    }
    const double joules = static_cast<double>(accumulator) * static_cast<double>(resolution) * 1e-6;
    return EnergyReading{.joules = accumulator_.update(joules),
                         .seconds = static_cast<double>(timestamp) * 1e-9};
  }

  [[nodiscard]] std::optional<double> powerLimit() const override { return limit_; }

  private:
  std::uint32_t index_;
  rsmi::DevEnergyCountGet* energy_;
  std::optional<double> limit_;
  Accumulator accumulator_;
};

class LevelZeroSensor : public Sensor {
  public:
  LevelZeroSensor(std::string name, zes::Handle domain, zes::PowerGetEnergyCounter* energy)
      : Sensor("levelzero", std::move(name), Part::Devices), domain_(domain), energy_(energy) {}

  std::optional<EnergyReading> read() override {
    zes::EnergyCounter counter{};
    if (energy_(domain_, &counter) != Success) {
      return {};
    }
    return EnergyReading{.joules = accumulator_.update(static_cast<double>(counter.energy) * 1e-6),
                         .seconds = static_cast<double>(counter.timestamp) * 1e-6};
  }

  private:
  zes::Handle domain_;
  zes::PowerGetEnergyCounter* energy_;
  Accumulator accumulator_;
};

//! @brief Keeps the sensor if it can be read now.
void keep(std::vector<std::unique_ptr<Sensor>>& sensors, std::unique_ptr<Sensor> sensor) {
  if (sensor->read().has_value()) {
    sensors.push_back(std::move(sensor));
  }
}

} // namespace

std::vector<std::unique_ptr<Sensor>> nvmlSensors() {
  static const auto Nvml = Library::open({"libnvidia-ml.so.1", "libnvidia-ml.so"});
  std::vector<std::unique_ptr<Sensor>> sensors;
  if (!Nvml.has_value()) {
    return sensors;
  }
  auto* init = Nvml->function<nvml::Init>("nvmlInit_v2");
  auto* count = Nvml->function<nvml::DeviceGetCount>("nvmlDeviceGetCount_v2");
  auto* handle = Nvml->function<nvml::DeviceGetHandleByIndex>("nvmlDeviceGetHandleByIndex_v2");
  auto* energy =
      Nvml->function<nvml::DeviceGetTotalEnergyConsumption>("nvmlDeviceGetTotalEnergyConsumption");
  auto* limit =
      Nvml->function<nvml::DeviceGetEnforcedPowerLimit>("nvmlDeviceGetEnforcedPowerLimit");
  auto* pci = Nvml->function<nvml::DeviceGetPciInfo>("nvmlDeviceGetPciInfo_v3");
  if (init == nullptr || count == nullptr || handle == nullptr || energy == nullptr ||
      init() != Success) {
    return sensors;
  }
  unsigned int devices = 0;
  if (count(&devices) != Success) {
    return sensors;
  }
  for (unsigned int index = 0; index < devices; ++index) {
    nvml::Device device = nullptr;
    if (handle(index, &device) != Success) {
      continue;
    }
    std::string name = "nvml/" + std::to_string(index);
    nvml::PciInfo info{};
    if (pci != nullptr && pci(device, &info) == Success) {
      name = "nvml/" + pciAddress(info.domain, info.bus, info.device, 0);
    }
    std::optional<double> watts;
    unsigned int milliwatts = 0;
    if (limit != nullptr && limit(device, &milliwatts) == Success) {
      watts = milliwatts * 1e-3;
    }
    keep(sensors, std::make_unique<NvmlSensor>(name, device, energy, watts));
  }
  return sensors;
}

std::vector<std::unique_ptr<Sensor>> amdSmiSensors() {
  static const auto AmdSmi = Library::open(rocmLibraryNames("libamd_smi.so"));
  std::vector<std::unique_ptr<Sensor>> sensors;
  if (!AmdSmi.has_value()) {
    return sensors;
  }
  auto* init = AmdSmi->function<amdsmi::Init>("amdsmi_init");
  auto* sockets = AmdSmi->function<amdsmi::GetSocketHandles>("amdsmi_get_socket_handles");
  auto* processors = AmdSmi->function<amdsmi::GetProcessorHandles>("amdsmi_get_processor_handles");
  auto* bdf = AmdSmi->function<amdsmi::GetGpuDeviceBdf>("amdsmi_get_gpu_device_bdf");
  auto* energy = AmdSmi->function<amdsmi::GetEnergyCount>("amdsmi_get_energy_count");
  auto* cap = AmdSmi->function<amdsmi::GetPowerCapInfo>("amdsmi_get_power_cap_info");
  if (init == nullptr || sockets == nullptr || processors == nullptr || energy == nullptr ||
      init(AmdSmiInitAmdGpus) != Success) {
    return sensors;
  }
  std::uint32_t socketCount = 0;
  if (sockets(&socketCount, nullptr) != Success) {
    return sensors;
  }
  std::vector<amdsmi::Handle> socketHandles(socketCount);
  if (socketCount == 0 || sockets(&socketCount, socketHandles.data()) != Success) {
    return sensors;
  }
  for (std::uint32_t socket = 0; socket < socketCount; ++socket) {
    std::uint32_t processorCount = 0;
    if (processors(socketHandles[socket], &processorCount, nullptr) != Success ||
        processorCount == 0) {
      continue;
    }
    std::vector<amdsmi::Handle> handles(processorCount);
    if (processors(socketHandles[socket], &processorCount, handles.data()) != Success) {
      continue;
    }
    for (std::uint32_t index = 0; index < processorCount; ++index) {
      std::string name = "amdsmi/" + std::to_string(socket) + "-" + std::to_string(index);
      std::uint64_t address = 0;
      if (bdf != nullptr && bdf(handles[index], &address) == Success) {
        // function: bits 0 to 2; device: 3 to 7; bus: 8 to 15; domain: 16 to 63
        name = "amdsmi/" + pciAddress(address >> 16U,
                                      (address >> 8U) & 0xffU,
                                      (address >> 3U) & 0x1fU,
                                      address & 0x7U);
      }
      std::optional<double> watts;
      amdsmi::PowerCapInfo info{};
      if (cap != nullptr && cap(handles[index], 0, &info) == Success && info.powerCap > 0) {
        watts = static_cast<double>(info.powerCap) * 1e-6;
      }
      keep(sensors, std::make_unique<AmdSmiSensor>(name, handles[index], energy, watts));
    }
  }
  return sensors;
}

std::vector<std::unique_ptr<Sensor>> rocmSmiSensors() {
  std::vector<std::string> names;
  for (const auto* version : {"librocm_smi64.so", "librocm_smi64.so.7", "librocm_smi64.so.6"}) {
    for (auto& name : rocmLibraryNames(version)) {
      names.push_back(std::move(name));
    }
  }
  static const auto RocmSmi = Library::open(names);
  std::vector<std::unique_ptr<Sensor>> sensors;
  if (!RocmSmi.has_value()) {
    return sensors;
  }
  auto* init = RocmSmi->function<rsmi::Init>("rsmi_init");
  auto* count = RocmSmi->function<rsmi::NumMonitorDevices>("rsmi_num_monitor_devices");
  auto* pci = RocmSmi->function<rsmi::DevPciIdGet>("rsmi_dev_pci_id_get");
  auto* energy = RocmSmi->function<rsmi::DevEnergyCountGet>("rsmi_dev_energy_count_get");
  auto* cap = RocmSmi->function<rsmi::DevPowerCapGet>("rsmi_dev_power_cap_get");
  if (init == nullptr || count == nullptr || energy == nullptr || init(0) != Success) {
    return sensors;
  }
  std::uint32_t devices = 0;
  if (count(&devices) != Success) {
    return sensors;
  }
  for (std::uint32_t index = 0; index < devices; ++index) {
    std::string name = "rocmsmi/" + std::to_string(index);
    std::uint64_t address = 0;
    if (pci != nullptr && pci(index, &address) == Success) {
      // domain: bits 32 to 63; bus: 8 to 15; device: 3 to 7; function: 0 to 2 (28 to 31 hold the
      // partition, which is no part of the address)
      name = "rocmsmi/" +
             pciAddress(
                 address >> 32U, (address >> 8U) & 0xffU, (address >> 3U) & 0x1fU, address & 0x7U);
    }
    std::optional<double> watts;
    std::uint64_t microwatts = 0;
    if (cap != nullptr && cap(index, 0, &microwatts) == Success && microwatts > 0) {
      watts = static_cast<double>(microwatts) * 1e-6;
    }
    keep(sensors, std::make_unique<RocmSmiSensor>(name, index, energy, watts));
  }
  return sensors;
}

std::vector<std::unique_ptr<Sensor>> levelZeroSensors() {
  static const auto LevelZero = Library::open({"libze_loader.so.1", "libze_loader.so"});
  std::vector<std::unique_ptr<Sensor>> sensors;
  if (!LevelZero.has_value()) {
    return sensors;
  }
  // Sysman of its own (zesInit) is there from Level Zero 1.5 on; before, it would have to be
  // enabled before the runtime initializes Level Zero, which SeisSol does not control
  auto* init = LevelZero->function<zes::Init>("zesInit");
  auto* drivers = LevelZero->function<zes::DriverGet>("zesDriverGet");
  auto* devices = LevelZero->function<zes::DeviceGet>("zesDeviceGet");
  auto* pci = LevelZero->function<zes::DevicePciGetProperties>("zesDevicePciGetProperties");
  auto* card = LevelZero->function<zes::DeviceGetCardPowerDomain>("zesDeviceGetCardPowerDomain");
  auto* domains = LevelZero->function<zes::DeviceEnumPowerDomains>("zesDeviceEnumPowerDomains");
  auto* properties = LevelZero->function<zes::PowerGetProperties>("zesPowerGetProperties");
  auto* energy = LevelZero->function<zes::PowerGetEnergyCounter>("zesPowerGetEnergyCounter");
  if (init == nullptr || drivers == nullptr || devices == nullptr || domains == nullptr ||
      properties == nullptr || energy == nullptr || init(0) != Success) {
    return sensors;
  }
  std::uint32_t driverCount = 0;
  if (drivers(&driverCount, nullptr) != Success || driverCount == 0) {
    return sensors;
  }
  std::vector<zes::Handle> driverHandles(driverCount);
  if (drivers(&driverCount, driverHandles.data()) != Success) {
    return sensors;
  }
  for (auto* driver : driverHandles) {
    std::uint32_t deviceCount = 0;
    if (devices(driver, &deviceCount, nullptr) != Success || deviceCount == 0) {
      continue;
    }
    std::vector<zes::Handle> deviceHandles(deviceCount);
    if (devices(driver, &deviceCount, deviceHandles.data()) != Success) {
      continue;
    }
    for (std::uint32_t index = 0; index < deviceCount; ++index) {
      auto* device = deviceHandles[index];
      std::string name = "levelzero/" + std::to_string(index);
      zes::PciProperties pciProperties{};
      pciProperties.stype = ZesStructurePciProperties;
      if (pci != nullptr && pci(device, &pciProperties) == Success) {
        const auto& address = pciProperties.address;
        name = "levelzero/" +
               pciAddress(address.domain, address.bus, address.device, address.function);
      }
      // the whole card if the driver has it; else the package, which holds the tiles
      zes::Handle cardDomain = nullptr;
      if (card != nullptr && card(device, &cardDomain) == Success && cardDomain != nullptr) {
        keep(sensors, std::make_unique<LevelZeroSensor>(name + "/card", cardDomain, energy));
        continue;
      }
      std::uint32_t domainCount = 0;
      if (domains(device, &domainCount, nullptr) != Success || domainCount == 0) {
        continue;
      }
      std::vector<zes::Handle> domainHandles(domainCount);
      if (domains(device, &domainCount, domainHandles.data()) != Success) {
        continue;
      }
      for (auto* domain : domainHandles) {
        zes::PowerProperties power{};
        power.stype = ZesStructurePowerProperties;
        if (properties(domain, &power) == Success && power.onSubdevice == 0) {
          keep(sensors, std::make_unique<LevelZeroSensor>(name + "/package", domain, energy));
        }
      }
    }
  }
  return sensors;
}

} // namespace seissol::monitoring::energy
