// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/Counters.h"
#include "Monitoring/Energy/Discovery.h"
#include "Monitoring/Energy/Library.h"
#include "Monitoring/Energy/Sensor.h"
#include "Monitoring/Energy/Sysfs.h"
#include "Monitoring/Energy/Vendors.h"
#include "Monitoring/EnergyMeter.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

namespace seissol::unit_test {

namespace energytest {
using monitoring::energy::Part;

//! @brief A directory that stands in for the root of a machine, with a /sys of its own: one per
//! process, as the test cases may run at the same time (ctest -j)
class FakeRoot {
  public:
  FakeRoot()
      : root_(std::filesystem::temp_directory_path() /
              ("seissol-energy-test-" + std::to_string(getpid()) + "-" +
               std::to_string(seissol::Mpi::mpi.rank()))) {
    std::filesystem::remove_all(root_);
  }
  ~FakeRoot() { std::filesystem::remove_all(root_); }

  FakeRoot(const FakeRoot&) = delete;
  FakeRoot(FakeRoot&&) = delete;
  auto operator=(const FakeRoot&) = delete;
  auto operator=(FakeRoot&&) = delete;

  void write(const std::string& path, const std::string& content) const {
    const auto file = root_ / path;
    std::filesystem::create_directories(file.parent_path());
    std::ofstream(file) << content << "\n";
  }

  [[nodiscard]] std::string str() const { return root_.string(); }

  private:
  std::filesystem::path root_;
};

inline const monitoring::energy::Sensor*
    find(const std::vector<std::unique_ptr<monitoring::energy::Sensor>>& sensors,
         const std::string& name) {
  const auto it =
      std::ranges::find_if(sensors, [&](const auto& sensor) { return sensor->name() == name; });
  return it == sensors.end() ? nullptr : it->get();
}
} // namespace energytest

using namespace energytest;

TEST_CASE("Energy accumulator" * doctest::test_suite("monitoring")) {
  SUBCASE("Counts from the first value") {
    monitoring::energy::Accumulator accumulator;
    CHECK(accumulator.update(100.0) == 0.0);
    CHECK(accumulator.update(130.0) == 30.0);
    CHECK(accumulator.update(130.0) == 30.0);
  }
  SUBCASE("A counter with a range runs over") {
    monitoring::energy::Accumulator accumulator(1000.0);
    accumulator.update(900.0);
    // 100 to the end of the range, 50 after it
    CHECK(accumulator.update(50.0) == doctest::Approx(150.0));
    CHECK(accumulator.update(80.0) == doctest::Approx(180.0));
  }
  SUBCASE("A counter without a range is reset") {
    monitoring::energy::Accumulator accumulator;
    accumulator.update(900.0);
    accumulator.update(950.0);
    // what it counted after the reset adds to what it counted before
    CHECK(accumulator.update(20.0) == doctest::Approx(70.0));
  }
}

TEST_CASE("The parts the work is put against" * doctest::test_suite("monitoring")) {
  // the names of the parts of the basis, for the parts measured
  const auto basis = [](const std::vector<Part>& measured) {
    const auto parts = monitoring::energy::basisParts(
        [&](Part part) { return std::ranges::find(measured, part) != measured.end(); });
    std::string names;
    for (const auto part : parts) {
      names += (names.empty() ? "" : "+") + std::string(monitoring::energy::partName(part));
    }
    return names;
  };
  SUBCASE("The node, where it is measured") {
    CHECK(basis({Part::Node, Part::Cpu, Part::Memory, Part::Accelerators, Part::Devices}) ==
          "node");
  }
  SUBCASE("Else the sum of the CPUs, the memory and the GPUs") {
    // NVIDIA Grace: its sockets and the GPUs
    CHECK(basis({Part::Cpu, Part::Devices}) == "cpu+devices");
    // RAPL
    CHECK(basis({Part::Memory, Part::Cpu}) == "cpu+memory");
    CHECK(basis({Part::Devices}) == "devices");
  }
  SUBCASE("The GPUs once: as their library counts them, else as the node does") {
    CHECK(basis({Part::Cpu, Part::Accelerators, Part::Devices}) == "cpu+devices");
    CHECK(basis({Part::Cpu, Part::Accelerators}) == "cpu+accelerators");
  }
  SUBCASE("None, where nothing is measured") { CHECK(basis({}).empty()); }
  SUBCASE("The names of the counters of the parts") {
    CHECK(std::string(monitoring::EnergyMeter::counterName(Part::Node)) == "nodeEnergy");
    CHECK(std::string(monitoring::EnergyMeter::counterName(Part::Cpu)) == "cpuEnergy");
    CHECK(std::string(monitoring::EnergyMeter::counterName(Part::Devices)) == "deviceEnergy");
  }
}

TEST_CASE("Energy from the files of HPE Cray power management" *
          doctest::test_suite("monitoring")) {
  REQUIRE(monitoring::energy::parseCrayPm("2416218307 J 1675950829406375 us").has_value());
  const auto value = monitoring::energy::parseCrayPm("2416218307 J 1675950829406375 us").value();
  CHECK(value.value == 2416218307.0);
  CHECK(value.unit == "J");
  CHECK(value.seconds.value() == doctest::Approx(1675950829.406375));
  CHECK_FALSE(monitoring::energy::parseCrayPm("0 W").value().seconds.has_value());
  CHECK_FALSE(monitoring::energy::parseCrayPm("unknown").has_value());

  const FakeRoot root;
  const std::string pm = "sys/cray/pm_counters/";
  root.write(pm + "energy", "1000 J 1000000000 us");
  root.write(pm + "cpu_energy", "300 J 1000000000 us");
  root.write(pm + "memory_energy", "100 J 1000000000 us");
  root.write(pm + "accel0_energy", "200 J 1000000000 us");
  root.write(pm + "accel1_energy", "200 J 1000000000 us");
  root.write(pm + "power", "500 W 1000000000 us");
  root.write(pm + "freshness", "12");
  root.write(pm + "power_cap", "0 W");

  auto sensors = monitoring::energy::crayPmSensors(root.str());
  REQUIRE(sensors.size() == 5);
  CHECK(find(sensors, "cray-pm/energy")->part() == Part::Node);
  CHECK(find(sensors, "cray-pm/cpu_energy")->part() == Part::Cpu);
  CHECK(find(sensors, "cray-pm/memory_energy")->part() == Part::Memory);
  CHECK(find(sensors, "cray-pm/accel1_energy")->part() == Part::Accelerators);
  CHECK_FALSE(monitoring::energy::crayPmPowerCap(root.str()).has_value());

  root.write(pm + "energy", "1450 J 1002000000 us");
  root.write(pm + "power_cap", "2500 W");
  for (auto& sensor : sensors) {
    if (sensor->name() == "cray-pm/energy") {
      const auto reading = sensor->read();
      REQUIRE(reading.has_value());
      // since the sensor was opened, on the clock of the node
      CHECK(reading.value().joules == 450.0);
      CHECK(reading.value().seconds.value() == doctest::Approx(1002.0));
    }
  }
  CHECK(monitoring::energy::crayPmPowerCap(root.str()).value() == 2500.0);
}

TEST_CASE("Energy from RAPL" * doctest::test_suite("monitoring")) {
  const FakeRoot root;
  const std::string powercap = "sys/class/powercap/";
  root.write(powercap + "intel-rapl:0/name", "package-0");
  root.write(powercap + "intel-rapl:0/energy_uj", "900000000");
  root.write(powercap + "intel-rapl:0/max_energy_range_uj", "1000000000");
  root.write(powercap + "intel-rapl:0:0/name", "dram");
  root.write(powercap + "intel-rapl:0:0/energy_uj", "5000000");
  root.write(powercap + "intel-rapl:0:1/name", "core");
  root.write(powercap + "intel-rapl:0:1/energy_uj", "7000000");
  // the same package once more, through MMIO
  root.write(powercap + "intel-rapl-mmio:0/name", "package-0");
  root.write(powercap + "intel-rapl-mmio:0/energy_uj", "900000000");

  auto sensors = monitoring::energy::powercapSensors(root.str());
  REQUIRE(sensors.size() == 2);
  auto* package = sensors[0]->part() == Part::Cpu ? sensors[0].get() : sensors[1].get();
  CHECK(package->name() == "powercap/intel-rapl:0/package-0");

  // it ran over its range of 1000 J: 100 J to the end of it, 50 J after
  root.write(powercap + "intel-rapl:0/energy_uj", "50000000");
  CHECK(package->read().value().joules == doctest::Approx(150.0));
}

TEST_CASE("The sensors of a node" * doctest::test_suite("monitoring")) {
  const FakeRoot root;
  root.write("sys/class/powercap/intel-rapl:0/name", "package-0");
  root.write("sys/class/powercap/intel-rapl:0/energy_uj", "1000");
  root.write("sys/class/powercap/intel-rapl:1/name", "package-1");
  root.write("sys/class/powercap/intel-rapl:1/energy_uj", "1000");
  root.write("sys/class/powercap/intel-rapl:0:0/name", "dram");
  root.write("sys/class/powercap/intel-rapl:0:0/energy_uj", "1000");
  // an AMD GPU and an ASPEED VGA controller of the board, which is no GPU to measure
  root.write("sys/bus/pci/devices/0000:c1:00.0/class", "0x038000");
  root.write("sys/bus/pci/devices/0000:c1:00.0/vendor", "0x1002");
  root.write("sys/bus/pci/devices/0000:02:00.0/class", "0x030000");
  root.write("sys/bus/pci/devices/0000:02:00.0/vendor", "0x1a03");
  root.write("sys/class/hwmon/hwmon2/name", "xe");
  root.write("sys/class/hwmon/hwmon2/energy1_input", "123456");

  SUBCASE("Without pm_counters, RAPL counts the CPUs and their memory") {
    const auto sensors = monitoring::energy::nodeSensors(root.str(), false);
    CHECK(sensors.size() == 3);
  }

  SUBCASE("With pm_counters, RAPL counts nothing they count") {
    root.write("sys/cray/pm_counters/energy", "1000 J 1000000000 us");
    root.write("sys/cray/pm_counters/cpu_energy", "300 J 1000000000 us");
    const auto sensors = monitoring::energy::nodeSensors(root.str(), false);
    REQUIRE(sensors.size() == 3);
    CHECK(find(sensors, "cray-pm/cpu_energy") != nullptr);
    CHECK(find(sensors, "powercap/intel-rapl:0/package-0") == nullptr);
    CHECK(find(sensors, "powercap/intel-rapl:0:0/dram") != nullptr);
  }

  SUBCASE("The GPUs of the node by their vendor") {
    const auto vendors = monitoring::energy::gpuVendors(root.str());
    REQUIRE(vendors.size() == 1);
    CHECK(*vendors.begin() == monitoring::energy::Vendor::Amd);
    CHECK(monitoring::energy::vendorName(monitoring::energy::Vendor::Amd) == "amd");
  }

  SUBCASE("The energy counters of the drivers of Intel GPUs") {
    const auto sensors = monitoring::energy::hwmonGpuSensors(root.str());
    REQUIRE(sensors.size() == 1);
    CHECK(sensors[0]->part() == Part::Devices);
    CHECK(sensors[0]->source() == "hwmon");
  }
}

TEST_CASE("The energy meter counts the nodes" * doctest::test_suite("monitoring")) {
  const FakeRoot root;
  const std::string pm = "sys/cray/pm_counters/";
  root.write(pm + "energy", "1000 J 1000000000 us");
  root.write(pm + "cpu_energy", "300 J 1000000000 us");
  root.write(pm + "accel0_energy", "200 J 1000000000 us");
  root.write(pm + "accel1_energy", "200 J 1000000000 us");
  root.write(pm + "power_cap", "2500 W");

  SUBCASE("Not enabled, it registers no counter") {
    monitoring::Counters counters;
    monitoring::EnergyMeter meter;
    meter.setup(counters, false, root.str());
    CHECK(counters.size() == 0);
    CHECK_FALSE(meter.describe().find("enabled")->asBool());
  }

  SUBCASE("The parts that the nodes measure") {
    monitoring::Counters counters;
    monitoring::EnergyMeter meter;
    meter.setup(counters, true, root.str());
    // the same on every rank, read on the first rank of every node
    REQUIRE(counters.size() == 3);
    CHECK(counters.info(0).name == "nodeEnergy");
    CHECK(counters.info(1).name == "cpuEnergy");
    CHECK(counters.info(2).name == "acceleratorEnergy");
    CHECK(counters.info(2).rateName == "acceleratorPower");
    CHECK(counters.info(0).perStage);
    CHECK(counters.info(0).scope == monitoring::CounterScope::Node);

    const auto from = counters.snapshot(0.0);
    root.write(pm + "energy", "1400 J 1002000000 us");
    root.write(pm + "accel0_energy", "250 J 1002000000 us");
    root.write(pm + "accel1_energy", "260 J 1002000000 us");
    const auto to = counters.snapshot(1.0);
    const auto totals = counters.reduce(from, to);
    // one node: the ranks share the files
    CHECK(totals.contributors[0] == 1);
    CHECK(totals.amounts[0] == doctest::Approx(400.0));
    // on the clock of the node: two seconds
    CHECK(totals.rates[0] == doctest::Approx(200.0));
    // two sensors updated together: on their clock, two seconds
    CHECK(totals.amounts[2] == doctest::Approx(110.0));
    CHECK(totals.rates[2] == doctest::Approx(55.0));

    // updated apart: on the clock of the rank, one second more
    root.write(pm + "accel0_energy", "300 J 1003000000 us");
    root.write(pm + "accel1_energy", "300 J 1003500000 us");
    const auto later = counters.snapshot(2.0);
    const auto apart = counters.reduce(to, later);
    CHECK(apart.amounts[2] == doctest::Approx(90.0));
    CHECK(apart.rates[2] == doctest::Approx(90.0));

    const auto description = meter.describe();
    const auto* parts = description.find("parts");
    REQUIRE(parts != nullptr);
    REQUIRE(parts->find("accelerators") != nullptr);
    CHECK(parts->find("accelerators")->find("counter")->asString() == "acceleratorEnergy");
    CHECK(parts->find("accelerators")->find("sources")->asArray().front().asString() == "cray-pm");
    CHECK(parts->find("accelerators")->find("nodes")->asDouble() == 1.0);
    CHECK(parts->find("accelerators")->find("sensorsPerNode")->find("max")->asDouble() == 2.0);
    CHECK(parts->find("memory") == nullptr);
    CHECK(description.find("nodePowerCap")->find("max")->asDouble() == 2500.0);
  }
}

TEST_CASE("Energy from the power meters of NVIDIA Grace" * doctest::test_suite("monitoring")) {
  using Clock = std::chrono::steady_clock;
  const auto seconds = [](Clock::duration duration) {
    return std::chrono::duration<double>(duration).count();
  };
  const FakeRoot root;
  const auto meter = [&](int index, const std::string& label, const std::string& microwatts) {
    const std::string hwmon = "sys/class/hwmon/hwmon" + std::to_string(index) + "/";
    root.write(hwmon + "name", "power_meter");
    if (!label.empty()) {
      root.write(hwmon + "power1_oem_info", label);
    }
    root.write(hwmon + "power1_average", microwatts);
    root.write(hwmon + "power1_average_interval", "10");
  };
  meter(1, "Grace Power Socket 0", "100000000");
  meter(2, "CPU Power Socket 0", "80000000");
  meter(3, "SysIO Power Socket 0", "10000000");
  meter(4, "Grace Power Socket 1", "50000000");
  // the meter of a server, which says not what it measures
  meter(5, "", "300000000");

  SUBCASE("The meter of the whole CPU of each socket, integrated") {
    const auto before = Clock::now();
    auto sensors = monitoring::energy::gracePowerSensors(root.str());
    const auto opened = Clock::now();
    REQUIRE(sensors.size() == 2);
    REQUIRE(find(sensors, "hwmon/hwmon1/power_meter/Grace Power Socket 0") != nullptr);
    REQUIRE(find(sensors, "hwmon/hwmon4/power_meter/Grace Power Socket 1") != nullptr);
    auto& socket = sensors[0]->name().ends_with("Socket 0") ? *sensors[0] : *sensors[1];
    CHECK(socket.part() == Part::Cpu);
    CHECK(socket.source() == "hwmon");

    // 100 W since it was opened, up to the read
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    const auto start = Clock::now();
    const auto first = socket.read();
    const auto end = Clock::now();
    REQUIRE(first.has_value());
    CHECK(first.value().joules >= 100.0 * seconds(start - opened) * (1 - 1e-9));
    CHECK(first.value().joules <= 100.0 * seconds(end - before) * (1 + 1e-9));

    // the power read last over the time since: a timeout of 10 s for the thread to read 200 W
    const auto file = std::filesystem::path(root.str()) / "sys/class/hwmon/hwmon1/power1_average";
    std::ofstream(file.string() + ".new") << "200000000\n";
    std::filesystem::rename(file.string() + ".new", file);
    const auto rate = [&] {
      const auto from = Clock::now();
      const auto joulesFrom = socket.read().value().joules;
      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      const auto to = Clock::now();
      return (socket.read().value().joules - joulesFrom) / seconds(to - from);
    };
    const auto deadline = Clock::now() + std::chrono::seconds(10);
    auto watts = rate();
    while (watts < 150.0 && Clock::now() < deadline) {
      watts = rate();
    }
    REQUIRE(watts >= 150.0);

    // a read of the meter that fails keeps the power it read last
    std::filesystem::remove(file);
    const auto fromStart = Clock::now();
    const auto from = socket.read().value().joules;
    const auto fromEnd = Clock::now();
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    const auto toStart = Clock::now();
    const auto to = socket.read().value().joules;
    const auto toEnd = Clock::now();
    CHECK(to - from >= 200.0 * seconds(toStart - fromEnd) * (1 - 1e-9));
    CHECK(to - from <= 200.0 * seconds(toEnd - fromStart) * (1 + 1e-9));
  }

  SUBCASE("For the CPUs, where neither pm_counters nor RAPL count them") {
    CHECK(monitoring::energy::nodeSensors(root.str(), false).size() == 2);
    root.write("sys/class/powercap/intel-rapl:0/name", "package-0");
    root.write("sys/class/powercap/intel-rapl:0/energy_uj", "1000");
    const auto sensors = monitoring::energy::nodeSensors(root.str(), false);
    REQUIRE(sensors.size() == 1);
    CHECK(sensors[0]->source() == "powercap");
  }

  SUBCASE("The energy meter counts them as the energy of the CPUs") {
    monitoring::Counters counters;
    monitoring::EnergyMeter meter;
    meter.setup(counters, true, root.str());
    REQUIRE(counters.size() == 1);
    CHECK(counters.info(0).name == "cpuEnergy");
    const auto description = meter.describe();
    const auto* cpu = description.find("parts")->find("cpu");
    REQUIRE(cpu != nullptr);
    CHECK(cpu->find("sources")->asArray().front().asString() == "hwmon");
    CHECK(cpu->find("sensorsPerNode")->find("max")->asDouble() == 2.0);
  }
}

TEST_CASE("Energy of GPUs without their libraries" * doctest::test_suite("monitoring")) {
  CHECK_FALSE(monitoring::energy::Library::open({"libseissol-not-there.so"}).has_value());
  const auto names = monitoring::energy::rocmLibraryNames("libamd_smi.so");
  CHECK(names.front() == "libamd_smi.so");
  CHECK(names.back() == "/opt/rocm/lib/libamd_smi.so");
  CHECK(monitoring::energy::partName(Part::Accelerators) == "accelerators");
}

} // namespace seissol::unit_test
