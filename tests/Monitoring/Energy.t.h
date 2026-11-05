// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/Energy/Discovery.h"
#include "Monitoring/Energy/Library.h"
#include "Monitoring/Energy/Sensor.h"
#include "Monitoring/Energy/Sysfs.h"
#include "Monitoring/Energy/Vendors.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace {
using monitoring::energy::Part;

//! @brief A directory that stands in for the root of a machine, with a /sys of its own.
class FakeRoot {
  public:
  FakeRoot()
      : root_(std::filesystem::temp_directory_path() /
              ("seissol-energy-test-" + std::to_string(seissol::Mpi::mpi.rank()))) {
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

const monitoring::energy::Sensor*
    find(const std::vector<std::unique_ptr<monitoring::energy::Sensor>>& sensors,
         const std::string& name) {
  const auto it = std::find_if(
      sensors.begin(), sensors.end(), [&](const auto& sensor) { return sensor->name() == name; });
  return it == sensors.end() ? nullptr : it->get();
}
} // namespace

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

TEST_CASE("Energy of GPUs without their libraries" * doctest::test_suite("monitoring")) {
  CHECK_FALSE(monitoring::energy::Library::open({"libseissol-not-there.so"}).has_value());
  const auto names = monitoring::energy::rocmLibraryNames("libamd_smi.so");
  CHECK(names.front() == "libamd_smi.so");
  CHECK(names.back() == "/opt/rocm/lib/libamd_smi.so");
  CHECK(monitoring::energy::partName(Part::Accelerators) == "accelerators");
}

} // namespace seissol::unit_test
