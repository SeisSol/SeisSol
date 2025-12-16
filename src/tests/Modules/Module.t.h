// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "tests/TestHelper.h"

#include <cstdint>
namespace seissol::unit_test {

class TestModule : public seissol::Module {
  public:
  ModuleHook state{ModuleHook::NullHook};
  double time{0};

  void preMPI() override { state = ModuleHook::PreMPI; }

  void postMPIInit() override { state = ModuleHook::PostMPIInit; }

  void preMesh() override { state = ModuleHook::PreMesh; }

  void postMesh() override { state = ModuleHook::PostMesh; }

  void preLtsInit() override { state = ModuleHook::PreLtsInit; }

  void postLtsInit() override { state = ModuleHook::PostLtsInit; }

  void preModel() override { state = ModuleHook::PreModel; }

  void postModel() override { state = ModuleHook::PostModel; }

  void simulationStart(std::optional<double> time) override {
    state = ModuleHook::SimulationStart;
    this->time = time.value_or(0);
  }

  void syncPoint(double time) override {
    state = ModuleHook::SynchronizationPoint;
    this->time = time;
  }

  void simulationEnd() override { state = ModuleHook::SimulationEnd; }

  void shutdown() override { state = ModuleHook::Shutdown; }

  TestModule(bool regAll, double interval) {
    setSyncInterval(interval);

    if (regAll) {
      seissol::Modules::registerHook(*this, ModuleHook::PreMPI);
      seissol::Modules::registerHook(*this, ModuleHook::PostMPIInit);
      seissol::Modules::registerHook(*this, ModuleHook::PreMesh);
    }

    seissol::Modules::registerHook(*this, ModuleHook::PostMesh);

    if (regAll) {
      seissol::Modules::registerHook(*this, ModuleHook::PreLtsInit);
      seissol::Modules::registerHook(*this, ModuleHook::PostLtsInit);
      seissol::Modules::registerHook(*this, ModuleHook::PreModel);
      seissol::Modules::registerHook(*this, ModuleHook::PostModel);
    }

    seissol::Modules::registerHook(*this, ModuleHook::SimulationStart);
    seissol::Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
    seissol::Modules::registerHook(*this, ModuleHook::SimulationEnd);
    seissol::Modules::registerHook(*this, ModuleHook::Shutdown);
  }
};

TEST_CASE("Module runthrough") {
  TestModule mod1(true, 2);
  TestModule mod2(false, 3);

  Modules::callHook<ModuleHook::PreMPI>();
  REQUIRE(mod1.state == ModuleHook::PreMPI);
  REQUIRE(mod2.state == ModuleHook::NullHook);

  Modules::callHook<ModuleHook::PostMPIInit>();
  REQUIRE(mod1.state == ModuleHook::PostMPIInit);
  REQUIRE(mod2.state == ModuleHook::NullHook);

  Modules::callHook<ModuleHook::PreMesh>();
  REQUIRE(mod1.state == ModuleHook::PreMesh);
  REQUIRE(mod2.state == ModuleHook::NullHook);

  Modules::callHook<ModuleHook::PostMesh>();
  REQUIRE(mod1.state == ModuleHook::PostMesh);
  REQUIRE(mod2.state == ModuleHook::PostMesh);

  Modules::callHook<ModuleHook::PreLtsInit>();
  REQUIRE(mod1.state == ModuleHook::PreLtsInit);
  REQUIRE(mod2.state == ModuleHook::PostMesh);

  Modules::callHook<ModuleHook::PostLtsInit>();
  REQUIRE(mod1.state == ModuleHook::PostLtsInit);
  REQUIRE(mod2.state == ModuleHook::PostMesh);

  Modules::callHook<ModuleHook::PreModel>();
  REQUIRE(mod1.state == ModuleHook::PreModel);
  REQUIRE(mod2.state == ModuleHook::PostMesh);

  Modules::callHook<ModuleHook::PostModel>();
  REQUIRE(mod1.state == ModuleHook::PostModel);
  REQUIRE(mod2.state == ModuleHook::PostMesh);

  const auto startTime = 1.4;

  Modules::callSimulationStartHook(startTime);
  REQUIRE(mod1.state == ModuleHook::SimulationStart);
  REQUIRE(mod2.state == ModuleHook::SimulationStart);

  // exact comparison is ok here
  REQUIRE(mod1.time == startTime);
  REQUIRE(mod2.time == startTime);

  std::vector<int32_t> expectedTimes{2, 3, 4, 6, 8, 9, 10, 12, 14, 15};

  double time = startTime;
  double mod1Time = startTime;
  double mod2Time = startTime;

  // first check; noone will be called
  time = Modules::callSyncHook(time, 1e-14, false);

  for (std::size_t i = 0; i < expectedTimes.size(); ++i) {
    REQUIRE(time == AbsApprox(expectedTimes[i]));

    if (expectedTimes[i] % 2 == 0) {
      mod1Time = expectedTimes[i];
    }
    if (expectedTimes[i] % 3 == 0) {
      mod2Time = expectedTimes[i];
    }

    time = Modules::callSyncHook(time, 1e-14, false);

    // exact comparison is ok here (or rather, _should_ be)

    REQUIRE(mod1.time == mod1Time);
    REQUIRE(mod2.time == mod2Time);

    if (mod1Time > startTime) {
      REQUIRE(mod1.state == ModuleHook::SynchronizationPoint);
    }
    if (mod2Time > startTime) {
      REQUIRE(mod2.state == ModuleHook::SynchronizationPoint);
    }
  }

  Modules::callHook<ModuleHook::SimulationEnd>();
  REQUIRE(mod1.state == ModuleHook::SimulationEnd);
  REQUIRE(mod2.state == ModuleHook::SimulationEnd);

  Modules::callHook<ModuleHook::Shutdown>();
  REQUIRE(mod1.state == ModuleHook::Shutdown);
  REQUIRE(mod2.state == ModuleHook::Shutdown);
}

} // namespace seissol::unit_test
