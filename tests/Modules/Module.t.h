// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "TestHelper.h"

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

class SyncOnlyModule : public seissol::Module {
  public:
  ModuleHook state{ModuleHook::NullHook};
  double time{0};
  explicit SyncOnlyModule(double interval) {
    setSyncInterval(interval);
    seissol::Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  }
  void syncPoint(double time) override {
    state = ModuleHook::SynchronizationPoint;
    this->time = time;
  }
};

TEST_CASE("Module runthrough" * doctest::test_suite("modules")) {
  TestModule mod1(true, 2);
  TestModule mod2(false, 3);
  SyncOnlyModule syncOnly(1000);
  SyncOnlyModule syncOnlyShort(2);

  Modules::callHook<ModuleHook::PreMPI>();
  CHECK(mod1.state == ModuleHook::PreMPI);
  CHECK(mod2.state == ModuleHook::NullHook);

  Modules::callHook<ModuleHook::PostMPIInit>();
  CHECK(mod1.state == ModuleHook::PostMPIInit);
  CHECK(mod2.state == ModuleHook::NullHook);

  Modules::callHook<ModuleHook::PreMesh>();
  CHECK(mod1.state == ModuleHook::PreMesh);
  CHECK(mod2.state == ModuleHook::NullHook);

  Modules::callHook<ModuleHook::PostMesh>();
  CHECK(mod1.state == ModuleHook::PostMesh);
  CHECK(mod2.state == ModuleHook::PostMesh);

  Modules::callHook<ModuleHook::PreLtsInit>();
  CHECK(mod1.state == ModuleHook::PreLtsInit);
  CHECK(mod2.state == ModuleHook::PostMesh);

  Modules::callHook<ModuleHook::PostLtsInit>();
  CHECK(mod1.state == ModuleHook::PostLtsInit);
  CHECK(mod2.state == ModuleHook::PostMesh);

  Modules::callHook<ModuleHook::PreModel>();
  CHECK(mod1.state == ModuleHook::PreModel);
  CHECK(mod2.state == ModuleHook::PostMesh);

  Modules::callHook<ModuleHook::PostModel>();
  CHECK(mod1.state == ModuleHook::PostModel);
  CHECK(mod2.state == ModuleHook::PostMesh);

  const auto startTime = 1.4;

  Modules::callSimulationStartHook(startTime);
  CHECK(mod1.state == ModuleHook::SimulationStart);
  CHECK(mod2.state == ModuleHook::SimulationStart);

  // exact comparison is ok here
  CHECK(mod1.time == startTime);
  CHECK(mod2.time == startTime);
  CHECK(mod1.state == ModuleHook::SimulationStart);
  CHECK(mod2.state == ModuleHook::SimulationStart);
  CHECK(syncOnly.state == ModuleHook::NullHook);
  CHECK(syncOnlyShort.state == ModuleHook::NullHook);

  // exact comparison is ok here
  CHECK(mod1.time == startTime);
  CHECK(mod2.time == startTime);
  CHECK(syncOnly.time == 0);
  CHECK(syncOnlyShort.time == 0);

  // Regression check: modules that register only for synchronization points
  // must still receive their initial schedule via setSimulationStartTime().
  // Without that, nextSyncPoint_ would remain at 0 and could cause a 0s->0s step.
  CHECK(syncOnly.potentialSyncPoint(startTime, 1e-14, false) == AbsApprox(1000));
  CHECK(syncOnlyShort.potentialSyncPoint(startTime, 1e-14, false) == AbsApprox(2.0));
  CHECK(syncOnly.state == ModuleHook::NullHook);
  CHECK(syncOnlyShort.state == ModuleHook::NullHook);
  CHECK(syncOnly.time == 0);
  CHECK(syncOnlyShort.time == 0);

  std::vector<int32_t> expectedTimes{2, 3, 4, 6, 8, 9, 10, 12, 14, 15};

  double time = startTime;
  double mod1Time = startTime;
  double mod2Time = startTime;
  double syncOnlyShortTime = 0;

  // first check; noone will be called
  time = Modules::callSyncHook(time, 1e-14, false);
  REQUIRE(syncOnly.state == ModuleHook::NullHook);
  REQUIRE(syncOnlyShort.state == ModuleHook::NullHook);
  REQUIRE(syncOnly.time == 0);
  REQUIRE(syncOnlyShort.time == 0);

  for (std::size_t i = 0; i < expectedTimes.size(); ++i) {
    CHECK(time == AbsApprox(expectedTimes[i]));

    if (expectedTimes[i] % 2 == 0) {
      mod1Time = expectedTimes[i];
    }
    if (expectedTimes[i] % 3 == 0) {
      mod2Time = expectedTimes[i];
    }
    if (expectedTimes[i] % 2 == 0) {
      syncOnlyShortTime = expectedTimes[i];
    }

    time = Modules::callSyncHook(time, 1e-14, false);

    // exact comparison is ok here (or rather, _should_ be)

    CHECK(mod1.time == mod1Time);
    CHECK(mod2.time == mod2Time);

    if (mod1Time > startTime) {
      CHECK(mod1.state == ModuleHook::SynchronizationPoint);
    }
    if (mod2Time > startTime) {
      CHECK(mod2.state == ModuleHook::SynchronizationPoint);
    }

    REQUIRE(syncOnly.state == ModuleHook::NullHook);
    REQUIRE(syncOnly.time == 0);

    REQUIRE(syncOnlyShort.time == AbsApprox(syncOnlyShortTime));
    if (syncOnlyShortTime > 0) {
      REQUIRE(syncOnlyShort.state == ModuleHook::SynchronizationPoint);
    } else {
      REQUIRE(syncOnlyShort.state == ModuleHook::NullHook);
    }
  }

  Modules::callHook<ModuleHook::SimulationEnd>();
  CHECK(mod1.state == ModuleHook::SimulationEnd);
  CHECK(mod2.state == ModuleHook::SimulationEnd);

  Modules::callHook<ModuleHook::Shutdown>();
  CHECK(mod1.state == ModuleHook::Shutdown);
  CHECK(mod2.state == ModuleHook::Shutdown);
}

} // namespace seissol::unit_test
