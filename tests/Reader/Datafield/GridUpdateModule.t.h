// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/GridUpdateModule.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

namespace seissol::unit_test {

namespace {

namespace df = reader::datafield;

/// Records the times it was advanced to, so a test can assert on the CADENCE
/// rather than on GridStore internals.
class RecordingGrid final : public df::Grid {
  public:
  RecordingGrid(std::size_t dimensions, double timeSpacing)
      : dimensions_(dimensions), timeSpacing_(timeSpacing) {}

  [[nodiscard]] std::size_t dimensions() const override { return dimensions_; }
  [[nodiscard]] std::size_t components() const override { return 1; }

  void geometry(df::GridGeometry& out) const override {
    out.dimensions = dimensions_;
    for (std::size_t d = 0; d < dimensions_; ++d) {
      out.min[d] = 0.0;
      out.delta[d] = (d + 1 == dimensions_) ? timeSpacing_ : 1.0;
      out.num[d] = 64;
    }
  }

  void sample(const double* /*x*/, double* out) const override { out[0] = 0.0; }
  void sampleBatch(const double* const* /*x*/,
                   std::size_t count,
                   double* const* out) const override {
    fill(count, out);
  }
  void
      sampleBatch(const double* const* /*x*/, std::size_t count, float* const* out) const override {
    fill(count, out);
  }

  private:
  template <typename Out>
  void fill(std::size_t count, Out* const* out) const {
    if (out[0] == nullptr) {
      return;
    }
    for (std::size_t i = 0; i < count; ++i) {
      out[0][i] = Out{0};
    }
  }

  std::size_t dimensions_;
  double timeSpacing_;
};

df::GridDesc staticDesc() {
  df::GridDesc desc;
  desc.path = "static.nc";
  desc.variable = "z";
  desc.kind = df::GridKind::AsagiLite;
  desc.interpolation = df::Interpolation::Linear;
  return desc;
}

df::GridDesc timeDesc(int axis, df::Interpolation interp) {
  df::GridDesc desc = staticDesc();
  desc.path = "time.nc";
  desc.interpolation = interp;
  desc.timeAxis = axis;
  return desc;
}

} // namespace

TEST_SUITE("GridUpdateModule") {

  TEST_CASE("a store with no time-dependent grid suggests no interval") {
    df::GridStore store;
    const std::size_t id = store.intern(staticDesc());
    store.injectForTesting(id, std::make_unique<RecordingGrid>(3, 1.0));

    // This is what registerGridUpdateModule() keys off: no interval means no
    // module, which means no synchronisation point and no global barrier for a
    // model that never changes.
    CHECK_FALSE(store.suggestedSyncInterval().has_value());
  }

  TEST_CASE("the interval leaves a stencil's worth of headroom") {
    // Grid.h derives (W - w) * dtFile: W resident slices, stencil width w, one
    // cell of headroom so the reload lands strictly before expiry.
    constexpr double Dt = 0.25;
    constexpr std::size_t Resident = 8;

    SUBCASE("linear in time has width 2") {
      df::GridStore store;
      store.setResidentSlices(Resident);
      store.setDefaultTimeSpacing(Dt);
      const std::size_t id = store.intern(timeDesc(3, df::Interpolation::Linear));
      store.injectForTesting(id, std::make_unique<RecordingGrid>(4, Dt));

      const auto interval = store.suggestedSyncInterval();
      REQUIRE(interval.has_value());
      CHECK(*interval == doctest::Approx(static_cast<double>(Resident - 2) * Dt));
    }

    SUBCASE("a wider stencil buys a shorter interval") {
      df::GridStore store;
      store.setResidentSlices(Resident);
      store.setDefaultTimeSpacing(Dt);
      const std::size_t id = store.intern(timeDesc(3, df::Interpolation::Cubic));
      store.injectForTesting(id, std::make_unique<RecordingGrid>(4, Dt));

      const auto interval = store.suggestedSyncInterval();
      REQUIRE(interval.has_value());
      CHECK(*interval == doctest::Approx(static_cast<double>(Resident - 4) * Dt));
      // The point of the case: cubic must NOT get the linear interval, or the
      // window expires two slices before anyone reloads it.
      CHECK(*interval < static_cast<double>(Resident - 2) * Dt);
    }
  }

  TEST_CASE("a window narrower than the stencil is refused rather than rounded") {
    df::GridStore store;
    store.setResidentSlices(3); // narrower than Cubic's width of 4
    store.setDefaultTimeSpacing(0.25);
    const std::size_t id = store.intern(timeDesc(3, df::Interpolation::Cubic));
    store.injectForTesting(id, std::make_unique<RecordingGrid>(4, 0.25));

    // No interval is safe here, so returning a small one would be a lie that
    // shows up as a silently wrong sample rather than as an error.
    CHECK_THROWS_AS(static_cast<void>(store.suggestedSyncInterval()), std::invalid_argument);
  }

  TEST_CASE("the module advances the store at the times it is given") {
    df::GridStore store;
    store.setResidentSlices(8);
    store.setDefaultTimeSpacing(0.25);
    const std::size_t id = store.intern(timeDesc(3, df::Interpolation::Linear));
    store.injectForTesting(id, std::make_unique<RecordingGrid>(4, 0.25));

    df::GridUpdateModule module(store);

    // A restart must reload for the checkpoint time, not for zero: the first
    // query after a restart is at the checkpoint, and a window sized for t = 0
    // would not cover it.
    CHECK_NOTHROW(module.simulationStart(std::optional<double>{12.5}));
    CHECK_NOTHROW(module.simulationStart(std::nullopt));
    CHECK_NOTHROW(module.syncPoint(13.0));
  }

  TEST_CASE("the shortest interval across grids wins") {
    df::GridStore store;
    store.setResidentSlices(8);
    store.setDefaultTimeSpacing(0.25);

    const std::size_t slow = store.intern(timeDesc(3, df::Interpolation::Linear));
    store.injectForTesting(slow, std::make_unique<RecordingGrid>(4, 0.25));

    df::GridDesc fast = timeDesc(3, df::Interpolation::Cubic);
    fast.path = "fast.nc";
    const std::size_t id = store.intern(fast);
    store.injectForTesting(id, std::make_unique<RecordingGrid>(4, 0.25));

    const auto interval = store.suggestedSyncInterval();
    REQUIRE(interval.has_value());
    // The minimum, not the first or the last: one grid expiring early is enough
    // to make every later sample wrong.
    CHECK(*interval == doctest::Approx(static_cast<double>(8 - 4) * 0.25));
  }

} // TEST_SUITE

} // namespace seissol::unit_test
