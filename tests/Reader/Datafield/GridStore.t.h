// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <memory>

namespace seissol::unit_test {

using namespace seissol::reader::datafield;

namespace {

/// Stands in for a backend. Nothing the store computes depends on where the
/// bytes come from, so none of this needs a file.
class StubGrid : public Grid {
  public:
  StubGrid(double dt, unsigned slices, std::size_t sliceBytes)
      : dt_(dt), slices_(slices), sliceBytes_(sliceBytes) {}

  [[nodiscard]] std::size_t dimensions() const override { return 3; }
  [[nodiscard]] std::size_t components() const override { return 1; }
  [[nodiscard]] std::size_t bytesPerTimeSlice() const override { return sliceBytes_; }
  void geometry(GridGeometry& out) const override {
    out = GridGeometry{};
    out.dimensions = 3;
    out.min[0] = 0.0;
    out.delta[0] = 100.0;
    out.num[0] = 11;
    out.min[1] = 0.0;
    out.delta[1] = 100.0;
    out.num[1] = 11;
    out.min[2] = 5.0;
    out.delta[2] = dt_;
    out.num[2] = slices_; // time
  }
  void sample(const double* /*x*/, double* out) const override { out[0] = 0.0; }
  void sampleBatch(const double* /*x*/, std::size_t count, double* out) const override {
    for (std::size_t i = 0; i < count; ++i) {
      out[i] = 0.0;
    }
  }

  private:
  double dt_;
  unsigned slices_;
  std::size_t sliceBytes_;
};

GridDesc makeDesc(Interpolation interp = Interpolation::Cubic) {
  GridDesc desc;
  desc.path = "model.nc";
  desc.variable = "rho";
  desc.interpolation = interp;
  desc.timeAxis = 2;
  return desc;
}

} // namespace

TEST_CASE("Datafield GridStore") {
  SUBCASE("descriptors deduplicate by value") {
    GridStore store;
    const auto a = store.intern(makeDesc());
    const auto b = store.intern(makeDesc());
    auto other = makeDesc();
    other.variable = "mu";
    const auto c = store.intern(other);
    REQUIRE(a == b);
    REQUIRE(a != c);
    REQUIRE(store.size() == 2);
  }

  SUBCASE("time spacing is the delta of the time axis") {
    GridStore store;
    const auto id = store.intern(makeDesc());
    store.injectForTesting(id, std::make_unique<StubGrid>(0.25, 400, 1024));
    const auto dt = store.timeSpacing(id);
    REQUIRE(dt.has_value());
    REQUIRE(*dt == doctest::Approx(0.25));
  }

  SUBCASE("bounds follow from the geometry") {
    GridStore store;
    const auto id = store.intern(makeDesc());
    store.injectForTesting(id, std::make_unique<StubGrid>(0.25, 40, 1024));
    double min[MaxGridDimensions];
    double max[MaxGridDimensions];
    store.get(id).bounds(min, max);
    REQUIRE(min[2] == doctest::Approx(5.0));
    REQUIRE(max[2] == doctest::Approx(5.0 + 39 * 0.25));
  }

  SUBCASE("window size follows the memory budget") {
    GridStore store;
    const auto id = store.intern(makeDesc());
    // 1 MiB slices, 64 MiB budget, one grid -> 64 slices.
    store.injectForTesting(id, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    store.setWindowMemoryBudget(64ULL << 20);
    REQUIRE(store.residentSlicesFor(id) == 64);

    // Halving the budget halves the window, and the interval with it.
    store.setWindowMemoryBudget(32ULL << 20);
    REQUIRE(store.residentSlicesFor(id) == 32);
  }

  SUBCASE("the budget is shared between grids") {
    GridStore store;
    const auto a = store.intern(makeDesc());
    auto second = makeDesc();
    second.variable = "mu";
    const auto b = store.intern(second);
    store.injectForTesting(a, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    store.injectForTesting(b, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    store.setWindowMemoryBudget(64ULL << 20);
    REQUIRE(store.residentSlicesFor(a) == 32);
    REQUIRE(store.residentSlicesFor(b) == 32);
  }

  SUBCASE("an expensive grid gets fewer slices than a cheap one") {
    GridStore store;
    const auto cheap = store.intern(makeDesc());
    auto expensiveDesc = makeDesc();
    expensiveDesc.variable = "volume";
    const auto expensive = store.intern(expensiveDesc);
    store.injectForTesting(cheap, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    store.injectForTesting(expensive, std::make_unique<StubGrid>(0.25, 400, 4ULL << 20));
    store.setWindowMemoryBudget(64ULL << 20);
    // Per-grid sizing is the point: a shared global count would force both to
    // the smaller of the two.
    REQUIRE(store.residentSlicesFor(cheap) == 32);
    REQUIRE(store.residentSlicesFor(expensive) == 8);
  }

  SUBCASE("the window never exceeds the file, nor falls below the scheme") {
    GridStore store;
    const auto id = store.intern(makeDesc());
    store.injectForTesting(id, std::make_unique<StubGrid>(0.25, 6, 1ULL << 20));
    store.setWindowMemoryBudget(1024ULL << 20);
    REQUIRE(store.residentSlicesFor(id) == 6); // capped by the file

    GridStore tiny;
    const auto tid = tiny.intern(makeDesc());
    tiny.injectForTesting(tid, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    tiny.setWindowMemoryBudget(1ULL << 10); // absurdly small
    // The floor wins over the budget: a window narrower than the stencil could
    // not be sampled at all, so exceeding the budget is the lesser evil — and
    // it is reported rather than done quietly.
    REQUIRE(tiny.residentSlicesFor(tid) == kernelOf(Interpolation::Cubic).width + 1);
  }

  SUBCASE("an explicit slice count overrides the budget") {
    GridStore store;
    const auto id = store.intern(makeDesc());
    store.injectForTesting(id, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    store.setWindowMemoryBudget(64ULL << 20);
    store.setResidentSlices(9);
    REQUIRE(store.residentSlicesFor(id) == 9);
    store.setResidentSlices(std::nullopt);
    REQUIRE(store.residentSlicesFor(id) == 64);
  }

  SUBCASE("sync interval is (W - w) * dt and scales with the budget") {
    GridStore store;
    const auto id = store.intern(makeDesc());
    store.injectForTesting(id, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    store.setWindowMemoryBudget(64ULL << 20); // -> 64 slices, cubic w = 4
    const auto interval = store.suggestedSyncInterval();
    REQUIRE(interval.has_value());
    REQUIRE(*interval == doctest::Approx((64 - 4) * 0.25));

    store.setWindowMemoryBudget(32ULL << 20);
    REQUIRE(*store.suggestedSyncInterval() == doctest::Approx((32 - 4) * 0.25));
  }

  SUBCASE("a cheaper scheme buys a longer interval at the same memory") {
    GridStore cubic;
    const auto cid = cubic.intern(makeDesc(Interpolation::Cubic));
    cubic.injectForTesting(cid, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    cubic.setWindowMemoryBudget(16ULL << 20);

    GridStore linear;
    const auto lid = linear.intern(makeDesc(Interpolation::Linear));
    linear.injectForTesting(lid, std::make_unique<StubGrid>(0.25, 400, 1ULL << 20));
    linear.setWindowMemoryBudget(16ULL << 20);

    REQUIRE(*linear.suggestedSyncInterval() > *cubic.suggestedSyncInterval());
  }

  SUBCASE("the tightest grid sets the interval") {
    GridStore store;
    const auto slow = store.intern(makeDesc());
    auto fastDesc = makeDesc();
    fastDesc.variable = "fast";
    const auto fast = store.intern(fastDesc);
    store.injectForTesting(slow, std::make_unique<StubGrid>(10.0, 400, 1ULL << 20));
    store.injectForTesting(fast, std::make_unique<StubGrid>(0.01, 400, 1ULL << 20));
    store.setWindowMemoryBudget(64ULL << 20); // 32 slices each
    REQUIRE(*store.suggestedSyncInterval() == doctest::Approx((32 - 4) * 0.01));
  }

  SUBCASE("a static grid has neither spacing nor interval") {
    GridStore store;
    auto desc = makeDesc();
    desc.timeAxis.reset();
    const auto id = store.intern(desc);
    store.injectForTesting(id, std::make_unique<StubGrid>(0.25, 40, 1024));
    REQUIRE_FALSE(store.timeSpacing(id).has_value());
    REQUIRE_FALSE(store.suggestedSyncInterval().has_value());
  }
}

} // namespace seissol::unit_test
