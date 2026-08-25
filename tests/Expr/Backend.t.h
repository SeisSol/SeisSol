// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"
#include "Reader/Scripting/DataTable.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace seissol::expr::test {

namespace {

namespace df = reader::datafield;
using reader::scripting::DataTable;
using reader::scripting::Direction;

// A grid backed by the real interpolation kernel rather than a constant. The
// point of these cases is the sampler adapter -- the pointer arrays, the
// null-skip, the f32 coordinate conversion -- and a fake that ignores its
// arguments would pass whatever the adapter did with them.
class LinearRampGrid final : public df::Grid {
  public:
  static constexpr unsigned Num = 8;
  static constexpr unsigned Comps = 2;

  LinearRampGrid() {
    data_.resize(static_cast<std::size_t>(Num) * Num * Num * Comps);
    std::size_t stride = Comps;
    for (unsigned d = 0; d < 3; ++d) {
      view_.min[d] = 0.0;
      view_.deltaInv[d] = 1.0;
      view_.num[d] = Num;
      view_.stride[d] = stride;
      stride *= Num;
    }
    view_.data = data_.data();
    view_.type = df::ElementType::Double;
    view_.dimensions = 3;
    view_.components = Comps;
    view_.componentStride = 1;

    for (unsigned k = 0; k < Num; ++k) {
      for (unsigned j = 0; j < Num; ++j) {
        for (unsigned i = 0; i < Num; ++i) {
          const std::size_t base = ((static_cast<std::size_t>(k) * Num + j) * Num + i) * Comps;
          const double value = value_(i, j, k);
          data_[base + 0] = value;
          data_[base + 1] = -value;
        }
      }
    }
  }

  // Trilinear in the coordinates, so linear interpolation reproduces it
  // exactly and the expected value is a closed form rather than a golden file.
  static double value_(double i, double j, double k) { return i + 10.0 * j + 100.0 * k; }

  [[nodiscard]] std::size_t dimensions() const override { return 3; }
  [[nodiscard]] std::size_t components() const override { return Comps; }

  void geometry(df::GridGeometry& out) const override {
    out.dimensions = 3;
    for (unsigned d = 0; d < 3; ++d) {
      out.min[d] = 0.0;
      out.delta[d] = 1.0;
      out.num[d] = Num;
    }
  }

  void sample(const double* x, double* out) const override {
    df::sample(view_, df::Interpolation::Linear, x, out);
  }

  void sampleBatch(const double* const* x, std::size_t count, double* const* out) const override {
    df::sampleBatch(view_, df::Interpolation::Linear, x, count, out, scratch_);
  }

  void sampleBatch(const double* const* x, std::size_t count, float* const* out) const override {
    df::sampleBatch(view_, df::Interpolation::Linear, x, count, out, scratch_);
  }

  private:
  std::vector<double> data_;
  df::ArrayView view_;
  mutable df::SampleScratch scratch_;
};

} // namespace

TEST_SUITE("ExprBackend") {

  TEST_CASE("a pointwise program runs end to end") {
    const Program program = compileSderiv("x * 2.0 + y", "out");
    constexpr std::size_t NumPoints = 8;

    std::vector<double> x(NumPoints);
    std::vector<double> y(NumPoints);
    std::vector<double> out(NumPoints, -1.0);
    for (std::size_t i = 0; i < NumPoints; ++i) {
      x[i] = static_cast<double>(i);
      y[i] = 100.0 * static_cast<double>(i);
    }

    DataTable table(NumPoints);
    table.bindView<double>("x", Direction::In, x.data());
    table.bindView<double>("y", Direction::In, y.data());
    table.bindView<double>("out", Direction::Out, out.data());

    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    const auto kernel = makeKernel(program, binding, store, {});
    REQUIRE(kernel != nullptr);
    CHECK(kernel->kind() == BackendKind::Interpreter);

    kernel->precompute(table);
    kernel->run(table);
    for (std::size_t i = 0; i < NumPoints; ++i) {
      CHECK(out[i] == doctest::Approx(2.0 * x[i] + y[i]));
    }
  }

  TEST_CASE("a grid lookup reaches the sampler and skips the other components") {
    const Program program =
        compileSderiv("grid m = \"asagi\", \"ramp.nc\", \"linear\", \"a\", \"b\"\n"
                      "m_b(x, y, z)",
                      "out");
    REQUIRE(program.grids().size() == 1);
    constexpr std::size_t NumPoints = 5;

    const std::vector<double> x = {1.5, 2.0, 3.25, 4.0, 5.5};
    const std::vector<double> y = {1.0, 1.0, 2.0, 2.0, 3.0};
    const std::vector<double> z = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> out(NumPoints, 0.0);

    DataTable table(NumPoints);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<double>("y", Direction::In, y.data());
    table.bindViewConst<double>("z", Direction::In, z.data());
    table.bindView<double>("out", Direction::Out, out.data());

    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    const std::size_t id = store.intern(program.grids()[0]);
    store.injectForTesting(id, std::make_unique<LinearRampGrid>());

    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);
    kernel->run(table);

    // Component 1 is the negated ramp; getting component 0's sign would mean
    // the null-skip pointed the destination at the wrong entry.
    for (std::size_t i = 0; i < NumPoints; ++i) {
      CHECK(out[i] == doctest::Approx(-LinearRampGrid::value_(x[i], y[i], z[i])));
    }
  }

  TEST_CASE("a grid lookup works from the f32 compute path") {
    Program program = compileSderiv("grid m = \"asagi\", \"ramp.nc\", \"linear\", \"a\", \"b\"\n"
                                    "m_a(x, y, z)",
                                    "out");
    program.setComputeType(ComputeType::F32);
    constexpr std::size_t NumPoints = 4;

    const std::vector<double> x = {1.5, 2.0, 3.25, 4.0};
    const std::vector<double> y = {1.0, 1.0, 2.0, 2.0};
    const std::vector<double> z = {0.0, 1.0, 2.0, 3.0};
    std::vector<float> out(NumPoints, 0.0F);

    DataTable table(NumPoints);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<double>("y", Direction::In, y.data());
    table.bindViewConst<double>("z", Direction::In, z.data());
    table.bindView<float>("out", Direction::Out, out.data());

    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    const std::size_t id = store.intern(program.grids()[0]);
    store.injectForTesting(id, std::make_unique<LinearRampGrid>());

    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);
    kernel->run(table);

    // Coordinates are converted to f64 for the lookup (Grid.h fixes them there),
    // so the tolerance is the f32 store at the end and nothing else.
    for (std::size_t i = 0; i < NumPoints; ++i) {
      CHECK(out[i] == doctest::Approx(LinearRampGrid::value_(x[i], y[i], z[i])).epsilon(1e-6));
    }
  }

  TEST_CASE("a partitioned point set lands back on the right points") {
    const Program program = compileSderiv("x + group", "out");
    constexpr std::size_t NumPoints = 6;

    const std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const std::vector<std::int32_t> group = {2, 1, 2, 1, 3, 1};
    std::vector<double> out(NumPoints, -1.0);

    DataTable table(NumPoints);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindViewConst<std::int32_t>("group", Direction::In, group.data());
    table.bindView<double>("out", Direction::Out, out.data());

    Binding binding = Binding::bind(program, table);
    REQUIRE(binding.groupRanges().size() == 3);
    REQUIRE_FALSE(binding.permutation().empty());

    df::GridStore store;
    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);
    kernel->run(table);

    for (std::size_t i = 0; i < NumPoints; ++i) {
      CHECK(out[i] == doctest::Approx(x[i] + static_cast<double>(group[i])));
    }
  }

  TEST_CASE("hoisting is applied and survives across run() calls") {
    const Program program = compileSderiv("sqrt(x) * y", "out");
    constexpr std::size_t NumPoints = 4;

    const std::vector<double> x = {4.0, 9.0, 16.0, 25.0};
    std::vector<double> y = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> out(NumPoints, 0.0);

    DataTable table(NumPoints);
    table.bindViewConst<double>("x", Direction::In, x.data());
    table.bindView<double>("y", Direction::In, y.data());
    table.bindView<double>("out", Direction::Out, out.data());

    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    BackendOptions options;
    options.lowering.invariantInputs = {"x"};
    options.lowering.hoistThreshold = 1;

    const auto kernel = makeKernel(program, binding, store, options);
    // Without this the case would pass even if hoisting silently did nothing,
    // which is exactly the failure it is meant to catch.
    REQUIRE(binding.persistentSlotCount() > 0);

    kernel->precompute(table);
    kernel->run(table);
    for (std::size_t i = 0; i < NumPoints; ++i) {
      CHECK(out[i] == doctest::Approx(std::sqrt(x[i]) * y[i]));
    }

    // A second call must reuse the hoisted values rather than recompute or lose
    // them -- the whole point of the split being across the timestep loop.
    for (std::size_t i = 0; i < NumPoints; ++i) {
      y[i] *= 10.0;
      out[i] = 0.0;
    }
    kernel->run(table);
    for (std::size_t i = 0; i < NumPoints; ++i) {
      CHECK(out[i] == doctest::Approx(std::sqrt(x[i]) * y[i]));
    }
  }

  TEST_CASE("an unimplemented backend falls back rather than returning null") {
    const Program program = compileSderiv("x + 1.0", "out");
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> out(2, 0.0);

    DataTable table(2);
    table.bindView<double>("x", Direction::In, x.data());
    table.bindView<double>("out", Direction::Out, out.data());

    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    BackendOptions options;
    // RtcCpu is implemented since Package 5, so the case now uses one that is
    // not. The property under test is the fallback, not which backend is
    // missing this week.
    options.preferred = BackendKind::RtcCuda;
    options.allowFallback = true;

    const auto kernel = makeKernel(program, binding, store, options);
    REQUIRE(kernel != nullptr);
    CHECK(kernel->kind() == BackendKind::Interpreter);
  }

  TEST_CASE("run(KernelArgs) agrees with run(table), at any granularity") {
    // The shape EasiBoundary::query needs: one binding, and only the bases move
    // per face. Bitwise agreement is the bar -- the two paths must not be two
    // numerical paths, or the element-wise call site would silently drift from
    // the dense one.
    const Program program = compileSderivModule(
        "def a = x*x + y*y\nout def u = a + z*group\nout def v = a - z\nout def w = a*z\n");
    constexpr std::size_t Nodes = 9;

    std::vector<double> face(3 * Nodes);
    std::vector<double> out(3 * Nodes, -1.0);
    for (std::size_t i = 0; i < face.size(); ++i) {
      face[i] = 0.1 * static_cast<double>(i + 1);
    }

    DataTable table(Nodes);
    table.bindViewConst<double>("x", Direction::In, face.data(), 3, 0);
    table.bindViewConst<double>("y", Direction::In, face.data(), 3, 1);
    table.bindViewConst<double>("z", Direction::In, face.data(), 3, 2);
    // bindConstant rather than bindComputed: the latter has no address
    // arithmetic behind it and would make the binding unaddressable.
    table.bindConstant<std::int32_t>("group", 1);
    table.bindView<double>("u", Direction::Out, out.data(), 3, 0);
    table.bindView<double>("v", Direction::Out, out.data(), 3, 1);
    table.bindView<double>("w", Direction::Out, out.data(), 3, 2);

    Binding binding = Binding::bind(program, table);
    REQUIRE(binding.addressable());

    df::GridStore store;
    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);
    kernel->run(table);
    const std::vector<double> reference = out;

    KernelArgs args{};
    const void* inputs[4] = {face.data(), face.data(), face.data(), nullptr};
    void* outputs[3] = {out.data(), out.data(), out.data()};
    args.inputs = inputs;
    args.inputCount = 4;
    args.outputs = outputs;
    args.outputCount = 3;

    SUBCASE("the whole range in one call") {
      std::fill(out.begin(), out.end(), -1.0);
      args.first = 0;
      args.count = Nodes;
      kernel->run(args);
      for (std::size_t i = 0; i < out.size(); ++i) {
        CHECK(out[i] == reference[i]);
      }
    }

    SUBCASE("one point per call") {
      std::fill(out.begin(), out.end(), -1.0);
      for (std::size_t n = 0; n < Nodes; ++n) {
        args.first = n;
        args.count = 1;
        kernel->run(args);
      }
      for (std::size_t i = 0; i < out.size(); ++i) {
        CHECK(out[i] == reference[i]);
      }
    }

    SUBCASE("a null base keeps the one the binding was built with") {
      // `group` is passed as null above, so it must still read the bound
      // constant rather than whatever inputs[3] happened to be.
      std::fill(out.begin(), out.end(), -1.0);
      args.first = 0;
      args.count = Nodes;
      kernel->run(args);
      CHECK(out[0] == reference[0]);
    }
  }

  TEST_CASE("a computed column makes a binding unaddressable") {
    const Program program = compileSderivModule("out def u = x + group\n");
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> out(2, 0.0);

    DataTable table(2);
    table.bindViewConst<double>("x", Direction::In, x.data());
    // The one bind form with no address arithmetic behind it.
    table.bindComputed("group", [](std::size_t) -> std::int32_t { return 1; });
    table.bindView<double>("u", Direction::Out, out.data());

    const Binding binding = Binding::bind(program, table);
    // Not an error -- run(table) still works. It is what makeKernel will check
    // before offering a device backend.
    CHECK_FALSE(binding.addressable());
  }

} // TEST_SUITE

} // namespace seissol::expr::test
