// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// The interpreter is the oracle for every other backend, so it needs an oracle
// of its own. The reference below is deliberately built the other way round from
// TileInterpreter: one point at a time, one value per node, no slot reuse, no
// stages, no tiles. Everything the tiled path can get wrong — a slot recycled
// one instruction too early, a value written in place over an operand still
// needed, a partial tail tile, a hoisted value read from the wrong point offset,
// a state slot updated sequentially instead of in parallel — shows up as a
// mismatch here, and the comparison is on bit patterns rather than a tolerance,
// because the two paths execute the same operations on the same inputs and any
// difference at all is a bug rather than rounding.

#include <doctest.h>

#include "Expr/Core.h"
#include "Expr/Expand.h"
#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"
#include "TestHelper.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace {

using DT = reader::scripting::DataType;

template <typename T>
bool sameBits(T a, T b) {
  return std::memcmp(&a, &b, sizeof(T)) == 0;
}

// Two grids with arrays of their own: one in two dimensions, sampled linearly, one in one
// dimension, sampled with the cubic stencil, two components each, values that no scheme reproduces
// exactly. The reference samples them with datafield::sampleBatch(), the interpreter with the
// instructions of the written-out lookups: the two have to agree bit for bit.
class TestGrids final : public GridArrays {
  public:
  TestGrids() {
    build(0, {9, 7}, -2.0, 0.5);
    build(1, {11}, -2.5, 0.375);
  }

  [[nodiscard]] const reader::datafield::ArrayView& array(GridId grid) const override {
    return views_.at(static_cast<std::size_t>(grid));
  }
  void report(GridId /*grid*/, CheckKind kind) const override {
    FAIL((kind == CheckKind::Window ? "no grid here has a window"
                                    : "no grid here fails beyond its volume"));
  }

  static reader::datafield::GridDesc desc(GridId grid) {
    reader::datafield::GridDesc desc;
    desc.path = grid == 0 ? "a.nc" : "b.nc";
    desc.interpolation = grid == 0 ? reader::datafield::Interpolation::Linear
                                   : reader::datafield::Interpolation::Cubic;
    return desc;
  }

  private:
  void build(std::size_t grid, const std::vector<unsigned>& num, double min, double delta) {
    auto& view = views_[grid];
    auto& data = data_[grid];
    view.type = reader::datafield::ElementType::Double;
    view.dimensions = static_cast<unsigned>(num.size());
    view.components = 2;
    view.componentStride = 1;
    std::size_t stride = 2;
    for (std::size_t d = 0; d < num.size(); ++d) {
      view.min[d] = min + 0.25 * static_cast<double>(d);
      view.deltaInv[d] = 1.0 / (delta + 0.125 * static_cast<double>(d));
      view.num[d] = num[d];
      view.stride[d] = stride;
      stride *= num[d];
    }
    data.resize(stride);
    for (std::size_t i = 0; i < stride; ++i) {
      data[i] = std::cos(0.61 * static_cast<double>(i) + static_cast<double>(grid)) * 3.0 +
                0.07 * static_cast<double>(i % 5);
    }
    view.data = data.data();
  }

  std::array<reader::datafield::ArrayView, 2> views_{};
  std::array<std::vector<double>, 2> data_{};
};

// Plain SoA arrays behind the Binding::gather / Binding::scatter contract, and the states
// slot-major over the points.
template <typename T>
class ArrayTileIo : public TileIo<T> {
  public:
  ArrayTileIo(const std::vector<T>& inputs,
              std::vector<T>& outputs,
              std::size_t numPoints,
              std::size_t numInputs,
              std::size_t numOutputs,
              std::vector<T>* state = nullptr)
      : inputs_(inputs), outputs_(outputs), numPoints_(numPoints), numInputs_(numInputs),
        numOutputs_(numOutputs), state_(state) {}

  void gather(std::size_t first, std::size_t count, T* dst) const override {
    for (std::size_t i = 0; i < numInputs_; ++i) {
      for (std::size_t l = 0; l < count; ++l) {
        dst[i * count + l] = inputs_[i * numPoints_ + first + l];
      }
    }
  }
  void scatter(std::size_t first, std::size_t count, const T* src) override {
    for (std::size_t o = 0; o < numOutputs_; ++o) {
      for (std::size_t l = 0; l < count; ++l) {
        outputs_[o * numPoints_ + first + l] = src[o * count + l];
      }
    }
  }
  void gatherState(std::size_t first, std::size_t count, T* dst) const override {
    for (std::size_t s = 0; s < state_->size() / numPoints_; ++s) {
      for (std::size_t l = 0; l < count; ++l) {
        dst[s * count + l] = (*state_)[s * numPoints_ + first + l];
      }
    }
  }
  void scatterState(std::size_t first, std::size_t count, const T* src) override {
    for (std::size_t s = 0; s < state_->size() / numPoints_; ++s) {
      for (std::size_t l = 0; l < count; ++l) {
        (*state_)[s * numPoints_ + first + l] = src[s * count + l];
      }
    }
  }

  private:
  const std::vector<T>& inputs_;
  std::vector<T>& outputs_;
  std::size_t numPoints_;
  std::size_t numInputs_;
  std::size_t numOutputs_;
  std::vector<T>* state_;
};

// Point-at-a-time evaluation of the DAG. Arena ids are topologically ordered, so
// one ascending sweep per point suffices and no recursion is involved.
template <typename T>
std::vector<T> referenceEvaluate(const Program& program,
                                 const std::vector<T>& inputs,
                                 std::size_t numPoints,
                                 const TestGrids* grids,
                                 std::vector<T>* state) {
  const Arena& arena = program.arena();
  std::vector<T> outputs(program.outputs().size() * numPoints, T(0));
  std::vector<T> nextState(state == nullptr ? 0 : state->size(), T(0));

  std::vector<int> inputOf(arena.size(), -1);
  std::vector<int> stateOf(arena.size(), -1);
  for (NodeId id = 0; id < static_cast<NodeId>(arena.size()); ++id) {
    if (arena[id].kind != Kind::Field) {
      continue;
    }
    const std::string& channel = arena.channelName(arena[id].ch);
    for (std::size_t i = 0; i < program.inputs().size(); ++i) {
      if (program.inputs()[i].name == channel) {
        inputOf[id] = static_cast<int>(i);
      }
    }
    for (std::size_t i = 0; i < program.state().size(); ++i) {
      if (program.state()[i].name == channel) {
        stateOf[id] = static_cast<int>(i);
      }
    }
  }

  std::vector<T> value(arena.size(), T(0));
  std::vector<NodeId> kids;
  for (std::size_t point = 0; point < numPoints; ++point) {
    for (NodeId id = 0; id < static_cast<NodeId>(arena.size()); ++id) {
      const Node& node = arena[id];
      switch (node.kind) {
      case Kind::Const:
        value[id] = static_cast<T>(node.value);
        break;
      case Kind::Field:
        if (stateOf[id] >= 0) {
          value[id] = (*state)[static_cast<std::size_t>(stateOf[id]) * numPoints + point];
        } else {
          value[id] = inputs[static_cast<std::size_t>(inputOf[id]) * numPoints + point];
        }
        break;
      case Kind::PW: {
        std::array<T, 3> args{};
        arena.children(id, kids);
        for (std::size_t k = 0; k < kids.size(); ++k) {
          args[k] = value[kids[k]];
        }
        value[id] = applyPw<T>(node.fn, args.data());
        break;
      }
      case Kind::Lookup: {
        REQUIRE(grids != nullptr);
        arena.children(id, kids);
        // the coordinates in binary64, the sample of the component, as Grid.h has it
        std::array<double, MaxLookupDimension> scalars{};
        std::array<const double*, MaxLookupDimension> pointers{};
        for (std::size_t k = 0; k < kids.size(); ++k) {
          scalars[k] = static_cast<double>(value[kids[k]]);
          pointers[k] = &scalars[k];
        }
        std::array<double, 2> samples{};
        std::array<double*, 2> outs{};
        outs[static_cast<std::size_t>(node.comp)] = &samples[static_cast<std::size_t>(node.comp)];
        reader::datafield::SampleScratch scratch;
        reader::datafield::sampleBatch(
            grids->array(node.grid),
            program.grids()[static_cast<std::size_t>(node.grid)].interpolation,
            pointers.data(),
            1,
            outs.data(),
            scratch);
        value[id] = static_cast<T>(samples[static_cast<std::size_t>(node.comp)]);
        break;
      }
      default:
        break;
      }
    }
    for (std::size_t o = 0; o < program.roots().size(); ++o) {
      outputs[o * numPoints + point] = value[program.roots()[o]];
    }
    for (std::size_t s = 0; s < program.state().size(); ++s) {
      nextState[s * numPoints + point] = value[program.state()[s].root];
    }
  }
  if (state != nullptr) {
    *state = nextState; // parallel assignment: every root saw the old values
  }
  return outputs;
}

struct Mismatch {
  bool any{false};
  std::size_t index{0};
};

template <typename T>
Mismatch firstDifference(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size()) {
    return {true, 0};
  }
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (!sameBits(a[i], b[i])) {
      return {true, i};
    }
  }
  return {};
}

// Builds a random but well-formed pointwise program over x, y, z, t and one grid.
Program randomProgram(std::mt19937& rng, int depth, bool withLookup, int outputCount) {
  static constexpr Fn Pool[] = {Fn::Add,  Fn::Sub, Fn::Mul, Fn::Div,   Fn::Min,   Fn::Max,
                                Fn::Sqrt, Fn::Abs, Fn::Exp, Fn::Log,   Fn::Sin,   Fn::Cos,
                                Fn::Tanh, Fn::Neg, Fn::Rcp, Fn::Sign,  Fn::Floor, Fn::Atan2,
                                Fn::Lt,   Fn::Le,  Fn::Pow, Fn::Select};
  Program program;
  for (const char* n : {"x", "y", "z", "t"}) {
    program.addInput(n, DT::F64);
  }
  Arena& arena = program.arena();
  const GridId plane = withLookup ? program.internGrid(TestGrids::desc(0)) : NoGrid;
  const GridId line = withLookup ? program.internGrid(TestGrids::desc(1)) : NoGrid;

  std::vector<NodeId> leaves{arena.field("x"),
                             arena.field("y"),
                             arena.field("z"),
                             arena.field("t"),
                             arena.konst(0.5),
                             arena.konst(-1.25),
                             arena.konst(2.0)};
  if (withLookup) {
    leaves.push_back(arena.lookup(plane, 0, {arena.field("x"), arena.field("y")}));
    leaves.push_back(arena.lookup(line, 1, {arena.field("z")}));
  }

  const auto build = [&](auto&& self, int d) -> NodeId {
    if (d == 0) {
      return leaves[rng() % leaves.size()];
    }
    const Fn f = Pool[rng() % (sizeof(Pool) / sizeof(Pool[0]))];
    const int n = arity(f);
    if (n == 1) {
      return arena.pw(f, self(self, d - 1));
    }
    if (n == 2) {
      return arena.pw(f, self(self, d - 1), self(self, d - 1));
    }
    return arena.pw(f, self(self, d - 1), self(self, d - 1), self(self, d - 1));
  };

  for (int o = 0; o < outputCount; ++o) {
    program.addOutput("u" + std::to_string(o), DT::F64, build(build, depth));
  }
  return program;
}

template <typename T>
std::vector<T> randomInputs(std::mt19937& rng, std::size_t numInputs, std::size_t numPoints) {
  std::uniform_real_distribution<double> dist(-2.0, 2.0);
  std::vector<T> values(numInputs * numPoints);
  for (T& v : values) {
    v = static_cast<T>(dist(rng));
  }
  return values;
}

} // namespace

TEST_SUITE("Expr::Interp") {

  TEST_CASE("tile width follows the live slot count") {
    // 16 KiB budget, f64: 8 slots -> 256 lanes, 32 slots -> 64 lanes.
    CHECK(chooseTileSize(8, ComputeType::F64, 16384) == 256);
    CHECK(chooseTileSize(32, ComputeType::F64, 16384) == 64);
    CHECK(chooseTileSize(8, ComputeType::F32, 16384) == 512);
    // Clamped at both ends, and always a whole number of vector lanes.
    CHECK(chooseTileSize(100000, ComputeType::F64, 16384) == MinTileSize);
    CHECK(chooseTileSize(1, ComputeType::F64, 1u << 30U) == MaxTileSize);
    CHECK(chooseTileSize(7, ComputeType::F64, 16384) % TileLaneGranularity == 0);
  }

  TEST_CASE("Mod is Lua's floored remainder") {
    // The results of Lua 5.4's `%` on floats, sign of zero included; fmod alone differs on every
    // pair with operands of opposite signs.
    constexpr double Inf = std::numeric_limits<double>::infinity();
    const std::array<std::array<double, 3>, 12> cases{{{-1.0, -3.0, -1.0},
                                                       {-1.0, 3.0, 2.0},
                                                       {1.0, -3.0, -2.0},
                                                       {5.5, -2.0, -0.5},
                                                       {-5.5, 2.0, 0.5},
                                                       {0.0, -2.0, 0.0},
                                                       {-0.0, 2.0, -0.0},
                                                       {3.0, Inf, 3.0},
                                                       {-3.0, Inf, Inf},
                                                       {3.0, -Inf, -Inf},
                                                       {-6.0, 3.0, -0.0},
                                                       {6.0, -3.0, 0.0}}};
    for (const auto& operands : cases) {
      const double a = operands[0];
      const double b = operands[1];
      const std::array<double, 2> args{a, b};
      const auto result = applyPw<double>(Fn::Mod, args.data());
      CAPTURE(a);
      CAPTURE(b);
      CHECK(unit_test::bitwiseEqual(result, operands[2]));
    }
  }

  TEST_CASE("matches the scalar reference bit for bit") {
    // with lookups: their instructions against datafield::sampleBatch()
    TestGrids grids;
    std::mt19937 rng(20260821);

    for (int trial = 0; trial < 40; ++trial) {
      const bool withLookup = (trial % 3) == 0;
      const Program program = randomProgram(
          rng, 3 + static_cast<int>(trial % 3), withLookup, 1 + static_cast<int>(trial % 3));
      validate(program);

      // Point counts that are not multiples of the tile width, so every run ends
      // on a partial tail tile.
      for (const std::size_t numPoints :
           {std::size_t{1}, std::size_t{7}, std::size_t{17}, std::size_t{100}, std::size_t{257}}) {
        const std::vector<double> inputs =
            randomInputs<double>(rng, program.inputs().size(), numPoints);
        const std::vector<double> expected =
            referenceEvaluate<double>(program, inputs, numPoints, &grids, nullptr);

        for (const std::size_t tile :
             {std::size_t{1}, std::size_t{3}, std::size_t{8}, std::size_t{64}}) {
          const LoweredProgram lowered = lower(program);
          std::vector<double> outputs(program.outputs().size() * numPoints, 0.0);
          ArrayTileIo<double> io(
              inputs, outputs, numPoints, program.inputs().size(), program.outputs().size());
          InterpreterOptions options;
          options.tileSize = tile;
          TileInterpreter<double> interp(program, lowered, &grids, options);
          std::vector<double> persistent(
              static_cast<std::size_t>(lowered.persistentSlotCount()) * numPoints, 0.0);
          interp.precompute(io, numPoints, persistent.data());
          interp.run(io, numPoints, persistent.data());

          const Mismatch diff = firstDifference(expected, outputs);
          INFO("trial=" << trial << " numPoints=" << numPoints << " tile=" << tile
                        << " index=" << diff.index);
          CHECK_FALSE(diff.any);
        }
      }
    }
  }

  // Hoisting must be a pure scheduling change. If it is not, the difference
  // shows up here as a bit mismatch rather than as a plausible-looking number.
  TEST_CASE("hoisting does not change results") {
    TestGrids grids;
    std::mt19937 rng(987654321);

    for (int trial = 0; trial < 30; ++trial) {
      const Program program = randomProgram(rng, 4, true, 2);
      validate(program);
      const std::size_t numPoints = 173;
      const std::vector<double> inputs =
          randomInputs<double>(rng, program.inputs().size(), numPoints);
      const std::vector<double> expected =
          referenceEvaluate<double>(program, inputs, numPoints, &grids, nullptr);

      LowerOptions options;
      options.invariantInputs = {"x", "y", "z"};
      options.invariantGrids = {0, 1};
      const LoweredProgram lowered = lower(program, options);
      REQUIRE(lowered.hasPrecompute());

      std::vector<double> outputs(program.outputs().size() * numPoints, 0.0);
      ArrayTileIo<double> io(
          inputs, outputs, numPoints, program.inputs().size(), program.outputs().size());
      InterpreterOptions interpOptions;
      interpOptions.tileSize = 24;
      TileInterpreter<double> interp(program, lowered, &grids, interpOptions);
      std::vector<double> persistent(
          static_cast<std::size_t>(lowered.persistentSlotCount()) * numPoints, 0.0);
      interp.precompute(io, numPoints, persistent.data());

      // Several calls: nothing may drift, and the second call must not need the
      // precompute to be repeated.
      for (int call = 0; call < 3; ++call) {
        interp.run(io, numPoints, persistent.data());
        const Mismatch diff = firstDifference(expected, outputs);
        INFO("trial=" << trial << " call=" << call << " index=" << diff.index);
        CHECK_FALSE(diff.any);
      }
    }
  }

  TEST_CASE("state accumulates across calls, in parallel assignment") {
    Program program;
    program.addInput("f", DT::F64);
    Arena& arena = program.arena();
    const NodeId acc = arena.field("acc");
    const NodeId prev = arena.field("prev");
    // Two mutually referring state slots: `acc` grows by f, `prev` records the
    // value acc had before this call. Sequential updates would give prev == acc.
    program.addState("acc", 1.0, arena.pw(Fn::Add, acc, arena.field("f")));
    program.addState("prev", 0.0, acc);
    program.addOutput("u", DT::F64, arena.pw(Fn::Mul, acc, prev));
    validate(program);

    const std::size_t numPoints = 37;
    std::mt19937 rng(4242);
    const std::vector<double> inputs = randomInputs<double>(rng, 1, numPoints);

    const LoweredProgram lowered = lower(program);
    REQUIRE(lowered.stateSlotCount() == 2);

    std::vector<double> outputs(numPoints, 0.0);
    std::vector<double> state(2 * numPoints, 0.0);
    initialiseState<double>(program, state.data(), numPoints);
    ArrayTileIo<double> io(inputs, outputs, numPoints, 1, 1, &state);
    InterpreterOptions options;
    options.tileSize = 5;
    TileInterpreter<double> interp(program, lowered, nullptr, options);

    // the states come through the tile source; nothing is hoisted
    REQUIRE(lowered.persistentSlotCount() == lowered.stateSlotCount());
    std::vector<double> persistent(1, 0.0);

    std::vector<double> referenceState(2 * numPoints, 0.0);
    initialiseState<double>(program, referenceState.data(), numPoints);

    for (int call = 0; call < 6; ++call) {
      interp.run(io, numPoints, persistent.data());
      const std::vector<double> expected =
          referenceEvaluate<double>(program, inputs, numPoints, nullptr, &referenceState);
      INFO("call=" << call);
      CHECK_FALSE(firstDifference(expected, outputs).any);
      CHECK_FALSE(firstDifference(referenceState, state).any);
    }

    // acc starts at 1 and gains f on every call; prev trails it by one.
    for (std::size_t p = 0; p < numPoints; ++p) {
      CHECK(state[p] == doctest::Approx(1.0 + 6.0 * inputs[p]));
      CHECK(state[numPoints + p] == doctest::Approx(1.0 + 5.0 * inputs[p]));
    }
  }

  TEST_CASE("tiling never crosses a partition boundary") {
    TestGrids grids;
    std::mt19937 rng(13579);
    const Program program = randomProgram(rng, 3, true, 1);
    validate(program);

    const std::size_t numPoints = 90;
    const std::vector<double> inputs =
        randomInputs<double>(rng, program.inputs().size(), numPoints);
    const std::vector<double> expected =
        referenceEvaluate<double>(program, inputs, numPoints, &grids, nullptr);

    const LoweredProgram lowered = lower(program);
    std::vector<double> outputs(numPoints, 0.0);
    ArrayTileIo<double> io(inputs, outputs, numPoints, program.inputs().size(), 1);
    InterpreterOptions options;
    options.tileSize = 16;
    TileInterpreter<double> interp(program, lowered, &grids, options);
    std::vector<double> persistent(
        static_cast<std::size_t>(lowered.persistentSlotCount()) * numPoints, 0.0);

    // Ranges deliberately misaligned against the tile width.
    const std::vector<PointRange> partitions{{0, 13}, {13, 50}, {50, 90}};
    interp.run(io, numPoints, persistent.data(), partitions);
    CHECK_FALSE(firstDifference(expected, outputs).any);
  }

  TEST_CASE("f32 programs run in f32") {
    Program program;
    program.setComputeType(ComputeType::F32);
    program.addInput("x", DT::F64);
    Arena& arena = program.arena();
    program.addOutput(
        "u", DT::F64, arena.pw(Fn::Mul, arena.field("x"), arena.pw(Fn::Sin, arena.field("x"))));
    validate(program);

    const std::size_t numPoints = 33;
    std::mt19937 rng(2468);
    const std::vector<float> inputs = randomInputs<float>(rng, 1, numPoints);
    const std::vector<float> expected =
        referenceEvaluate<float>(program, inputs, numPoints, nullptr, nullptr);

    const LoweredProgram lowered = lower(program);
    std::vector<float> outputs(numPoints, 0.0F);
    ArrayTileIo<float> io(inputs, outputs, numPoints, 1, 1);
    TileInterpreter<float> interp(program, lowered, nullptr, InterpreterOptions{});
    std::vector<float> persistent(
        static_cast<std::size_t>(lowered.persistentSlotCount()) * numPoints, 0.0F);
    interp.run(io, numPoints, persistent.data());
    CHECK_FALSE(firstDifference(expected, outputs).any);
  }

  TEST_CASE("a lookup at a time the window of a grid does not hold is reported") {
    // A grid over x and the time (axis 1) whose resident window holds four slices of its file, at
    // the times 10 to 11.5, with `before` slices of the file before it and `after` after it. A time
    // beyond the window where the file goes on is an error of the update (CheckKind::Window); one
    // beyond the file is clamped onto it, or, where the grid fails there, outside it.
    class WindowGrid final : public GridArrays {
  public:
      WindowGrid(std::size_t before, std::size_t after) {
        view_.type = reader::datafield::ElementType::Double;
        view_.dimensions = 2;
        view_.components = 1;
        view_.componentStride = 1;
        view_.min[0] = 0.0;
        view_.deltaInv[0] = 1.0;
        view_.num[0] = 3;
        view_.stride[0] = 1;
        view_.min[1] = 10.0;
        view_.deltaInv[1] = 2.0;
        view_.num[1] = 4;
        view_.stride[1] = 3;
        view_.before[1] = before;
        view_.after[1] = after;
        for (std::size_t i = 0; i < data_.size(); ++i) {
          data_[i] = static_cast<double>(i);
        }
        view_.data = data_.data();
      }
      [[nodiscard]] const reader::datafield::ArrayView& array(GridId /*grid*/) const override {
        return view_;
      }
      void report(GridId /*grid*/, CheckKind kind) const override { reports.push_back(kind); }

      mutable std::vector<CheckKind> reports;

  private:
      reader::datafield::ArrayView view_{};
      std::array<double, 12> data_{};
    };

    const auto sample = [](reader::datafield::BoundaryMode boundary,
                           std::size_t before,
                           std::size_t after,
                           double time,
                           std::vector<CheckKind>& reports,
                           reader::datafield::Interpolation interpolation =
                               reader::datafield::Interpolation::Linear) {
      Program program;
      program.addInput("x", DT::F64);
      program.addInput("t", DT::F64);
      Arena& arena = program.arena();
      reader::datafield::GridDesc desc;
      desc.path = "paddle.nc";
      desc.interpolation = interpolation;
      desc.boundary = boundary;
      desc.timeAxis = 1;
      const GridId grid = program.internGrid(desc);
      program.addOutput("u", DT::F64, arena.lookup(grid, 0, {arena.field("x"), arena.field("t")}));
      validate(program);

      const WindowGrid grids(before, after);
      const std::vector<double> inputs = {1.0, time};
      std::vector<double> outputs(1, 0.0);
      ArrayTileIo<double> io(inputs, outputs, 1, 2, 1);
      const LoweredProgram lowered = lower(program);
      TileInterpreter<double> interp(program, lowered, &grids, InterpreterOptions{});
      std::vector<double> persistent(
          std::max<std::size_t>(static_cast<std::size_t>(lowered.persistentSlotCount()), 1), 0.0);
      interp.run(io, 1, persistent.data());
      reports = grids.reports;
      return outputs[0];
    };

    using reader::datafield::BoundaryMode;
    std::vector<CheckKind> reports;
    // whether the checks reported just `kinds`
    const auto reported = [&reports](std::vector<CheckKind> kinds) { return reports == kinds; };
    for (const auto boundary : {BoundaryMode::Clamp, BoundaryMode::Fail}) {
      CAPTURE(static_cast<int>(boundary));
      const bool fails = boundary == BoundaryMode::Fail;
      // within the window: the sample, linear in time between the slices 1 and 2 at x = 1
      CHECK(sample(boundary, 2, 2, 10.75, reports) == 0.5 * 4.0 + 0.5 * 7.0);
      CHECK(reports.empty());
      // on its last node; however little beyond it, where the file goes on, it is reported
      CHECK(sample(boundary, 2, 2, 11.5, reports) == 10.0);
      CHECK(reports.empty());
      sample(boundary, 2, 2, 11.5 + 1e-9, reports);
      CHECK(reported({CheckKind::Window}));
      sample(boundary, 2, 2, 11.75, reports);
      CHECK(reported({CheckKind::Window}));
      sample(boundary, 2, 2, 9.75, reports);
      CHECK(reported({CheckKind::Window}));
      sample(boundary, 0, 2, 11.75, reports);
      CHECK(reported({CheckKind::Window}));
      // in its first cell as well, which the times a window is placed for never reach: it keeps a
      // slice to spare before them (TimeWindowSlack), and the check costs no binary64 operation
      sample(boundary, 2, 2, 10.0, reports);
      CHECK(reported({CheckKind::Window}));
      sample(boundary, 2, 2, 10.25, reports);
      CHECK(reported({CheckKind::Window}));
      CHECK(sample(boundary, 2, 2, 10.5, reports) == 4.0);
      CHECK(reports.empty());
      // beyond the file: clamped, or outside a grid that fails there
      const std::vector<CheckKind> beyond =
          fails ? std::vector<CheckKind>{CheckKind::Outside} : std::vector<CheckKind>{};
      CHECK(sample(boundary, 2, 0, 11.75, reports) == 10.0);
      CHECK(reported(beyond));
      CHECK(sample(boundary, 0, 2, 9.75, reports) == 1.0);
      CHECK(reported(beyond));

      // the cubic stencil reaches a node back and one on: the window serves the cell between its
      // second and its third node (10.5 to 11, the third itself reported as well), and from its
      // ends on where the file ends there
      const auto cubic = reader::datafield::Interpolation::Cubic;
      sample(boundary, 2, 2, 10.5, reports, cubic);
      CHECK(reports.empty());
      sample(boundary, 2, 2, 10.75, reports, cubic);
      CHECK(reports.empty());
      sample(boundary, 2, 2, 11.0, reports, cubic);
      CHECK(reported({CheckKind::Window}));
      sample(boundary, 2, 2, 11.25, reports, cubic);
      CHECK(reported({CheckKind::Window}));
      sample(boundary, 2, 2, 10.25, reports, cubic);
      CHECK(reported({CheckKind::Window}));
      sample(boundary, 0, 0, 10.25, reports, cubic);
      CHECK(reports.empty());
      sample(boundary, 0, 0, 11.25, reports, cubic);
      CHECK(reports.empty());
      // the nearest node: up to the last node of the window, not half a cell beyond it
      const auto nearest = reader::datafield::Interpolation::Nearest;
      CHECK(sample(boundary, 2, 2, 11.5, reports, nearest) == 10.0);
      CHECK(reports.empty());
      sample(boundary, 2, 2, 11.6, reports, nearest);
      CHECK(reported({CheckKind::Window}));
    }
  }

  TEST_CASE("a point that does not meet a check is reported, with its inputs") {
    // where its condition is zero; a NaN is true, as every value other than zero
    class Reporter final : public GridArrays {
  public:
      [[nodiscard]] const reader::datafield::ArrayView& array(GridId /*grid*/) const override {
        throw std::logic_error("no grid here");
      }
      void report(GridId /*grid*/, CheckKind /*kind*/) const override {
        FAIL("no grid here is checked");
      }
      void failed(std::int32_t check,
                  std::size_t point,
                  const std::vector<double>& inputs) const override {
        failures.push_back({check, point, inputs});
      }
      struct Failure {
        std::int32_t check;
        std::size_t point;
        std::vector<double> inputs;
      };
      mutable std::vector<Failure> failures;
    };

    Program program;
    program.addInput("x", DT::F64);
    program.addInput("y", DT::F64);
    Arena& arena = program.arena();
    const NodeId x = arena.field("x");
    const NodeId y = arena.field("y");
    program.addOutput("u", DT::F64, arena.pw(Fn::Add, x, y));
    program.addCheck(arena.pw(Fn::Le, x, arena.konst(1.0)), "x at most 1");
    program.addCheck(arena.pw(Fn::Sub, y, arena.konst(3.0)), "y not 3");
    validate(program);

    constexpr double Nan = std::numeric_limits<double>::quiet_NaN();
    const std::size_t numPoints = 5;
    // input-major: x of every point, then y
    const std::vector<double> inputs = {0.5, 2.0, 0.7, 1.0, Nan, 1.0, 3.0, 3.0, Nan, 2.0};
    std::vector<double> outputs(numPoints, 0.0);
    ArrayTileIo<double> io(inputs, outputs, numPoints, 2, 1);
    const LoweredProgram lowered = lower(program);
    const Reporter reporter;
    InterpreterOptions options;
    options.tileSize = 2;
    TileInterpreter<double> interp(program, lowered, &reporter, options);
    std::vector<double> persistent(1, 0.0);
    interp.run(io, numPoints, persistent.data());

    // per tile and check, the first point that fails it: in the tile of points 0 and 1, point 1
    // fails both (x = 2, y = 3); in that of points 2 and 3, point 2 the second (y = 3), while
    // y - 3 = NaN at point 3 is true; and point 4 the first, le(NaN, 1) being 0
    REQUIRE(reporter.failures.size() == 4);
    CHECK(reporter.failures[0].check == 0);
    CHECK(reporter.failures[0].point == 1);
    CHECK(reporter.failures[0].inputs == std::vector<double>{2.0, 3.0});
    CHECK(reporter.failures[1].check == 1);
    CHECK(reporter.failures[1].point == 1);
    CHECK(reporter.failures[2].check == 1);
    CHECK(reporter.failures[2].point == 2);
    CHECK(reporter.failures[2].inputs == std::vector<double>{0.7, 3.0});
    CHECK(reporter.failures[3].check == 0);
    CHECK(reporter.failures[3].point == 4);
  }

  TEST_CASE("a perturbed evaluation moves products, and marks the points whose selects flip") {
    // u = select(lt(x * y, 1), 1, 2): at x = 1 + 2^-27, y = 1 - 2^-27 the product rounds to 1, a
    // product moved down by an ulp is below it; at x = 2, y = 3 no ulp decides
    Program program;
    program.addInput("x", DT::F64);
    program.addInput("y", DT::F64);
    Arena& arena = program.arena();
    const NodeId product = arena.pw(Fn::Mul, arena.field("x"), arena.field("y"));
    program.addOutput("u",
                      DT::F64,
                      arena.pw(Fn::Select,
                               arena.pw(Fn::Lt, product, arena.konst(1.0)),
                               arena.konst(1.0),
                               arena.konst(2.0)));
    program.addOutput("p", DT::F64, product);
    validate(program);
    const std::vector<double> inputs = {
        1.0 + std::ldexp(1.0, -27), 2.0, 1.0 - std::ldexp(1.0, -27), 3.0};
    const LoweredProgram lowered = lower(program);
    const auto evaluate = [&](const Perturbation& perturbation,
                              std::vector<std::uint64_t>& decisions) {
      std::vector<double> outputs(4, 0.0);
      ArrayTileIo<double> io(inputs, outputs, 2, 2, 2);
      InterpreterOptions options;
      options.perturbation = perturbation;
      options.decisions = &decisions;
      TileInterpreter<double> interp(program, lowered, nullptr, options);
      std::vector<double> persistent(1, 0.0);
      interp.run(io, 2, persistent.data());
      return outputs;
    };
    std::vector<std::uint64_t> reference(2, 0);
    const auto exact = evaluate(Perturbation{}, reference);
    CHECK(exact[0] == 2.0);
    CHECK(exact[2] == 1.0);
    bool flipped = false;
    for (std::uint64_t seed = 1; seed <= 8; ++seed) {
      std::vector<std::uint64_t> decisions(2, 0);
      const auto moved = evaluate(Perturbation{true, seed, 1, 4}, decisions);
      // the product, an ulp off, the same way for the same seed
      CHECK(std::fabs(moved[2] - 1.0) <= std::ldexp(1.0, -52));
      CHECK(moved[2] != 1.0);
      CHECK(std::fabs(moved[3] - 6.0) <= std::ldexp(1.0, -50));
      std::vector<std::uint64_t> again(2, 0);
      CHECK(evaluate(Perturbation{true, seed, 1, 4}, again)[2] == moved[2]);
      // the select of the second point never flips; that of the first where the product is below 1
      CHECK(decisions[1] == reference[1]);
      CHECK((decisions[0] != reference[0]) == (moved[2] < 1.0));
      flipped = flipped || decisions[0] != reference[0];
    }
    CHECK(flipped);
  }

  TEST_CASE("a select whose value takes the slot of its condition records what it decided") {
    // w = select(lt(x * y, 1), v, u): the value of the select takes the slot of its condition, and
    // the decision is the condition's, not the value's (which is not zero either way). At x = 1 +
    // 2^-27, y = 1 - 2^-27 the product rounds to 1; a product moved down by an ulp is below it.
    const Program program = loadCore("#!sderiv 1\nin x, y, z\ndef p = x * y\ndef u = p + z\n"
                                     "def s = sin(x)\ndef v = s * y\ndef c = lt(p, 1)\n"
                                     "def w = select(c, v, u)\nout u = u\nout v = v\nout w = w\n",
                                     "");
    const std::vector<double> inputs = {
        1.0 + std::ldexp(1.0, -27), 1.0 - std::ldexp(1.0, -27), 0.25};
    const LoweredProgram lowered = lower(program);
    REQUIRE(lowered.dump().find("t0:f64 = pw select t0") != std::string::npos);
    const auto evaluate = [&](const Perturbation& perturbation,
                              std::vector<std::uint64_t>& decisions,
                              std::map<std::size_t, std::vector<std::uint8_t>>* branches) {
      std::vector<double> outputs(3, 0.0);
      ArrayTileIo<double> io(inputs, outputs, 1, 3, 3);
      InterpreterOptions options;
      options.perturbation = perturbation;
      options.decisions = &decisions;
      options.branches = branches;
      TileInterpreter<double> interp(program, lowered, nullptr, options);
      std::vector<double> persistent(1, 0.0);
      interp.run(io, 1, persistent.data());
      return outputs;
    };
    std::vector<std::uint64_t> reference(1, 0);
    std::map<std::size_t, std::vector<std::uint8_t>> branches;
    const auto exact = evaluate(Perturbation{}, reference, &branches);
    CHECK(exact[2] == exact[0]);
    REQUIRE(branches.size() == 1);
    CHECK(branches.begin()->second == std::vector<std::uint8_t>{0});
    bool flipped = false;
    for (std::uint64_t run = 1; run <= 8; ++run) {
      std::vector<std::uint64_t> decisions(1, 0);
      const auto moved = evaluate(Perturbation{true, splitMix64(run), 1, 4}, decisions, nullptr);
      const bool other = moved[2] != moved[0];
      CHECK(other == (decisions[0] != reference[0]));
      flipped = flipped || other;
    }
    CHECK(flipped);
  }

  TEST_CASE("the interpreter computes as IEEE 754 says") {
    // no -ffast-math: NaN, fmax, division, contraction, infinities, signed zero
    CHECK(interpreterFault().empty());
  }

  TEST_CASE("a program with a lookup but no grid arrays is rejected, not ignored") {
    Program program;
    program.addInput("x", DT::F64);
    Arena& arena = program.arena();
    const GridId grid = program.internGrid(reader::datafield::GridDesc{});
    program.addOutput("u", DT::F64, arena.lookup(grid, 0, {arena.field("x")}));
    validate(program);

    const std::size_t numPoints = 4;
    const std::vector<double> inputs(numPoints, 1.0);
    std::vector<double> outputs(numPoints, 0.0);
    ArrayTileIo<double> io(inputs, outputs, numPoints, 1, 1);
    const LoweredProgram lowered = lower(program);
    TileInterpreter<double> interp(program, lowered, nullptr, InterpreterOptions{});
    std::vector<double> persistent(1, 0.0);
    CHECK_THROWS_AS(interp.run(io, numPoints, persistent.data()), std::invalid_argument);
  }
}

} // namespace seissol::expr::test
