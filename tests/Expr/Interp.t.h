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

#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"

#include <array>
#include <cstring>
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

// A pure elementwise sampler: whatever it computes, it must compute the same
// thing for count = 1 and count = 256, so the batched and scalar paths agree.
//
// It accumulates one dimension at a time, over the destination buffer, because
// that is how a separable interpolation stencil is written and it is what makes
// `dst` aliasing a coordinate array unsafe. The interpreter therefore must not
// give a Lookup a destination slot recycled from one of its own coordinate
// operands; this sampler is what turns that rule into a failing test rather than
// a comment.
class ToySampler : public GridSampler {
  public:
  void sampleBatch(GridId grid,
                   std::int32_t component,
                   const double* const* coords,
                   std::int32_t dimension,
                   std::size_t count,
                   double* dst) const override {
    apply(grid, component, coords, dimension, count, dst);
  }
  void sampleBatch(GridId grid,
                   std::int32_t component,
                   const float* const* coords,
                   std::int32_t dimension,
                   std::size_t count,
                   float* dst) const override {
    apply(grid, component, coords, dimension, count, dst);
  }

  private:
  template <typename T>
  static void apply(GridId grid,
                    std::int32_t component,
                    const T* const* coords,
                    std::int32_t dimension,
                    std::size_t count,
                    T* dst) {
    for (std::size_t l = 0; l < count; ++l) {
      dst[l] = static_cast<T>(grid) * T(0.25) + static_cast<T>(component) * T(0.125) + coords[0][l];
    }
    for (std::int32_t d = 1; d < dimension; ++d) {
      for (std::size_t l = 0; l < count; ++l) {
        dst[l] = dst[l] + coords[d][l] * static_cast<T>(d + 1);
      }
    }
  }
};

// Plain SoA arrays behind the Binding::gather / Binding::scatter contract.
template <typename T>
class ArrayTileIo : public TileIo<T> {
  public:
  ArrayTileIo(const std::vector<T>& inputs,
              std::vector<T>& outputs,
              std::size_t numPoints,
              std::size_t numInputs,
              std::size_t numOutputs)
      : inputs_(inputs), outputs_(outputs), numPoints_(numPoints), numInputs_(numInputs),
        numOutputs_(numOutputs) {}

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

  private:
  const std::vector<T>& inputs_;
  std::vector<T>& outputs_;
  std::size_t numPoints_;
  std::size_t numInputs_;
  std::size_t numOutputs_;
};

// Point-at-a-time evaluation of the DAG. Arena ids are topologically ordered, so
// one ascending sweep per point suffices and no recursion is involved.
template <typename T>
std::vector<T> referenceEvaluate(const Program& program,
                                 const std::vector<T>& inputs,
                                 std::size_t numPoints,
                                 const GridSampler* sampler,
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
        REQUIRE(sampler != nullptr);
        arena.children(id, kids);
        std::array<T, MaxLookupDimension> scalars{};
        std::array<const T*, MaxLookupDimension> pointers{};
        for (std::size_t k = 0; k < kids.size(); ++k) {
          scalars[k] = value[kids[k]];
          pointers[k] = &scalars[k];
        }
        T out = T(0);
        sampler->sampleBatch(
            node.grid, node.comp, pointers.data(), static_cast<std::int32_t>(kids.size()), 1, &out);
        value[id] = out;
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
  const GridId grid = withLookup ? program.internGrid(datafield::GridDesc{}) : NoGrid;

  std::vector<NodeId> leaves{arena.field("x"),
                             arena.field("y"),
                             arena.field("z"),
                             arena.field("t"),
                             arena.konst(0.5),
                             arena.konst(-1.25),
                             arena.konst(2.0)};
  if (withLookup) {
    leaves.push_back(arena.lookup(grid, 0, {arena.field("x"), arena.field("y")}));
    leaves.push_back(arena.lookup(grid, 1, {arena.field("z")}));
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

  TEST_CASE("matches the scalar reference bit for bit") {
    ToySampler sampler;
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
            referenceEvaluate<double>(program, inputs, numPoints, &sampler, nullptr);

        for (const std::size_t tile :
             {std::size_t{1}, std::size_t{3}, std::size_t{8}, std::size_t{64}}) {
          const LoweredProgram lowered = lower(program);
          std::vector<double> outputs(program.outputs().size() * numPoints, 0.0);
          ArrayTileIo<double> io(
              inputs, outputs, numPoints, program.inputs().size(), program.outputs().size());
          InterpreterOptions options;
          options.tileSize = tile;
          TileInterpreter<double> interp(program, lowered, &sampler, options);
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
    ToySampler sampler;
    std::mt19937 rng(987654321);

    for (int trial = 0; trial < 30; ++trial) {
      const Program program = randomProgram(rng, 4, true, 2);
      validate(program);
      const std::size_t numPoints = 173;
      const std::vector<double> inputs =
          randomInputs<double>(rng, program.inputs().size(), numPoints);
      const std::vector<double> expected =
          referenceEvaluate<double>(program, inputs, numPoints, &sampler, nullptr);

      LowerOptions options;
      options.invariantInputs = {"x", "y", "z"};
      options.invariantGrids = {0};
      const LoweredProgram lowered = lower(program, options);
      REQUIRE(lowered.hasPrecompute());

      std::vector<double> outputs(program.outputs().size() * numPoints, 0.0);
      ArrayTileIo<double> io(
          inputs, outputs, numPoints, program.inputs().size(), program.outputs().size());
      InterpreterOptions interpOptions;
      interpOptions.tileSize = 24;
      TileInterpreter<double> interp(program, lowered, &sampler, interpOptions);
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
    ArrayTileIo<double> io(inputs, outputs, numPoints, 1, 1);
    InterpreterOptions options;
    options.tileSize = 5;
    TileInterpreter<double> interp(program, lowered, nullptr, options);

    std::vector<double> persistent(
        static_cast<std::size_t>(lowered.persistentSlotCount()) * numPoints, 0.0);
    initialiseState<double>(program, persistent.data(), numPoints);

    std::vector<double> referenceState(2 * numPoints, 0.0);
    initialiseState<double>(program, referenceState.data(), numPoints);

    for (int call = 0; call < 6; ++call) {
      interp.run(io, numPoints, persistent.data());
      const std::vector<double> expected =
          referenceEvaluate<double>(program, inputs, numPoints, nullptr, &referenceState);
      INFO("call=" << call);
      CHECK_FALSE(firstDifference(expected, outputs).any);
      CHECK_FALSE(firstDifference(referenceState, persistent).any);
    }

    // acc starts at 1 and gains f on every call; prev trails it by one.
    for (std::size_t p = 0; p < numPoints; ++p) {
      CHECK(persistent[p] == doctest::Approx(1.0 + 6.0 * inputs[p]));
      CHECK(persistent[numPoints + p] == doctest::Approx(1.0 + 5.0 * inputs[p]));
    }
  }

  TEST_CASE("tiling never crosses a partition boundary") {
    ToySampler sampler;
    std::mt19937 rng(13579);
    const Program program = randomProgram(rng, 3, true, 1);
    validate(program);

    const std::size_t numPoints = 90;
    const std::vector<double> inputs =
        randomInputs<double>(rng, program.inputs().size(), numPoints);
    const std::vector<double> expected =
        referenceEvaluate<double>(program, inputs, numPoints, &sampler, nullptr);

    const LoweredProgram lowered = lower(program);
    std::vector<double> outputs(numPoints, 0.0);
    ArrayTileIo<double> io(inputs, outputs, numPoints, program.inputs().size(), 1);
    InterpreterOptions options;
    options.tileSize = 16;
    TileInterpreter<double> interp(program, lowered, &sampler, options);
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

  TEST_CASE("a program with a lookup but no sampler is rejected, not ignored") {
    Program program;
    program.addInput("x", DT::F64);
    Arena& arena = program.arena();
    const GridId grid = program.internGrid(datafield::GridDesc{});
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
