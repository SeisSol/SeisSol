// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Branches.h"
#include "Expr/Core.h"
#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/RtcCpu.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/InlineGrid.h"
#include "Reader/Scripting/DataTable.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace {

namespace branches {

using DT = reader::scripting::DataType;
constexpr double Nan = std::numeric_limits<double>::quiet_NaN();
constexpr double Inf = std::numeric_limits<double>::infinity();

/// What `roots` need of the arena of `program`, in ascending order.
std::vector<NodeId> needed(const Program& program, const std::vector<NodeId>& roots) {
  const Arena& arena = program.arena();
  std::vector<char> seen(arena.size(), 0);
  std::vector<NodeId> stack(roots.begin(), roots.end());
  std::vector<NodeId> kids;
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (seen[static_cast<std::size_t>(id)] != 0) {
      continue;
    }
    seen[static_cast<std::size_t>(id)] = 1;
    arena.children(id, kids);
    stack.insert(stack.end(), kids.begin(), kids.end());
  }
  std::vector<NodeId> values;
  for (std::size_t id = 0; id < seen.size(); ++id) {
    if (seen[id] != 0) {
      values.push_back(static_cast<NodeId>(id));
    }
  }
  return values;
}

/// The inline grids of a program, as a kernel reads their arrays, and the checks they report.
class InlineGrids final : public GridArrays {
  public:
  explicit InlineGrids(const Program& program) {
    for (const auto& desc : program.grids()) {
      grids_.push_back(std::make_unique<reader::datafield::InlineGrid>(
          desc.data, std::max<std::size_t>(desc.componentNames.size(), 1)));
    }
  }
  [[nodiscard]] const reader::datafield::ArrayView& array(GridId grid) const override {
    return *grids_.at(static_cast<std::size_t>(grid))->arrayView();
  }
  void report(GridId /*grid*/, CheckKind kind) const override { reports.push_back(kind); }

  mutable std::vector<CheckKind> reports;

  private:
  std::vector<std::unique_ptr<reader::datafield::InlineGrid>> grids_;
};

/// Inputs and outputs of the interpreter, [input][point] and [output][point].
class Columns final : public TileIo<double> {
  public:
  Columns(const std::vector<double>& inputs,
          std::vector<double>& outputs,
          std::size_t points,
          std::size_t numInputs,
          std::size_t numOutputs)
      : inputs_(inputs), outputs_(outputs), points_(points), numInputs_(numInputs),
        numOutputs_(numOutputs) {}
  void gather(std::size_t first, std::size_t count, double* dst) const override {
    for (std::size_t i = 0; i < numInputs_; ++i) {
      for (std::size_t l = 0; l < count; ++l) {
        dst[i * count + l] = inputs_[i * points_ + first + l];
      }
    }
  }
  void scatter(std::size_t first, std::size_t count, const double* src) override {
    for (std::size_t o = 0; o < numOutputs_; ++o) {
      for (std::size_t l = 0; l < count; ++l) {
        outputs_[o * points_ + first + l] = src[o * count + l];
      }
    }
  }

  private:
  const std::vector<double>& inputs_;
  std::vector<double>& outputs_;
  std::size_t points_;
  std::size_t numInputs_;
  std::size_t numOutputs_;
};

/// The outputs of `program` by the interpreter at `points` points of `inputs` ([input][point]),
/// lowered with the threshold `threshold` for the branches, in tiles of `tile` points; the checks
/// the grids reported in `reports`.
std::vector<double> interpret(const Program& program,
                              const std::vector<double>& inputs,
                              std::size_t points,
                              std::int32_t threshold,
                              std::size_t tile,
                              std::vector<CheckKind>* reports = nullptr) {
  LowerOptions options;
  options.branchThreshold = threshold;
  const LoweredProgram lowered = lower(program, options);
  const InlineGrids grids(program);
  std::vector<double> outputs(program.outputs().size() * points, -1.0);
  Columns io(inputs, outputs, points, program.inputs().size(), program.outputs().size());
  InterpreterOptions interpreterOptions;
  interpreterOptions.tileSize = tile;
  TileInterpreter<double> interpreter(program, lowered, &grids, interpreterOptions);
  std::vector<double> persistent(
      std::max<std::size_t>(static_cast<std::size_t>(lowered.persistentSlotCount()) * points, 1),
      0.0);
  interpreter.precompute(io, points, persistent.data());
  interpreter.run(io, points, persistent.data());
  if (reports != nullptr) {
    *reports = grids.reports;
  }
  return outputs;
}

/// The outputs of `program` at `points` points of `inputs` on `backend` in the shape `shape`, with
/// the lowering's default branches; whether that backend ran.
bool evaluate(const Program& program,
              const std::vector<double>& inputs,
              std::size_t points,
              BackendKind backend,
              KernelShape shape,
              std::vector<double>& outputs) {
  outputs.assign(program.outputs().size() * points, -1.0);
  reader::scripting::DataTable table(points);
  for (std::size_t i = 0; i < program.inputs().size(); ++i) {
    table.bindViewConst<double>(
        program.inputs()[i].name, reader::scripting::Direction::In, inputs.data() + i * points);
  }
  for (std::size_t o = 0; o < program.outputs().size(); ++o) {
    table.bindView<double>(
        program.outputs()[o].name, reader::scripting::Direction::Out, outputs.data() + o * points);
  }
  Binding binding = Binding::bind(program, table);
  reader::datafield::GridStore store;
  BackendOptions options;
  options.quiet = true;
  options.preferred = backend;
  options.shape = shape;
  const auto kernel = makeKernel(program, binding, store, options);
  kernel->precompute(table);
  kernel->run(table);
  return kernel->kind() == backend;
}

/// The same bits, or both NaN: which NaN an operation of two NaNs gives depends on the order of its
/// operands in the machine code, as the backends compare (compareAt()).
bool sameBits(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (std::size_t k = 0; k < a.size(); ++k) {
    if (std::memcmp(&a[k], &b[k], sizeof(double)) != 0 && !(std::isnan(a[k]) && std::isnan(b[k]))) {
      return false;
    }
  }
  return true;
}

/// A grid of x and z by its knots, of two components, linear on both axes, which fails beyond
/// them; and a program that samples it in branches of selects, nested, with a value that weighs
/// enough in the other branch, and a select of each side alone.
const char* const Dispatch =
    "#!sderiv 1\nin x, z\n"
    "dataset g (kind = \"inline\", components = (\"v\", \"w\"), knots = ((0, 1, 3), (10, 20)), "
    "values = (((0, 5), (1, 6), (3, 7)), ((100, 50), (101, 51), (103, 53))), "
    "interpolation = \"linear\", boundary = \"%s\")\n"
    "def c = lt(x, 2)\n"
    "def d = lt(z, 15)\n"
    "def a = g.v(x, z)\n"
    "def b = exp(x / 10) * sin(z) + log(abs(z) + 1) * cos(x)\n"
    "def e = g.w(x, z)\n"
    "out u = select(c, select(d, a, a * e + 1), b)\n"
    "out v = select(c, a + 2, x)\n"
    "out w = select(lt(z, 30), x * z, sqrt(abs(x)) * exp(z / 100) + log(abs(x) + 2))\n";

Program dispatch(const char* boundary) {
  std::string text(Dispatch);
  text.replace(text.find("%s"), 2, boundary);
  return loadCore(text, "");
}

/// Points around the grid of Dispatch: in either branch, inside and beyond the knots, NaN and the
/// infinities; and runs of points that take the same branches, so that tiles of a few take them
/// all alike or not.
void dispatchPoints(std::vector<double>& inputs, std::size_t& points) {
  std::vector<double> x;
  std::vector<double> z;
  for (int i = 0; i < 24; ++i) {
    x.push_back(0.1 + 0.05 * i);
    z.push_back(11.0 + 0.3 * i);
  }
  for (int i = 0; i < 24; ++i) {
    x.push_back(2.5 + 0.1 * i);
    z.push_back(12.0 + 0.7 * i);
  }
  const std::vector<double> specials = {Nan, Inf, -Inf, 0.0, -0.0, 1.0, 2.0, 3.0, 1.999, 2.001};
  for (std::size_t i = 0; i < specials.size(); ++i) {
    for (std::size_t j = 0; j < specials.size(); j += 3) {
      x.push_back(specials[i]);
      z.push_back(specials[j] + (j % 2 == 0 ? 15.0 : 40.0));
    }
  }
  for (int i = 0; i < 40; ++i) {
    // alternating branches
    x.push_back(i % 2 == 0 ? 0.7 + 0.01 * i : 2.7 + 0.02 * i);
    z.push_back(i % 3 == 0 ? 12.5 + 0.1 * i : 17.5 - 0.05 * i);
  }
  points = x.size();
  inputs = x;
  inputs.insert(inputs.end(), z.begin(), z.end());
}

} // namespace branches

} // namespace

TEST_SUITE("expr branches") {
  using namespace branches;

  TEST_CASE("a value one branch of a select uses lies in that branch") {
    Program program;
    program.addInput("x", DT::F64);
    Arena& a = program.arena();
    const NodeId x = a.field("x");
    const NodeId c = a.pw(Fn::Lt, x, a.konst(0.0));
    const NodeId e = a.pw(Fn::Exp, x);
    const NodeId l = a.pw(Fn::Log, x);
    const NodeId shared = a.pw(Fn::Sqrt, x);
    const NodeId first = a.pw(Fn::Add, e, shared);
    const NodeId second = a.pw(Fn::Add, l, shared);
    const NodeId s = a.pw(Fn::Select, c, first, second);
    // another select of the same condition: it shares the branches
    const NodeId twice = a.pw(Fn::Mul, e, a.konst(2.0));
    const NodeId t = a.pw(Fn::Select, c, twice, a.konst(1.0));
    // and one of another condition, of what the first branch has
    const NodeId d = a.pw(Fn::Lt, x, a.konst(-1.0));
    const NodeId inner = a.pw(Fn::Sin, e);
    const NodeId r = a.pw(Fn::Select, c, a.pw(Fn::Select, d, inner, e), a.konst(0.0));
    const std::vector<NodeId> roots = {s, t, r};
    const Regions regions(a, needed(program, roots), roots);

    CHECK(regions.of(x) == 0);
    CHECK(regions.of(c) == 0);
    CHECK(regions.of(shared) == 0);
    const std::int32_t taken = regions.of(e);
    REQUIRE(taken > 0);
    CHECK(regions[taken].condition == c);
    CHECK(regions[taken].side);
    CHECK(regions[taken].parent == 0);
    CHECK(regions.of(first) == taken);
    CHECK(regions.of(twice) == taken);
    CHECK(regions.of(d) == taken);
    const std::int32_t other = regions.of(l);
    REQUIRE(other > 0);
    CHECK(regions[other].condition == c);
    CHECK_FALSE(regions[other].side);
    CHECK(regions.of(second) == other);
    const std::int32_t nested = regions.of(inner);
    REQUIRE(nested > 0);
    CHECK(regions[nested].parent == taken);
    CHECK(regions[nested].condition == d);
    CHECK(regions.common(nested, other) == 0);
    CHECK(regions.common(nested, taken) == taken);
  }

  TEST_CASE("a heavy branch of a select becomes a branch of the code, a light one does not") {
    Program program;
    program.addInput("x", DT::F64);
    Arena& a = program.arena();
    const NodeId x = a.field("x");
    const NodeId heavy = a.pw(Fn::Add,
                              a.pw(Fn::Mul, a.pw(Fn::Exp, x), a.pw(Fn::Sin, x)),
                              a.pw(Fn::Log, a.pw(Fn::Abs, x)));
    const NodeId light = a.pw(Fn::Add, x, a.konst(1.0));
    program.addOutput("u", DT::F64, a.pw(Fn::Select, a.pw(Fn::Lt, x, a.konst(0.0)), heavy, light));
    validate(program);

    const LoweredProgram branched = lower(program);
    REQUIRE(branched.run().branches.size() == 1);
    const Branch& branch = branched.run().branches[0];
    CHECK(branch.side);
    CHECK(branch.end - branch.begin == 6);
    // the condition before the branch, the select after it
    const auto& code = branched.run().code;
    CHECK(code[static_cast<std::size_t>(branch.end)].fn == Fn::Select);
    CHECK(branched.dump().find("if t") != std::string::npos);

    LowerOptions none;
    none.branchThreshold = std::numeric_limits<std::int32_t>::max();
    const LoweredProgram straight = lower(program, none);
    CHECK(straight.run().branches.empty());
    CHECK(straight.run().code.size() == code.size());
    // the options shape the kernel
    CHECK(none.fingerprint() != LowerOptions{}.fingerprint());
  }

  TEST_CASE("a dispatch over layers starts each branch once, and keeps as few values alive") {
    // a grid per layer, each sampled in its layer: a branch per layer and one per rest, and none in
    // the lookups, whose selects choose between the cases at the ends of an axis
    std::string text = "#!sderiv 1\nin x, z\n";
    constexpr int Layers = 3;
    std::string value = "0";
    for (int l = 0; l < Layers; ++l) {
      const std::string base = std::to_string(10 * l);
      text += "dataset g" + std::to_string(l) +
              " (kind = \"inline\", components = (\"v\"), knots = ((0, 1, 3), (10, 20)), values = "
              "(((" +
              base +
              "), (1), (3)), ((100), (101), (103))), interpolation = \"linear\", "
              "boundary = \"clamp\")\n";
    }
    for (int l = Layers - 1; l >= 0; --l) {
      value = "select(lt(z, " + std::to_string(12 + 3 * l) + "), g" + std::to_string(l) +
              ".v(x, z) * exp(x) + sin(z), " + value + ")";
    }
    text += "out u = " + value + "\n";
    const Program program = loadCore(text, "");
    REQUIRE(program.grids().size() == static_cast<std::size_t>(Layers));
    const LoweredProgram branched = lower(program);
    LowerOptions none;
    none.branchThreshold = std::numeric_limits<std::int32_t>::max();
    const LoweredProgram straight = lower(program, none);
    CHECK(branched.run().branches.size() == 2 * Layers - 1);
    CHECK(branched.run().slotCount <= straight.run().slotCount);
    CHECK(branched.run().code.size() == straight.run().code.size());
  }

  TEST_CASE("branches give the bits of computing every value, in every tile") {
    const Program program = dispatch("clamp");
    std::vector<double> inputs;
    std::size_t points = 0;
    dispatchPoints(inputs, points);
    REQUIRE_FALSE(lower(program).run().branches.empty());
    const auto everything =
        interpret(program, inputs, points, std::numeric_limits<std::int32_t>::max(), 64);
    for (const std::size_t tile : {1UL, 3UL, 8UL, 64UL, 256UL}) {
      CAPTURE(tile);
      CHECK(sameBits(interpret(program, inputs, points, DefaultBranchThreshold, tile), everything));
      CHECK(sameBits(interpret(program, inputs, points, 0, tile), everything));
    }

    // the compiled CPU kernel, one point after the other and in blocks of points
    for (const KernelShape shape : {KernelShape::Point, KernelShape::Tile}) {
      CAPTURE(static_cast<int>(shape));
      std::vector<double> compiled;
      if (!evaluate(program, inputs, points, BackendKind::RtcCpu, shape, compiled)) {
        WARN_MESSAGE(false, "no usable C++ compiler; the compiled CPU kernel was not run");
        continue;
      }
      CHECK(sameBits(compiled, everything));
    }
  }

  TEST_CASE("the compiled CPU kernel runs a branch in an if") {
    const Program program = dispatch("clamp");
    const LoweredProgram lowered = lower(program);
    const std::string point = emitCpuSource(program, lowered, KernelShape::Point);
    CHECK(point.find("if (s") != std::string::npos);
    CHECK(point.find("} else {") != std::string::npos);
    const std::string tile = emitCpuSource(program, lowered, KernelShape::Tile);
    CHECK(tile.find("reduction(|:b0)") != std::string::npos);
    CHECK(tile.find("[32] = {0};") != std::string::npos);
  }

  TEST_CASE("a sample in a branch not taken neither fails nor counts beyond its grid") {
    const Program program = dispatch("fail");
    // one point at a time, and points of either branch in one tile, each with and without branches
    const auto reported =
        [&](const std::vector<double>& x, const std::vector<double>& z, std::int32_t threshold) {
          std::vector<double> inputs = x;
          inputs.insert(inputs.end(), z.begin(), z.end());
          std::vector<CheckKind> reports;
          interpret(program, inputs, x.size(), threshold, 64, &reports);
          return !reports.empty();
        };
    for (const std::int32_t threshold :
         {DefaultBranchThreshold, std::numeric_limits<std::int32_t>::max()}) {
      CAPTURE(threshold);
      // inside, and taken
      CHECK_FALSE(reported({1.0}, {15.0}, threshold));
      // beyond the knots of z, where the branch is taken
      CHECK(reported({1.0}, {25.0}, threshold));
      CHECK(reported({1.0}, {Nan}, threshold));
      // and where it is not: x >= 2, or NaN, which `lt` takes as false
      CHECK_FALSE(reported({3.0}, {25.0}, threshold));
      CHECK_FALSE(reported({Nan}, {25.0}, threshold));
      CHECK_FALSE(reported({Inf}, {-Inf}, threshold));
      // a tile of both: the branch runs for the point that does not take it as well
      CHECK_FALSE(reported({3.0, 1.0}, {25.0, 15.0}, threshold));
      CHECK(reported({3.0, 1.0}, {15.0, 25.0}, threshold));
    }

    // the compiled CPU kernel does not fail at a point beyond the grid that does not use the sample
    for (const KernelShape shape : {KernelShape::Point, KernelShape::Tile}) {
      CAPTURE(static_cast<int>(shape));
      const std::vector<double> inputs = {3.0, 4.0, 1.0, 25.0, -7.0, 15.0};
      std::vector<double> outputs;
      if (!evaluate(program, inputs, 3, BackendKind::RtcCpu, shape, outputs)) {
        WARN_MESSAGE(false, "no usable C++ compiler; the compiled CPU kernel was not run");
        continue;
      }
      // u and v at (1, 15): the sample 51 and the other component, 28.5
      CHECK(outputs[2] == 51.0 * 28.5 + 1.0);
      CHECK(outputs[5] == 53.0);
    }

    // nor does the evaluation record count it
    reader::datafield::GridStore store;
    const std::vector<double> inputs = {3.0, 1.0, 3.0, 25.0, 25.0, 15.0};
    const auto clamps = clampsAt(program, store, inputs, 3);
    REQUIRE(clamps.size() == 1);
    REQUIRE(clamps[0].axes.size() == 2);
    CHECK(clamps[0].axes[1][0] == 0);
    CHECK(clamps[0].axes[1][1] == GridClamps::Above);
    CHECK(clamps[0].axes[1][2] == 0);
  }

  TEST_CASE("a sample that branches of two conditions use counts where either is taken") {
    const Program program = loadCore(
        "#!sderiv 1\nin x, z\n"
        "dataset g (kind = \"inline\", components = (\"v\"), knots = ((0, 1, 3), (10, 20)), "
        "values = (((0), (1), (3)), ((100), (101), (103))), interpolation = \"linear\", "
        "boundary = \"fail\")\n"
        "def s = g.v(x, z)\n"
        "out p = select(lt(x, 0.5), s * exp(x) + sin(s), 0)\n"
        "out q = select(lt(z, 12), s * log(abs(z) + 1) + cos(s), 1)\n",
        "");
    const auto reported = [&](double x, double z) {
      std::vector<CheckKind> reports;
      interpret(program, {x, z}, 1, DefaultBranchThreshold, 64, &reports);
      return !reports.empty();
    };
    CHECK_FALSE(reported(1.0, 25.0));
    CHECK(reported(0.2, 25.0));
    CHECK(reported(1.0, 5.0));
    CHECK(reported(-1.0, 5.0));
    CHECK_FALSE(reported(1.0, 15.0));
    CHECK_FALSE(reported(Nan, 30.0));
  }

  TEST_CASE("the checks of a program and those of its samples hold without branches as well") {
    // the selects that restrict a check are those of the program: computed everywhere, the check
    // finds the same
    const Program program = dispatch("fail");
    std::vector<CheckKind> with;
    std::vector<CheckKind> without;
    const std::vector<double> inputs = {3.0, 1.0, Nan, 25.0, 15.0, 25.0};
    interpret(program, inputs, 3, DefaultBranchThreshold, 1, &with);
    interpret(program, inputs, 3, std::numeric_limits<std::int32_t>::max(), 1, &without);
    CHECK(with.empty());
    CHECK(without.empty());
  }
}

} // namespace seissol::expr::test
