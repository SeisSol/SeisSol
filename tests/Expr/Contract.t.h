// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Cost.h"
#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/Rewrite.h"
#include "Expr/RtcCpu.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace df = reader::datafield;
using reader::scripting::DataTable;
using reader::scripting::DataType;
using reader::scripting::Direction;

// Points per cell, coefficients per cell, and the padded leading dimension of the matrix.
constexpr std::size_t Rows = 3;
constexpr std::size_t Cols = 4;
constexpr std::size_t Ld = 5;
// The cells of the point set, and the cells in memory they are gathered from.
constexpr std::size_t NumCells = 5;
constexpr std::size_t MemoryCells = 7;
constexpr std::size_t NumPoints = NumCells * Rows;
// Coefficients are stored per cell as [quantity][padded coefficient], two quantities.
constexpr std::size_t PaddedCols = 6;
constexpr std::size_t CellStride = 2 * PaddedCols;

/// Deterministic, not nice values: a contraction that silently reads the wrong coefficient or
/// sums in another order has to show up.
inline double value(std::size_t i, double scale) {
  return scale * std::sin(1.3 * static_cast<double>(i) + 0.7) + 1e-3 * static_cast<double>(i % 7);
}

struct Fixture {
  std::vector<double> matrix = std::vector<double>(Ld * Cols);
  std::vector<double> otherMatrix = std::vector<double>(Ld * Cols);
  std::vector<double> dofs = std::vector<double>(MemoryCells * CellStride);
  // a gathered subset of the memory cells, out of order
  std::vector<std::uint32_t> cellIndex = {6, 2, 0, 5, 3};
  std::vector<double> t = std::vector<double>(NumPoints);

  Fixture() {
    for (std::size_t i = 0; i < matrix.size(); ++i) {
      matrix[i] = value(i, 1.0);
      otherMatrix[i] = value(i + 100, 2.0);
    }
    for (std::size_t i = 0; i < dofs.size(); ++i) {
      dofs[i] = value(i + 1000, 10.0);
    }
    for (std::size_t i = 0; i < t.size(); ++i) {
      t[i] = value(i + 2000, 3.0);
    }
  }

  void bindBlocks(DataTable& table) const {
    table.bindBlock<double>("v", dofs.data(), Cols, CellStride, 1, cellIndex.data());
    table.bindBlock<double>("w", dofs.data() + PaddedCols, Cols, CellStride, 1, cellIndex.data());
    table.bindMatrix<double>("proj", matrix.data(), Rows, Cols, Ld);
  }

  /// The reference value of the contraction of block `quantity` (0: v, 1: w) at every point.
  [[nodiscard]] std::vector<double> projected(std::size_t quantity, const double* m) const {
    ContractOperands operands;
    operands.matrix = m;
    operands.matrixType = DataType::F64;
    operands.rows = Rows;
    operands.cols = Cols;
    operands.leadingDimension = Ld;
    operands.block = dofs.data() + quantity * PaddedCols;
    operands.blockType = DataType::F64;
    operands.cellStride = CellStride * sizeof(double);
    operands.modeStride = sizeof(double);
    operands.cellIndex = cellIndex.data();
    std::vector<double> out(NumPoints);
    contractLanes<double>(operands, nullptr, 0, NumPoints, out.data());
    return out;
  }
};

/// The pointwise model: it reads v and w as ordinary channels.
const char* const Model = "def a = 2.0 * v - w * t\n"
                          "out def u = a * a + sqrt(abs(v))\n"
                          "out def r = select(lt(v, w), v, w) + t\n";

inline Program contracted() {
  Program program = compileSderivModule(Model);
  const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
  const auto v = program.internBlock("v", Cols);
  const auto w = program.internBlock("w", Cols);
  substituteByContraction(program, matrix, {{"v", v}, {"w", w}});
  return program;
}

inline std::vector<double> evaluate(const Program& program,
                                    const Fixture& fixture,
                                    BackendKind backend,
                                    BackendKind* used = nullptr) {
  const std::size_t outputs = program.outputs().size();
  std::vector<double> out(outputs * NumPoints, -1.0);
  DataTable table(NumPoints);
  table.bindViewConst<double>("t", Direction::In, fixture.t.data());
  fixture.bindBlocks(table);
  for (std::size_t i = 0; i < outputs; ++i) {
    table.bindView<double>(program.outputs()[i].name, Direction::Out, &out[i * NumPoints]);
  }
  Binding binding = Binding::bind(program, table);
  df::GridStore store;
  BackendOptions options;
  options.preferred = backend;
  const auto kernel = makeKernel(program, binding, store, options);
  kernel->precompute(table);
  kernel->run(table);
  if (used != nullptr) {
    *used = kernel->kind();
  }
  return out;
}

inline bool sameBits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

TEST_SUITE("ExprContract") {

  TEST_CASE("a contraction is a leaf, interned by matrix and block") {
    Program program;
    const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
    const auto v = program.internBlock("v", Cols);
    const auto w = program.internBlock("w", Cols);
    Arena& arena = program.arena();
    const NodeId a = arena.contract(matrix, v);
    CHECK(arena.contract(matrix, v) == a);
    CHECK(arena.contract(matrix, w) != a);
    CHECK(arena.children(a).empty());
    CHECK(program.pointsPerCell() == Rows);
  }

  TEST_CASE("matrices and blocks are declared by name, with one form each") {
    Program program;
    const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
    CHECK(program.internMatrix("proj", MatrixShape{Rows, Cols, Ld}) == matrix);
    CHECK_THROWS_AS(program.internMatrix("proj", MatrixShape{Rows, Cols, Ld + 1}),
                    std::invalid_argument);
    // one point set, one cell structure
    CHECK_THROWS_AS(program.internMatrix("other", MatrixShape{Rows + 1, Cols, Ld + 1}),
                    std::invalid_argument);
    CHECK_THROWS_AS(program.internMatrix("narrow", MatrixShape{Rows, Cols, Rows - 1}),
                    std::invalid_argument);
    const auto v = program.internBlock("v", Cols);
    CHECK(program.internBlock("v", Cols) == v);
    CHECK_THROWS_AS(program.internBlock("v", Cols + 1), std::invalid_argument);
  }

  TEST_CASE("validation checks the block against the columns of the matrix") {
    Program program = compileSderivModule("out def u = v\n");
    const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
    const auto shortBlock = program.internBlock("v", Cols - 1);
    CHECK_THROWS_AS(substituteByContraction(program, matrix, {{"v", shortBlock}}),
                    std::invalid_argument);
  }

  TEST_CASE("the rewrite replaces the channel, drops the input and changes the identity") {
    const Program plain = compileSderivModule(Model);
    const Program program = contracted();

    std::vector<std::string> inputs;
    for (const auto& input : program.inputs()) {
      inputs.push_back(input.name);
    }
    CHECK(inputs == std::vector<std::string>{"t"});
    CHECK(program.fingerprint() != plain.fingerprint());

    std::size_t contractions = 0;
    for (std::size_t i = 0; i < program.arena().size(); ++i) {
      contractions += program.arena()[static_cast<NodeId>(i)].kind == Kind::Contract ? 1 : 0;
    }
    // v is read three times and w twice, but each is one node
    CHECK(contractions == 2);

    // a name the program does not read is not an error
    Program again = program;
    const auto matrix = again.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
    const auto unused = again.internBlock("unused", Cols);
    CHECK_NOTHROW(substituteByContraction(again, matrix, {{"nothere", unused}}));
    CHECK(again.fingerprint() != program.fingerprint()); // one more declared block
  }

  TEST_CASE("the interpreter contracts as the reference does, through a gathered cell index") {
    const Fixture fixture;
    Program program = compileSderivModule("out def u = v\nout def r = w\n");
    const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
    const auto v = program.internBlock("v", Cols);
    const auto w = program.internBlock("w", Cols);
    substituteByContraction(program, matrix, {{"v", v}, {"w", w}});

    const auto out = evaluate(program, fixture, BackendKind::Interpreter);
    const auto expectedV = fixture.projected(0, fixture.matrix.data());
    const auto expectedW = fixture.projected(1, fixture.matrix.data());
    for (std::size_t p = 0; p < NumPoints; ++p) {
      CAPTURE(p);
      CHECK(out[p] == expectedV[p]);
      CHECK(out[NumPoints + p] == expectedW[p]);
    }

    // and the reference is the plain dot product of the gathered cell
    for (std::size_t p = 0; p < NumPoints; ++p) {
      const std::size_t cell = fixture.cellIndex[p / Rows];
      double sum = 0.0;
      for (std::size_t m = 0; m < Cols; ++m) {
        sum += fixture.matrix[p % Rows + m * Ld] * fixture.dofs[cell * CellStride + m];
      }
      CHECK(expectedV[p] == doctest::Approx(sum).epsilon(1e-14));
    }
  }

  TEST_CASE("a program against projected columns and against coefficients agree bit for bit") {
    const Fixture fixture;
    const Program plain = compileSderivModule(Model);
    const Program program = contracted();

    // the projected columns, computed by the reference contraction
    const auto v = fixture.projected(0, fixture.matrix.data());
    const auto w = fixture.projected(1, fixture.matrix.data());
    std::vector<double> viaColumns(plain.outputs().size() * NumPoints, -1.0);
    {
      DataTable table(NumPoints);
      table.bindViewConst<double>("t", Direction::In, fixture.t.data());
      table.bindViewConst<double>("v", Direction::In, v.data());
      table.bindViewConst<double>("w", Direction::In, w.data());
      for (std::size_t i = 0; i < plain.outputs().size(); ++i) {
        table.bindView<double>(plain.outputs()[i].name, Direction::Out, &viaColumns[i * NumPoints]);
      }
      Binding binding = Binding::bind(plain, table);
      df::GridStore store;
      const auto kernel = makeKernel(plain, binding, store, {});
      kernel->precompute(table);
      kernel->run(table);
    }

    const auto viaBlocks = evaluate(program, fixture, BackendKind::Interpreter);
    CHECK(sameBits(viaColumns, viaBlocks));
  }

  TEST_CASE("the compiled CPU kernel contracts bit for bit as the interpreter does") {
    const Fixture fixture;
    const Program program = contracted();
    const auto interpreted = evaluate(program, fixture, BackendKind::Interpreter);
    BackendKind used = BackendKind::Interpreter;
    const auto compiled = evaluate(program, fixture, BackendKind::RtcCpu, &used);
    if (used != BackendKind::RtcCpu) {
      WARN_MESSAGE(false, "no usable C++ compiler; the compiled CPU kernel was not executed");
      return;
    }
    CHECK(sameBits(interpreted, compiled));
  }

  TEST_CASE("only the matrix base moves between the subcells") {
    const Fixture fixture;
    const Program program = contracted();

    for (const auto backend : {BackendKind::Interpreter, BackendKind::RtcCpu}) {
      CAPTURE(name(backend));
      const std::size_t outputs = program.outputs().size();
      std::vector<double> out(outputs * NumPoints, -1.0);
      DataTable table(NumPoints);
      table.bindViewConst<double>("t", Direction::In, fixture.t.data());
      fixture.bindBlocks(table);
      for (std::size_t i = 0; i < outputs; ++i) {
        table.bindView<double>(program.outputs()[i].name, Direction::Out, &out[i * NumPoints]);
      }
      Binding binding = Binding::bind(program, table);
      df::GridStore store;
      BackendOptions options;
      options.preferred = backend;
      const auto kernel = makeKernel(program, binding, store, options);
      kernel->precompute(table);

      // a second subcell: the same program and binding, another matrix
      std::vector<double> moved(outputs * NumPoints, -1.0);
      std::vector<void*> bases(outputs);
      for (std::size_t i = 0; i < outputs; ++i) {
        bases[i] = &moved[i * NumPoints];
      }
      const void* matrices[] = {fixture.otherMatrix.data()};
      KernelArgs args{};
      args.outputs = bases.data();
      args.outputCount = outputs;
      args.matrices = matrices;
      args.matrixCount = 1;
      args.first = 0;
      args.count = NumPoints;
      kernel->run(args);

      // against a binding to that matrix outright
      Fixture other = fixture;
      other.matrix = fixture.otherMatrix;
      const auto expected = evaluate(program, other, BackendKind::Interpreter);
      CHECK(sameBits(moved, expected));
    }
  }

  TEST_CASE("a reordered point set contracts at the points it was given") {
    const Fixture fixture;
    Program program = compileSderivModule("out def u = v + group\n");
    const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
    const auto v = program.internBlock("v", Cols);
    substituteByContraction(program, matrix, {{"v", v}});

    // groups that force a permutation
    std::vector<std::int32_t> group(NumPoints);
    for (std::size_t p = 0; p < NumPoints; ++p) {
      group[p] = static_cast<std::int32_t>((p * 7) % 3);
    }
    for (const auto backend : {BackendKind::Interpreter, BackendKind::RtcCpu}) {
      CAPTURE(name(backend));
      std::vector<double> out(NumPoints, -1.0);
      DataTable table(NumPoints);
      table.bindViewConst<std::int32_t>("group", Direction::In, group.data());
      fixture.bindBlocks(table);
      table.bindView<double>("u", Direction::Out, out.data());
      Binding binding = Binding::bind(program, table);
      REQUIRE_FALSE(binding.permutation().empty());
      df::GridStore store;
      BackendOptions options;
      options.preferred = backend;
      const auto kernel = makeKernel(program, binding, store, options);
      kernel->precompute(table);
      kernel->run(table);

      const auto expected = fixture.projected(0, fixture.matrix.data());
      for (std::size_t p = 0; p < NumPoints; ++p) {
        CAPTURE(p);
        CHECK(out[p] == expected[p] + static_cast<double>(group[p]));
      }
    }
  }

  TEST_CASE("binding checks the bound blocks and matrices") {
    const Fixture fixture;
    const Program program = contracted();
    std::vector<double> u(NumPoints);
    std::vector<double> r(NumPoints);

    const auto table = [&](std::size_t points, std::size_t length, std::size_t ld) {
      DataTable result(points);
      result.bindViewConst<double>("t", Direction::In, fixture.t.data());
      result.bindView<double>("u", Direction::Out, u.data());
      result.bindView<double>("r", Direction::Out, r.data());
      result.bindBlock<double>("v", fixture.dofs.data(), length, CellStride);
      result.bindBlock<double>("w", fixture.dofs.data() + PaddedCols, length, CellStride);
      result.bindMatrix<double>("proj", fixture.matrix.data(), Rows, Cols, ld);
      return result;
    };

    CHECK_NOTHROW(Binding::bind(program, table(NumPoints, Cols, Ld)));
    // more coefficients than read are fine, fewer are not
    CHECK_NOTHROW(Binding::bind(program, table(NumPoints, PaddedCols, Ld)));
    CHECK_THROWS_AS(Binding::bind(program, table(NumPoints, Cols - 1, Ld)), std::invalid_argument);
    // the matrix has to have the form the program was built for
    CHECK_THROWS_AS(Binding::bind(program, table(NumPoints, Cols, Ld + 1)), std::invalid_argument);
    // and the points have to form whole cells
    CHECK_THROWS_AS(Binding::bind(program, table(NumPoints - 1, Cols, Ld)), std::invalid_argument);

    DataTable missing(NumPoints);
    missing.bindViewConst<double>("t", Direction::In, fixture.t.data());
    missing.bindView<double>("u", Direction::Out, u.data());
    missing.bindView<double>("r", Direction::Out, r.data());
    CHECK_THROWS_AS(Binding::bind(program, missing), std::invalid_argument);
  }

  TEST_CASE("a channel can be replaced by an expression over new channels") {
    // the physical derivative of v by the chain rule over reference derivatives and the rows of
    // the inverse Jacobian -- the form the derived output offers dx_<quantity> in
    Program program = compileSderivModule("out def u = dx_v * 2.0 + w\n");
    std::map<std::string, ChannelBuilder> builders;
    builders.emplace("dx_v", [](Arena& arena) {
      NodeId sum = NoNode;
      for (int k = 0; k < 3; ++k) {
        // one statement per channel: the order of evaluation of call arguments is unspecified
        const NodeId jinv = arena.field("jinv" + std::to_string(k) + "0");
        const NodeId derivative = arena.field("v_r" + std::to_string(k));
        const NodeId term = arena.pw(Fn::Mul, jinv, derivative);
        sum = sum == NoNode ? term : arena.pw(Fn::Add, sum, term);
      }
      return sum;
    });
    substituteChannels(program, builders);

    std::vector<std::string> inputs;
    for (const auto& input : program.inputs()) {
      inputs.push_back(input.name);
    }
    // w stays, dx_v goes, the new channels come in creation order
    CHECK(inputs ==
          std::vector<std::string>{"w", "jinv00", "v_r0", "jinv10", "v_r1", "jinv20", "v_r2"});

    // and evaluates as the expression it was replaced by
    std::vector<double> values(inputs.size());
    DataTable table(1);
    for (std::size_t i = 0; i < inputs.size(); ++i) {
      values[i] = 0.5 + static_cast<double>(i);
      table.bindViewConst<double>(inputs[i], Direction::In, &values[i]);
    }
    double u = -1.0;
    table.bindView<double>("u", Direction::Out, &u);
    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);
    kernel->run(table);
    const double dxv = values[1] * values[2] + values[3] * values[4] + values[5] * values[6];
    CHECK(u == doctest::Approx(dxv * 2.0 + values[0]));
  }

  TEST_CASE("a contraction costs one multiply-add and two loads per coefficient") {
    Program program = compileSderivModule("out def u = v\n");
    const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
    const auto v = program.internBlock("v", Cols);
    substituteByContraction(program, matrix, {{"v", v}});
    const auto costs = cost(program, lower(program), program.computeType());
    CHECK(costs.run.fma == Cols);
    CHECK(costs.run.loads == 2 * Cols);
    CHECK(costs.run.operations() == Cols);
  }

  TEST_CASE("the emitted CPU source runs the reference loop") {
    const Program program = contracted();
    const std::string source = emitCpuSource(program, lower(program));
    CHECK(source.find("for (unsigned long c_m = 0; c_m < 4; ++c_m)") != std::string::npos);
    CHECK(source.find("c_matrix[c_m * 5]") != std::string::npos);
    CHECK(source.find("pointIndex != 0 ? pointIndex[l] : first + l") != std::string::npos);
  }

  TEST_CASE("two programs merge into one, which computes what they share once") {
    const Program first = compileSderivModule("out def a = v * 2.0 + x\n");
    const Program second = compileSderivModule(
        "state s = 0.0\nout def s = s + v * 2.0\nout def b = sin(v * 2.0) + y\n");
    const Program merged = mergePrograms(first, second);

    std::vector<std::string> inputs;
    for (const auto& input : merged.inputs()) {
      inputs.push_back(input.name);
    }
    CHECK(inputs == std::vector<std::string>{"v", "x", "y"});
    REQUIRE(merged.outputs().size() == 3);
    CHECK(merged.outputs()[0].name == "a");
    CHECK(merged.outputs()[1].name == "s");
    CHECK(merged.outputs()[2].name == "b");
    REQUIRE(merged.state().size() == 1);
    CHECK(merged.state()[0].name == "s");
    // v * 2.0 is one node, which both read
    const auto& arena = merged.arena();
    const NodeId shared = arena[merged.roots()[0]].a;
    CHECK(arena[shared].fn == Fn::Mul);
    CHECK(arena[arena[merged.roots()[2]].a].a == shared);

    // the merged program computes what both compute, also the state over calls
    constexpr std::size_t Points = 4;
    const std::vector<double> v = {0.5, -1.0, 2.0, 3.25};
    const std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
    const std::vector<double> y = {-1.0, 0.0, 1.0, 2.0};
    const auto run = [&](const Program& program, std::vector<std::vector<double>>& values) {
      DataTable table(Points);
      table.bindViewConst<double>("v", Direction::In, v.data());
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindViewConst<double>("y", Direction::In, y.data());
      values.assign(program.outputs().size(), std::vector<double>(Points, -1.0));
      for (std::size_t i = 0; i < program.outputs().size(); ++i) {
        table.bindView<double>(program.outputs()[i].name, Direction::Out, values[i].data());
      }
      Binding binding = Binding::bind(program, table);
      df::GridStore store;
      const auto kernel = makeKernel(program, binding, store, {});
      kernel->precompute(table);
      kernel->run(table);
      kernel->run(table);
    };
    std::vector<std::vector<double>> ofFirst;
    std::vector<std::vector<double>> ofSecond;
    std::vector<std::vector<double>> ofMerged;
    run(first, ofFirst);
    run(second, ofSecond);
    run(merged, ofMerged);
    CHECK(ofMerged[0] == ofFirst[0]);
    CHECK(ofMerged[1] == ofSecond[0]);
    CHECK(ofMerged[2] == ofSecond[1]);
    // the state after two calls
    CHECK(ofMerged[1][1] == -4.0);
  }

  TEST_CASE("programs that both have a name do not merge") {
    const Program first = compileSderivModule("out def a = v + 1.0\n");
    CHECK_THROWS_AS(static_cast<void>(mergePrograms(first, first)), std::invalid_argument);
    // a state of one is an input of the other
    const Program reader = compileSderivModule("out def c = s + 1.0\n");
    const Program state = compileSderivModule("state s = 0.0\nout def s = s + v\n");
    CHECK_THROWS_AS(static_cast<void>(mergePrograms(reader, state)), std::invalid_argument);
    CHECK_THROWS_AS(static_cast<void>(mergePrograms(state, reader)), std::invalid_argument);
    CHECK_THROWS_AS(static_cast<void>(mergePrograms(state, state)), std::invalid_argument);
  }

  TEST_CASE("a program restricted to what is there computes the rest as before") {
    const Program program = compileSderivModule(
        "state acc = 0.0\nout def acc = acc + Qp * 2.0\nout def a = v * 2.0 + x\n"
        "out def b = Qp + v\nout def pprime = -v\nout def c = acc + 1.0\n");
    Restriction restriction;
    const Program restricted = restrictProgram(
        program,
        [](const std::string& name) { return name != "Qp"; },
        [](const std::string& name) { return name != "pprime"; },
        &restriction);

    // a reads only what is there; b reads Qp directly, acc and c through the state
    REQUIRE(restricted.outputs().size() == 1);
    CHECK(restricted.outputs()[0].name == "a");
    CHECK(restricted.state().empty());
    std::vector<std::string> inputs;
    for (const auto& input : restricted.inputs()) {
      inputs.push_back(input.name);
    }
    CHECK(inputs == std::vector<std::string>{"v", "x"});
    CHECK(restriction.unwanted == std::vector<std::string>{"pprime"});
    REQUIRE(restriction.unavailable.size() == 3);
    CHECK(restriction.unavailable[0] == std::pair<std::string, std::string>{"acc", "Qp"});
    CHECK(restriction.unavailable[1] == std::pair<std::string, std::string>{"b", "Qp"});
    CHECK(restriction.unavailable[2] == std::pair<std::string, std::string>{"c", "Qp"});

    // with everything there, the state stays, and the outputs compute as before
    const Program full = restrictProgram(
        program,
        [](const std::string&) { return true; },
        [](const std::string& name) { return name != "a"; });
    REQUIRE(full.outputs().size() == 4);
    REQUIRE(full.state().size() == 1);
    constexpr std::size_t Points = 3;
    const std::vector<double> v = {0.5, -1.0, 2.0};
    const std::vector<double> x = {1.0, 2.0, 3.0};
    const std::vector<double> qp = {10.0, 20.0, 30.0};
    const auto run = [&](const Program& p, std::vector<std::vector<double>>& values) {
      DataTable table(Points);
      table.bindViewConst<double>("v", Direction::In, v.data());
      table.bindViewConst<double>("x", Direction::In, x.data());
      table.bindViewConst<double>("Qp", Direction::In, qp.data());
      values.assign(p.outputs().size(), std::vector<double>(Points, -1.0));
      for (std::size_t i = 0; i < p.outputs().size(); ++i) {
        table.bindView<double>(p.outputs()[i].name, Direction::Out, values[i].data());
      }
      Binding binding = Binding::bind(p, table);
      df::GridStore store;
      const auto kernel = makeKernel(p, binding, store, {});
      kernel->precompute(table);
      kernel->run(table);
      kernel->run(table);
    };
    std::vector<std::vector<double>> ofProgram;
    std::vector<std::vector<double>> ofRestricted;
    std::vector<std::vector<double>> ofFull;
    run(program, ofProgram);
    run(restricted, ofRestricted);
    run(full, ofFull);
    CHECK(ofRestricted[0] == ofProgram[1]);
    CHECK(ofFull[0] == ofProgram[0]);
    CHECK(ofFull[1] == ofProgram[2]);
    CHECK(ofFull[2] == ofProgram[3]);
    CHECK(ofFull[3] == ofProgram[4]);
  }

} // TEST_SUITE

} // namespace seissol::expr::test
