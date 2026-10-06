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

namespace {

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
double value(std::size_t i, double scale) {
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

Program contracted() {
  Program program = compileSderivModule(Model);
  const auto matrix = program.internMatrix("proj", MatrixShape{Rows, Cols, Ld});
  const auto v = program.internBlock("v", Cols);
  const auto w = program.internBlock("w", Cols);
  substituteByContraction(program, matrix, {{"v", v}, {"w", w}});
  return program;
}

std::vector<double> evaluate(const Program& program,
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

bool sameBits(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

} // namespace

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

} // TEST_SUITE

} // namespace seissol::expr::test
