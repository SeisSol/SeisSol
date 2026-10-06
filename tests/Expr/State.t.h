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
#include "Reader/Scripting/DataTable.h"
#include "TestHelper.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace df = reader::datafield;
using reader::scripting::DataTable;
using reader::scripting::Direction;

inline bool hasInput(const Program& program, const std::string& name) {
  return std::any_of(program.inputs().begin(), program.inputs().end(), [&](const VarSpec& v) {
    return v.name == name;
  });
}

TEST_SUITE("ExprState") {

  TEST_CASE("a state is declared, updated by its definition and not an input") {
    const Program program = compileSderivModule("state pgv = 0.0\n"
                                                "out def pgv = max(pgv, sqrt(v1*v1 + v2*v2))\n");
    REQUIRE(program.state().size() == 1);
    CHECK(program.state()[0].name == "pgv");
    CHECK(program.state()[0].initial == 0.0);
    REQUIRE(program.outputs().size() == 1);
    CHECK(program.outputs()[0].name == "pgv");
    // the same name is state and output: one root for both
    CHECK(program.roots()[0] == program.state()[0].root);
    CHECK_FALSE(hasInput(program, "pgv"));
    CHECK(hasInput(program, "v1"));
  }

  TEST_CASE("a state updated by a plain definition is no output, and its reads are inputs") {
    const Program program = compileSderivModule("state last = -1.0\n"
                                                "out def fresh = lt(last, window)\n"
                                                "def last = window2\n");
    REQUIRE(program.state().size() == 1);
    CHECK(program.state()[0].initial == -1.0);
    REQUIRE(program.outputs().size() == 1);
    CHECK(program.outputs()[0].name == "fresh");
    CHECK(hasInput(program, "window"));
    // read by the update only
    CHECK(hasInput(program, "window2"));
  }

  TEST_CASE("state declarations are checked") {
    const auto rejects = [](const std::string& source, const std::string& fragment) {
      CAPTURE(source);
      bool thrown = false;
      try {
        static_cast<void>(compileSderivModule(source));
      } catch (const SderivError& error) {
        thrown = true;
        const std::string message = error.what();
        CAPTURE(message);
        CHECK(message.find(fragment) != std::string::npos);
      }
      CHECK(thrown);
    };
    rejects("state a = 0.0\nout def b = a\n", "never updated");
    rejects("state a = 0.0\nstate a = 1.0\nout def a = a\n", "duplicate state");
    rejects("state a = 0.0\ndef a(x) = x\nout def b = a(1.0)\n", "cannot take parameters");
    rejects("state sqrt = 0.0\nout def sqrt = 1.0\n", "name of a constant, function");
  }

  TEST_CASE("a maximum over time is a state per point") {
    const Program program = compileSderivModule(
        "state pgv = 0.0\nout def pgv = max(pgv, sqrt(v1*v1 + v2*v2 + v3*v3))\n");
    constexpr std::size_t NumPoints = 6;
    constexpr std::size_t Calls = 5;

    for (const auto backend : {BackendKind::Interpreter, BackendKind::RtcCpu}) {
      CAPTURE(name(backend));
      std::vector<double> v1(NumPoints);
      std::vector<double> v2(NumPoints);
      std::vector<double> v3(NumPoints);
      std::vector<double> pgv(NumPoints, -1.0);
      DataTable table(NumPoints);
      table.bindViewConst<double>("v1", Direction::In, v1.data());
      table.bindViewConst<double>("v2", Direction::In, v2.data());
      table.bindViewConst<double>("v3", Direction::In, v3.data());
      table.bindView<double>("pgv", Direction::Out, pgv.data());
      Binding binding = Binding::bind(program, table);
      df::GridStore store;
      BackendOptions options;
      options.preferred = backend;
      const auto kernel = makeKernel(program, binding, store, options);
      kernel->precompute(table);

      std::vector<double> reference(NumPoints, 0.0);
      for (std::size_t call = 0; call < Calls; ++call) {
        for (std::size_t p = 0; p < NumPoints; ++p) {
          const double phase = 0.9 * static_cast<double>(call) + 0.4 * static_cast<double>(p);
          v1[p] = std::sin(phase);
          v2[p] = 0.5 * std::cos(2.0 * phase);
          v3[p] = 0.1 * static_cast<double>(call) - 0.2;
          reference[p] =
              std::fmax(reference[p], std::sqrt(v1[p] * v1[p] + v2[p] * v2[p] + v3[p] * v3[p]));
        }
        kernel->run(table);
        for (std::size_t p = 0; p < NumPoints; ++p) {
          CAPTURE(call);
          CAPTURE(p);
          CHECK(pgv[p] == doctest::Approx(reference[p]).epsilon(1e-15));
        }
      }
    }
  }

  TEST_CASE("a window maximum restarts itself through a stride-0 counter") {
    // The counter is a variable of the consumer, bound once as a stride-0 view: the next call
    // sees its new value without a rebind. The restart rule lives in the model.
    const Program program = compileSderivModule("state last = -1.0\n"
                                                "state peak = 0.0\n"
                                                "def fresh = lt(last, window)\n"
                                                "out def peak = select(fresh, v, max(peak, v))\n"
                                                "out def last = window\n");
    double window = 0.0;
    double v = 0.0;
    double peak = -1.0;
    double last = -1.0;
    DataTable table(1);
    table.bindViewConst<double>("window", Direction::In, &window, 0);
    table.bindViewConst<double>("v", Direction::In, &v, 0);
    table.bindView<double>("peak", Direction::Out, &peak);
    table.bindView<double>("last", Direction::Out, &last);
    Binding binding = Binding::bind(program, table);
    df::GridStore store;
    const auto kernel = makeKernel(program, binding, store, {});
    kernel->precompute(table);

    const auto step = [&](double value) {
      v = value;
      kernel->run(table);
      return peak;
    };
    CHECK(step(3.0) == 3.0);
    CHECK(step(7.0) == 7.0);
    CHECK(step(2.0) == 7.0);
    window = 1.0;
    CHECK(step(1.0) == 1.0); // restarted, not 7
    CHECK(step(4.0) == 4.0);
    CHECK(last == 1.0);
  }

  TEST_CASE("a state the table keeps lies where the table put it") {
    // Two states, one of them updated by a plain definition, kept per cell of a gathered subset
    // of two layers whose bases move per call -- [state][point][simulation] in a cell, as a
    // solver's storage holds them -- against the same program with the states in the Binding.
    const Program program = compileSderivModule("state peak = -1.0\n"
                                                "state total = 0.5\n"
                                                "out def peak = max(peak, v)\n"
                                                "out def ratio = total / (1.0 + abs(peak))\n"
                                                "def total = total + v * w\n");
    constexpr std::size_t Rows = 3;
    constexpr std::size_t States = 2;
    constexpr std::size_t Simulations = 2;
    constexpr std::size_t Simulation = 1;
    constexpr std::size_t CellValues = States * Rows * Simulations;
    constexpr std::size_t LayerCells = 5;
    // the point set: three cells of the first layer, then two of the second
    const std::vector<std::uint32_t> cellIndex = {4, 1, 3, 0, 2};
    constexpr std::size_t FirstLayerCells = 3;
    constexpr std::size_t NumPoints = 5 * Rows;
    constexpr double Unused = 1234.5;

    std::vector<double> v(NumPoints);
    std::vector<double> w(NumPoints);
    const auto setInputs = [&](std::size_t call) {
      for (std::size_t p = 0; p < NumPoints; ++p) {
        v[p] = std::sin(0.7 * static_cast<double>(call) + 0.3 * static_cast<double>(p));
        w[p] = 0.25 * static_cast<double>((p + call) % 4) - 0.4;
      }
    };

    for (const auto backend : {BackendKind::Interpreter, BackendKind::RtcCpu}) {
      CAPTURE(name(backend));
      // the reference: the states in the Binding
      std::vector<double> peakReference(NumPoints);
      std::vector<double> ratioReference(NumPoints);
      DataTable referenceTable(NumPoints);
      referenceTable.bindViewConst<double>("v", Direction::In, v.data());
      referenceTable.bindViewConst<double>("w", Direction::In, w.data());
      referenceTable.bindView<double>("peak", Direction::Out, peakReference.data());
      referenceTable.bindView<double>("ratio", Direction::Out, ratioReference.data());
      Binding referenceBinding = Binding::bind(program, referenceTable);
      CHECK_FALSE(referenceBinding.stateKept(0));
      df::GridStore store;
      BackendOptions interpreted;
      interpreted.preferred = BackendKind::Interpreter;
      const auto reference = makeKernel(program, referenceBinding, store, interpreted);
      reference->precompute(referenceTable);

      // two layers of cells; the states of the simulation not evaluated stay as they are
      std::vector<std::vector<double>> layers(2,
                                              std::vector<double>(LayerCells * CellValues, Unused));
      for (auto& layer : layers) {
        for (std::size_t cell = 0; cell < LayerCells; ++cell) {
          for (std::size_t s = 0; s < States; ++s) {
            for (std::size_t row = 0; row < Rows; ++row) {
              layer[cell * CellValues + (s * Rows + row) * Simulations + Simulation] =
                  program.state()[s].initial;
            }
          }
        }
      }
      const auto stateBase = [&](std::size_t layer, std::size_t s) {
        return layers[layer].data() + s * Rows * Simulations + Simulation;
      };

      std::vector<double> peak(NumPoints);
      std::vector<double> ratio(NumPoints);
      DataTable table(NumPoints);
      table.bindViewConst<double>("v", Direction::In, v.data());
      table.bindViewConst<double>("w", Direction::In, w.data());
      table.bindView<double>("peak", Direction::Out, peak.data());
      table.bindView<double>("ratio", Direction::Out, ratio.data());
      for (std::size_t s = 0; s < States; ++s) {
        table.bindState<double>(program.state()[s].name,
                                stateBase(0, s),
                                Rows,
                                CellValues,
                                Simulations,
                                cellIndex.data());
      }
      Binding binding = Binding::bind(program, table);
      CHECK(binding.stateKept(0));
      CHECK(binding.stateKept(1));
      BackendOptions options;
      options.preferred = backend;
      options.quiet = true;
      const auto kernel = makeKernel(program, binding, store, options);
      kernel->precompute(table);

      for (std::size_t call = 0; call < 4; ++call) {
        CAPTURE(call);
        setInputs(call);
        reference->run(referenceTable);
        // a call per layer, with the state bases of its layer
        for (std::size_t layer = 0; layer < 2; ++layer) {
          std::array<void*, States> bases{stateBase(layer, 0), stateBase(layer, 1)};
          KernelArgs args;
          args.states = bases.data();
          args.stateCount = bases.size();
          args.first = layer == 0 ? 0 : FirstLayerCells * Rows;
          args.count = (layer == 0 ? FirstLayerCells : 5 - FirstLayerCells) * Rows;
          kernel->run(args);
        }
        CHECK(unit_test::bitwiseEqual(peak.data(), peakReference.data(), NumPoints));
        CHECK(unit_test::bitwiseEqual(ratio.data(), ratioReference.data(), NumPoints));
      }

      // the states lie where the table put them, as the Binding keeps its own
      for (std::size_t point = 0; point < NumPoints; ++point) {
        const std::size_t element = point / Rows;
        const std::size_t layer = element < FirstLayerCells ? 0 : 1;
        const std::size_t cell = cellIndex[element];
        for (std::size_t s = 0; s < States; ++s) {
          const double kept = stateBase(layer, s)[cell * CellValues + (point % Rows) * Simulations];
          const double own = static_cast<const double*>(referenceBinding.stateBase(s))[point];
          CHECK(unit_test::bitwiseEqual(kept, own));
        }
      }
      // and nothing else was touched
      for (std::size_t layer = 0; layer < 2; ++layer) {
        for (std::size_t i = 0; i < layers[layer].size(); ++i) {
          if (i % Simulations != Simulation) {
            CHECK(layers[layer][i] == Unused);
          }
        }
        // the cells outside the point set keep their initial state
        const std::size_t untouched = layer == 0 ? 0 : 1;
        CHECK(layers[layer][untouched * CellValues + Simulation] == program.state()[0].initial);
      }
    }
  }

  TEST_CASE("a state kept by the table is checked") {
    const Program program = compileSderivModule("state peak = 0.0\nout def peak = max(peak, v)\n");
    std::vector<double> v(4);
    std::vector<double> peak(4);
    std::vector<double> storage(8);
    std::vector<float> narrow(8);
    const auto bind = [&](const std::function<void(DataTable&)>& keep) {
      DataTable table(4);
      table.bindViewConst<double>("v", Direction::In, v.data());
      table.bindView<double>("peak", Direction::Out, peak.data());
      keep(table);
      return Binding::bind(program, table);
    };
    CHECK_NOTHROW(
        bind([&](DataTable& table) { table.bindState<double>("peak", storage.data(), 2, 2, 1); }));
    // in another type than the program computes in
    CHECK_THROWS_AS(
        bind([&](DataTable& table) { table.bindState<float>("peak", narrow.data(), 2, 2, 1); }),
        std::invalid_argument);
    // per cell of points that do not divide the point set
    CHECK_THROWS_AS(
        bind([&](DataTable& table) { table.bindState<double>("peak", storage.data(), 3, 3, 1); }),
        std::invalid_argument);
    // twice
    CHECK_THROWS_AS(bind([&](DataTable& table) {
                      table.bindState<double>("peak", storage.data(), 2, 2, 1);
                      table.bindState<double>("peak", storage.data(), 2, 2, 1);
                    }),
                    std::invalid_argument);
    // a state the program does not declare is not an error, like an extra column
    CHECK_NOTHROW(
        bind([&](DataTable& table) { table.bindState<double>("other", storage.data(), 2, 2, 1); }));
  }

} // TEST_SUITE

} // namespace seissol::expr::test
