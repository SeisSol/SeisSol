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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace {

namespace df = reader::datafield;
using reader::scripting::DataTable;
using reader::scripting::Direction;

bool hasInput(const Program& program, const std::string& name) {
  return std::any_of(program.inputs().begin(), program.inputs().end(), [&](const VarSpec& v) {
    return v.name == name;
  });
}

} // namespace

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

} // TEST_SUITE

} // namespace seissol::expr::test
