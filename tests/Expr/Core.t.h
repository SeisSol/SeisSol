// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// The core format of sderiv (Expr/Core.h): what the loader reads and rejects, with the line and the
// column, and the normal form the printer writes, which reads back into the same program.

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Core.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace seissol::unit_test {

using namespace seissol::expr;

namespace {

/// The outputs of `program` at one point, with the inputs `inputs` (by name), by the interpreter,
/// in the order of program.outputs().
std::vector<double> evaluateAt(const Program& program,
                               const std::vector<std::pair<std::string, double>>& inputs) {
  reader::scripting::DataTable table(1);
  std::vector<double> in(inputs.size());
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    in[i] = inputs[i].second;
    table.bindViewConst<double>(inputs[i].first, reader::scripting::Direction::In, &in[i]);
  }
  std::vector<double> out(program.outputs().size(), -1.0);
  for (std::size_t o = 0; o < out.size(); ++o) {
    table.bindView<double>(program.outputs()[o].name, reader::scripting::Direction::Out, &out[o]);
  }
  auto binding = Binding::bind(program, table);
  reader::datafield::GridStore store;
  BackendOptions options;
  options.quiet = true;
  const auto kernel = makeKernel(program, binding, store, options);
  kernel->precompute(table);
  kernel->run(table);
  return out;
}

bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

/// The error loading `text` gives: its line, its column and its message.
struct Failure {
  int line{0};
  int column{0};
  std::string message;
};

Failure failureOf(const std::string& text) {
  try {
    static_cast<void>(loadCore(text, "model.sderiv"));
  } catch (const CoreError& error) {
    return Failure{error.where().line, error.where().column, error.what()};
  }
  FAIL("the core text loaded: " << text);
  return {};
}

/// The material of the specification, without its check, and with its grid on one line, as the
/// printer writes it.
constexpr auto Material = R"(#!sderiv 1
meta generator = "easi 1.8.0"
meta source = "material.yaml", "sha256:9f86d08188"
in x, y, z
grid m (kind = "asagi", path = "vmodel.nc", variable = "data", components = ("rho", "mu", "lambda"), interpolation = "linear", boundary = "clamp")
def t0 = lt(z, -2000)
def t1 = m.rho(x, y, z)
def t2 = -0.5 * z
def t3 = 2600 + t2
def t4 = select(t0, t1, t3)
out rho = t4
)";

/// The derived output of the specification: the peak ground velocity.
constexpr auto PeakVelocity = R"(#!sderiv 1
in v1, v2, v3
state pgv = 0
def t0 = v1 * v1
def t1 = v2 * v2
def t2 = t0 + t1
def t3 = v3 * v3
def t4 = t2 + t3
def t5 = sqrt(t4)
def t6 = max(pgv, t5)
next pgv = t6
out pgv = t6
)";

} // namespace

TEST_SUITE("ExprCore") {

  TEST_CASE("the examples of the specification load, and print as they are written") {
    for (const char* text : {Material, PeakVelocity}) {
      CAPTURE(text);
      const Program program = loadCore(text, "");
      // the normal form is what the specification writes
      CHECK(printCore(program) == text);
    }
    const Program material = loadCore(Material, "/data/run/material.sderiv");
    REQUIRE(material.grids().size() == 1);
    // a relative path is relative to the directory of the core file
    CHECK(material.grids()[0].path == "/data/run/vmodel.nc");
    CHECK(material.grids()[0].variable == "data");
    REQUIRE(material.meta().size() == 2);
    CHECK(material.meta()[1].values[1].text == "sha256:9f86d08188");
  }

  TEST_CASE("a program prints, loads back the same, and prints the same again") {
    // from the full language: names, nesting, functions, constants, states and types
    const Program program = compileSderivModule(R"(
grid m (kind = "asagi", path = "/data/vmodel.nc", variable = "data",
        components = ("rho", "vs"), interpolation = "cubic", boundary = "fail", time = 2)
def f(p) = p * p + 1.5
def two = 2
state f64 sum = 0.0
def sum = f64(m_rho(x, y, t) * f(z)) + sum
out def u = real(sum) - pi * select(lt(x, y), m_vs(x, y, t), -0.25)
out def w = exp(-x) / two
)");
    const std::string text = printCore(program);
    CAPTURE(text);
    const Program loaded = loadCore(text, "");
    CHECK(printCore(loaded) == text);
    CHECK(text.find("grid m (kind = \"asagi\", path = \"/data/vmodel.nc\", variable = \"data\", "
                    "components = (\"rho\", \"vs\"), interpolation = \"cubic\", boundary = "
                    "\"fail\", time = 2)\n") != std::string::npos);
    CHECK(text.find("state f64 sum = 0\n") != std::string::npos);

    // and computes the same bits
    const Program plain = compileSderivModule("out def u = (x - y) * exp(-x) / (z + 0.1)\n"
                                              "out def v = select(le(x, 0.3), pow(y, 2), -z)\n");
    const Program back = loadCore(printCore(plain), "");
    const auto expected = evaluateAt(plain, {{"x", 0.25}, {"y", -1.5}, {"z", 3.0}});
    const auto actual = evaluateAt(back, {{"x", 0.25}, {"y", -1.5}, {"z", 3.0}});
    REQUIRE(expected.size() == actual.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(sameBits(expected[i], actual[i]));
    }
  }

  TEST_CASE("numbers are read with correct rounding, and print in their shortest form") {
    const auto valueOf = [](const std::string& number) {
      const Program program = loadCore("#!sderiv 1\nout u = " + number + "\n", "");
      return program.arena()[program.roots()[0]].value;
    };
    CHECK(valueOf(".5") == 0.5);
    CHECK(valueOf("2.") == 2.0);
    CHECK(valueOf("1e-3") == 1e-3);
    CHECK(valueOf("2.5E+3") == 2500.0);
    CHECK(valueOf("0x1.8p+1") == 3.0);
    CHECK(valueOf("0X1P-1") == 0.5);
    CHECK(valueOf("0.1") == 0.1);
    CHECK(valueOf("-2.5") == -2.5);
    CHECK(std::isinf(valueOf("inf")));
    CHECK(valueOf("-inf") < 0);
    CHECK(std::isnan(valueOf("nan")));
    // beyond the range: an infinity or zero, not an error
    CHECK(valueOf("1e999") == std::numeric_limits<double>::infinity());
    CHECK(valueOf("-1e999") == -std::numeric_limits<double>::infinity());
    CHECK(valueOf("1e-999") == 0.0);
    CHECK(valueOf("0x1p-1080") == 0.0);
    CHECK(valueOf("4.9e-324") == std::numeric_limits<double>::denorm_min());
    // negation of a number is a number; of anything else an operation
    const Program negation = loadCore("#!sderiv 1\nin x\nout u = -(2.5)\nout v = -x\n", "");
    CHECK(negation.arena()[negation.roots()[0]].kind == Kind::Const);
    CHECK(negation.arena()[negation.roots()[1]].kind == Kind::PW);

    const auto printed = [](double value) {
      Program program;
      program.addOutput("u", reader::scripting::DataType::F64, program.arena().konst(value));
      const std::string text = printCore(program);
      return text.substr(text.find("out u = ") + 8, text.size() - text.find("out u = ") - 9);
    };
    CHECK(printed(0.1) == "0.1");
    CHECK(printed(1e23) == "1e+23");
    CHECK(printed(-0.0) == "-0");
    CHECK(printed(std::numeric_limits<double>::infinity()) == "inf");
    CHECK(printed(-std::numeric_limits<double>::infinity()) == "-inf");
    CHECK(printed(std::numeric_limits<double>::quiet_NaN()) == "nan");
    for (const double value : {0.1, 1.0 / 3.0, 5e-324, -1e308, 123456789.0}) {
      CHECK(sameBits(valueOf(printed(value)), value));
    }
  }

  TEST_CASE("names: escaped where they are reserved, and operations are not reserved") {
    const Program program =
        loadCore("#!sderiv 1\nin `in`, `next`, max, `a b`\nout u = max(max, `in`) + `a b`\n"
                 "out `out` = `next`\n",
                 "");
    REQUIRE(program.inputs().size() == 4);
    CHECK(program.inputs()[0].name == "in");
    const std::string text = printCore(program);
    CHECK(text.find("in `a b`, `in`, max, `next`\n") != std::string::npos);
    CHECK(text.find("def t0 = max(max, `in`)\n") != std::string::npos);
    CHECK(text.find("out `out` = `next`\n") != std::string::npos);
    CHECK(printCore(loadCore(text, "")) == text);
    // the reserved words of the full language are escaped as well
    Program full;
    full.addInput("let", reader::scripting::DataType::F64);
    full.addOutput("u", reader::scripting::DataType::F64, full.arena().field("let"));
    CHECK(printCore(full).find("in `let`") != std::string::npos);
  }

  TEST_CASE("every operation of the program has a spelling, and reads back") {
    for (int f = 0; f <= static_cast<int>(Fn::Select); ++f) {
      const auto fn = static_cast<Fn>(f);
      CAPTURE(name(fn));
      Program program;
      program.addInput("x", reader::scripting::DataType::F64);
      program.addInput("y", reader::scripting::DataType::F64);
      program.addInput("z", reader::scripting::DataType::F64);
      Arena& arena = program.arena();
      const NodeId x = arena.field("x");
      const NodeId y = arena.field("y");
      const NodeId z = arena.field("z");
      const NodeId value = arity(fn) == 1   ? arena.pw(fn, x)
                           : arity(fn) == 2 ? arena.pw(fn, x, y)
                                            : arena.pw(fn, x, y, z);
      program.addOutput("u", reader::scripting::DataType::F64, value);
      const std::string text = printCore(program);
      const Program loaded = loadCore(text, "");
      const Node& node = loaded.arena()[loaded.roots()[0]];
      CHECK(node.kind == Kind::PW);
      CHECK(node.fn == fn);
    }
  }

  TEST_CASE("the types of values: conversions and states of a type of their own") {
    const std::string text = R"(#!sderiv 1
in t
state f64 sum = 0
def t0 = f64(t)
def t1 = sum + t0
def t2 = real(t1)
def t3 = f32(t2)
next sum = t1
out s = t3
out u = t2
)";
    const Program program = loadCore(text, "");
    REQUIRE(program.state().size() == 1);
    CHECK(program.state()[0].type == ValueType::F64);
    CHECK(printCore(program) == text);
  }

  TEST_CASE("an error names the file, the line and the column") {
    struct Case {
      std::string text;
      int line;
      int column;
      std::string fragment;
    };
    const std::vector<Case> cases = {
        {"in x\n", 1, 1, "starts with the header"},
        {"\xEF\xBB\xBF#!sderiv 1\n", 1, 1, "byte order mark"},
        {"#!sderiv\n", 1, 9, "names the version"},
        {"#!sderiv 2\n", 1, 10, "needs version 2"},
        {"#!sderiv 1\nin x\nout u = x + y\n", 3, 13, "`y` is not declared"},
        {"#!sderiv 1\nin x\nout u = x\nout v = u\n", 4, 9, "an output cannot be read"},
        {"#!sderiv 1\nin x, x\n", 2, 7, "declared twice"},
        {"#!sderiv 1\nin x\ndef x = 1\n", 3, 5, "declared twice"},
        {"#!sderiv 1\nstate s = 0\n", 2, 7, "has no next value"},
        {"#!sderiv 1\nin x\nnext x = 1\n", 3, 6, "is none"},
        {"#!sderiv 1\nin x\nout u = foo(x)\n", 3, 9, "no operation `foo`"},
        {"#!sderiv 1\nin x\nout u = max(x)\n", 3, 9, "takes 2 arguments, not 1"},
        {"#!sderiv 1\nin x\nout u = x\nmeta a = \"b\"\n", 4, 1, "directly after the header"},
        {"#!sderiv 1\nin x\nout u = x $ 1\n", 3, 11, "unexpected character `$`"},
        {"#!sderiv 1\nin x\nout u = 2x\n", 3, 10, "followed by `x`"},
        {"#!sderiv 1\nin x\nout u = 0x1.8\n", 3, 9, "binary exponent"},
        {"#!sderiv 1\nin x\nout u = \"x\n", 3, 9, "ends with \" on its line"},
        {"#!sderiv 1\nin def\n", 2, 4, "reserved word"},
        {"#!sderiv 1\ngrid m (kind = \"asagi\", path = \"a.nc\", variable = \"z\", "
         "components = (\"rho\"), interpolation = \"linear\")\n",
         2,
         6,
         "gives no boundary"},
        {"#!sderiv 1\ngrid m (kind = \"asagi\", colour = \"red\")\n",
         2,
         25,
         "no attribute `colour`"},
        {"#!sderiv 1\ngrid m (kind = \"netcdf\", path = \"a.nc\", variable = \"z\", components = "
         "(\"rho\"), interpolation = \"linear\", boundary = \"clamp\")\n",
         2,
         6,
         "\"asagi\" or \"scec\""},
        {"#!sderiv 1\nin x\ngrid m (kind = \"asagi\", path = \"a.nc\", variable = \"z\", "
         "components "
         "= (\"rho\"), interpolation = \"linear\", boundary = \"clamp\")\nout u = m.vp(x)\n",
         4,
         11,
         "no component `vp`"},
        {"#!sderiv 1\nin x\nout u = x\nout u = x\n", 4, 5, "given twice"},
        {"#!sderiv 1\nin x\n# ä comment\nout u = ä\n", 4, 9, "unexpected character"},
    };
    for (const auto& c : cases) {
      CAPTURE(c.text);
      const Failure failure = failureOf(c.text);
      CAPTURE(failure.message);
      CHECK(failure.line == c.line);
      CHECK(failure.column == c.column);
      CHECK(failure.message.find(c.fragment) != std::string::npos);
      CHECK(failure.message.rfind("model.sderiv:" + std::to_string(c.line) + ":", 0) == 0);
    }
  }

  TEST_CASE("deep nesting is an error, and long chains take no stack") {
    std::string deep = "#!sderiv 1\nin x\nout u = ";
    for (int i = 0; i < 5000; ++i) {
      deep += "(";
    }
    deep += "x";
    for (int i = 0; i < 5000; ++i) {
      deep += ")";
    }
    CHECK(failureOf(deep + "\n").message.find("nests deeper") != std::string::npos);
    // 256 levels load
    std::string nested = "#!sderiv 1\nin x\nout u = ";
    for (int i = 0; i < 256; ++i) {
      nested += "sqrt(";
    }
    nested += "x";
    for (int i = 0; i < 256; ++i) {
      nested += ")";
    }
    CHECK_NOTHROW(static_cast<void>(loadCore(nested + "\n", "")));

    // a chain of definitions, as deep as it is long, loads and prints
    constexpr int Chain = 100000;
    std::string chain = "#!sderiv 1\nin x\ndef t0 = x + 1\n";
    for (int i = 1; i < Chain; ++i) {
      chain += "def t" + std::to_string(i) + " = t" + std::to_string(i - 1) + " + 1\n";
    }
    chain += "out u = t" + std::to_string(Chain - 1) + "\n";
    const Program program = loadCore(chain, "");
    CHECK(printCore(program) == chain);
  }

  TEST_CASE("the printer refuses what the format cannot express") {
    Program program;
    const MatrixId matrix = program.internMatrix("basis", MatrixShape{1, 4, 1});
    const BlockId block = program.internBlock("dofs", 4);
    program.addOutput(
        "u", reader::scripting::DataType::F64, program.arena().contract(matrix, block));
    CHECK_THROWS_AS(static_cast<void>(printCore(program)), std::invalid_argument);
  }

  TEST_CASE("checks: in the order given, after the next values, before the outputs") {
    // the material of the specification with its check
    const std::string text = R"(#!sderiv 1
in x, y, z
grid m (kind = "asagi", path = "vmodel.nc", variable = "data", components = ("rho", "mu", "lambda"), interpolation = "linear", boundary = "clamp")
state s = 0
def t0 = lt(z, -2000)
def t1 = m.rho(x, y, z)
def t2 = -0.5 * z
def t3 = 2600 + t2
def t4 = select(t0, t1, t3)
def t5 = s + 1
def t6 = le(z, 0)
next s = t5
check t6, "no material above the surface"
check x, "x is \"zero\""
out rho = t4
)";
    const Program program = loadCore(text, "");
    REQUIRE(program.checks().size() == 2);
    CHECK(program.checks()[0].message == "no material above the surface");
    CHECK(program.checks()[1].message == "x is \"zero\"");
    CHECK(printCore(program) == text);
    // a check is part of what the kernel computes
    const Program without = loadCore("#!sderiv 1\nin x\nout u = x\n", "");
    const Program with = loadCore("#!sderiv 1\nin x\ncheck x, \"a\"\nout u = x\n", "");
    CHECK(without.fingerprint() != with.fingerprint());

    // the full language has them as well
    const Program full = compileSderivModule("check lt(x, 2), \"x below 2\"\nout def u = x\n");
    REQUIRE(full.checks().size() == 1);
    CHECK(printCore(full).find("def t0 = lt(x, 2)\ncheck t0, \"x below 2\"\n") !=
          std::string::npos);
  }

  TEST_CASE("a program with checks runs in the interpreter, whichever backend is asked for") {
    const Program program =
        loadCore("#!sderiv 1\nin x\ncheck le(x, 1), \"x at most 1\"\nout u = 2 * x\n", "");
    std::vector<double> x = {0.5, 1.0};
    std::vector<double> u(2, 0.0);
    reader::scripting::DataTable table(2);
    table.bindViewConst<double>("x", reader::scripting::Direction::In, x.data());
    table.bindView<double>("u", reader::scripting::Direction::Out, u.data());
    auto binding = Binding::bind(program, table);
    reader::datafield::GridStore store;
    BackendOptions options;
    options.quiet = true;
    options.preferred = BackendKind::RtcCpu;
    const auto kernel = makeKernel(program, binding, store, options);
    CHECK(kernel->kind() == BackendKind::Interpreter);
    kernel->precompute(table);
    kernel->run(table);
    // every point meets the check
    CHECK(u[0] == 1.0);
    CHECK(u[1] == 2.0);
  }

  TEST_CASE("a module in the core format reads through the sderiv entry point") {
    SderivOptions options;
    options.file = "/data/model.sderiv";
    const Program program = compileSderivModule(PeakVelocity, options);
    REQUIRE(program.outputs().size() == 1);
    CHECK(program.outputs()[0].name == "pgv");
    REQUIRE(program.state().size() == 1);
    CHECK(program.inputs().size() == 3);
  }
}

} // namespace seissol::unit_test
