// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// Sources inline, positive / negative / differential.
//
// The negative half is short, and that is a property of this frontend rather
// than an omission: sderiv has no control flow, so there is no `if`, no
// truthiness, no iteration order and no side effects to translate (compare
// LuaFrontend.h). What remains to test is a grammar and a set of name checks.

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace seissol::unit_test {

using namespace seissol::expr;

namespace {

constexpr auto GridAndDef = R"(
# density and Lame parameters from a velocity model
grid m (kind = "asagi", path = "model.nc", components = ("rho", "vp", "vs"),
        interpolation = "linear")

def mu(px, py, pz) = m.rho(px, py, pz) * m.vs(px, py, pz) ** 2

mu(x, y, z)
)";

constexpr auto TwoGridsInterleaved = R"(
grid vel (kind = "asagi", path = "vel.nc", components = ("vp", "vs"), interpolation = "linear")
def half(a) = a / 2
grid fault (kind = "scec", path = "fault.txt", components = ("mu_s"), interpolation = "nearest")
half(vel.vp(x, y, z)) + fault.mu_s(x, y)
)";

constexpr auto SelectOverLookup = R"(
grid m (kind = "asagi", path = "model.nc", components = ("rho"), interpolation = "linear")
select(lt(z, -1000), m.rho(x, y, z), 2700)
)";

constexpr auto EscapedPath =
    R"(grid m (kind = "asagi", path = "C:\\models\\a\"b.nc", components = ("rho"), interpolation = "linear")
m.rho(x, y, z))";

constexpr auto SharedSubexpression = R"(
let r = sqrt(x*x + y*y) in r + r * 2
)";

constexpr auto StringInExpression =
    R"(grid m (kind = "asagi", path = "model.nc", components = ("rho"), interpolation = "linear")
m.rho(x, y, z) + "oops")";

constexpr auto UnknownKind =
    R"(grid m (kind = "hdf5", path = "m.nc", components = ("rho"), interpolation = "linear")
m.rho(x, y, z))";

constexpr auto UnknownInterpolation =
    R"(grid m (kind = "asagi", path = "model.nc", components = ("rho"), interpolation = "model.nc")
m.rho(x, y, z))";

constexpr auto NoComponents =
    R"(grid m (kind = "asagi", path = "model.nc", interpolation = "linear")
x)";

constexpr auto DuplicateComponent =
    R"(grid m (kind = "asagi", path = "model.nc", components = ("rho", "rho"), interpolation = "linear")
m.rho(x, y, z))";

constexpr auto GridNamedLikeDefinition =
    R"(grid m (kind = "asagi", path = "model.nc", components = ("rho"), interpolation = "linear")
def m(a) = a
m(x))";

constexpr auto NoSuchComponent =
    R"(grid m (kind = "asagi", path = "model.nc", components = ("rho"), interpolation = "linear")
m.vp(x, y, z))";

constexpr auto NoSuchGrid = R"(n.rho(x, y, z))";

constexpr auto GridAsValue =
    R"(grid m (kind = "asagi", path = "model.nc", components = ("rho"), interpolation = "linear")
m + 1)";

constexpr auto ComponentWithoutCoordinates =
    R"(grid m (kind = "asagi", path = "model.nc", components = ("rho"), interpolation = "linear")
m.rho + 1)";

constexpr auto UnterminatedString = R"(grid m (kind = "asagi", path = "model.nc)
x)";

constexpr auto BareNameForString =
    R"(grid m (kind = asagi, path = "model.nc", components = ("rho"), interpolation = "linear")
m.rho(x, y, z))";

constexpr auto WrongCoordinateCount =
    R"(grid m (kind = "asagi", path = "model.nc", components = ("rho"), interpolation = "linear")
m.rho(x, y, z, t, x, y, z))";

constexpr auto PositionalGrid = R"(grid m = "asagi", "model.nc", "linear", "rho"
m.rho(x, y, z))";

Program mustCompile(const std::string& source) {
  return compileSderiv(source, /*outputName=*/"out");
}

void mustReject(const std::string& source, const std::string& fragment) {
  // The message is part of the contract here: these are user-facing syntax
  // errors in a parameter file, and "parse error" alone sends someone to the
  // wrong line. (Checked by hand: doctest::Contains needs a newer doctest than
  // the bundled one.)
  bool thrown = false;
  try {
    static_cast<void>(compileSderiv(source, "out"));
  } catch (const SderivError& error) {
    thrown = true;
    const std::string message = error.what();
    CAPTURE(message);
    CAPTURE(fragment);
    CHECK(message.find(fragment) != std::string::npos);
  }
  REQUIRE(thrown);
}

} // namespace

TEST_SUITE("SderivFrontend") {

  TEST_CASE("a grid declaration registers one callable per component") {
    const auto program = mustCompile(GridAndDef);
    REQUIRE(program.grids().size() == 1);
    int lookups = 0;
    for (std::size_t i = 0; i < program.arena().size(); ++i) {
      if (program.arena()[static_cast<NodeId>(i)].kind == Kind::Lookup) {
        ++lookups;
      }
    }
    // m.rho and m.vs are read; m.vp is declared but unused and must not appear.
    CHECK(lookups == 2);
  }

  TEST_CASE("the full language writes where its definitions, states and outputs come from") {
    SderivOptions options;
    options.file = "model.sderiv";
    const auto program = compileSderivModule("def a = v * v\n"
                                             "state peak = 0.0\n"
                                             "out def peak = max(peak, a)\n"
                                             "check le(0, v), \"v is not negative\"\n"
                                             "out def b = sqrt(a) + 1\n",
                                             options);
    const std::string printed = printCore(program);
    CAPTURE(printed);
    // the definition of a, the state, the outputs and the condition of the check, each where its
    // construct starts: at `def`, `state`, the `def` of `out def`, and the condition
    CHECK(printed.find("meta line = \"t0\", \"model.sderiv\", 1, 1\n") != std::string::npos);
    CHECK(printed.find("meta line = \"peak\", \"model.sderiv\", 2, 1\n") != std::string::npos);
    CHECK(printed.find("meta line = \"peak\", \"model.sderiv\", 3, 5\n") != std::string::npos);
    CHECK(printed.find("meta line = \"b\", \"model.sderiv\", 5, 5\n") != std::string::npos);
    CHECK(printed.find("\"model.sderiv\", 4, 7\n") != std::string::npos);
    const auto b = std::find_if(program.outputs().begin(),
                                program.outputs().end(),
                                [](const auto& output) { return output.name == "b"; });
    REQUIRE(b != program.outputs().end());
    const auto positions =
        program.sourcesOfOutput(static_cast<std::size_t>(b - program.outputs().begin()));
    REQUIRE(!positions.empty());
    CHECK(positions.front().text() == "model.sderiv:5:5");
    // the normal form keeps them
    CHECK(printCore(loadCore(printed, "")) == printed);
  }

  TEST_CASE("a grid of layers gives its knots and values, and a file its axes") {
    const auto program = mustCompile(
        "grid layers (kind = \"inline\", components = (\"rho\", \"mu\"),\n"
        "             knots = ((-3000, -2000, -1000)),\n"
        "             values = ((3000, 40), (2700, 20), (2600, 10)),\n"
        "             interpolation = \"next\", boundary = \"clamp\")\n"
        "grid m (kind = \"asagi\", path = \"model.nc\", variable = \"data\", components = "
        "(\"vs\"),\n"
        "        axes = (\"uniform\", \"knots\"), interpolation = (\"cubic\", \"previous\"),\n"
        "        boundary = (\"clamp\", \"fail\"))\n"
        "layers.rho(z) + layers.mu(z) + m.vs(x, z)");
    REQUIRE(program.grids().size() == 2);
    const auto& layers = program.grids()[0];
    CHECK(layers.kind == reader::datafield::GridKind::Inline);
    REQUIRE(layers.data != nullptr);
    REQUIRE(layers.data->knots.size() == 1);
    CHECK(layers.data->knots[0] == std::vector<double>{-3000, -2000, -1000});
    CHECK(layers.data->values == std::vector<double>{3000, 40, 2700, 20, 2600, 10});
    CHECK(layers.interpolationOf(0) == reader::datafield::Interpolation::Next);
    CHECK(layers.axisOf(0) == reader::datafield::AxisKind::Knots);
    const auto& m = program.grids()[1];
    CHECK(m.axisOf(0) == reader::datafield::AxisKind::Uniform);
    CHECK(m.axisOf(1) == reader::datafield::AxisKind::Knots);
    CHECK(m.interpolationOf(0) == reader::datafield::Interpolation::Cubic);
    CHECK(m.interpolationOf(1) == reader::datafield::Interpolation::Previous);
    CHECK(m.boundaryOf(0) == reader::datafield::BoundaryMode::Clamp);
    CHECK(m.boundaryOf(1) == reader::datafield::BoundaryMode::Fail);

    // an axis of knots is not cubic; a file reads its knots
    mustReject("grid m (kind = \"asagi\", path = \"model.nc\", variable = \"data\", components = "
               "(\"vs\"), axes = (\"knots\"), interpolation = \"cubic\", boundary = \"clamp\")\n"
               "m.vs(z)",
               "not cubic");
    mustReject("grid m (kind = \"asagi\", path = \"model.nc\", variable = \"data\", components = "
               "(\"vs\"), knots = ((1, 2)), interpolation = \"linear\", boundary = \"clamp\")\n"
               "m.vs(z)",
               "a grid of the kind \"asagi\" has no attribute `knots`");
  }

  TEST_CASE("the values of an inline grid nest as its axes and components") {
    // two axes: a row per knot of axis 1, a node per knot of axis 0, the components of a node
    const auto program = mustCompile("grid t (kind = \"inline\", components = (\"a\", \"b\"),\n"
                                     "        knots = ((0, 1, 2), (10, 20)),\n"
                                     "        values = (((1, 2), (3, 4), (5, 6)),\n"
                                     "                  ((7, 8), (9, 10), (11, 12))),\n"
                                     "        interpolation = \"linear\")\n"
                                     "t.a(x, z) + t.b(x, z)");
    REQUIRE(program.grids().size() == 1);
    REQUIRE(program.grids()[0].data != nullptr);
    CHECK(program.grids()[0].data->values ==
          std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});

    const std::string head = "grid t (kind = \"inline\", components = (\"a\", \"b\"), knots = "
                             "((0, 1, 2), (10, 20)), interpolation = \"linear\", values = ";
    // the flat form of before, a row too short, a node without its second component, a number
    // where a node is, and a tuple where a number is
    mustReject(head + "(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))\nt.a(x, z)",
               "expected the 2 nodes along axis 1, not 12");
    mustReject(head + "(((1, 2), (3, 4)), ((7, 8), (9, 10), (11, 12))))\nt.a(x, z)",
               "expected the 3 nodes along axis 0, not 2");
    mustReject(head + "(((1, 2), (3, 4), (5)), ((7, 8), (9, 10), (11, 12))))\nt.a(x, z)",
               "expected the 2 components of a node (a, b), not 1");
    mustReject(head + "(((1, 2), (3, 4), 5), ((7, 8), (9, 10), (11, 12))))\nt.a(x, z)",
               "expected a tuple of the 2 components of a node (a, b)");
    mustReject(head + "(((1, 2), (3, 4), (5, (6))), ((7, 8), (9, 10), (11, 12))))\nt.a(x, z)",
               "a value of a node is a number, not a tuple");
  }

  TEST_CASE("a grid gives the attributes of its kind, and the defaults fill the rest") {
    // an ASAGI file reads the variable z where it names none, and every grid clamps by default
    const auto plain = mustCompile(GridAndDef);
    CHECK(plain.grids()[0].variable == "z");
    CHECK(plain.grids()[0].boundary == reader::datafield::BoundaryMode::Clamp);
    // a SCEC file and a table have no variable, and a table its coordinates
    const auto fault = mustCompile(TwoGridsInterleaved);
    CHECK(fault.grids()[1].kind == reader::datafield::GridKind::Scec);
    CHECK(fault.grids()[1].variable.empty());
    CHECK(fault.grids()[1].columns == std::vector<std::string>{"mu_s"});
    const auto table =
        mustCompile("grid t (kind = \"table\", path = \"profile.csv\", coordinates = (\"depth\"),\n"
                    "        components = (\"rho\", \"vs\"), interpolation = \"linear\")\n"
                    "t.vs(z)");
    CHECK(table.grids()[0].kind == reader::datafield::GridKind::Table);
    CHECK(table.grids()[0].coordinates == std::vector<std::string>{"depth"});
    CHECK(table.grids()[0].columns == std::vector<std::string>{"rho", "vs"});
    CHECK(table.grids()[0].axisOf(0) == reader::datafield::AxisKind::Knots);

    mustReject("grid f (kind = \"scec\", path = \"f.txt\", variable = \"z\", components = "
               "(\"mu_s\"), interpolation = \"linear\")\nf.mu_s(x, z)",
               "a grid of the kind \"scec\" has no attribute `variable`");
    mustReject("grid f (kind = \"scec\", path = \"f.txt\", components = (\"friction\"), "
               "interpolation = \"linear\")\nf.friction(x, z)",
               "a SCEC file has no component \"friction\"; it has s_normal, s_strike");
    mustReject("grid t (kind = \"table\", path = \"t.csv\", components = (\"v\"), interpolation "
               "= \"linear\")\nt.v(z)",
               "gives no coordinates");
    mustReject("grid t (kind = \"table\", path = \"t.csv\", coordinates = (\"z\"), components = "
               "(\"z\"), interpolation = \"linear\")\nt.z(z)",
               "the column \"z\" of a table is a coordinate and a component");
    mustReject("grid t (kind = \"table\", path = \"t.csv\", coordinates = (\"z\"), components = "
               "(\"v\"), interpolation = \"cubic\")\nt.v(z)",
               "not cubic");
    mustReject("grid l (kind = \"inline\", path = \"l.txt\", components = (\"v\"), knots = ((1, "
               "2)), values = ((1), (2)), interpolation = \"linear\")\nl.v(z)",
               "a grid of the kind \"inline\" has no attribute `path`");
  }

  TEST_CASE("a grid declared with attributes states its boundary and its time axis") {
    const auto program = mustCompile(
        "grid m (kind = \"asagi\", path = \"model.nc\", variable = \"data\",\n"
        "        components = (\"rho\", \"mu\"), interpolation = \"cubic\", boundary = \"fail\",\n"
        "        time = 3)\n"
        "m.mu(x, y, z, t)");
    REQUIRE(program.grids().size() == 1);
    const auto& desc = program.grids()[0];
    CHECK(desc.path == "model.nc");
    CHECK(desc.variable == "data");
    CHECK(desc.interpolation == reader::datafield::Interpolation::Cubic);
    CHECK(desc.boundary == reader::datafield::BoundaryMode::Fail);
    REQUIRE(desc.timeAxis.has_value());
    CHECK(*desc.timeAxis == 3);
    bool sampled = false;
    for (std::size_t i = 0; i < program.arena().size(); ++i) {
      const auto& node = program.arena()[static_cast<NodeId>(i)];
      sampled = sampled || (node.kind == Kind::Lookup && node.comp == 1 && node.argCount == 4);
    }
    CHECK(sampled);

    // without them, it clamps and has no time axis
    const auto plain = mustCompile(GridAndDef);
    CHECK(plain.grids()[0].boundary == reader::datafield::BoundaryMode::Clamp);
    CHECK_FALSE(plain.grids()[0].timeAxis.has_value());

    // one component may stand alone; the time axis is optional
    const auto single = mustCompile(
        "grid g (kind = \"asagi\", path = \"g.nc\", variable = \"z\", components = \"v\",\n"
        "        interpolation = \"linear\", boundary = \"clamp\")\ng.v(x, y)");
    CHECK_FALSE(single.grids()[0].timeAxis.has_value());
  }

  TEST_CASE("a grid declared with attributes gives all of them, once, and only those") {
    const std::string head = "grid m (kind = \"asagi\", path = \"m.nc\", variable = \"z\", "
                             "components = (\"rho\")";
    mustReject(head + ")\nm.rho(x, y, z)", "gives no interpolation");
    mustReject(head + ", interpolation = \"linear\", boundary = \"wrap\")\nm.rho(x, y, z)",
               "unknown boundary \"wrap\"");
    mustReject(head + ", interpolation = \"linear\", colour = \"red\")\nm.rho(x, y, z)",
               "unknown attribute `colour`");
    mustReject(head + ", interpolation = \"linear\", boundary = \"clamp\", boundary = "
                      "\"fail\")\nm.rho(x, y, z)",
               "`boundary` is given twice");
    mustReject(head + ", interpolation = \"linear\", time = 1.5)\nm.rho(x, y, z)",
               "the time is an axis");
    mustReject(head + ", interpolation = \"linear\", time = 3)\nm.rho(x, y, z)",
               "axis 3 is the time");
  }

  TEST_CASE("the input signature is inferred with no declaration at all") {
    // Every unresolved name that is not a const, builtin, def, let or grid
    // component is an input channel. Unlike the Lua path this needs no
    // parameter-name trick and has no counterpart to output_parameters.
    const auto program = mustCompile(GridAndDef);
    std::vector<std::string> inputs;
    for (const auto& in : program.inputs()) {
      inputs.push_back(in.name);
    }
    CHECK(inputs == std::vector<std::string>{"x", "y", "z"});
  }

  TEST_CASE("grid and def declarations interleave in any order") {
    const auto program = mustCompile(TwoGridsInterleaved);
    CHECK(program.grids().size() == 2);
  }

  TEST_CASE("a two-dimensional grid is sampled with two coordinates") {
    const auto program = mustCompile(TwoGridsInterleaved);
    bool sawTwo = false;
    for (std::size_t i = 0; i < program.arena().size(); ++i) {
      const auto& node = program.arena()[static_cast<NodeId>(i)];
      if (node.kind == Kind::Lookup && node.argCount == 2) {
        sawTwo = true;
      }
    }
    CHECK(sawTwo);
  }

  TEST_CASE("select and the comparison builtins spell the same as on the Lua side") {
    const auto program = mustCompile(SelectOverLookup);
    const auto& root = program.arena()[program.roots().at(0)];
    REQUIRE(root.kind == Kind::PW);
    CHECK(root.fn == Fn::Select);
    CHECK(program.arena()[root.a].fn == Fn::Lt);
  }

  TEST_CASE("comparisons and logic as operators, and if ... then ... else, are the functions") {
    // the same programs as with the functions of the core format
    const auto same = [](const std::string& operators, const std::string& functions) {
      CAPTURE(operators);
      CHECK(compileSderivModule(operators).fingerprint() ==
            compileSderivModule(functions).fingerprint());
    };
    same("out def r = if x < 0 then a else b\n", "out def r = select(lt(x, 0), a, b)\n");
    same("out def r = x <= y\n", "out def r = le(x, y)\n");
    // gt and ge as in the full language: lt and le of the swapped operands, which the inputs read
    // in the order of the source
    same("out def r = x > y\n", "out def r = gt(x, y)\n");
    same("out def r = x >= y\n", "out def r = ge(x, y)\n");
    same("out def r = x == y\n", "out def r = eq(x, y)\n");
    same("out def r = x != y\n", "out def r = eq(eq(x, y), 0)\n");
    same("out def r = not x\n", "out def r = eq(x, 0)\n");
    // a truth value as it is, another number as one
    same("out def r = x < 0 and y < 0\n", "out def r = select(lt(x, 0), lt(y, 0), 0)\n");
    same("out def r = x or y\n", "out def r = select(x, 1, lor(y, 0))\n");
    // the precedences: not over the comparisons, and over or, the sums over the comparisons; the
    // last branch of an if as far as it goes
    same("out def r = not x < y\n", "out def r = eq(lt(x, y), 0)\n");
    same("out def r = a or b and c\n", "out def r = select(a, 1, select(b, lor(c, 0), 0))\n");
    same("out def r = x + 1 < 2 * y\n", "out def r = lt(x + 1, 2 * y)\n");
    same("out def r = if a then 1 else if b then 2 else 3 + x\n",
         "out def r = select(a, 1, select(b, 2, 3 + x))\n");
    same("out def r = 1 + (if a then 1 else 2)\n", "out def r = 1 + select(a, 1, 2)\n");
    // field by field
    same("out if z < 0 then {a = 1, b = x} else {a = 2, b = y}\n",
         "out def a = select(lt(z, 0), 1, 2)\nout def b = select(lt(z, 0), x, y)\n");
  }

  TEST_CASE("the operators give 1 and 0, and a NaN is true") {
    const auto program = compileSderivModule("out def less = x < y\n"
                                             "out def greater = x > y\n"
                                             "out def differ = x != y\n"
                                             "out def both = x and y\n"
                                             "out def either = x or y\n"
                                             "out def neither = not x\n"
                                             "out def pick = if x then 2 else 3\n");
    reader::datafield::GridStore store;
    const double nan = std::nan("");
    // x, y at four points: (1, 2), (0, 0), (NaN, 1), (2.5, 0)
    const std::vector<double> inputs{1.0, 0.0, nan, 2.5, 2.0, 0.0, 1.0, 0.0};
    const auto values = interpretValues(program, store, inputs, 4);
    REQUIRE(values.size() == 7 * 4);
    const auto row = [&](std::size_t output) {
      return std::vector<double>(values.begin() + static_cast<std::ptrdiff_t>(output * 4),
                                 values.begin() + static_cast<std::ptrdiff_t>(output * 4 + 4));
    };
    CHECK(row(0) == std::vector<double>{1, 0, 0, 0});
    CHECK(row(1) == std::vector<double>{0, 0, 0, 1});
    CHECK(row(2) == std::vector<double>{1, 0, 1, 1});
    CHECK(row(3) == std::vector<double>{1, 0, 1, 0});
    CHECK(row(4) == std::vector<double>{1, 0, 1, 1});
    CHECK(row(5) == std::vector<double>{0, 1, 0, 0});
    CHECK(row(6) == std::vector<double>{2, 3, 2, 2});
  }

  TEST_CASE("comparisons do not chain, and the negation is not") {
    mustReject("x < y < z", "comparisons do not chain");
    mustReject("x == y != z", "comparisons do not chain");
    mustReject("!x", "the negation is `not`");
    mustReject("if x then 1", "expected `else`");
    mustReject("1 + if x then 1 else 2", "unexpected");
  }

  TEST_CASE("backslashes and quotes survive a file path") {
    const auto program = mustCompile(EscapedPath);
    REQUIRE(program.grids().size() == 1);
    CHECK(program.grids().at(0).path == R"(C:\models\a"b.nc)");
  }

  TEST_CASE("let bindings share rather than duplicate") {
    const auto program = mustCompile(SharedSubexpression);
    int sqrtCount = 0;
    for (std::size_t i = 0; i < program.arena().size(); ++i) {
      const auto& node = program.arena()[static_cast<NodeId>(i)];
      if (node.kind == Kind::PW && node.fn == Fn::Sqrt) {
        ++sqrtCount;
      }
    }
    CHECK(sqrtCount == 1);
  }

  // ------------------------------------------------- records and arrays ---

  TEST_CASE("a record holds named values, which operations take field by field") {
    // the fields of shallow in another order: matched by name, in the order of deep
    const auto records = compileSderivModule("def deep = {rho = 2700, mu = 3e10}\n"
                                             "def shallow = {mu = 1e10 + z, rho = 2600 - 0.5 * z}\n"
                                             "out select(lt(z, -2000), deep, shallow)\n");
    const auto numbers =
        compileSderivModule("out def rho = select(lt(z, -2000), 2700, 2600 - 0.5 * z)\n"
                            "out def mu = select(lt(z, -2000), 3e10, 1e10 + z)\n");
    REQUIRE(records.outputs().size() == 2);
    CHECK(records.outputs()[0].name == "rho");
    CHECK(records.outputs()[1].name == "mu");
    CHECK(records.fingerprint() == numbers.fingerprint());
    // a number with every field, and a field of a record
    const auto scaled = compileSderivModule("def p = {a = x, b = 2 * x}\n"
                                            "out 0.5 * p + {b = 1, a = 2}\n"
                                            "out def c = p.b\n");
    reader::datafield::GridStore store;
    const auto values = interpretValues(scaled, store, {3.0}, 1);
    CHECK(values == std::vector<double>{3.5, 4.0, 6.0});
  }

  TEST_CASE("a function takes and gives records, and a grid samples all of its components") {
    const auto records = compileSderivModule(
        "grid m (kind = \"asagi\", path = \"model.nc\", components = (\"rho\", \"vp\", \"vs\"),\n"
        "        interpolation = \"linear\")\n"
        "def lame(p) = {rho = p.rho, mu = p.rho * p.vs ** 2,\n"
        "               lambda = p.rho * (p.vp ** 2 - 2 * p.vs ** 2)}\n"
        "out lame(m(x, y, z))\n");
    const auto numbers = compileSderivModule(
        "grid m (kind = \"asagi\", path = \"model.nc\", components = (\"rho\", \"vp\", \"vs\"),\n"
        "        interpolation = \"linear\")\n"
        "out def rho = m.rho(x, y, z)\n"
        "out def mu = m.rho(x, y, z) * m.vs(x, y, z) ** 2\n"
        "out def lambda = m.rho(x, y, z) * (m.vp(x, y, z) ** 2 - 2 * m.vs(x, y, z) ** 2)\n");
    CHECK(records.fingerprint() == numbers.fingerprint());
    int lookups = 0;
    for (std::size_t i = 0; i < records.arena().size(); ++i) {
      lookups += records.arena()[static_cast<NodeId>(i)].kind == Kind::Lookup ? 1 : 0;
    }
    CHECK(lookups == 3);
  }

  TEST_CASE("arrays and matrices, their elements and the intrinsics of them") {
    const auto program = compileSderivModule(
        "def n = [x, y, z]\n"
        "def s = [[1, 2, 3], [2, 4, 5], [3, 5, 6]]\n"
        "out def a = dot(matmul(s, n), n)\n"
        "out def b = norm(cross(n, [0, 0, 1]))\n"
        "out def c = det(s) + trace(s)\n"
        "out def d = transpose(s)[0][2] + s[2][0] + len(n) + sum(range(4)) + max([x, y, z]) +\n"
        "            min(n)\n"
        "out def e = matmul(n, s)[1]\n"
        "out def f = matmul(eye(3), n)[3 - 1] + (n * 2)[0]\n"
        "def twice(v) = 2 * v\n"
        "out def h = sum(twice(range(3)))\n");
    reader::datafield::GridStore store;
    const auto values = interpretValues(program, store, {1.0, 2.0, 3.0}, 1);
    REQUIRE(values.size() == 7);
    CHECK(values[0] == 157.0);
    CHECK(values[1] == std::sqrt(5.0));
    CHECK(values[2] == 10.0);
    CHECK(values[3] == 19.0);
    CHECK(values[4] == 25.0);
    CHECK(values[5] == 5.0);
    CHECK(values[6] == 6.0);
    // a definition may take the name of an intrinsic, which it hides
    const auto hidden = compileSderivModule("def sum(a, b) = a - b\nout def r = sum(x, 1)\n");
    CHECK(interpretValues(hidden, store, {5.0}, 1) == std::vector<double>{4.0});
  }

  TEST_CASE("records and arrays that do not fit, and what is not a number, are errors") {
    const auto reject = [](const std::string& source, const std::string& fragment) {
      bool thrown = false;
      try {
        static_cast<void>(compileSderivModule(source));
      } catch (const SderivError& error) {
        thrown = true;
        const std::string message = error.what();
        CAPTURE(message);
        CAPTURE(fragment);
        CHECK(message.find(fragment) != std::string::npos);
      }
      CAPTURE(source);
      CHECK(thrown);
    };
    reject("out {a = 1} + {b = 2}\n",
           "the operands of `+` do not fit: a record of a and a record of b");
    reject("out def r = sum([1, 2] + [1, 2, 3])\n",
           "the operands of `+` do not fit: an array of 2 and an array of 3");
    reject("out def r = x.rho\n", "`.rho` is a field of a record, not of a number");
    reject("out def r = {a = 1}.b\n", "a record of a has no field `b`");
    reject("out def r = [1, 2][2]\n", "the index is a whole number from 0 to 1, not 2");
    reject("out def r = [1, 2][x]\n", "the index is a whole number the program fixes");
    reject("out def r = {a = 1}\n", "the output `r` is a number, not a record of a");
    reject("out [1, 2]\n", "`out` of an expression gives the fields of a record as outputs");
    reject("out {a = [1, 2]}\n", "the output `a` is a number, not an array of 2");
    reject("out def r = {a = 1, a = 2}.a\n", "the record has the field `a` twice");
    reject("out def r = matmul([1, 2], [3, 4])\n", "their product is dot(a, b)");
    reject("out def r = det([[1, 2], [3, 4], [5, 6]])\n", "takes a square matrix, not 3 by 2");
    reject("out def r = det(eye(4))\n", "takes a matrix of 1 by 1 to 3 by 3");
    reject("out def r = cross([1, 2], [3, 4])\n", "takes two vectors of 3");
    reject("out def r = sum(x)\n", "`sum` takes an array, not a number");
    reject("out def r = len(range(x))[0]\n", "a whole number the program fixes");
    reject("check {a = 1}, \"never\"\nout def r = x\n",
           "the condition of a check is a number, not a record of a");
  }

  TEST_CASE("a type names the form of a value, which a value of it takes") {
    registerRecordType("elastic", {"rho", "mu", "lambda"});
    // a record with more fields is of the type, as its fields of the type
    const auto typed = compileSderivModule(
        "type velocities = {rho, vp, vs}\n"
        "def lame(p : velocities) : elastic =\n"
        "  {rho = p.rho, mu = p.rho * p.vs ** 2, lambda = p.rho * (p.vp ** 2 - 2 * p.vs ** 2)}\n"
        "def v = {qp = 100, rho = 2700 + x, vp = 6000, vs = 3464}\n"
        "out lame(v) : elastic\n");
    REQUIRE(typed.outputs().size() == 3);
    CHECK(typed.outputs()[0].name == "rho");
    CHECK(typed.outputs()[1].name == "mu");
    CHECK(typed.outputs()[2].name == "lambda");
    const auto untyped =
        compileSderivModule("out def rho = 2700 + x\n"
                            "out def mu = (2700 + x) * 3464 ** 2\n"
                            "out def lambda = (2700 + x) * (6000 ** 2 - 2 * 3464 ** 2)\n");
    CHECK(typed.fingerprint() == untyped.fingerprint());
    // the material of the point: material.rho reads the input rho, material : t those of t
    const auto material = compileSderivModule("type shear = {rho, mu}\n"
                                              "def m = material : shear\n"
                                              "out def vs = sqrt(m.mu / m.rho)\n"
                                              "out def r = material.rho * 2\n");
    std::set<std::string> inputs;
    for (const auto& input : material.inputs()) {
      inputs.insert(input.name);
    }
    CHECK(inputs == std::set<std::string>{"rho", "mu"});
    // arrays of a length, and records of records
    const auto nested = compileSderivModule("type side = {material : elastic, normal : number[3]}\n"
                                            "def s : side = {normal = [x, 0, 1], material =\n"
                                            "  {lambda = 2, mu = 1, rho = 3}, extra = 5}\n"
                                            "out def a = s.material.mu + s.normal[0]\n");
    reader::datafield::GridStore store;
    CHECK(interpretValues(nested, store, {4.0}, 1) == std::vector<double>{5.0});
  }

  TEST_CASE("a value that is no value of its type, and a type that is none, are errors") {
    registerRecordType("elastic", {"rho", "mu", "lambda"});
    const auto reject = [](const std::string& source, const std::string& fragment) {
      bool thrown = false;
      try {
        static_cast<void>(compileSderivModule(source));
      } catch (const SderivError& error) {
        thrown = true;
        const std::string message = error.what();
        CAPTURE(message);
        CAPTURE(fragment);
        CHECK(message.find(fragment) != std::string::npos);
      }
      CAPTURE(source);
      CHECK(thrown);
    };
    reject("out {rho = 1, mu = 2} : elastic\n",
           "the value is no `elastic`: it has no field `lambda`");
    reject("def f(p : elastic) = p.rho\nout def r = f({rho = 1})\n",
           "the argument p of f is no `elastic`: it has no field `mu`");
    reject("def f(a) : number = [a]\nout def r = f(1)\n",
           "the value of f is no `number`: it is an array of 1");
    reject("out def r = sum([1, 2] : number[3])\n", "is no `number[3]`: it is an array of 2");
    reject("type t = {a, b}\nout def r = ({a = 1, b = [2]} : t).a\n",
           "the value, field b, is no `number`: it is an array of 1");
    reject("out def r = (x : vector)[0]\n", "unknown type `vector`");
    reject("type elastic = {rho}\nout def r = x\n", "the type `elastic` is one of SeisSol");
    reject("type t = {a}\ntype t = {b}\nout def r = x\n", "duplicate type `t`");
    reject("type a = {b : c}\ntype c = {d : a}\nout def r = x\n", "is defined in terms of itself");
    reject("type t = {a, a}\nout def r = x\n", "the type has the field `a` twice");
    reject("out def r = material\n", "`material` is the material of the point");
    reject("out def r = (material : number[3])[0]\n", "`material` is a record of the parameters");
  }

  TEST_CASE("@, comprehensions, and states of arrays and records") {
    const auto program = compileSderivModule(
        "def m = [[1, 2], [3, 4]]\n"
        "def v = [x, y]\n"
        "out def a = (m @ v)[1] + (v @ m)[0]\n"
        "out def b = sum([i * i for i in range(4)])\n"
        "out def c = [[r[j] for r in m] for j in range(2)][1][0]\n"
        "state peak = {v = [0, 0], n = 0}\n"
        "def peak = {v = [max(peak.v[0], x), max(peak.v[1], y)], n = peak.n + 1}\n"
        "out def p = peak.v[1] + peak.n\n");
    std::vector<std::string> states;
    for (const auto& state : program.state()) {
      states.push_back(state.name);
    }
    CHECK(states == std::vector<std::string>{"peak.v[0]", "peak.v[1]", "peak.n"});
    reader::datafield::GridStore store;
    const std::vector<double> stateValues = {5.0, 6.0, 2.0};
    const auto values = interpretValues(program, store, {1.0, 2.0}, 1, nullptr, &stateValues);
    // m @ v = [5, 11], v @ m = [7, 10]; 0 + 1 + 4 + 9; the transpose of m, [[1, 3], [2, 4]]; and
    // the states of the previous call; then the next states
    REQUIRE(values.size() == 7);
    CHECK(values[0] == 18.0);
    CHECK(values[1] == 14.0);
    CHECK(values[2] == 2.0);
    CHECK(values[3] == 8.0);
    CHECK(values[4] == 5.0);
    CHECK(values[5] == 6.0);
    CHECK(values[6] == 3.0);
    // the next value has the shape of the initial one
    bool thrown = false;
    try {
      static_cast<void>(compileSderivModule("state s = [0, 0]\ndef s = [x]\nout def r = s[0]\n"));
    } catch (const SderivError& error) {
      thrown = true;
      CHECK(std::string(error.what())
                .find("the next value of the state `s` has the shape of "
                      "its initial value") != std::string::npos);
    }
    CHECK(thrown);
  }

  TEST_CASE("a module imports modules of the library and files, by qualified names") {
    registerRecordType("elastic", {"rho", "mu", "lambda"});
    const auto directory = std::filesystem::temp_directory_path() /
                           ("seissol-test-import-" + std::to_string(::getpid()));
    std::filesystem::create_directories(directory / "sub");
    const auto write = [&](const std::string& name, const std::string& text) {
      std::ofstream(directory / name) << text;
    };
    // a file, which imports another, relative to itself; its own names unqualified
    write("lib.sderiv",
          "import \"sub/base.sderiv\" as base\n"
          "type pair = {a, b}\n"
          "def half = 0.5\n"
          "def mean(p : pair) = (p.a + p.b) * half\n"
          "def scaled(p : pair) : pair = {a = base.twice(p.a), b = p.b}\n");
    write("sub/base.sderiv", "def twice(v) = 2 * v\n");
    SderivOptions options;
    options.file = (directory / "model.sderiv").string();
    const auto program =
        compileSderivModule("import \"lib.sderiv\" as lib\n"
                            "import motion\n"
                            "def p : lib.pair = {a = x, b = 3, c = 4}\n"
                            "out def m = lib.mean(lib.scaled(p))\n"
                            "out def h = lib.half\n"
                            "out def r = max(motion.rotated(x, 0, motion.angles(4)))\n",
                            options);
    reader::datafield::GridStore store;
    // mean of 2 x and 3, at x = 2.5; and |2.5 cos t| at 0, 45, 90 and 135 degrees
    CHECK(interpretValues(program, store, {2.5}, 1) == std::vector<double>{4.0, 0.5, 2.5});

    const auto reject =
        [&](const std::string& library, const std::string& source, const std::string& fragment) {
          write("bad.sderiv", library);
          bool thrown = false;
          try {
            static_cast<void>(compileSderivModule(source, options));
          } catch (const SderivError& error) {
            thrown = true;
            const std::string message = error.what();
            CAPTURE(message);
            CAPTURE(fragment);
            CHECK(message.find(fragment) != std::string::npos);
          }
          CAPTURE(source);
          CHECK(thrown);
        };
    const std::string uses = "import \"bad.sderiv\" as bad\nout def r = bad.f(1)\n";
    reject("def f(v) = v\nout def g = 1\n",
           uses,
           "a library declares definitions, types and imports, not `out`");
    // an error in a module names its file and its line
    reject("def f(v) = v + z\n",
           uses,
           "bad.sderiv:1:16: resolve: a definition of a library reads its parameters, the "
           "definitions of its module and constants, and `z` is none of them");
    reject("import \"bad.sderiv\" as again\ndef f(v) = v\n", uses, "imports itself");
    reject(
        "def f(v) = v\n", "import tensors\nout def r = 1\n", "the library has no module `tensors`");
    reject("def f(v) = v\n",
           "import \"bad.sderiv\" as bad\nimport motion as bad\nout def r = bad.f(1)\n",
           "the module is imported as `bad` twice");
    reject("def f(v) = v\n",
           "import \"bad.sderiv\" as bad\nimport motion\ndef motion = 1\nout def r = bad.f(1)\n",
           "the module is imported as `motion`, the name of a definition");
    reject("def f(v) = v\n",
           "import \"bad.sderiv\" as bad\nout def r = bad.g(1)\n",
           "the module `bad` has no definition `g`");
    reject(
        "def f(v) = v\n", "import \"bad.sderiv\" as bad\nout def r = bad\n", "`bad` is a module");
    std::filesystem::remove_all(directory);
  }

  TEST_CASE("the library: rocks, rotations, and the motion along every angle") {
    registerRecordType("elastic", {"rho", "mu", "lambda"});
    // Gassmann: lambda + alpha^2 M, the density of solid and fluid
    const auto rock = compileSderivModule(
        "import rock\n"
        "out rock.gassmann({bulk_solid = 3.6e10, bulk_fluid = 2.2e9, porosity = phi, rho = 2650,\n"
        "                   rho_fluid = 1000, mu = 1e10, lambda = 5e9}) : elastic\n");
    reader::datafield::GridStore store;
    const double phi = 0.2;
    const double matrix = 5e9 + 2.0 / 3.0 * 1e10;
    const double alpha = (1.0 - matrix / 3.6e10) * (1.0 - matrix / 3.6e10);
    const double inverse = phi / 2.2e9 + (1.0 - phi) / 3.6e10 - matrix / (3.6e10 * 3.6e10);
    const auto rocks = interpretValues(rock, store, {phi}, 1);
    REQUIRE(rocks.size() == 3);
    CHECK(rocks[0] == doctest::Approx(phi * 1000 + (1 - phi) * 2650).epsilon(1e-15));
    CHECK(rocks[1] == 1e10);
    CHECK(rocks[2] == doctest::Approx(5e9 + alpha / inverse).epsilon(1e-15));
    // the rotation of a point source: by strike 90, dip 90 and rake 0 the axes are permuted, and
    // a rotated tensor keeps its trace
    const auto rotation = compileSderivModule(
        "import rotation\n"
        "def r = rotation.matrix(90, 90, 0)\n"
        "def m = rotation.rotate([[1, 0, 0], [0, 2, 0], [0, 0, 3]], rotation.matrix(s, 40, 10))\n"
        "out def r20 = r[2][0]\n"
        "out def tr = trace(m)\n");
    const auto rotated = interpretValues(rotation, store, {25.0}, 1);
    CHECK(rotated[0] == 1.0);
    CHECK(rotated[1] == doctest::Approx(6.0).epsilon(1e-14));
    // RotD100 over two calls, its state an array of a state per angle
    const auto motion = compileSderivModule("import motion as m\n"
                                            "def angles = m.angles(180)\n"
                                            "state peak = [0 for t in angles]\n"
                                            "def now = max(peak, m.rotated(v1, v2, angles))\n"
                                            "def peak = now\n"
                                            "out def rotd100 = max(now)\n");
    CHECK(motion.state().size() == 180);
    std::vector<double> states(180, 0.0);
    const auto first = interpretValues(motion, store, {3.0, 4.0}, 1, nullptr, &states);
    // |3 cos t + 4 sin t| is largest, 5, at t = atan(4 / 3) = 53.13 degrees: the angle 53
    CHECK(first[0] == doctest::Approx(5.0).epsilon(1e-3));
    CHECK(first[0] <= 5.0);
  }

  // ------------------------------------------------------------ negative ---

  TEST_CASE("strings are confined to grid declarations") {
    mustReject(StringInExpression, "only allowed in a `grid` declaration");
  }

  TEST_CASE("a bare name where a string belongs names the expected kind") {
    // The obvious first mistake with this syntax, so the message says what was
    // wanted rather than "unexpected token".
    mustReject(BareNameForString, "expected a quoted string");
  }

  TEST_CASE("an unknown grid kind is rejected at parse time") {
    mustReject(UnknownKind, "unknown kind");
  }

  TEST_CASE("an interpolation outside the vocabulary is caught before load") {
    mustReject(UnknownInterpolation, "unknown interpolation \"model.nc\"");
  }

  TEST_CASE("a grid must declare at least one component") {
    mustReject(NoComponents, "gives no components");
  }

  TEST_CASE("a duplicate component is rejected") {
    mustReject(DuplicateComponent, "duplicate component");
  }

  TEST_CASE("a grid shares the namespace of the definitions") {
    mustReject(GridNamedLikeDefinition, "has the name of a definition");
  }

  TEST_CASE("a sample names a grid and one of its components, and has coordinates") {
    mustReject(NoSuchComponent, "the grid `m` has no component `vp`; it has rho");
    mustReject(NoSuchGrid, "`n` is no grid");
    mustReject(GridAsValue, "`m` is a grid; a component of it is sampled as m.c(x, y, z)");
    mustReject(ComponentWithoutCoordinates, "sampled with coordinates: m.rho(x, y, z)");
  }

  TEST_CASE("a grid gives its attributes in parentheses") {
    mustReject(PositionalGrid, "its attributes follow in parentheses");
  }

  TEST_CASE("an unterminated string points at its opening quote") {
    mustReject(UnterminatedString, "unterminated string");
  }

  TEST_CASE("the coordinate count is checked against the IR limit") {
    mustReject(WrongCoordinateCount, "1..6");
  }

  // -------------------------------------------------------- cross-frontend -

  TEST_CASE("the same model traces to the same fingerprint from Lua and sderiv") {
    // The point of one shared IR: two frontends, one Program. If this ever
    // fails it is a lowering difference between the frontends, not a modelling
    // difference, because both sources are written to be the same expression.
    const auto sderiv = mustCompile(R"(
grid m (kind = "asagi", path = "model.nc", variable = "z", components = ("rho"),
        interpolation = "linear", boundary = "clamp")
select(lt(z, -1000), m.rho(x, y, z), 2700)
)");
    const auto lua = mustCompile(SelectOverLookup);
    CHECK(sderiv.fingerprint() == lua.fingerprint());
  }

  // ------------------------------------------------------ out def ----------

  TEST_CASE("a module names its own outputs and they can read each other") {
    const auto program = compileSderivModule(R"SD(
# a small material model
grid m (kind = "asagi", path = "model.nc", variable = "data", components = ("rho", "vp", "vs"),
        interpolation = "linear")

def vs = m.vs(x, y, z)

out def rho    = m.rho(x, y, z)
out def mu     = rho * vs * vs
out def lambda = rho * m.vp(x, y, z)**2.0 - 2.0 * mu
)SD");

    REQUIRE(program.outputs().size() == 3);
    CHECK(program.outputs()[0].name == "rho");
    CHECK(program.outputs()[1].name == "mu");
    CHECK(program.outputs()[2].name == "lambda");
    CHECK(program.grids().size() == 1);
    CHECK(program.roots().size() == 3);
  }

  TEST_CASE("an output read by a later output is interned, not duplicated") {
    // The reason the outputs live in ONE module rather than three sources: the
    // shared subexpression must be one node, or the module form buys nothing
    // over three calls to compileSderiv.
    const auto shared = compileSderivModule("out def a = sqrt(x)\nout def b = a + a\n");
    const auto separate =
        compileSderivModule("out def a = sqrt(x)\nout def b = sqrt(x) + sqrt(x)\n");
    CHECK(shared.arena().size() == separate.arena().size());
  }

  TEST_CASE("comments and blank lines are ignored") {
    const auto withComments = compileSderivModule(R"SD(
# leading comment
out def a = x + 1.0   # trailing comment
#
out def b = y
)SD");
    const auto without = compileSderivModule("out def a = x + 1.0\nout def b = y\n");
    CHECK(withComments.fingerprint() == without.fingerprint());
  }

  TEST_CASE("out def is refused where it would not mean anything") {
    // A function has no call site to take arguments from, so an exported one
    // cannot be evaluated per point.
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("out def f(a) = a * 2.0\n")),
                    SderivError);

    // No outputs at all: the module produces nothing, which is a configuration
    // error rather than an empty program.
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("def a = x\n")), SderivError);

    // Both forms at once would leave it open which one is "the" output.
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("out def a = x\nx + 1.0\n")),
                    SderivError);

    // A duplicate plain def is an error, not last-wins.
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("def a = x\ndef a = y\nout def b = a\n")),
                    SderivError);
  }

  TEST_CASE("a definition may not take the name of a constant or function") {
    // a read of `g` would be the constant, and the definition ignored
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("def g = 2.0 * x\nout def u = g\n")),
                    SderivError);
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("def exp = x\nout def u = exp\n")),
                    SderivError);
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("out def pi = x\n")), SderivError);
  }

  TEST_CASE("a definition that reads itself is an error, not a recursion") {
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("out def a = a + 1.0\n")), SderivError);
    CHECK_THROWS_AS(
        static_cast<void>(compileSderivModule("def a = b * 2.0\ndef b = a\nout def c = a\n")),
        SderivError);
    CHECK_THROWS_AS(
        static_cast<void>(compileSderivModule("def f(p) = f(p) + 1.0\nout def c = f(x)\n")),
        SderivError);
    // a function applied to itself is not a cycle
    CHECK_NOTHROW(
        static_cast<void>(compileSderivModule("def f(p) = 2.0 * p\nout def c = f(f(x))\n")));
  }

  TEST_CASE("an input the consumer names reads the input next to a definition of its name") {
    // the inner state of a boundary onto its ghost state, as one parallel assignment
    SderivOptions options;
    options.inputs = {"v1", "s_xy"};
    const auto program = compileSderivModule(
        "out def v1 = -v1\nout def s_xy = s_xy + 0.5 * v1\nout def s_yy = 3.0\n", options);
    std::set<std::string> inputs;
    for (const auto& input : program.inputs()) {
      inputs.insert(input.name);
    }
    CHECK(inputs == std::set<std::string>{"v1", "s_xy"});
    REQUIRE(program.outputs().size() == 3);
    CHECK(channelsRead(program, {program.roots()[1]}) == std::set<std::string>{"v1", "s_xy"});

    // evaluated under other names for the outputs, as a table has one column per name
    Program renamed;
    renamed.arena() = program.arena();
    for (const auto& input : program.inputs()) {
      renamed.addInput(input.name, input.type);
    }
    for (std::size_t i = 0; i < program.outputs().size(); ++i) {
      renamed.addOutput(
          "ghost_" + program.outputs()[i].name, program.outputs()[i].type, program.roots()[i]);
    }
    const double v1 = 2.0;
    const double sxy = 3.0;
    std::vector<double> ghost(3, -1.0);
    reader::scripting::DataTable table(1);
    table.bindViewConst<double>("v1", reader::scripting::Direction::In, &v1);
    table.bindViewConst<double>("s_xy", reader::scripting::Direction::In, &sxy);
    for (std::size_t i = 0; i < ghost.size(); ++i) {
      table.bindView<double>(
          renamed.outputs()[i].name, reader::scripting::Direction::Out, &ghost[i]);
    }
    auto binding = Binding::bind(renamed, table);
    reader::datafield::GridStore store;
    const auto kernel = makeKernel(renamed, binding, store, BackendOptions{});
    kernel->precompute(table);
    kernel->run(table);
    CHECK(ghost[0] == -2.0);
    // the input v1, not the ghost value
    CHECK(ghost[1] == 4.0);
    CHECK(ghost[2] == 3.0);

    // without the option, the same module reads itself
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule("out def v1 = -v1\n")), SderivError);
  }

  TEST_CASE("every definition is lowered once, however often it is read") {
    // every definition reads the previous one three times: lowered again at every read, the chain
    // would take 3^60 steps
    std::string source = "def v0 = x\n";
    for (int i = 1; i <= 60; ++i) {
      const auto previous = "v" + std::to_string(i - 1);
      source += "def v" + std::to_string(i) + " = " + previous + " * " + previous + " - " +
                previous + "\n";
    }
    source += "out def u = v60\n";
    const auto program = compileSderivModule(source);
    REQUIRE(program.outputs().size() == 1);
    // two operations per definition, the input and nothing else
    CHECK(program.arena().size() <= 3 * 60 + 2);

    // as does a function, per argument
    std::string calls = "def f(p) = p * p + 1.0\ndef w0 = x\n";
    for (int i = 1; i <= 40; ++i) {
      const auto previous = "w" + std::to_string(i - 1);
      calls += "def w" + std::to_string(i) + " = f(" + previous + ") + f(" + previous + ")\n";
    }
    calls += "out def u = w40\n";
    CHECK_NOTHROW(static_cast<void>(compileSderivModule(calls)));
  }

  TEST_CASE("long chains of definitions and of terms take no stack") {
    constexpr int Chain = 30000;
    std::string chain = "def d0 = x\n";
    for (int i = 1; i <= Chain; ++i) {
      chain += "def d" + std::to_string(i) + " = d" + std::to_string(i - 1) + " + 1.0\n";
    }
    chain += "out def u = d" + std::to_string(Chain) + "\n";
    std::string sum = "out def u = x";
    for (int i = 0; i < 100000; ++i) {
      sum += " + 1.0";
    }
    sum += "\n";
    for (const auto& [source, expected] :
         {std::pair{chain, Chain + 0.5}, std::pair{sum, 100000.5}}) {
      const auto program = compileSderivModule(source);
      const double x = 0.5;
      double u = 0;
      reader::scripting::DataTable table(1);
      table.bindViewConst<double>("x", reader::scripting::Direction::In, &x);
      table.bindView<double>("u", reader::scripting::Direction::Out, &u);
      auto binding = Binding::bind(program, table);
      reader::datafield::GridStore store;
      const auto kernel = makeKernel(program, binding, store, BackendOptions{});
      kernel->precompute(table);
      kernel->run(table);
      CHECK(u == expected);
    }
  }

  TEST_CASE("an expression nested too deep is an error, not a crash") {
    const auto nested = [](int depth) {
      return "out def u = " + std::string(static_cast<std::size_t>(depth), '(') + "x" +
             std::string(static_cast<std::size_t>(depth), ')') + "\n";
    };
    CHECK_NOTHROW(static_cast<void>(compileSderivModule(nested(256))));
    CHECK_THROWS_AS(static_cast<void>(compileSderivModule(nested(100000))), SderivError);
    CHECK_THROWS_AS(
        static_cast<void>(compileSderivModule("out def u = " + std::string(100000, '-') + "x\n")),
        SderivError);
  }

  TEST_CASE("the inputs are in the order the source reads them first") {
    // a definition is read where it is named, before what follows
    const auto program =
        compileSderivModule("def v = y + t\nout def u = x + v * z\nout def w = gt(s, x)\n");
    std::vector<std::string> inputs;
    for (const auto& input : program.inputs()) {
      inputs.push_back(input.name);
    }
    CHECK(inputs == std::vector<std::string>{"x", "y", "t", "z", "s"});
  }

  TEST_CASE("the two module forms stay separate") {
    // The trailing-expression form still works and is still named by the caller.
    const auto named = compileSderiv("x * 2.0", "rho");
    REQUIRE(named.outputs().size() == 1);
    CHECK(named.outputs()[0].name == "rho");

    // Naming an output from outside a module that names its own is a mismatch,
    // not something to resolve silently in either direction.
    CHECK_THROWS_AS(static_cast<void>(compileSderiv("out def rho = x\n", "sigma")), SderivError);
  }
}

} // namespace seissol::unit_test
