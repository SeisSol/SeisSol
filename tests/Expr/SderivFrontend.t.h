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
#include <set>
#include <string>
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
