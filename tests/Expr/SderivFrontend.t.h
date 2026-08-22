// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// Same shape as LuaTracer.t.h: sources inline, positive / negative /
// differential.
//
// The negative half is shorter than the Lua one and that is the headline
// property of this frontend rather than an omission. sderiv has no control
// flow, so there is no `if` to detect, no probe ladder, no truthiness, no
// iteration order and no side effects -- the whole four-net apparatus in
// LuaTracer.h has nothing to do here. What remains to test is a grammar and a
// set of name checks.

#include <doctest.h>

#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Scripting/DataTable.h"

#include <string>
#include <vector>

namespace seissol::unit_test {

using namespace seissol::expr;

namespace {

constexpr auto GridAndDef = R"(
# density and Lame parameters from a velocity model
grid m = "asagi", "model.nc", "linear", "rho", "vp", "vs"

def mu(px, py, pz) = m_rho(px, py, pz) * m_vs(px, py, pz) ** 2

mu(x, y, z)
)";

constexpr auto TwoGridsInterleaved = R"(
grid vel = "asagi", "vel.nc", "linear", "vp", "vs"
def half(a) = a / 2
grid topo = "scec", "topo.bin", "nearest", "elevation"
half(vel_vp(x, y, z)) + topo_elevation(x, y)
)";

constexpr auto SelectOverLookup = R"(
grid m = "asagi", "model.nc", "linear", "rho"
select(lt(z, -1000), m_rho(x, y, z), 2700)
)";

constexpr auto EscapedPath = R"(grid m = "asagi", "C:\\models\\a\"b.nc", "linear", "rho"
m_rho(x, y, z))";

constexpr auto SharedSubexpression = R"(
let r = sqrt(x*x + y*y) in r + r * 2
)";

constexpr auto StringInExpression = R"(grid m = "asagi", "model.nc", "linear", "rho"
m_rho(x, y, z) + "oops")";

constexpr auto UnknownKind = R"(grid m = "hdf5", "m.nc", "linear", "rho"
m_rho(x, y, z))";

constexpr auto SwappedFileAndInterpolation = R"(grid m = "asagi", "linear", "model.nc", "rho"
m_rho(x, y, z))";

constexpr auto NoComponents = R"(grid m = "asagi", "model.nc", "linear"
x)";

constexpr auto DuplicateComponent = R"(grid m = "asagi", "model.nc", "linear", "rho", "rho"
m_rho(x, y, z))";

constexpr auto CollidesWithDef = R"(grid m = "asagi", "model.nc", "linear", "rho"
def m_rho(a) = a
m_rho(x))";

constexpr auto GridNameCollidesWithComponent = R"(grid m = "asagi", "a.nc", "linear", "rho"
grid m_rho = "asagi", "b.nc", "linear", "q"
x)";

constexpr auto UnterminatedString = R"(grid m = "asagi", "model.nc, "linear", "rho"
x)";

constexpr auto BareNameForString = R"(grid m = asagi, "model.nc", "linear", "rho"
m_rho(x, y, z))";

constexpr auto WrongCoordinateCount = R"(grid m = "asagi", "model.nc", "linear", "rho"
m_rho(x, y, z, t, x, y, z))";

Program mustCompile(const std::string& source) {
  return compileSderiv(source, /*outputName=*/"out");
}

void mustReject(const std::string& source, const std::string& fragment) {
  // The message is part of the contract here: these are user-facing syntax
  // errors in a parameter file, and "parse error" alone sends someone to the
  // wrong line.
  REQUIRE_THROWS_WITH_AS(
      compileSderiv(source, "out"), doctest::Contains(fragment.c_str()), SderivError);
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
    // m_rho and m_vs are read; m_vp is declared but unused and must not appear.
    CHECK(lookups == 2);
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
    // GridDesc::file became GridDesc::path in Package 4, when the frontends
    // were ported onto Grid.h's contract.
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

  TEST_CASE("swapping the file and the interpolation is caught before load") {
    // Both slots are strings, so only the closed vocabulary of the
    // interpolation slot makes this detectable without touching the disk.
    mustReject(SwappedFileAndInterpolation, "swapped");
  }

  TEST_CASE("a grid must declare at least one component") {
    mustReject(NoComponents, "no components");
  }

  TEST_CASE("a duplicate component is rejected") {
    mustReject(DuplicateComponent, "duplicate component");
  }

  TEST_CASE("a generated component name may not collide with a def") {
    mustReject(CollidesWithDef, "collides");
  }

  TEST_CASE("a grid may not be named like another grid's component function") {
    // Not formally ambiguous, since a grid name is never an expression -- but
    // `m_rho(...)` would be unreadable. This is what underscore flattening
    // costs, and the price of not adding a `.` token.
    mustReject(GridNameCollidesWithComponent, "same name as a component function");
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
grid m = "asagi", "model.nc", "linear", "rho"
select(lt(z, -1000), m_rho(x, y, z), 2700)
)");
    const auto lua = mustCompile(SelectOverLookup);
    CHECK(sderiv.fingerprint() == lua.fingerprint());
  }

  // ------------------------------------------------------ out def ----------

  TEST_CASE("a module names its own outputs and they can read each other") {
    const auto program = compileSderivModule(R"SD(
# a small material model
grid m = "asagi", "model.nc", "data", "linear", "rho", "vp", "vs"

def vs = m_vs(x, y, z)

out def rho    = m_rho(x, y, z)
out def mu     = rho * vs * vs
out def lambda = rho * m_vp(x, y, z)**2.0 - 2.0 * mu
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
    CHECK_THROWS_AS(compileSderivModule("out def f(a) = a * 2.0\n"), SderivError);

    // No outputs at all: the module produces nothing, which is a configuration
    // error rather than an empty program.
    CHECK_THROWS_AS(compileSderivModule("def a = x\n"), SderivError);

    // Both forms at once would leave it open which one is "the" output.
    CHECK_THROWS_AS(compileSderivModule("out def a = x\nx + 1.0\n"), SderivError);

    // A duplicate plain def used to be taken last-wins, silently.
    CHECK_THROWS_AS(compileSderivModule("def a = x\ndef a = y\nout def b = a\n"), SderivError);
  }

  TEST_CASE("the two module forms stay separate") {
    // The trailing-expression form still works and is still named by the caller.
    const auto named = compileSderiv("x * 2.0", "rho");
    REQUIRE(named.outputs().size() == 1);
    CHECK(named.outputs()[0].name == "rho");

    // Naming an output from outside a module that names its own is a mismatch,
    // not something to resolve silently in either direction.
    CHECK_THROWS_AS(compileSderiv("out def rho = x\n", "sigma"), SderivError);
  }
}

} // namespace seissol::unit_test
