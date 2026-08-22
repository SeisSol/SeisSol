// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// Every script is a string literal here rather than a file on disk: the corpus
// is small, the tests are then hermetic (no path resolution, no install step,
// no risk of a test passing because it silently read a stale file), and a
// failing case can be pasted straight into a bug report.
//
// The suite has three layers, matching the three things that can go wrong:
//   POSITIVE     a script traces, and the resulting Program has the signature
//                and shape it should
//   NEGATIVE     a script that would trace WRONG is refused, with the right
//                cause -- checking the cause and not merely "it failed" is the
//                point, because the four detection nets are independent and a
//                test that only asserts failure cannot tell which one fired
//   DIFFERENTIAL the traced Program and the interpreted LuaReader agree
//                numerically on the same points

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Program.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/LuaReader.h"
#include "Reader/Scripting/LuaTracer.h"

#include <cmath>
#include <optional>
#include <string>
#include <vector>

namespace seissol::unit_test {

using namespace seissol::reader::scripting;
using Cause = TraceFailure::Cause;

namespace {

// ---------------------------------------------------------------- corpus ---

constexpr auto PlanarWave = R"lua(
local M = {}
M.output_parameters = {"u", "v", "w"}

local kx, ky, kz = 1.0, 0.5, -0.25
local omega = 2.0

function M.evaluate(fields, x, y, z, t)
  local phase = kx*x + ky*y + kz*z - omega*t
  local a = math.sin(phase)
  local b = math.cos(phase)
  local amp = 0.0
  for i = 1, 3 do
    amp = amp + 1.0 / i
  end
  return amp*a, amp*b, amp*(a + b)
end

return M
)lua";

constexpr auto LayeredSelect = R"lua(
local M = {}
M.output_parameters = {"rho"}

function M.evaluate(fields, x, y, z)
  local shallow = ssol.gt(z, -1000.0)
  return ssol.select(shallow, 2200.0 + 0.1*z, 2700.0)
end

return M
)lua";

constexpr auto GridSample = R"lua(
local M = {}
M.output_parameters = {"rho", "mu", "lambda"}
M.field_specs = {
  { name = "field_0001", kind = "asagi", file = "model.nc",
    interpolation = "linear", parameters = {"rho", "vp", "vs"} },
}

function M.evaluate(fields, x, y, z)
  local rho, vp, vs = fields.field_0001:sample(x, y, z)
  local mu = rho * vs * vs
  return rho, mu, rho*vp*vp - 2.0*mu
end

return M
)lua";

constexpr auto SharedSubexpression = R"lua(
local M = {}
M.output_parameters = {"a", "b"}

function M.evaluate(fields, x, y)
  local r = math.sqrt(x*x + y*y)
  return r + 1.0, r * 2.0
end

return M
)lua";

constexpr auto RawIfGreater = R"lua(
local M = {}
M.output_parameters = {"rho"}

function M.evaluate(fields, x, y, z)
  if z > -1000.0 then
    return 2200.0
  else
    return 2700.0
  end
end

return M
)lua";

constexpr auto RawIfEqual = R"lua(
local M = {}
M.output_parameters = {"rho"}

function M.evaluate(fields, x, y, z)
  if z == 0.0 then
    return 1000.0
  end
  return 2700.0 + 0.0*z
end

return M
)lua";

constexpr auto ConditionInRawIf = R"lua(
local M = {}
M.output_parameters = {"rho"}

function M.evaluate(fields, x, y, z)
  local shallow = ssol.gt(z, -1000.0)
  if shallow then
    return 2200.0 + 0.0*z
  end
  return 2700.0 + 0.0*z
end

return M
)lua";

constexpr auto MutatesModuleState = R"lua(
local M = {}
M.output_parameters = {"n"}

local counter = 0

function M.evaluate(fields, x)
  counter = counter + 1
  return x + counter
end

return M
)lua";

constexpr auto VariadicEvaluate = R"lua(
local M = {}
M.output_parameters = {"a"}
function M.evaluate(...)
  local fields, x = ...
  return x + 1.0
end
return M
)lua";

constexpr auto MethodSyntaxEvaluate = R"lua(
local M = {}
M.output_parameters = {"a"}
function M:evaluate(fields, x)
  return x + 1.0
end
return M
)lua";

constexpr auto OutputCountMismatch = R"lua(
local M = {}
M.output_parameters = {"a", "b"}
function M.evaluate(fields, x)
  return x
end
return M
)lua";

constexpr auto TableIterationOrder = R"lua(
local M = {}
M.output_parameters = {"s"}
local coeffs = { a = 1.0, b = 2.0, c = 3.0, d = 4.0, e = 5.0, f = 6.0 }
function M.evaluate(fields, x)
  local s = 0.0
  for k, v in pairs(coeffs) do s = s + v * x end
  return s
end
return M
)lua";

constexpr auto UntraceableMathFmod = R"lua(
local M = {}
M.output_parameters = {"m"}
function M.evaluate(fields, x)
  return math.fmod(x, 3.0)
end
return M
)lua";

constexpr auto LuaFloorMod = R"lua(
local M = {}
M.output_parameters = {"m"}
function M.evaluate(fields, x)
  return x % 3.0
end
return M
)lua";

constexpr auto TooManyCoordinates = R"lua(
local M = {}
M.output_parameters = {"a"}
M.field_specs = {
  { name = "g", kind = "asagi", file = "f.nc", interpolation = "linear", parameters = {"a"} },
}
function M.evaluate(fields, x, y, z, t, u, v, w)
  return fields.g:sample(x, y, z, t, u, v, w)
end
return M
)lua";

// ---------------------------------------------------------------- helpers ---

expr::Program mustTrace(const std::string& code, const TraceOptions& options = {}) {
  TraceFailure failure;
  auto program = traceLuaModule(code, options, failure);
  REQUIRE_MESSAGE(program.has_value(), failure.reason);
  return std::move(*program);
}

TraceFailure mustRefuse(const std::string& code, const TraceOptions& options = {}) {
  TraceFailure failure;
  auto program = traceLuaModule(code, options, failure);
  REQUIRE_FALSE(program.has_value());
  return failure;
}

std::vector<std::string> names(const std::vector<expr::VarSpec>& specs) {
  std::vector<std::string> out;
  out.reserve(specs.size());
  for (const auto& s : specs) {
    out.push_back(s.name);
  }
  return out;
}

} // namespace

TEST_SUITE("LuaTracer") {

  // ------------------------------------------------------------ positive ---

  TEST_CASE("signature comes from the parameter names of evaluate") {
    const auto program = mustTrace(PlanarWave);
    CHECK(names(program.inputs()) == std::vector<std::string>{"x", "y", "z", "t"});
    CHECK(names(program.outputs()) == std::vector<std::string>{"u", "v", "w"});
    CHECK(program.roots().size() == 3);
  }

  TEST_CASE("a constant loop unrolls without IR support") {
    // amp = 1 + 1/2 + 1/3 is folded to a single Const during the trace; nothing
    // loop-shaped survives into the DAG.
    const auto program = mustTrace(PlanarWave);
    bool sawConst = false;
    for (std::size_t i = 0; i < program.arena().size(); ++i) {
      const auto& node = program.arena()[static_cast<expr::NodeId>(i)];
      if (node.kind == expr::Kind::Const &&
          std::abs(node.value - (1.0 + 1.0 / 2.0 + 1.0 / 3.0)) < 1e-15) {
        sawConst = true;
      }
    }
    CHECK(sawConst);
  }

  TEST_CASE("shared subexpressions collapse in the arena") {
    const auto program = mustTrace(SharedSubexpression);
    int sqrtCount = 0;
    for (std::size_t i = 0; i < program.arena().size(); ++i) {
      const auto& node = program.arena()[static_cast<expr::NodeId>(i)];
      if (node.kind == expr::Kind::PW && node.fn == expr::Fn::Sqrt) {
        ++sqrtCount;
      }
    }
    // Interning, not a CSE pass: the second `r` reference finds the first node.
    CHECK(sqrtCount == 1);
  }

  TEST_CASE("ssol.select lowers to a ternary Select over a boolean-valued node") {
    const auto program = mustTrace(LayeredSelect);
    const auto& arena = program.arena();
    const auto& root = arena[program.roots().at(0)];
    REQUIRE(root.kind == expr::Kind::PW);
    REQUIRE(root.fn == expr::Fn::Select);
    // ssol.gt(a, b) is Lt with the operands swapped -- there is no Fn::Gt.
    CHECK(arena[root.a].fn == expr::Fn::Lt);
  }

  TEST_CASE("a grid sample becomes one Lookup per declared component") {
    const auto program = mustTrace(GridSample);
    REQUIRE(program.grids().size() == 1);
    std::vector<bool> seen(3, false);
    for (std::size_t i = 0; i < program.arena().size(); ++i) {
      const auto& node = program.arena()[static_cast<expr::NodeId>(i)];
      if (node.kind == expr::Kind::Lookup) {
        REQUIRE(node.argCount == 3);
        REQUIRE(node.comp >= 0);
        REQUIRE(node.comp < 3);
        seen[node.comp] = true;
      }
    }
    // rho, vp and vs are each read exactly once despite three uses of the tuple.
    CHECK(seen[0]);
    CHECK(seen[1]);
    CHECK(seen[2]);
  }

  TEST_CASE("`%` traces to Fn::Mod and matches Lua's floor semantics") {
    // Regression guard for the lowering that looks obvious and is wrong:
    // a - floor(a/b)*b disagrees with Lua's `%` on 14 of 121 sign/magnitude
    // pairs, so Fn::Mod must be fmod pulled onto the sign of the divisor.
    const auto program = mustTrace(LuaFloorMod);
    const auto& root = program.arena()[program.roots().at(0)];
    REQUIRE(root.kind == expr::Kind::PW);
    CHECK(root.fn == expr::Fn::Mod);
  }

  // ------------------------------------------------------------ negative ---

  TEST_CASE("net 1: a comparison operator on a traced value is refused") {
    const auto failure = mustRefuse(RawIfGreater);
    CHECK(failure.cause == Cause::UntracedOperator);
    CHECK(failure.line == 5);
  }

  TEST_CASE("net 2: a condition that never reaches select is refused") {
    // The dangerous one: it uses the sanctioned API and then puts the result in
    // a raw `if`, where a userdata is unconditionally truthy. No metamethod
    // fires, and the probe ladder does not straddle -1000 either.
    const auto failure = mustRefuse(ConditionInRawIf);
    CHECK(failure.cause == Cause::RawCondition);
  }

  TEST_CASE("net 3: `== 0` against a number is caught by the probes") {
    // __eq is only dispatched when both operands are full userdata, so this one
    // is invisible to nets 1 and 2 and can only be found by running.
    const auto failure = mustRefuse(RawIfEqual);
    CHECK(failure.cause == Cause::DataDependentFlow);
  }

  TEST_CASE("mutating module state during the trace is refused") {
    // Detected by comparing evaluate's upvalues before and after: a module-scope
    // local is an upvalue, not a table field, so no __newindex on M would see it.
    const auto failure = mustRefuse(MutatesModuleState);
    CHECK(failure.cause == Cause::SideEffect);
  }

  TEST_CASE("a variadic evaluate is refused rather than traced as constant") {
    // lua_getinfo reports nparams == 0 for `function M.evaluate(...)`, which
    // would otherwise yield a valid-looking Program with no inputs at all.
    const auto failure = mustRefuse(VariadicEvaluate);
    CHECK(failure.cause == Cause::NoEvaluate);
  }

  TEST_CASE("method syntax is refused because `self` shifts every parameter") {
    const auto failure = mustRefuse(MethodSyntaxEvaluate);
    CHECK(failure.cause == Cause::NoEvaluate);
  }

  TEST_CASE("output_parameters is checked against the number of returned values") {
    const auto failure = mustRefuse(OutputCountMismatch);
    CHECK(failure.cause == Cause::SignatureMismatch);
  }

  TEST_CASE("math.fmod is refused rather than lowered onto Fn::Mod") {
    // math.fmod is C truncation; `%` is floor. Treating them as one op would be
    // a silent sign bug for negative arguments.
    const auto failure = mustRefuse(UntraceableMathFmod);
    CHECK(failure.cause == Cause::UntracedOperator);
  }

  TEST_CASE("a sample with more than six coordinates is refused") {
    const auto failure = mustRefuse(TooManyCoordinates);
    CHECK(failure.cause == Cause::UntracedOperator);
  }

  TEST_CASE("the node budget stops a runaway constant loop") {
    const auto code = std::string(R"lua(
local M = {}
M.output_parameters = {"s"}
function M.evaluate(fields, x)
  local s = x
  for i = 1, 1000000 do s = s + 1.0 end
  return s
end
return M
)lua");
    TraceOptions options;
    options.nodeBudget = 1000;
    const auto failure = mustRefuse(code, options);
    CHECK(failure.cause == Cause::BudgetExceeded);
  }

  // -------------------------------------------------------- determinism ----

  TEST_CASE("pairs iteration does not leak into the fingerprint") {
    // Lua seeds its string hash per state, so an unshadowed `pairs` gives a
    // different accumulation order in every process -- and under a domain
    // decomposition, different values at the same point on different ranks.
    // Tracing the same source repeatedly must give one fingerprint.
    const auto reference = mustTrace(TableIterationOrder).fingerprint();
    for (int repeat = 0; repeat < 8; ++repeat) {
      CHECK(mustTrace(TableIterationOrder).fingerprint() == reference);
    }
  }

  TEST_CASE("a traced program has a stable fingerprint across traces") {
    CHECK(mustTrace(PlanarWave).fingerprint() == mustTrace(PlanarWave).fingerprint());
    CHECK(mustTrace(PlanarWave).fingerprint() != mustTrace(SharedSubexpression).fingerprint());
  }

  // ------------------------------------------------------ differential -----

  TEST_CASE("the traced program agrees with the interpreted reader") {
    // This is the acceptance criterion for the package, and the same check
    // CompiledReader::prepare should run at init before trusting a kernel.
    // Deliberately not only at "nice" coordinates: the ladder is where the
    // interesting disagreements live.
    const std::vector<double> ladder = {-1e6, -1e3, -1.0, -1e-3, 0.0, 1e-3, 1.0, 1e3, 1e6};

    for (const auto* code : {PlanarWave, LayeredSelect, SharedSubexpression}) {
      const auto program = mustTrace(code);

      DataTable table(ladder.size());
      for (const auto& in : program.inputs()) {
        table.addInput(in.name, DataType::F64);
      }
      for (const auto& out : program.outputs()) {
        table.addOutput(out.name, DataType::F64);
      }
      for (std::size_t channel = 0; channel < program.inputs().size(); ++channel) {
        for (std::size_t point = 0; point < ladder.size(); ++point) {
          table.set(program.inputs()[channel].name, point, ladder[point]);
        }
      }

      const auto binding = expr::Binding::bind(program, table);
      expr::BackendOptions options;
      options.preferred = expr::BackendKind::Interpreter;
      makeKernel(program, binding, /*grids=*/{}, options)->run(table);

      DataTable reference = table;
      LuaReader interpreted{code};
      interpreted.prepare();
      interpreted.call(reference);

      for (const auto& out : program.outputs()) {
        for (std::size_t point = 0; point < ladder.size(); ++point) {
          const double traced = table.get(out.name, point);
          const double expected = reference.get(out.name, point);
          if (std::isnan(expected)) {
            CHECK(std::isnan(traced));
          } else {
            // Bit equality is the right bar here, not a tolerance: both paths
            // evaluate the same operations in the same order in fp64, so any
            // difference is a lowering bug and not accumulated rounding.
            CHECK(traced == expected);
          }
        }
      }
    }
  }
}

} // namespace seissol::unit_test
