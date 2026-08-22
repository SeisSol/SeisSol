// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_LUATRACER_H_
#define SEISSOL_SRC_READER_SCRIPTING_LUATRACER_H_

// Traces a Lua model module into a seissol::expr::Program by running
// `M.evaluate` ONCE with symbolic values that build IR nodes instead of
// computing. No parser: the Lua VM does the parsing, the metatable does the
// rest. Constant loops unroll for free, and common subexpressions collapse in
// the Arena's interning rather than in a pass.
//
// The module convention is unchanged from LuaReader, so the same script feeds
// the traced and the interpreted path:
//     M.evaluate(fields, <inputs...>) -> <outputs...>
//     M.output_parameters, M.field_specs, M.version, M.source_file
//
// WHAT CHANGES: M.input_parameters is no longer read. The input signature comes
// from the PARAMETER NAMES of M.evaluate, via lua_getlocal(L, nullptr, i) on the
// function object. Measured behaviour (5.4.6), which is why this needs three
// guards rather than a bare read:
//   - locals are NOT leaked; exactly the parameters come back
//   - `function M:evaluate(...)` inserts `self` as parameter 1, shifting every
//     name by one -> require names[0] == "fields"
//   - a variadic `evaluate(...)` reports ZERO parameters, which would trace as
//     a constant program with no inputs -> reject isvararg outright
//   - a stripped chunk yields no names while lua_getinfo still reports the
//     count -> require names.size() == nparams. luaL_loadbufferx(mode="t")
//     already refuses binary chunks, so this is belt and braces.
// That removes one of the two hand-maintained declarations whose
// runtime string matching has repeatedly resolved to the wrong column. The
// other one, M.output_parameters, CANNOT be inferred: Lua returns values
// positionally and there is no name to recover. It stays, but the tracer knows
// the number of returned values and so a count mismatch is now a load-time
// error instead of a silent misbinding.
//
// ==== WHY TRACING NEEDS FOUR NETS, NOT ONE ==================================
//
// A traced value is a userdata. Lua's VM will happily reduce control flow on
// one to a constant without ever consulting our metatable, so the interesting
// work is detecting that rather than performing it. Four independent nets, in
// increasing cost and decreasing precision:
//
//   1. METAMETHOD REFUSAL. __lt/__le/__eq exist but the VM casts their result
//      to a boolean, so `if x > 0` would pick one branch silently. They are
//      wired to error() instead. Conditions are built explicitly with
//      ssol.lt/le/gt/ge/eq/land/lor/lnot and consumed by ssol.select.
//      Also refused, because they have no symbolic meaning and their natural
//      failure mode is confusing: // .. # & | ~ << >> and indexing/calling.
//      LIMIT: __eq is only dispatched when BOTH operands are full userdata and
//      they are not primitively equal. `x == 0` against a number never reaches
//      the metatable — it is silently false. Net 3 covers that.
//
//   2. UNCONSUMED-CONDITION CHECK. A node built by ssol.lt/le/... that is never
//      passed to select/land/lor/lnot and is not an output root was almost
//      certainly used as a raw `if` condition — which the symbolic run cannot
//      observe, since every userdata is truthy. Deterministic and cheap, and it
//      is the only net that catches a threshold outside the probe ladder.
//
//   3. PROBE RUNS. The module is re-loaded with the REAL math library and
//      concrete numbers, under lua_sethook(LUA_MASKLINE), once per rung of a
//      geometric ladder (±1e9 … ±1e-3, 0). Two runs whose visited-line
//      sequences differ prove a data-dependent branch, whatever mechanism
//      produced it. The ladder, not the count, is what matters: an earlier
//      version used four "interesting" vectors and missed `if z > -1000`
//      because every draw sat above the threshold. LIMIT: a threshold beyond
//      the ends of the ladder is not straddled and is not caught here.
//
//   4. DIFFERENTIAL CHECK (CompiledReader::prepare, Package 4). The traced
//      program is evaluated against the interpreted LuaReader on a sample of
//      the REAL point set. This is the only net that uses the actual data
//      distribution, and it is what should gate the fallback to the
//      interpreted reader. Probabilistic, and deliberately the last line
//      rather than the first.
//
// ==== DETERMINISM ===========================================================
//
// Lua seeds its string hash per state from an address and the clock, so `pairs`
// visits a table in a different order in every process. A script that
// accumulates over a table then produces a different summation order — hence a
// different fingerprint AND different rounding — on every rank, which under a
// domain decomposition means two ranks computing different values at the same
// point. `pairs` is therefore replaced by a sorted-key version in the trace
// environment (numbers by value first, then strings). os, io, debug, require,
// load and math.random are absent from the environment for the same reason: a
// traced program must not be able to bake in a per-run value.
//
// ==== SEMANTICS PINNED HERE =================================================
//
//   `%`         -> Fn::Mod, defined as fmod(a,b) pulled onto the sign of b.
//                  NOT a - floor(a/b)*b: over 121 sign/magnitude pairs the
//                  naive form disagrees with Lua's `%` on 14 of them, which
//                  would show up as an unexplained interpreter/kernel drift.
//   math.fmod   -> refused. It is C truncation semantics, a different function
//                  from `%`, and has no IR op.
//   expm1, log1p, hypot, cbrt, erfc, gamma, copysign, …
//               -> refused rather than lowered. exp(x)-1 and log(1+x) lose
//                  exactly the precision those functions exist for, so a
//                  lowered version would disagree with the interpreted reader
//                  precisely where it matters.
//   conditions  -> 0.0/1.0 in the compute type, so land/lor are Min/Max and
//                  lnot is 1-c, all branchless. Select tests c != 0.
//                  Fn::And/Or must be Min/Max in the interpreter and every
//                  backend for this to hold; that is an assumption on Ir.h,
//                  not something this header can enforce.
//
// ==== NOT DONE HERE =========================================================
//
// Two-pass branch enumeration (record the comparison, return a concrete
// boolean, re-trace with the opposite branch forced, merge with select) is the
// additive next step and does not break scripts written against this
// interface. It subsumes what an AST rewrite would buy — and Lua 5.4 exposes no
// AST to rewrite, since lparser.c emits bytecode in one pass — but it can only
// force branches that go through our comparison functions, so nets 2 and 3
// remain necessary either way.

#include "Expr/Program.h"
#include "Reader/Scripting/LuaReader.h"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace seissol::reader::scripting {

// Grids come from M.field_specs. `kind` selects the loader ("asagi", "scec",
// …) but every kind is a regular grid as far as the IR is concerned: one
// Kind::Lookup per requested component. NOTE that this makes FieldSpec-
// ::parameters load-bearing for every kind, including SCEC, where LuaReader.h
// currently documents it as ignored: the tracer must know the component count
// before it can emit the Lookup nodes. Making it required turns it into a check
// the loader can validate its own component list against, rather than metadata
// that is silently unused.

struct TraceOptions {
  // Ceiling on Arena size, so a constant loop with an absurd trip count fails
  // with a diagnosis instead of exhausting memory.
  std::size_t nodeBudget{200000};
  // Run the probe pass. Off only for tests that need to inspect a DAG the
  // probes would reject.
  bool probe{true};
  // Script parameters folded into the DAG as Const. See the note below.
  std::vector<std::pair<std::string, double>> parameters;
};

// Why a diagnostic type rather than logError: a failed trace is NOT an error at
// this level. ReaderBuilder is expected to fall back to the interpreted reader
// and logWarning with `reason` plus the note that the interpreted path is two
// to three orders of magnitude slower per point. logError aborts, and aborting
// on an untraceable-but-valid script would be a regression against today's
// behaviour.
struct TraceFailure {
  enum class Cause : std::uint8_t {
    LoadError,         // the chunk does not compile or does not return a table
    NoEvaluate,        // no M.evaluate, or its first parameter is not `fields`
    UntracedOperator,  // a refused metamethod or an untraceable math function
    RawCondition,      // net 2: a condition that never reached select
    DataDependentFlow, // net 3: probe line traces diverged
    SideEffect,        // evaluate mutated module state during the trace
    SignatureMismatch, // output_parameters vs. the number of returned values
    BudgetExceeded
  };
  Cause cause{Cause::LoadError};
  std::string reason; // already formatted for logWarning, includes the line
  int line{-1};       // in the user's script; the lmathx prelude is offset out
};

// Trace `code` into a Program. Returns std::nullopt on any of the causes above
// and fills `failure`; never throws for a script-level problem.
[[nodiscard]] std::optional<expr::Program>
    traceLuaModule(const std::string& code, const TraceOptions& options, TraceFailure& failure);

// --- script parameters -------------------------------------------------------
// OPEN DECISION FROM THE HANDOVER, resolved as follows: the choice does not
// belong at trace time. Tracing is the expensive and fragile step and should
// happen once per script, not once per parametrisation. So parameters named in
// TraceOptions are folded to Const during the trace (Feature B: the built-in
// ICs, constant per run, few variants — a specialised kernel each), and a
// script that wants one kernel across parametrisations simply declares them as
// ordinary inputs and leaves TraceOptions::parameters empty.
//
// Nothing in Ir.h changes for this: an unfolded parameter is a Kind::Field on a
// channel whose Binding is broadcast rather than per-point. That does require
// Binding to grow a stride-0 column; it is the only place the two forms differ,
// and it is a smaller change than a new node kind.

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_LUATRACER_H_
