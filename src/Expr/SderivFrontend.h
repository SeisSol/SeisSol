// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_SDERIVFRONTEND_H_
#define SEISSOL_SRC_EXPR_SDERIVFRONTEND_H_

// The sderiv expression language as a second frontend onto the shared IR,
// alongside the Lua tracer. Same Program out, so every backend, the binding
// layer and the kernel cache are shared.
//
// WHY A SECOND LANGUAGE. It has no control flow, and that is the feature. The
// Lua tracer needs four independent detection nets because the Lua VM will
// happily reduce a branch on a traced value to a constant without consulting
// anything we control; here there is no `if` to detect, no probe ladder, no
// truthiness, no table iteration order and no side effects. What is left is a
// grammar and a set of name checks. For a new model this is the language to
// reach for; Lua stays for the models that already exist and for the ones that
// genuinely need a general-purpose language.
//
// Two things also fall out for free. The input signature is exact with no
// declaration at all -- every unresolved name that is not a constant, builtin,
// def, let or grid component IS an input channel, so there is no counterpart to
// M.input_parameters and no lua_getlocal trick. And `select`, `lt`, `le`, `eq`,
// `land`, `lor`, `lnot` are ordinary builtins spelled exactly as on the Lua
// side, so a model can be transliterated between the two without renaming.
//
// GRAMMAR (deltas against the derived-output frontend marked +):
//
//   + program    := { definition | spread | griddecl | statedecl | check } ( expr | ) EOF
//   + definition := [ 'out' ] 'def' NAME [ '(' NAME { ',' NAME } ')' ] '=' expr
//   + spread     := 'out' expr          // the fields of a record, an output each
//   + griddecl   := 'grid' NAME '(' attribute { ',' attribute } ')'
//   + attribute  := NAME '=' ( STR | NUM | '(' value { ',' value } ')' )
//                    // those of the core format; an ASAGI file's variable and the boundary have
//                    // defaults ("z", "clamp")
//   + statedecl  := 'state' NAME '=' [ '-' ] NUM
//   + check      := 'check' expr ',' STR    // a condition every point meets (Program::checks())
//   + comment    := '#' { any } NEWLINE
//
//     expr       := 'let' NAME '=' expr 'in' expr | add
//     add        := mul { ('+'|'-') mul }
//     mul        := unary { ('*'|'/') unary }
//     unary      := '-' unary | power
//     power      := postfix [ '**' unary ]           // right-associative
//   + postfix    := primary { '.' NAME | '[' expr ']' }
//     primary    := NUM | NAME [ '(' args ')' ] | NAME '.' NAME '(' args ')' | '(' expr ')'
//                 | '{' NAME '=' expr { ',' NAME '=' expr } '}' | '[' expr { ',' expr } ']'
//
// RECORDS AND ARRAYS. A value is a number, a record of named values ({rho = 2600, mu = 1e10}),
// or an array of values ([x, y, z]); a matrix is an array of arrays. `p.rho` is a field, `a[0]` an
// element, its index a whole number the program fixes. They exist only here: the program has
// their numbers, each its own, so the core format and the backends see numbers only. Operators,
// the builtins and the conversions take them leaf by leaf: a number with every number of a
// record or an array, records of the same fields field by field (by name, in the order of the
// first), arrays of one length element by element -- so select(c, deep, shallow) of two records
// of a material selects every parameter, and a definition with parameters takes any of them, as
// its body is substituted at the call. `m(x, y, z)` samples every component of the grid m, a
// record of their names. `out e` gives the fields of the record e as outputs of their names. The
// intrinsics of arrays (Intrinsics): sum, min and max of one array; dot, norm and cross of
// vectors; matmul of matrices and vectors, transpose, trace and det (to 3 by 3) of matrices;
// eye(n), len(a) and range(n). A definition may take the name of an intrinsic, which it hides.
//
// `out` IS A MODIFIER ON def, not a sibling keyword, and the grammar says so
// with one production rather than two that happen to look alike. An exported
// definition binds exactly like an ordinary one -- same symbol table, same
// visibility, same reuse -- so a later output can read an earlier one by name:
//
//     def vs     = m.vs(x, y, z)          # internal
//     out def rho    = m.rho(x, y, z)
//     out def mu     = rho * vs * vs      # reads the output above
//
// The alternative -- a separate `export rho, mu` list -- was rejected because it
// is a hand-maintained parallel list to the definitions, which is the same
// defect class as Lua's M.output_parameters against the actual return values.
// A marker on the definition cannot drift.
//
// `out def NAME(params)` is refused: an output is a value per point, not a
// function, and there is no call site for it to take arguments from.
//
// STATE is a value per point that survives from one call to the next
// (Program::state()). `state NAME = initial` declares it, and the definition of
// the same name gives its next value; every reference to NAME -- including the
// one in its own definition -- reads the value of the previous call, so all
// states update as one parallel assignment. Made an output with `out def`, the
// same name is state AND output: the running value is written along with every
// call, and the consumer reads it at its own cadence. A maximum over time:
//
//     state pgv = 0.0
//     out def pgv = max(pgv, sqrt(v1*v1 + v2*v2 + v3*v3))
//
// and one that restarts whenever a window counter advances:
//
//     state last = -1.0
//     state peak = 0.0
//     def fresh = lt(last, window)
//     out def peak = select(fresh, v, max(peak, v))
//     def last = window
//
// INPUTS A CONSUMER NAMES (SderivOptions::inputs) read the input of that name wherever they are
// referenced, as a state reads its previous value, also where a definition of the same name is in
// scope; that definition only gives the output of the name. So a module may map a state onto the
// next one name by name, the inner state of a boundary onto its ghost state, say, as one parallel
// assignment:
//
//     out def v1 = -v1
//     out def s_xy = s_xy + 0.1 * v1      # v1 is the input here too
//
// Otherwise a definition that reads itself, directly or through others, is an error.
//
// A module is EITHER a trailing expression, whose name the caller supplies via
// compileSderiv(source, name), OR one or more `out def`s -- never both, because
// then it would be open which one is "the" output.
//
// STRINGS ARE CONFINED TO griddecl. They are lexed everywhere but primary() has
// no string case, so a string in an expression is a parse error that names the
// restriction. That is what keeps `asagi`, `linear` and the component names out
// of the keyword set: `grid` and `out` are the only words this frontend adds.
//
// A GRID COMPONENT IS SAMPLED AS A CALL. `grid m (kind = "asagi", path = "model.nc", components =
// ("rho", "vp", "vs"), interpolation = "linear")` makes m.rho, m.vp and m.vs callable with one to
// six coordinates, each call one Kind::Lookup, as in the core format. There is deliberately no
// `sample(field, ...)` special form: a special form would not compose -- no `let g = m in
// sample(g, ...)`, no grid through a def -- whereas a plain call does. A dotted name is no name a
// declaration can have, so a component collides with nothing; the grid itself shares the namespace
// of the definitions.

#include "Expr/Program.h"

#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr {

// Carries a byte offset into the source so the caller can point at the right
// place in a parameter file rather than saying "syntax error".
class SderivError : public std::runtime_error {
  public:
  SderivError(const std::string& stage, const std::string& what, int position)
      : std::runtime_error(stage + ": " + what), position_(position) {}
  [[nodiscard]] int position() const { return position_; }

  private:
  int position_;
};

// One source per output. The caller supplies the output names because the
// language has no place to put them: a source is one expression, and the
// consumer's DataTable already knows which columns it wants filled.
struct SderivOutput {
  std::string name;
  std::string source;
  reader::scripting::DataType type{reader::scripting::DataType::F64};
};

// Compile one or more sources into a single Program. All outputs share one
// Arena, so a subexpression common to two outputs is interned once and the
// backends see the sharing.
//
// Throws SderivError on any lexing, parsing, name-resolution or grid-check
// failure. Unlike the Lua path there is no fallback to fall back TO, so a
// failure here is a hard configuration error and the caller should logError.
[[nodiscard]] Program compileSderiv(const std::vector<SderivOutput>& outputs);

/// Compile ONE source that names its own outputs with `out def`.
///
/// The form a `.sderiv` file on disk takes, and the reason it exists: the
/// overload above takes one source PER OUTPUT, so a three-quantity material
/// model would repeat its grid declarations and helper definitions three times.
/// Here the prologue is written once and every output can read the ones above
/// it by name.
///
/// Throws SderivError if the source has no `out def`, or if it ends in a bare
/// expression instead.
[[nodiscard]] Program compileSderivModule(const std::string& source);

/// What a consumer tells the compiler of a module beyond its source.
struct SderivOptions {
  /// Names that read the input of that name wherever they are referenced (see above).
  std::set<std::string> inputs;
  /// The file the source comes from, if any: a source in the core format (Core.h) names it in its
  /// errors, and its relative grid paths are relative to its directory.
  std::string file;
};

/// The same, with `options`. A source that starts with the header `#!sderiv` is in the core format
/// (Core.h), which declares its inputs itself; `options.inputs` does not apply to it.
[[nodiscard]] Program compileSderivModule(const std::string& source, const SderivOptions& options);

// Convenience for the single-output case and for tests.
[[nodiscard]] Program compileSderiv(const std::string& source, const std::string& outputName);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_SDERIVFRONTEND_H_
