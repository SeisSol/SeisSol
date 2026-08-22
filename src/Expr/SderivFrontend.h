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
//   + program    := { definition | griddecl } ( expr | ) EOF
//   + definition := [ 'out' ] 'def' NAME [ '(' NAME { ',' NAME } ')' ] '=' expr
//   + griddecl   := 'grid' NAME '=' STR ',' STR [ ',' STR ] ',' STR { ',' STR }
//                    //     name     kind  file  [variable]  interp  components...
//   + comment    := '#' { any } NEWLINE
//
//     expr       := 'let' NAME '=' expr 'in' expr | add
//     add        := mul { ('+'|'-') mul }
//     mul        := unary { ('*'|'/') unary }
//     unary      := '-' unary | power
//     power      := primary [ '**' unary ]           // right-associative
//     primary    := NUM | NAME [ '(' args ')' ] | '(' expr ')'
//
// `out` IS A MODIFIER ON def, not a sibling keyword, and the grammar says so
// with one production rather than two that happen to look alike. An exported
// definition binds exactly like an ordinary one -- same symbol table, same
// visibility, same reuse -- so a later output can read an earlier one by name:
//
//     def vs     = m_vs(x, y, z)          # internal
//     out def rho    = m_rho(x, y, z)
//     out def mu     = rho * vs * vs      # reads the output above
//
// The alternative -- a separate `export rho, mu` list -- was rejected because it
// is a hand-maintained parallel list to the definitions, which is the same
// defect class as Lua's M.output_parameters against the actual return values,
// and the reason the tracer exists. A marker on the definition cannot drift.
//
// `out def NAME(params)` is refused: an output is a value per point, not a
// function, and there is no call site for it to take arguments from.
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
// GRID COMPONENTS ARE ORDINARY FUNCTIONS. `grid m = "asagi", "model.nc",
// "linear", "rho", "vp", "vs"` registers m_rho, m_vp and m_vs, each callable
// with one to six coordinates and each lowering to one Kind::Lookup. There is
// deliberately no `sample(field, ...)` special form: a special form would not
// compose -- no `let g = m in sample(g, ...)`, no grid through a def -- whereas
// a plain call does. The price of flattening with an underscore is that a grid
// named `m_rho` alongside a grid `m` with a component `rho` is unreadable even
// though it is not formally ambiguous, so that combination is rejected.
//
// NOT SUPPORTED HERE, on purpose: dx, the temporal reducers and Sample. Those
// are derived-output nodes and Program::validate rejects them on the scripting
// path anyway; the catalogs below simply do not contain them.

#include "Expr/Program.h"

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

// Convenience for the single-output case and for tests.
[[nodiscard]] Program compileSderiv(const std::string& source, const std::string& outputName);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_SDERIVFRONTEND_H_
