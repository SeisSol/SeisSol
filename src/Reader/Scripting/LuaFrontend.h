// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_LUAFRONTEND_H_
#define SEISSOL_SRC_READER_SCRIPTING_LUAFRONTEND_H_

// The Lua frontend: a script of Lua 5.4 (LuaSyntax.h), evaluated once for all points at once into
// a program (expr::Program), without a Lua VM. The program computes what Lua computes, bit for bit,
// with the math library easi and SeisSol's Lua reader load (Lua's own, with lmathx's functions
// over it).
//
// A value that depends on the point is a node of the program; everything else -- numbers the
// script fixes, strings, tables, functions -- the frontend computes itself, as Lua does. So:
//   - an `if` on a value of the point runs both branches, and merges what they leave with
//     `select`; an early `return` or `break` keeps the points it stops apart until the end of the
//     function or the loop, so that a chain of them becomes a chain of selects; `c and x or y`
//     becomes a select as well;
//   - a loop whose bounds the script fixes is unrolled, `for`, `while` and `repeat` alike, and a
//     generic `for` over a table whose keys it fixes;
//   - a function is inlined at its call; its closures, varargs and multiple results are Lua's;
//   - a number is an integer or a float of Lua, as the script makes it: `7 // 2` is 3, and the
//     zero of an integer has no sign.
// Refused, with the line: a loop whose end depends on the point, `goto`, recursion, a table index
// that depends on the point, a string that depends on the point, metatables, and the libraries
// beyond math (print and io.write do nothing). An error of Lua at some points becomes a check of
// the program; at every point, an error of the translation.
//
// Two conventions: the function `f(x)` of easi's !LuaMap, and the module of a `lua:` file:
//     M.evaluate(fields, <inputs...>) -> <outputs...>
// The inputs are the names of the parameters after `fields` (or M.input_parameters, in their
// order). The outputs are the returned values in the order of M.output_parameters, or one returned
// table that names them (all its names, or those of M.output_parameters). M.state = {name =
// initial, ...} declares what the program keeps from one call to the next: a parameter of the name
// of a state reads its value of the previous call, the returned value of that name is its next
// one. M.field_specs declares grids, {name, kind, file, variable, interpolation, parameters}, which
// `fields.<name>:sample(<coordinates>)` samples, one value per parameter. `ssol` holds the
// vocabulary of the earlier tracer (lt, le, gt, ge, eq, land, lor, lnot, select) and the builtins
// of sderiv that Lua spells differently (sign, mod, atan2, pow, pi, g).

#include "Expr/Program.h"
#include "Reader/Scripting/LuaSyntax.h"

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace seissol::reader::scripting::lua {

/// A script, or a part of one, that is Lua but outside the subset the frontend translates.
class LuaUnsupported : public LuaError {
  public:
  using LuaError::LuaError;
};

struct LuaOptions {
  /// math with lmathx's functions over Lua's own, as easi and SeisSol load it where they are built
  /// with it: floor, ceil and round of C, which give floats, and fmin, fmax, trunc, log2 and the
  /// like. Without, math is Lua's own: its floor and ceil give integers.
  bool lmathx{true};
  /// The most nodes the program may have, and the most steps (statements and iterations) the
  /// evaluation may take: a loop of absurd length fails with a diagnosis.
  std::size_t nodeBudget{std::size_t{1} << 22U};
  std::size_t stepBudget{std::size_t{1} << 24U};
};

/// Where a script is, for the messages of its checks: the file, the line of the file the first
/// line of the script is on, and how a position is written: `file@line` as easi writes it, or
/// `file:line` as Lua does.
struct LuaOrigin {
  std::string file;
  int firstLine{1};
  char separator{':'};
  [[nodiscard]] std::string at(int line) const;
};

/// A check of a translation: a condition that holds at every point at which Lua does not fail,
/// and the message where it does not.
struct LuaCheck {
  expr::NodeId holds{expr::NoNode};
  std::string message;
};

/// The values of an easi !LuaMap, and its checks, built into the program they were made for.
struct LuaMapTranslation {
  std::map<std::string, expr::NodeId> values;
  std::vector<LuaCheck> checks;
  /// The lines of the `return` statements that give the values.
  std::vector<int> lines;
};

/// The function `f` of an easi !LuaMap: the values of `returns` at the dimensions `in` (by name, in
/// easi's order, which is that of the names), for the points at which `reach` holds (NoNode: every
/// point). As easi has it: `x[1]`, `x["name"]` and `x.name` read a dimension, `f` returns a table
/// with the values. Throws LuaError (syntax, or a failure at every point) and LuaUnsupported.
[[nodiscard]] LuaMapTranslation
    translateLuaMap(expr::Program& program,
                    const std::string& code,
                    const std::vector<std::pair<std::string, expr::NodeId>>& in,
                    const std::set<std::string>& returns,
                    expr::NodeId reach,
                    const LuaOrigin& origin,
                    const LuaOptions& options = {});

/// The program of a Lua model module (a `lua:` file), as above. Throws LuaError and LuaUnsupported.
[[nodiscard]] expr::Program translateLuaModule(const std::string& code,
                                               const LuaOrigin& origin,
                                               const LuaOptions& options = {});

} // namespace seissol::reader::scripting::lua

#endif // SEISSOL_SRC_READER_SCRIPTING_LUAFRONTEND_H_
