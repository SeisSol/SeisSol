// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_CORE_H_
#define SEISSOL_SRC_EXPR_CORE_H_

// The core format of sderiv: the text form of a Program, as docs/sderiv-spec.rst specifies it.
// Every construct of the format maps to one construct of the program, so that a program prints and
// loads back without loss. It is what frontends emit and what SeisSol reads; the full language
// (SderivFrontend.h) is the syntax for writing by hand.
//
//   #!sderiv 1
//   meta generator = "easi 1.8.0"
//   in x, y, z
//   grid m (kind = "asagi", path = "vmodel.nc", variable = "data",
//           components = ("rho", "mu", "lambda"),
//           interpolation = "linear", boundary = "clamp")
//   def t0 = lt(z, -2000)
//   def t1 = m.rho(x, y, z)
//   def t2 = select(t0, t1, 2600 - 0.5 * z)
//   out rho = t2
//
// Beyond the draft of the specification, the format has the types of values of the program
// (ValueType), without which a program that converts a value would not print: the conversions
// f64(x), f32(x) and real(x) as operations, and `state f64 s = 0` for a state of a type of its own.
// And its grids are of the kinds of the reader (reader::datafield::gridAttributesOf()): besides
// ASAGI and SCEC files a table (`coordinates = ("x", "z")`) and an inline grid, whose values nest
// as its axes, the last outermost, down to the components of a node:
//
//   grid layers (kind = "inline",
//                components = ("rho", "mu"),
//                knots = ((-3000, -2000, -1000)),
//                values = ((3000, 40),
//                          (2700, 20),
//                          (2600, 10)),
//                interpolation = "previous", boundary = "clamp")
//
// A loader declares every name before it is used, and once: inputs, states, definitions and grids
// share one namespace, outputs have their own. There are no implicit inputs and no constants; a
// number is read with correct rounding, and no operation on numbers is evaluated, except negation.
// Every error names the file, the line and the column (CoreError).

#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"

#include <cstddef>
#include <string>

namespace seissol::expr {

/// The newest version of the core format this loader reads; it reads every older one.
inline constexpr int CoreVersion = 1;

/// Where in a core file: line and column count from 1, the column in characters.
struct CorePosition {
  std::string file;
  int line{0};
  int column{0};
};

/// An error in a core file, at its position. A SderivError, so that what catches the errors of the
/// full language catches these as well; its what() starts with "file:line:column:".
class CoreError : public SderivError {
  public:
  CoreError(const CorePosition& position, const std::string& what);
  [[nodiscard]] const CorePosition& where() const { return where_; }

  private:
  CorePosition where_;
};

/// The value of the spelling of a number of sderiv: decimal (`2`, `2.5`, `.5`, `2.`, `1e-3`) or
/// hexadecimal floating point (`0x1.8p+1`), with correct rounding and without the locale; beyond
/// the range of binary64 an infinity or zero. Throws std::invalid_argument for anything else.
[[nodiscard]] double numberOf(const std::string& spelling);

/// The position of the byte `offset` of `text`: its line and its column in characters.
[[nodiscard]] CorePosition
    positionOf(const std::string& text, std::size_t offset, const std::string& file);

/// Whether `text` is core text, i.e. starts with the header `#!sderiv`. A file without it is read
/// as the full language.
[[nodiscard]] bool isCoreText(const std::string& text);

/// The program of the core text `text`. `file` names it in errors, and a relative grid path is
/// relative to its directory. Throws CoreError.
[[nodiscard]] Program loadCore(const std::string& text, const std::string& file = "");

/// `program` in the normal form of the core format: one text per program, its inputs, grids, states
/// and outputs in the order of the program, which loadCore() reads back into the same program, of
/// the same fingerprint. Throws std::invalid_argument for a program the format cannot
/// express: one with contractions, which a consumer introduces when it binds a program, or one
/// the lowering has made (written-out lookups, whole numbers, truth values).
[[nodiscard]] std::string printCore(const Program& program);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_CORE_H_
