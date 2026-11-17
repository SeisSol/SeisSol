// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_IMPALATOLUA_H_
#define SEISSOL_SRC_READER_SCRIPTING_IMPALATOLUA_H_

// The Impala of a !FunctionMap as Lua, with the meaning ImpalaJIT gives it (easi up to 1.4) and
// easi's conversion into Lua (convertImpalaToLua()) gives it again from its fix on: + - and * /
// left-associative, `a - b + c` as `(a - b) + c`; the inputs local copies, which an assignment
// changes; every other variable local to its function; the constants in binary64; the functions
// those of C's math library (min, max and abs as fmin, fmax and fabs). easi runs that Lua as it
// runs a !LuaMap, and the walker hands it to the Lua frontend (LuaFrontend.h) as it hands it a
// !LuaMap, so that Lua's rules decide the rest at one place. The Lua says what easi's says, though
// not in its words (easi's looks the functions up through a metatable, which the frontend does not
// read); the code of every entry is on the lines of the file it is on, so that a message of the
// frontend names the line of the file.
//
// easi 1.5 to 1.7 grouped the chains to the right (`a - (b + c)`), and took the inputs, the
// variables and the constants otherwise; impalaToLua() names the chains whose meaning that
// changed (ImpalaChain).

#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::reader::scripting {

/// An entry of the map of a !FunctionMap: the parameter, its Impala, and the line of the file its
/// text starts on. The lines of the text are those of the file where `literal` (a literal block
/// scalar, `|`); else the text is placed on its first line.
struct ImpalaEntry {
  std::string name;
  std::string code;
  int line{1};
  bool literal{false};
};

/// A chain whose meaning the grouping of easi 1.5 to 1.7 changed, as in `a - b + c`: the entry,
/// the line, the operator, and what follows it.
struct ImpalaChain {
  std::string name;
  int line{1};
  std::string op;
  std::string following;
};

/// The Lua of a !FunctionMap: a global function `f` of a table of the inputs that returns a table
/// of the parameters, as an !LuaMap has it; its first line is the line `firstLine` of the file.
struct ImpalaLua {
  std::string code;
  int firstLine{1};
  std::vector<ImpalaChain> chains;
};

/// Impala that easi's conversion does not read, at a line of the file.
class ImpalaError : public std::runtime_error {
  public:
  ImpalaError(int line, const std::string& what) : std::runtime_error(what), line_(line) {}
  [[nodiscard]] int line() const { return line_; }

  private:
  int line_;
};

/// The Lua of `entries` (in the order of the file) for the inputs `inputs`, of a name the last
/// entry, as easi reads them. Throws ImpalaError.
[[nodiscard]] ImpalaLua impalaToLua(const std::vector<ImpalaEntry>& entries,
                                    const std::set<std::string>& inputs);

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_IMPALATOLUA_H_
