// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_IMPALATOLUA_H_
#define SEISSOL_SRC_READER_SCRIPTING_IMPALATOLUA_H_

// easi's conversion of the Impala of a !FunctionMap into Lua (convertImpalaToLua(), easi >= 1.5),
// token for token. easi runs that Lua as it runs a !LuaMap, and the walker hands it to the Lua
// frontend (LuaFrontend.h) as it hands it a !LuaMap: what the Impala means -- the integers of Lua,
// the functions of Lua's math with lmathx, the variables, which are globals of Lua -- is decided
// at one place, by Lua's rules.
//
// The grammar is easi's, the grouping of its chains included: easi reads `a - b + c` as
// `a - (b + c)` and `a / b * c` as `a / (b * c)`, and so does the Lua it makes. The Lua differs
// from easi's in its line breaks only: the code of every entry is on the lines of the file it is
// on, so that a message of the frontend names the line of the file.

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

/// A chain whose meaning easi's grouping changes, as in `a - b + c`: the entry, the line, the
/// operator, and what follows it.
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

/// The Lua easi makes of `entries` (in the order of the file) for the inputs `inputs`. The
/// parameters are returned in the order easi returns them, that of its hash map. Throws
/// ImpalaError.
[[nodiscard]] ImpalaLua impalaToLua(const std::vector<ImpalaEntry>& entries,
                                    const std::set<std::string>& inputs);

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_IMPALATOLUA_H_
