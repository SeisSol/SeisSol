// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// The Lua frontend (Reader/Scripting/LuaFrontend.h): the parser, the programs it makes and their
// agreement with Lua bit for bit -- the interpreted reader runs Lua 5.4 with lmathx, as easi and
// SeisSol do --, its checks where Lua fails at some points, and what it refuses, with the line.

#include <doctest.h>

#include "Expr/Interp.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Scripting/CompiledReader.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/LuaFrontend.h"
#include "Reader/Scripting/LuaSyntax.h"
#include "TestHelper.h"

#ifdef USE_LUA
#include "Reader/Scripting/LuaReader.h"
#endif

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace luafrontendtest {

namespace sc = reader::scripting;
namespace lua = reader::scripting::lua;

inline expr::Program translated(const std::string& code, const lua::LuaOptions& options = {}) {
  return lua::translateLuaModule(code, lua::LuaOrigin{"model.lua", 1, ':'}, options);
}

/// The program at the points of `columns` (its inputs, by name), with the interpreter: the
/// outputs, [output][point], and the points that fail a check, with the check.
struct Evaluation {
  std::vector<std::vector<double>> outputs;
  std::vector<std::pair<std::int32_t, std::size_t>> failures;
  [[nodiscard]] bool fails(std::size_t point) const {
    for (const auto& [check, at] : failures) {
      if (at == point) {
        return true;
      }
    }
    return false;
  }
};

inline Evaluation evaluate(const expr::Program& program,
                           const std::map<std::string, std::vector<double>>& columns,
                           std::size_t count) {
  Evaluation evaluation;
  evaluation.outputs.assign(program.outputs().size(), std::vector<double>(count, -1.0));
  class Io final : public expr::TileIo<double> {
public:
    Io(const expr::Program& program,
       const std::map<std::string, std::vector<double>>& columns,
       Evaluation& evaluation)
        : program_(program), columns_(columns), evaluation_(evaluation) {}
    void gather(std::size_t first, std::size_t count, double* dst) const override {
      for (std::size_t i = 0; i < program_.inputs().size(); ++i) {
        const auto& column = columns_.at(program_.inputs()[i].name);
        for (std::size_t l = 0; l < count; ++l) {
          dst[i * count + l] = column[first + l];
        }
      }
    }
    void scatter(std::size_t first, std::size_t count, const double* src) override {
      for (std::size_t o = 0; o < evaluation_.outputs.size(); ++o) {
        for (std::size_t l = 0; l < count; ++l) {
          evaluation_.outputs[o][first + l] = src[o * count + l];
        }
      }
    }

private:
    const expr::Program& program_;
    const std::map<std::string, std::vector<double>>& columns_;
    Evaluation& evaluation_;
  };
  class Arrays final : public expr::GridArrays {
public:
    explicit Arrays(Evaluation& evaluation) : evaluation_(evaluation) {}
    [[nodiscard]] const reader::datafield::ArrayView& array(expr::GridId /*grid*/) const override {
      throw std::logic_error("no grids here");
    }
    void report(expr::GridId /*grid*/, expr::CheckKind /*kind*/) const override {}
    void failed(std::int32_t check,
                std::size_t point,
                const std::vector<double>& /*inputs*/) const override {
      evaluation_.failures.emplace_back(check, point);
    }

private:
    Evaluation& evaluation_;
  };
  const expr::LoweredProgram lowered = expr::lower(program);
  Io io(program, columns, evaluation);
  const Arrays arrays(evaluation);
  expr::InterpreterOptions options;
  options.tileSize = 1;
  expr::TileInterpreter<double> interpreter(program, lowered, &arrays, options);
  std::vector<double> persistent(
      std::max<std::size_t>(static_cast<std::size_t>(lowered.persistentSlotCount()), 1), 0.0);
  interpreter.run(io, count, persistent.data());
  return evaluation;
}

/// The values of a ladder that crosses the thresholds of the corpus, both zeros, and some others.
inline std::vector<double> ladder() {
  return {-1e6, -4500.0, -1000.0, -999.5, -2.5, -1.0, -0.0, 0.0, 1e-3, 0.75,
          1.0,  2.5,     3.0,     10.0,   1e3,  1e6,  -7.0, 6.0, -0.5, 0.5};
}

/// Columns of `count` points for the inputs `names`: the ladder, turned for each.
inline std::map<std::string, std::vector<double>>
    columnsFor(const std::vector<std::string>& names) {
  const auto values = ladder();
  std::map<std::string, std::vector<double>> columns;
  for (std::size_t i = 0; i < names.size(); ++i) {
    std::vector<double> column(values.size() * 3);
    for (std::size_t p = 0; p < column.size(); ++p) {
      column[p] = values[(p * (i + 1) + i * 7 + p / values.size()) % values.size()];
    }
    columns[names[i]] = std::move(column);
  }
  return columns;
}

inline bool sameValue(double a, double b) {
  return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof(double)) == 0;
}

#ifdef USE_LUA
/// Where the program of `code` and Lua first give other bits, at the ladder, over `calls` calls
/// (for states); empty where they agree.
inline std::string disagreement(const std::string& code, int calls = 1) {
  const expr::Program program = translated(code);
  std::vector<std::string> inputs;
  for (const auto& input : program.inputs()) {
    inputs.push_back(input.name);
  }
  const auto columns = columnsFor(inputs);
  const std::size_t count = ladder().size() * 3;
  std::vector<std::vector<double>> compiled(program.outputs().size(), std::vector<double>(count));
  std::vector<std::vector<double>> interpreted(program.outputs().size(),
                                               std::vector<double>(count));
  sc::DataTable compiledTable(count);
  sc::DataTable referenceTable(count);
  for (const auto& [name, column] : columns) {
    compiledTable.bindViewConst<double>(name, sc::Direction::In, column.data());
    referenceTable.bindViewConst<double>(name, sc::Direction::In, column.data());
  }
  for (std::size_t o = 0; o < program.outputs().size(); ++o) {
    compiledTable.bindView<double>(
        program.outputs()[o].name, sc::Direction::Out, compiled[o].data());
    referenceTable.bindView<double>(
        program.outputs()[o].name, sc::Direction::Out, interpreted[o].data());
  }
  sc::CompiledReaderOptions options;
  options.backend.preferred = expr::BackendKind::Interpreter;
  options.differentialSamples = 0;
  sc::CompiledReader reader(program, nullptr, options);
  sc::LuaReader lua(code);
  reader.prepare(compiledTable);
  lua.prepare(referenceTable);
  for (int call = 0; call < calls; ++call) {
    reader.call(compiledTable);
    lua.call(referenceTable);
    for (std::size_t o = 0; o < program.outputs().size(); ++o) {
      for (std::size_t p = 0; p < count; ++p) {
        if (sameValue(compiled[o][p], interpreted[o][p])) {
          continue;
        }
        std::ostringstream text;
        text << std::hexfloat << "call " << call << ", " << program.outputs()[o].name << " at";
        for (const auto& [name, column] : columns) {
          text << " " << name << " = " << column[p];
        }
        text << ": " << compiled[o][p] << " for Lua's " << interpreted[o][p];
        return text.str();
      }
    }
  }
  return "";
}

/// That the program of `code` gives Lua's bits.
inline void agreesWithLua(const std::string& code, int calls = 1) {
  const std::string found = disagreement(code, calls);
  INFO(found);
  CHECK(found.empty());
}

/// Random modules of the subset, for SEISSOL_LUA_FUZZ: branches on the point and early returns,
/// loops with breaks (for, ipairs, while, repeat), functions and closures over upvalues, tables
/// written under conditions, <const> locals, integers and floats of Lua, lmathx -- and nothing that
/// fails in Lua, which would stop the interpreted reader.
class RandomModule {
  public:
  explicit RandomModule(std::uint64_t seed) : random_(seed) {}

  std::string make() {
    std::ostringstream out;
    out << "local M = {}\nM.output_parameters = {'a', 'b', 'c'}\n";
    out << "local k = {" << literal() << ", " << literal() << ", " << literal() << "}\n";
    scopes_ = {Scope{{"k[1]", "k[2]", "k[3]"}, {"k"}, {}}};
    out << "local function h(p, q)\n";
    scopes_.push_back(Scope{{"p", "q"}, {}, {}});
    arity_ = 1;
    block(out, 2, 2, false);
    out << "  return " << expression(2) << "\nend\n";
    scopes_.pop_back();
    out << "function M.evaluate(fields, x, y, z)\n";
    scopes_.push_back(Scope{{"x", "y", "z"}, {}, {}});
    canCall_ = true;
    arity_ = 3;
    block(out, 2, 3, false);
    out << "  return " << expression(2) << ", " << expression(2) << ", " << expression(2)
        << "\nend\nreturn M\n";
    return out.str();
  }

  private:
  struct Scope {
    std::vector<std::string> numbers;
    std::vector<std::string> tables;
    std::vector<std::string> functions;
  };

  int pick(int n) { return std::uniform_int_distribution<int>(0, n - 1)(random_); }

  template <typename T>
  const T& among(const std::vector<T>& items) {
    return items[static_cast<std::size_t>(pick(static_cast<int>(items.size())))];
  }

  std::vector<std::string> visible(std::vector<std::string> Scope::* kind) const {
    std::vector<std::string> out;
    for (const auto& scope : scopes_) {
      out.insert(out.end(), (scope.*kind).begin(), (scope.*kind).end());
    }
    return out;
  }

  std::string literal() {
    static const std::vector<std::string> Literals = {
        "0", "1", "2", "-3", "7", "0.0", "-0.0", "0.5", "2.5", "1e3", "0.1", "-1.5", "1000", "3"};
    return among(Literals);
  }

  std::string atom() {
    if (pick(4) == 0) {
      return literal();
    }
    return among(visible(&Scope::numbers));
  }

  /// A divisor of `//` and `%`: none is an integer zero, which fails in Lua.
  std::string divisor(int depth) {
    switch (pick(4)) {
    case 0:
      return "2";
    case 1:
      return "-3";
    case 2:
      return "1.5";
    default:
      return "(" + expression(depth - 1) + " + 0.5)";
    }
  }

  std::string expression(int depth) {
    if (depth <= 0) {
      return atom();
    }
    static const std::vector<std::string> Functions = {"math.floor",
                                                       "math.ceil",
                                                       "math.abs",
                                                       "math.sin",
                                                       "math.cos",
                                                       "math.trunc",
                                                       "math.round",
                                                       "math.tanh",
                                                       "math.deg",
                                                       "math.rad",
                                                       "math.exp"};
    switch (pick(16)) {
    case 0:
    case 1:
    case 2:
    case 3: {
      static const std::vector<std::string> Ops = {"+", "-", "*", "/", "+", "*", "-"};
      return "(" + expression(depth - 1) + " " + among(Ops) + " " + expression(depth - 1) + ")";
    }
    case 4:
      return "(" + expression(depth - 1) + " // " + divisor(depth) + ")";
    case 5:
      return "(" + expression(depth - 1) + " % " + divisor(depth) + ")";
    case 6:
      return "(" + expression(depth - 1) + " ^ " + (pick(2) == 0 ? "2" : expression(depth - 1)) +
             ")";
    case 7:
      return "(- " + expression(depth - 1) + ")";
    case 8: {
      const std::string& f = among(Functions);
      return f + "(" +
             (f == "math.exp" ? "(" + expression(depth - 1) + ") / 1e3" : expression(depth - 1)) +
             ")";
    }
    case 9:
      return "(" + condition(depth - 1) + " and " + expression(depth - 1) + " or " +
             expression(depth - 1) + ")";
    case 10:
      return std::string(pick(2) == 0 ? "math.max(" : "math.min(") + expression(depth - 1) + ", " +
             expression(depth - 1) + (pick(2) == 0 ? ", " + literal() : std::string()) + ")";
    case 11:
      return std::string(pick(2) == 0 ? "math.fmax(" : "math.fmin(") + expression(depth - 1) +
             ", " + expression(depth - 1) + ")";
    case 12:
      if (canCall_) {
        return "h(" + expression(depth - 1) + ", " + expression(depth - 1) + ")";
      }
      return "math.atan(" + expression(depth - 1) + ", " + expression(depth - 1) + ")";
    case 13: {
      const auto functions = visible(&Scope::functions);
      if (!functions.empty()) {
        return among(functions) + "(" + expression(depth - 1) + ")";
      }
      return "(" + expression(depth - 1) + " - " + (pick(2) == 0 ? "0" : "-2") + ")";
    }
    case 14:
      return "(" + expression(depth - 1) + " - " +
             among(std::vector<std::string>{"0", "(1 - 1)", "-0", "3", "0.0"}) + ")";
    default:
      return "math.log(math.abs(" + expression(depth - 1) + ") + 1" +
             (pick(2) == 0 ? std::string(", ") + (pick(2) == 0 ? "10" : "2") : std::string()) + ")";
    }
  }

  std::string condition(int depth) {
    switch (pick(5)) {
    case 0:
      return "(" + condition(depth - 1 < 0 ? 0 : depth - 1) + " and " + comparison(depth) + ")";
    case 1:
      return "(" + comparison(depth) + " or " + comparison(depth) + ")";
    case 2:
      return "not " + comparison(depth);
    default:
      return comparison(depth);
    }
  }

  std::string comparison(int depth) {
    static const std::vector<std::string> Comparisons = {"<", "<=", ">", ">=", "==", "~="};
    return "(" + expression(depth > 0 ? depth - 1 : 0) + " " + among(Comparisons) + " " +
           expression(depth > 0 ? depth - 1 : 0) + ")";
  }

  static void indent(std::ostream& out, int level) {
    for (int i = 0; i < level; ++i) {
      out << "  ";
    }
  }

  void block(std::ostream& out, int level, int depth, bool inLoop) {
    const int count = 1 + pick(4);
    for (int i = 0; i < count; ++i) {
      statement(out, level, depth, inLoop);
    }
  }

  /// A block of its own scope.
  void inner(std::ostream& out, int level, int depth, bool inLoop, Scope scope = {}) {
    scopes_.push_back(std::move(scope));
    block(out, level, depth, inLoop);
    scopes_.pop_back();
  }

  void statement(std::ostream& out, int level, int depth, bool inLoop) {
    const int kind = depth <= 0 ? pick(2) : pick(13);
    indent(out, level);
    const std::string n = std::to_string(next_++);
    switch (kind) {
    case 0:
      out << "local v" << n << " = " << expression(2) << "\n";
      scopes_.back().numbers.push_back("v" + n);
      return;
    case 1: {
      std::vector<std::string> assignable;
      for (const auto& name : visible(&Scope::numbers)) {
        if (name.front() == 'v' || name.front() == 't') {
          assignable.push_back(name);
        }
      }
      if (assignable.empty()) {
        out << "local v" << n << " = " << expression(1) << "\n";
        scopes_.back().numbers.push_back("v" + n);
        return;
      }
      out << among(assignable) << " = " << expression(2) << "\n";
      return;
    }
    case 2:
    case 3:
      out << "if " << condition(1) << " then\n";
      inner(out, level + 1, depth - 1, inLoop);
      if (pick(2) == 0) {
        indent(out, level);
        out << "elseif " << condition(1) << " then\n";
        inner(out, level + 1, depth - 1, inLoop);
      }
      if (pick(2) == 0) {
        indent(out, level);
        out << "else\n";
        inner(out, level + 1, depth - 1, inLoop);
      }
      indent(out, level);
      out << "end\n";
      return;
    case 4:
      out << "for i" << n << " = " << (pick(2) == 0 ? "1, 3" : "4, 1, -1.5") << " do\n";
      inner(out, level + 1, depth - 1, true, Scope{{"i" + n}, {}, {}});
      indent(out, level);
      out << "end\n";
      return;
    case 5:
      if (inLoop) {
        out << "if " << condition(1) << " then break end\n";
        return;
      }
      [[fallthrough]];
    case 6:
      out << "if " << condition(1) << " then return " << expression(1);
      for (int i = 1; i < arity_; ++i) {
        out << ", " << expression(1);
      }
      out << " end\n";
      return;
    case 7:
      out << "local t" << n << " = {" << expression(1) << ", " << expression(1) << ", "
          << expression(1) << "}\n";
      scopes_.back().tables.push_back("t" + n);
      for (const char* i : {"[1]", "[2]", "[3]"}) {
        scopes_.back().numbers.push_back("t" + n + i);
      }
      return;
    case 8: {
      const auto tables = visible(&Scope::tables);
      out << "for _, e" << n << " in ipairs(" << among(tables) << ") do\n";
      inner(out, level + 1, depth - 1, true, Scope{{"e" + n}, {}, {}});
      indent(out, level);
      out << "end\n";
      return;
    }
    case 9:
      out << "local c" << n << " <const> = "
          << among(std::vector<std::string>{"0", "2", "-1", "1 - 1", "0.5", "3 // 2"}) << "\n";
      scopes_.back().numbers.push_back("c" + n);
      return;
    case 10:
      out << "local w" << n << " = 0\n";
      indent(out, level);
      out << "while w" << n << " < 3 do\n";
      indent(out, level + 1);
      out << "w" << n << " = w" << n << " + 1\n";
      inner(out, level + 1, depth - 1, true);
      indent(out, level);
      out << "end\n";
      scopes_.back().numbers.push_back("w" + n);
      return;
    case 11:
      out << "local r" << n << " = 0\n";
      indent(out, level);
      out << "repeat\n";
      indent(out, level + 1);
      out << "r" << n << " = r" << n << " + 1\n";
      inner(out, level + 1, depth - 1, true);
      indent(out, level);
      out << "until r" << n << " >= 2\n";
      scopes_.back().numbers.push_back("r" + n);
      return;
    default: {
      // a closure over the locals so far, which later assignments change
      out << "local function f" << n << "(a" << n << ")\n";
      const int arity = arity_;
      arity_ = 1;
      inner(out, level + 1, depth - 1, false, Scope{{"a" + n}, {}, {}});
      arity_ = arity;
      indent(out, level + 1);
      scopes_.push_back(Scope{{"a" + n}, {}, {}});
      out << "return " << expression(2) << "\n";
      scopes_.pop_back();
      indent(out, level);
      out << "end\n";
      scopes_.back().functions.push_back("f" + n);
      return;
    }
    }
  }

  std::mt19937_64 random_;
  std::vector<Scope> scopes_;
  int next_{0};
  int arity_{1};
  bool canCall_{false};
};
#endif

/// The refusal of `code`, which has to be one.
inline lua::LuaUnsupported refusalOf(const std::string& code) {
  try {
    static_cast<void>(translated(code));
  } catch (const lua::LuaUnsupported& error) {
    return error;
  }
  FAIL("the script was translated");
  return lua::LuaUnsupported(0, "");
}

// --------------------------------------------------------------- corpus ---

constexpr auto Branches = R"lua(
local M = {}
M.output_parameters = {"rho", "mu"}
function M.evaluate(fields, x, y, z)
  if z > -1000 then
    return 2200 + 0.1 * z, 1e9
  elseif z > -5000 then
    if x < 0 then
      return 2500, 2e9
    end
    return 2600, 3e9 + y
  end
  return 2700, 4e9
end
return M
)lua";

constexpr auto LayerSearch = R"lua(
local M = {}
M.output_parameters = {"vs", "layer", "below"}
local tops = {10.0, 2.5, -1.0, -999.5}
local speeds = {1000, 2000.5, 3000, 4500}
function M.evaluate(fields, x, y, z)
  local vs, layer = 5000, #tops + 1
  for i = 1, #tops do
    if z > tops[i] then
      vs = speeds[i]
      layer = i
      break
    end
  end
  local below = 0
  for i, top in ipairs(tops) do
    if z <= top then below = below + 1 end
  end
  return vs, layer, below
end
return M
)lua";

constexpr auto Integers = R"lua(
local M = {}
M.output_parameters = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"}
local zero <const> = 0
function M.evaluate(fields, x, y, z)
  local m = math.max(x, 0)
  local a = m * -1
  local b = 7 // 2 + x
  local c = -7 % 3 + 2^-1
  local d = (x > 0 and 0 or 1) * -1
  local e = x ^ 2 - x * x
  local f = 2 ^ 0.5 + (x // 1) + (x % 1.5) + (y // -2.5)
  local g = -(math.min(y, 0)) + 3 // 2
  local h = 0x10 + 1e2 + 3 / 2 + (5 % -3) + (-5 // 2)
  -- Lua computes `a - k` for a small integer constant k as `a + -k`: at x = -0.0, 0.0
  local i, j, k, l = x - 0, x - zero, x - (1 - 1), x - 0.0
  local m = (-0.0) - 0 + (y - -0) * 0
  return a, b, c, d, e, f, g, h, i, j, k, l, m
end
return M
)lua";

constexpr auto Functions = R"lua(
local M = {}
M.output_parameters = {"s", "n", "p", "q"}
local function sum(...)
  local total = 0
  for i = 1, select("#", ...) do
    total = total + select(i, ...)
  end
  return total, select("#", ...)
end
local function counter(start)
  local n = start
  return function(step)
    n = n + step
    return n
  end
end
local x0 = 1
local function shifted(v) return v + x0 end
local x0 = 100
function M.evaluate(fields, x, y, z)
  local c = counter(x)
  c(1)
  c(y)
  local s, n = sum(x, y, z, 1.5)
  local t = {x, y, z}
  local p = 1
  for _, v in ipairs(t) do
    p = p * (v + 2)
  end
  local a, b = (function() return 1, 2 end)()
  return s + c(0), n, p, shifted(a + b) + x0
end
return M
)lua";

constexpr auto Tables = R"lua(
local M = {}
M.output_parameters = {"value", "count", "length"}
local layer = {rho = 2500, scale = function(self, f) return self.rho * f end}
local grid = {}
for i = 1, 3 do
  grid[i] = {top = -i * 1000, rho = 2000 + 100 * i}
end
table.insert(grid, {top = -1e9, rho = 3000})
function M.evaluate(fields, x, y, z)
  local count = 0
  for k, v in pairs({a = 1, b = 2, c = 3}) do
    count = count + v
  end
  local rho = grid[#grid].rho
  for _, entry in ipairs(grid) do
    if z > entry.top then
      rho = entry.rho
      break
    end
  end
  local t = {}
  t[1], t[2] = x, y
  t.sum = t[1] + t[2]
  return layer:scale(x) + rho + t.sum, count, #grid + #t
end
return M
)lua";

constexpr auto Logic = R"lua(
local M = {}
M.output_parameters = {"v", "u", "w", "k"}
function M.evaluate(fields, x, y, z)
  local v = (x > 0 and y > 0) and x * y or -1
  local w
  if x ~= y then
    w = x - y
  end
  local u = w or 0
  local flag = not (z < 0)
  local k = 0
  if x > 1 then
    k = 1
  elseif x > 0 then
    k = 2
  else
    if y > 0 then
      k = 3
    end
  end
  return v, u, flag and 1 or 0, k
end
return M
)lua";

constexpr auto Loops = R"lua(
local M = {}
M.output_parameters = {"acc", "found", "steps"}
function M.evaluate(fields, x, y, z)
  local i, acc = 0, 0
  while i < 5 do
    i = i + 1
    acc = acc + x / i
  end
  local j = 0
  repeat
    j = j + 1
    acc = acc * 0.5 + j
  until j >= 3
  local found = -1
  for k = 10, 1, -3 do
    if y < k then
      found = k
    end
  end
  local steps = 0
  for t = 0.0, 1.0, 0.25 do
    steps = steps + t
  end
  return acc, found, steps
end
return M
)lua";

constexpr auto MathLibrary = R"lua(
local M = {}
M.output_parameters = {"floor", "ceil", "fmax", "fmin", "min", "max", "abs", "sqrt", "exp",
  "log10", "log2", "log3", "atan2", "atan", "trunc", "round", "rad", "deg", "pow", "cube",
  "tanh", "erf", "nan", "fdim", "huge", "pi"}
function M.evaluate(fields, x, y, z)
  return math.floor(x), math.ceil(y), math.fmax(x, y), math.fmin(x, z), math.min(x, y, z),
    math.max(x, y, z), math.abs(x), math.sqrt(math.abs(y)), math.exp(-math.abs(z) / 1000),
    math.log(math.abs(x) + 1, 10), math.log(math.abs(x) + 1, 2), math.log(math.abs(x) + 2, 3),
    math.atan(y, x), math.atan(x), math.trunc(y), math.round(z), math.rad(x), math.deg(y),
    math.pow(x, 3), x ^ 3, math.tanh(x / 1000), math.erf(y / 1000),
    math.isnan(x / x) and 1 or 0, math.fdim(x, y), math.huge, math.pi
end
return M
)lua";

constexpr auto Ssol = R"lua(
local M = {}
M.output_parameters = {"s", "sign", "mod", "atan2", "pow"}
function M.evaluate(fields, x, y, z)
  local c = ssol.land(ssol.gt(z, -1000), ssol.lt(x, 0))
  return ssol.select(c, 1.5, 2.5), ssol.sign(x), ssol.mod(x, 3), ssol.atan2(y, x), ssol.pow(x, 2)
end
return M
)lua";

constexpr auto Named = R"lua(
local M = {}
function M.evaluate(fields, x, y)
  local out = {a = x + y}
  if x > 0 then
    out.b = x * 2
  else
    out.b = -y
  end
  return out
end
return M
)lua";

constexpr auto Peak = R"lua(
local M = {}
M.state = { peak = 0.0 }
function M.evaluate(fields, v, peak)
  if math.abs(v) > peak then
    peak = math.abs(v)
  end
  return { peak = peak }
end
return M
)lua";

} // namespace luafrontendtest

using namespace luafrontendtest;

TEST_CASE("lua frontend: the parser reads Lua 5.4" * doctest::test_suite("reader")) {
  constexpr auto Everything = R"lua(#!/usr/bin/env lua
-- a comment
--[==[ a long
comment ]==]
local a <const>, b = 1, [[long
string]]
local t = {1, 2; x = 3, ["y"] = 4, [5] = 6, f = function(...) return ... end}
t.x, t[1] = t[1], t.x
local s = 'quote \' "\n\t\x41\65\z
           end' .. "dq" .. #"abc"
goto_ = 0x1p4 + 0xFF + 3e-2 + .5 + 1. + 7 // 2 % 3 ^ 2 ^ 0.5
local c = a < b and a <= b or a > b and a >= b or a == b or a ~= b
local d = 1 << 2 >> 1 & 3 | 4 ~ 5
local e = -~1 + not nil and 1
do local f = 1 end
while false do break end
repeat local g = 1 until g == 1
for i = 1, 3, 1 do end
for k, v in pairs(t) do end
function t.m(self) end
function t:n(x) return self, x end
local function h() return end
if a then elseif b then else end
::label::
print(t:n(1))
print "x"
print {1}
return t
)lua";
  const auto chunk = lua::parseChunk(Everything);
  REQUIRE(chunk != nullptr);
  CHECK(chunk->block.statements.size() == 22);
  CHECK(chunk->block.statements.back().kind == lua::Statement::Kind::Return);

  // the grouping of the operators, as Lua's priorities have it
  std::vector<std::shared_ptr<lua::FunctionBody>> parsed;
  const auto expressionOf = [&](const std::string& code) {
    parsed.push_back(lua::parseChunk("return " + code));
    return parsed.back()->block.statements.front().values.front().get();
  };
  const auto* power = expressionOf("2 ^ 3 ^ 2");
  CHECK(power->op == lua::Expression::Op::Pow);
  CHECK(power->b->op == lua::Expression::Op::Pow);
  const auto* negated = expressionOf("-x ^ 2");
  CHECK(negated->kind == lua::Expression::Kind::Unary);
  CHECK(negated->a->op == lua::Expression::Op::Pow);
  const auto* chain = expressionOf("a - b + c");
  CHECK(chain->op == lua::Expression::Op::Add);
  CHECK(chain->a->op == lua::Expression::Op::Sub);
  const auto* concat = expressionOf("a .. b .. c");
  CHECK(concat->b->op == lua::Expression::Op::Concat);
  const auto* logic = expressionOf("a or b and c");
  CHECK(logic->op == lua::Expression::Op::Or);
  CHECK(logic->b->op == lua::Expression::Op::And);
  const auto* integer = expressionOf("9007199254740993");
  CHECK(integer->isInteger);
  CHECK(integer->integer == 9007199254740993LL);
  CHECK_FALSE(expressionOf("9223372036854775808")->isInteger);
  CHECK(expressionOf("0xffffffffffffffff")->integer == -1);

  const auto lineOf = [](const std::string& code) {
    try {
      static_cast<void>(lua::parseChunk(code));
    } catch (const lua::LuaError& error) {
      return error.line();
    }
    return -1;
  };
  CHECK(lineOf("local x = 1\nlocal y = = 2") == 2);
  CHECK(lineOf("local s = 'open\nx = 1") == 1);
  CHECK(lineOf("x = 1 +") == 1);
  CHECK(lineOf("f(\n1,\n") == 3);
  CHECK(lineOf("x.y") == 1);
  CHECK(lineOf("if x then\n") == 2);
}

TEST_CASE("lua frontend: an if on a value of the point becomes a select" *
          doctest::test_suite("reader")) {
  const auto program = translated(Branches);
  REQUIRE(program.outputs().size() == 2);
  const auto columns = columnsFor({"x", "y", "z"});
  const auto evaluation = evaluate(program, columns, columns.at("x").size());
  for (std::size_t p = 0; p < columns.at("x").size(); ++p) {
    const double x = columns.at("x")[p];
    const double y = columns.at("y")[p];
    const double z = columns.at("z")[p];
    const double rho = z > -1000 ? 2200 + 0.1 * z : (z > -5000 ? (x < 0 ? 2500 : 2600) : 2700);
    const double mu = z > -1000 ? 1e9 : (z > -5000 ? (x < 0 ? 2e9 : 3e9 + y) : 4e9);
    CHECK(sameValue(evaluation.outputs[0][p], rho));
    CHECK(sameValue(evaluation.outputs[1][p], mu));
  }
  CHECK(evaluation.failures.empty());
}

#ifdef USE_LUA
TEST_CASE("lua frontend: programs give Lua's bits" * doctest::test_suite("reader")) {
  for (const auto* code : {Branches,
                           LayerSearch,
                           Integers,
                           Functions,
                           Tables,
                           Logic,
                           Loops,
                           MathLibrary,
                           Ssol,
                           Named}) {
    CAPTURE(code);
    agreesWithLua(code);
  }
}

TEST_CASE("lua frontend: a state carries over calls as in Lua" * doctest::test_suite("reader")) {
  agreesWithLua(Peak, 3);
}

TEST_CASE("lua frontend: random modules give Lua's bits (SEISSOL_LUA_FUZZ)" *
          doctest::test_suite("reader")) {
  // opt-in: SEISSOL_LUA_FUZZ modules, made from SEISSOL_LUA_FUZZ_SEED on (0 without)
  const char* count = std::getenv("SEISSOL_LUA_FUZZ");
  if (count == nullptr) {
    return;
  }
  const char* seed = std::getenv("SEISSOL_LUA_FUZZ_SEED");
  const std::uint64_t first = seed == nullptr ? 0 : std::stoull(seed);
  const std::uint64_t modules = std::stoull(count);
  std::size_t agreeing = 0;
  for (std::uint64_t i = first; i < first + modules; ++i) {
    RandomModule random(i);
    const std::string code = random.make();
    std::string found;
    try {
      found = disagreement(code);
    } catch (const std::exception& error) {
      found = std::string("not translated: ") + error.what();
    }
    if (!found.empty()) {
      FAIL_CHECK("module " << i << ": " << found << "\n" << code);
      continue;
    }
    ++agreeing;
  }
  MESSAGE(agreeing << " of " << modules << " random modules agree with Lua");
}
#endif

TEST_CASE("lua frontend: a layer search becomes a chain of selects" *
          doctest::test_suite("reader")) {
  // an early return or break keeps its points apart until the end: one select per layer, and no
  // select on the flow
  const auto program = translated(LayerSearch);
  std::size_t selects = 0;
  for (std::size_t id = 0; id < program.arena().size(); ++id) {
    const auto& node = program.arena()[static_cast<expr::NodeId>(id)];
    if (node.kind == expr::Kind::PW && node.fn == expr::Fn::Select) {
      ++selects;
    }
  }
  // vs and layer: four selects each; below: four
  CHECK(selects <= 12);
}

TEST_CASE("lua frontend: where Lua fails at some points, a check fails there" *
          doctest::test_suite("reader")) {
  constexpr auto Failing = R"lua(
local M = {}
M.output_parameters = {"v"}
function M.evaluate(fields, x)
  if x < -100 then
    error("too deep")
  end
  assert(x < 100, "too high")
  local t = {}
  if x > 1 then
    t.v = x
  end
  return t.v + 1
end
return M
)lua";
  const auto program = translated(Failing);
  CHECK(program.checks().size() == 3);
  const std::vector<double> x = {-1000.0, -100.0, 0.0, 1.0, 2.0, 99.0, 100.0, 1e3};
  const auto evaluation = evaluate(program, {{"x", x}}, x.size());
  for (std::size_t p = 0; p < x.size(); ++p) {
    CAPTURE(x[p]);
    CHECK(evaluation.fails(p) == !(x[p] > 1.0 && x[p] < 100.0));
    if (!evaluation.fails(p)) {
      CHECK(evaluation.outputs[0][p] == x[p] + 1);
    }
  }
  bool named = false;
  for (const auto& check : program.checks()) {
    named = named || check.message.find("model.lua:6: too deep") != std::string::npos;
  }
  CHECK(named);
}

TEST_CASE("lua frontend: what is outside the subset is refused with its line" *
          doctest::test_suite("reader")) {
  const auto module = [](const std::string& body) {
    return "local M = {}\nM.output_parameters = {'v'}\nfunction M.evaluate(fields, x)\n" + body +
           "\nend\nreturn M\n";
  };
  // line 4 is the first of the body
  struct Case {
    std::string body;
    int line;
    std::string says;
  };
  const std::vector<Case> cases = {
      {"local v = 0\nwhile v < x do v = v + 1 end\nreturn v", 5, "loop"},
      {"for i = 1, x do end\nreturn 1", 4, "bounds"},
      {"local t = {1, 2}\nreturn t[x]", 5, "index"},
      {"return tostring(x)", 4, "tostring"},
      {"return setmetatable({}, {}).v", 4, "setmetatable"},
      {"return math.random()", 4, "random"},
      {"local function f(n) if n > 0 then return f(n - 1) end return 0 end\nreturn f(2)",
       4,
       "recursion"},
      {"return math.cbrt(x)", 4, "cbrt"},
      {"return x .. 'm'", 4, "string"},
      {"local f = x > 0 and math.sin or math.cos\nreturn f(x)", 5, "function"},
      {"goto done\n::done::\nreturn 1", 4, "goto"},
  };
  for (const auto& c : cases) {
    CAPTURE(c.body);
    const auto error = refusalOf(module(c.body));
    CHECK(error.line() == c.line);
    CHECK(error.message().find(c.says) != std::string::npos);
  }

  // f keeps a value from one point to the next: the global `count`
  constexpr auto Leaks = R"lua(
local M = {}
M.output_parameters = {"v"}
function M.evaluate(fields, x)
  count = (count or 0) + 1
  return x + count
end
return M
)lua";
  CHECK(refusalOf(Leaks).message().find("count") != std::string::npos);
  // a global written before it is read keeps nothing
  constexpr auto Scratch = R"lua(
local M = {}
M.output_parameters = {"v"}
function M.evaluate(fields, x)
  r = {}
  r.v = x * 2
  return r.v
end
return M
)lua";
  CHECK_NOTHROW(static_cast<void>(translated(Scratch)));
}

TEST_CASE("lua frontend: an error at every point is one of the translation" *
          doctest::test_suite("reader")) {
  constexpr auto Always = R"lua(
local M = {}
M.output_parameters = {"v"}
function M.evaluate(fields, x)
  local t
  return t.v
end
return M
)lua";
  try {
    static_cast<void>(translated(Always));
    FAIL("translated");
  } catch (const lua::LuaUnsupported& error) {
    FAIL("refused: " << error.what());
  } catch (const lua::LuaError& error) {
    CHECK(error.line() == 6);
    CHECK(error.message().find("attempt to index a nil value") != std::string::npos);
  }
}

} // namespace seissol::unit_test
