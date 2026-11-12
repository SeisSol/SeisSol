// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "LuaFrontend.h"

#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/LuaSyntax.h"

#include <algorithm>
#include <bit>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace seissol::reader::scripting::lua {

std::string LuaOrigin::at(int line) const {
  return file + separator + std::to_string(firstLine + std::max(line, 1) - 1);
}

namespace {

using expr::Fn;
using expr::NodeId;
using expr::NoNode;
using Op = Expression::Op;
using ExprKind = Expression::Kind;
using StatKind = Statement::Kind;

constexpr double TwoTo63 = 9223372036854775808.0;
constexpr double Pi = 3.141592653589793238462643383279502884;
constexpr double Nan = std::numeric_limits<double>::quiet_NaN();
constexpr double Inf = std::numeric_limits<double>::infinity();
constexpr std::size_t MaxCallDepth = 200;

// ================================================================= numbers ==

/// A number of Lua that the script fixes: an integer or a float, Lua 5.4's two subtypes.
struct Constant {
  bool isInt{false};
  std::int64_t i{0};
  double f{0.0};

  static Constant integer(std::int64_t value) {
    Constant c;
    c.isInt = true;
    c.i = value;
    return c;
  }
  static Constant real(double value) {
    Constant c;
    c.f = value;
    return c;
  }
  /// The float of the number, as Lua converts an integer (cast_num).
  [[nodiscard]] double asFloat() const { return isInt ? static_cast<double>(i) : f; }
  [[nodiscard]] bool same(const Constant& other) const {
    return isInt == other.isInt &&
           (isInt ? i == other.i
                  : std::bit_cast<std::uint64_t>(f) == std::bit_cast<std::uint64_t>(other.f));
  }
};

std::int64_t wrapped(std::uint64_t value) { return static_cast<std::int64_t>(value); }

/// A float that is an integer in Lua's range, as that integer (luaV_flttointns, F2Ieq).
std::optional<std::int64_t> exactInteger(double value) {
  if (!(value >= -TwoTo63 && value < TwoTo63) || std::floor(value) != value) {
    return std::nullopt;
  }
  return static_cast<std::int64_t>(value);
}

std::optional<std::int64_t> exactInteger(const Constant& c) {
  return c.isInt ? std::optional<std::int64_t>(c.i) : exactInteger(c.f);
}

/// The order of an integer and a float, exactly, as Lua compares them: -1, 0, 1, or none for NaN.
std::optional<int> compareIntFloat(std::int64_t i, double f) {
  if (std::isnan(f)) {
    return std::nullopt;
  }
  if (f >= TwoTo63) {
    return -1;
  }
  if (f < -TwoTo63) {
    return 1;
  }
  const double below = std::floor(f);
  const auto whole = static_cast<std::int64_t>(below);
  if (i < whole) {
    return -1;
  }
  if (i > whole) {
    return 1;
  }
  return below == f ? 0 : -1;
}

std::optional<int> compareConstants(const Constant& a, const Constant& b) {
  if (a.isInt && b.isInt) {
    return a.i < b.i ? -1 : (a.i > b.i ? 1 : 0);
  }
  if (!a.isInt && !b.isInt) {
    if (std::isnan(a.f) || std::isnan(b.f)) {
      return std::nullopt;
    }
    return a.f < b.f ? -1 : (a.f > b.f ? 1 : 0);
  }
  if (a.isInt) {
    return compareIntFloat(a.i, b.f);
  }
  const auto order = compareIntFloat(b.i, a.f);
  if (!order.has_value()) {
    return std::nullopt;
  }
  return -*order;
}

/// Lua's `//` and `%` of integers (luaV_idiv, luaV_mod); none for a zero divisor.
std::optional<std::int64_t> floorDivide(std::int64_t m, std::int64_t n) {
  if (n == 0) {
    return std::nullopt;
  }
  if (n == -1) {
    return wrapped(0U - static_cast<std::uint64_t>(m));
  }
  std::int64_t q = m / n;
  if ((m ^ n) < 0 && m % n != 0) {
    q -= 1;
  }
  return q;
}

std::optional<std::int64_t> floorModulo(std::int64_t m, std::int64_t n) {
  if (n == 0) {
    return std::nullopt;
  }
  if (n == -1) {
    return 0;
  }
  std::int64_t r = m % n;
  if (r != 0 && (r ^ n) < 0) {
    r += n;
  }
  return r;
}

std::int64_t shiftLeft(std::int64_t x, std::int64_t y) {
  if (y < 0) {
    if (y <= -64) {
      return 0;
    }
    return wrapped(static_cast<std::uint64_t>(x) >> static_cast<unsigned>(-y));
  }
  if (y >= 64) {
    return 0;
  }
  return wrapped(static_cast<std::uint64_t>(x) << static_cast<unsigned>(y));
}

/// An operation of the program on numbers the script fixes, as the program computes it.
double apply(Fn fn, double x, double y = 0.0) {
  const double args[2] = {x, y};
  return expr::applyPw<double>(fn, args);
}

/// A number as Lua writes it (tostring, `..`).
std::string numberText(const Constant& c) {
  char buffer[64];
  if (c.isInt) {
    std::snprintf(buffer, sizeof buffer, "%lld", static_cast<long long>(c.i));
    return buffer;
  }
  std::snprintf(buffer, sizeof buffer, "%.14g", c.f);
  std::string text = buffer;
  if (text.find_first_not_of("-0123456789") == std::string::npos) {
    text += ".0";
  }
  return text;
}

/// A string as Lua reads a number in it (lua_stringtonumber), if it is one.
std::optional<Constant> stringToNumber(const std::string& text) {
  const char* const space = " \f\n\r\t\v";
  const auto begin = text.find_first_not_of(space);
  if (begin == std::string::npos) {
    return std::nullopt;
  }
  const auto end = text.find_last_not_of(space) + 1;
  const std::string s = text.substr(begin, end - begin);
  std::size_t at = 0;
  bool negative = false;
  if (s[at] == '-' || s[at] == '+') {
    negative = s[at] == '-';
    ++at;
  }
  // an integer
  if (s.compare(at, 2, "0x") == 0 || s.compare(at, 2, "0X") == 0) {
    std::uint64_t value = 0;
    std::size_t digits = 0;
    std::size_t k = at + 2;
    for (; k < s.size() && std::isxdigit(static_cast<unsigned char>(s[k])) != 0; ++k, ++digits) {
      const int digit = std::isdigit(static_cast<unsigned char>(s[k])) != 0
                            ? s[k] - '0'
                            : std::tolower(static_cast<unsigned char>(s[k])) - 'a' + 10;
      value = value * 16 + static_cast<std::uint64_t>(digit);
    }
    if (digits > 0 && k == s.size()) {
      return Constant::integer(wrapped(negative ? 0U - value : value));
    }
  } else {
    std::uint64_t value = 0;
    std::size_t digits = 0;
    std::size_t k = at;
    bool overflow = false;
    constexpr std::uint64_t Limit = std::uint64_t{1} << 63U;
    for (; k < s.size() && std::isdigit(static_cast<unsigned char>(s[k])) != 0; ++k, ++digits) {
      const auto digit = static_cast<std::uint64_t>(s[k] - '0');
      if (value > (Limit - digit) / 10) {
        overflow = true;
      }
      value = value * 10 + digit;
    }
    if (digits > 0 && k == s.size() && !overflow && (negative || value < Limit)) {
      return Constant::integer(wrapped(negative ? 0U - value : value));
    }
  }
  // a float
  if (s.find_first_of("nN") != std::string::npos) {
    return std::nullopt;
  }
  char* stop = nullptr;
  const double value = std::strtod(s.c_str(), &stop);
  if (stop == s.c_str() || stop != s.c_str() + s.size()) {
    return std::nullopt;
  }
  return Constant::real(value);
}

// ================================================================== truths ==

/// A truth of Lua at a point: known, or a node of the program that is 1 where it holds and 0
/// where it does not.
struct Truth {
  NodeId node{NoNode};
  bool value{true};
  [[nodiscard]] bool known() const { return node == NoNode; }
  [[nodiscard]] bool isTrue() const { return node == NoNode && value; }
  [[nodiscard]] bool isFalse() const { return node == NoNode && !value; }
  [[nodiscard]] bool operator==(const Truth& other) const {
    return node == other.node && (node != NoNode || value == other.value);
  }
};

constexpr Truth True{NoNode, true};
constexpr Truth False{NoNode, false};
Truth known(bool value) { return value ? True : False; }

// ================================================================== values ==

/// A number of Lua at a point: one the script fixes, or a node of the program with its value in
/// binary64, and where it is an integer of Lua.
struct Number {
  bool known{true};
  Constant k;
  NodeId node{NoNode};
  Truth intness{False};

  static Number of(const Constant& c) {
    Number n;
    n.k = c;
    return n;
  }
  static Number at(NodeId node, Truth intness = False) {
    Number n;
    n.known = false;
    n.node = node;
    n.intness = intness;
    return n;
  }
};

struct Scope;

/// A function of the script, and the locals it sees: those of `scope` declared before it.
struct Closure {
  std::shared_ptr<FunctionBody> body;
  std::shared_ptr<Scope> scope;
  std::size_t visible{0};
  int id{0};
};

/// The local variables of a block: names and their cells, in the order of their declarations; of
/// the parent's, it sees those declared before it.
struct Scope {
  std::shared_ptr<Scope> parent;
  std::size_t parentVisible{0};
  std::vector<std::pair<std::string, int>> names;
  std::vector<bool> constants;
  /// the value of a <const> local the compiler knows (luaK_exp2const)
  std::vector<std::optional<Constant>> compileTime;
  bool captured{false};
};

using ScopePtr = std::shared_ptr<Scope>;

ScopePtr childOf(const ScopePtr& parent) {
  auto scope = std::make_shared<Scope>();
  scope->parent = parent;
  scope->parentVisible = parent != nullptr ? parent->names.size() : 0;
  return scope;
}

enum class Type : std::uint8_t { Nil, Boolean, Number, String, Table, Function, Builtin, Stale };

const char* typeName(Type type) {
  switch (type) {
  case Type::Nil:
    return "nil";
  case Type::Boolean:
    return "boolean";
  case Type::Number:
    return "number";
  case Type::String:
    return "string";
  case Type::Table:
    return "table";
  case Type::Function:
  case Type::Builtin:
    return "function";
  case Type::Stale:
    return "stale";
  }
  return "?";
}

/// One of the values a value of Lua may be at a point, and where (`guard`, among the points the
/// value is defined at).
struct Alternative {
  Truth guard{True};
  Type type{Type::Nil};
  Truth boolean{False};
  Number number;
  /// String: the text; Builtin: its name; Stale: what it is
  std::string text;
  /// Table: its location; Function: the closure; Builtin: its datum; Stale: the line that read it
  int id{-1};
  std::shared_ptr<const Closure> closure;

  [[nodiscard]] bool sameIdentity(const Alternative& other) const {
    if (type != other.type) {
      return false;
    }
    switch (type) {
    case Type::String:
      return text == other.text;
    case Type::Table:
    case Type::Function:
      return id == other.id;
    case Type::Builtin:
      return text == other.text && id == other.id;
    default:
      return true;
    }
  }
};

/// A value of Lua at a point: its alternatives, whose guards part the points. At most one
/// alternative is nil, a boolean, a number or stale; the others differ in their identity.
struct Value {
  std::vector<Alternative> alternatives;
  [[nodiscard]] bool single() const { return alternatives.size() == 1; }
  [[nodiscard]] const Alternative& only() const { return alternatives.front(); }
  [[nodiscard]] bool is(Type type) const { return single() && only().type == type; }
};

Value atom(Alternative alternative) {
  alternative.guard = True;
  Value value;
  value.alternatives.push_back(std::move(alternative));
  return value;
}
Value nil() { return atom(Alternative{}); }
Value boolean(Truth truth) {
  Alternative a;
  a.type = Type::Boolean;
  a.boolean = truth;
  return atom(std::move(a));
}
Value number(const Number& n) {
  Alternative a;
  a.type = Type::Number;
  a.number = n;
  return atom(std::move(a));
}
Value number(const Constant& c) { return number(Number::of(c)); }
Value string(std::string text) {
  Alternative a;
  a.type = Type::String;
  a.text = std::move(text);
  return atom(std::move(a));
}
Value table(int id) {
  Alternative a;
  a.type = Type::Table;
  a.id = id;
  return atom(std::move(a));
}
Value function(std::shared_ptr<const Closure> closure) {
  Alternative a;
  a.type = Type::Function;
  a.id = closure->id;
  a.closure = std::move(closure);
  return atom(std::move(a));
}
Value builtin(std::string name, int datum = 0) {
  Alternative a;
  a.type = Type::Builtin;
  a.text = std::move(name);
  a.id = datum;
  return atom(std::move(a));
}
Value stale(std::string what) {
  Alternative a;
  a.type = Type::Stale;
  a.text = std::move(what);
  return atom(std::move(a));
}

bool sameNumber(const Number& a, const Number& b) {
  if (a.known != b.known) {
    return false;
  }
  return a.known ? a.k.same(b.k) : a.node == b.node && a.intness == b.intness;
}

bool identical(const Value& a, const Value& b) {
  if (a.alternatives.size() != b.alternatives.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.alternatives.size(); ++i) {
    const auto& x = a.alternatives[i];
    const auto& y = b.alternatives[i];
    if (!x.sameIdentity(y) || !(x.guard == y.guard) || !(x.boolean == y.boolean) ||
        (x.type == Type::Number && !sameNumber(x.number, y.number)) ||
        (x.type == Type::Stale && x.text != y.text)) {
      return false;
    }
  }
  return true;
}

bool definitelyNil(const Value& v) { return v.is(Type::Nil); }

// ================================================================== tables ==

/// A key of a table: a value the script fixes, normalized as Lua does (a float that is an integer
/// is that integer).
struct Key {
  enum class Kind : std::uint8_t { Number, String, Boolean, Table, Function, Builtin };
  Kind kind{Kind::Number};
  Constant number;
  std::string text;
  int id{0};

  static Key integer(std::int64_t value) {
    Key key;
    key.number = Constant::integer(value);
    return key;
  }
  static Key string(std::string text) {
    Key key;
    key.kind = Kind::String;
    key.text = std::move(text);
    return key;
  }
};

bool operator<(const Key& a, const Key& b) {
  if (a.kind != b.kind) {
    return a.kind < b.kind;
  }
  switch (a.kind) {
  case Key::Kind::Number: {
    const auto order = compareConstants(a.number, b.number);
    return order.has_value() && *order < 0;
  }
  case Key::Kind::String:
    return a.text < b.text;
  case Key::Kind::Builtin:
    return a.text != b.text ? a.text < b.text : a.id < b.id;
  default:
    return a.id < b.id;
  }
}

std::string keyText(const Key& key) {
  switch (key.kind) {
  case Key::Kind::Number:
    return numberText(key.number);
  case Key::Kind::String:
    return key.text;
  case Key::Kind::Boolean:
    return key.id != 0 ? "true" : "false";
  default:
    return "?";
  }
}

struct Table {
  std::map<Key, Value> fields;
};

// =================================================================== state ==

/// The cells of the variables and the tables, by location.
struct Store {
  std::map<int, std::shared_ptr<Value>> cells;
  std::map<int, std::shared_ptr<Table>> tables;
};

/// How the points of the current function stand, among those the store is of: still running
/// (`live`), returned (or failed), or left the innermost loop; and what they returned.
struct Flow {
  Truth live{True};
  Truth returned{False};
  Truth broken{False};
  bool hasResults{false};
  std::vector<Value> results;
  int version{0};
};

struct Pending;

/// Where the evaluation stands: the store and the flow of the points it runs (`domain`), what is
/// known on the way there (`facts`), and the states of points that stopped on the way (`pendings`).
struct State {
  Store store;
  Flow flow;
  Truth domain{True};
  std::vector<std::pair<NodeId, bool>> facts;
  std::vector<std::shared_ptr<const Pending>> pendings;
};

/// The state of the points of one side of an `if` that stopped there (returned, broke or failed),
/// kept apart while the other side goes on, until the end of the loop or the function.
struct Pending {
  Truth condition;
  bool stoppedOnTrue{true};
  Store store;
  Flow flow;
  Truth domain;
  std::size_t facts{0};
  int version{0};
};

/// A location that existed before the call of the entry function: a cell, or a field of a table.
struct Location {
  int cell{-1};
  int table{-1};
  Key key;
  [[nodiscard]] bool operator<(const Location& other) const {
    if (cell != other.cell) {
      return cell < other.cell;
    }
    if (table != other.table) {
      return table < other.table;
    }
    return key < other.key;
  }
};

// =============================================================== evaluator ==

class Evaluator {
  public:
  Evaluator(expr::Program& program,
            const LuaOptions& options,
            const LuaOrigin& origin,
            NodeId reach,
            const std::set<Location>* stale,
            bool module)
      : program_(program), arena_(program.arena()), options_(options), origin_(origin),
        stale_(stale), module_(module) {
    baseReach_ = reach == NoNode ? True : truthOf(reach);
    makeGlobals();
  }

  [[nodiscard]] const std::set<Location>& leaks() const { return leaks_; }
  [[nodiscard]] std::vector<LuaCheck> checks() const { return checks_; }
  [[nodiscard]] const std::vector<int>& returnLines() const { return returnLines_; }

  /// Runs the chunk, and gives what it returns.
  std::vector<Value> runChunk(const std::shared_ptr<FunctionBody>& chunk) {
    auto root = std::make_shared<Scope>();
    root->captured = true;
    declare(*root, "_ENV", table(globals_), false);
    auto closure = std::make_shared<Closure>();
    closure->body = chunk;
    closure->scope = root;
    closure->visible = root->names.size();
    closure->id = nextClosure_++;
    return callClosure(*closure, {}, 0);
  }

  // --------------------------------------------------------------- easi ---

  LuaMapTranslation luaMap(const std::shared_ptr<FunctionBody>& chunk,
                           const std::vector<std::pair<std::string, NodeId>>& in,
                           const std::set<std::string>& returns) {
    static_cast<void>(runChunk(chunk));
    // lua_getglobal(L, "f"), at every point
    const Value f = readField(globals_, Key::string("f"));
    const int line = f.is(Type::Function) ? f.only().closure->body->line : 0;
    beginEntry();
    const int x = newTable();
    for (std::size_t i = 0; i < in.size(); ++i) {
      const Value value = number(Number::at(in[i].second));
      rawSet(x, Key::integer(static_cast<std::int64_t>(i + 1)), value);
      rawSet(x, Key::string(in[i].first), value);
    }
    const auto results = call(f, {table(x)}, line, " (global 'f')");
    endEntry();
    const Value result = results.empty() ? nil() : results.front();
    LuaMapTranslation translation;
    for (const auto& name : returns) {
      const Value field = outputField(result, name, line, "f does not return a table");
      translation.values[name] =
          outputNode(field, "Tried key " + name + " but it did not return a number.", line);
    }
    translation.checks = checks_;
    translation.lines = returnLines_;
    return translation;
  }

  // ------------------------------------------------------------- module ---

  void module(const std::shared_ptr<FunctionBody>& chunk) {
    const auto returned = runChunk(chunk);
    if (returned.empty() || !returned.front().is(Type::Table)) {
      throw LuaError(0, "the module does not return a table");
    }
    const int m = returned.front().only().id;
    const auto declared = [&](const char* name) {
      return !definitelyNil(readField(m, Key::string(name)));
    };
    const bool declaresInputs = declared("input_parameters");
    const bool declaresOutputs = declared("output_parameters");
    const auto declaredInputs = stringList(readField(m, Key::string("input_parameters")));
    auto outputNames = stringList(readField(m, Key::string("output_parameters")));

    // M.state: the states, sorted by name, with the values they start from
    std::vector<std::pair<std::string, double>> states;
    const Value stateTable = readField(m, Key::string("state"));
    if (stateTable.is(Type::Table)) {
      for (const auto& [key, value] : tableAt(stateTable.only().id).fields) {
        if (key.kind != Key::Kind::String || !value.is(Type::Number) ||
            !value.only().number.known) {
          throw LuaError(0, "`state` maps names to the numbers they start from");
        }
        states.emplace_back(key.text, value.only().number.k.asFloat());
      }
    } else if (!definitelyNil(stateTable)) {
      throw LuaError(0, "`state` has to be a table of initial values");
    }
    const auto isState = [&](const std::string& name) {
      return std::any_of(
          states.begin(), states.end(), [&](const auto& entry) { return entry.first == name; });
    };

    // M.field_specs: the grids, and the object of each in `fields`
    const int fields = newTable();
    const Value specs = readField(m, Key::string("field_specs"));
    if (specs.is(Type::Table)) {
      const auto count = border(specs.only().id, 0);
      for (std::int64_t i = 1; i <= count; ++i) {
        gridOf(readField(specs.only().id, Key::integer(i)), fields);
      }
    }

    const Value evaluate = readField(m, Key::string("evaluate"));
    if (!evaluate.is(Type::Function)) {
      throw LuaError(0, "the module has no `evaluate` function");
    }
    const auto& body = *evaluate.only().closure->body;
    const int line = body.line;
    if (body.name.find(':') != std::string::npos && !declaresInputs) {
      throw LuaUnsupported(line,
                           "`evaluate` is a method: `function M:evaluate(...)` puts `self` before "
                           "`fields`, and every input one parameter later");
    }
    if (declaresInputs && !states.empty()) {
      throw LuaUnsupported(line,
                           "a module with state names its inputs by the parameters of `evaluate`; "
                           "remove `input_parameters`");
    }
    if (!declaresInputs && body.vararg) {
      throw LuaUnsupported(line,
                           "`evaluate` is variadic; its inputs are the names of its parameters, or "
                           "those of `input_parameters`");
    }

    // the arguments after `fields`: the declared inputs, else the parameters, where those of the
    // name of a state read its value of the previous call
    std::vector<std::string> inputs;
    std::vector<Value> arguments{table(fields)};
    const auto names =
        declaresInputs
            ? declaredInputs
            : std::vector<std::string>(body.parameters.begin() +
                                           std::min<std::size_t>(1, body.parameters.size()),
                                       body.parameters.end());
    for (const auto& name : names) {
      arguments.push_back(number(Number::at(arena_.field(name))));
      if (!isState(name)) {
        inputs.push_back(name);
      }
    }

    beginEntry();
    const auto results = call(evaluate, std::move(arguments), line, " (field 'evaluate')");
    endEntry();

    // the outputs: of a returned table by their names, else in the order of output_parameters
    std::vector<NodeId> roots;
    std::vector<NodeId> stateRoots;
    const bool named = results.size() == 1 &&
                       std::any_of(results.front().alternatives.begin(),
                                   results.front().alternatives.end(),
                                   [](const Alternative& a) { return a.type == Type::Table; });
    const std::string wrongType = "Incorrect Lua data return type.";
    if (named) {
      const Value& result = results.front();
      std::set<std::string> keys;
      for (const auto& alternative : result.alternatives) {
        if (alternative.type != Type::Table) {
          throw LuaUnsupported(line,
                               "`evaluate` returns a table at some points and a " +
                                   std::string(typeName(alternative.type)) + " at others");
        }
        for (const auto& [key, value] : tableAt(alternative.id).fields) {
          if (key.kind != Key::Kind::String) {
            throw LuaUnsupported(line, "`evaluate` returns a table with a key that is no name");
          }
          keys.insert(key.text);
        }
      }
      if (!declaresOutputs) {
        // every name the table returns, states included, in their order
        outputNames.assign(keys.begin(), keys.end());
      }
      for (const auto& key : keys) {
        if (std::find(outputNames.begin(), outputNames.end(), key) == outputNames.end() &&
            !isState(key)) {
          throw LuaUnsupported(
              line, "`evaluate` returns `" + key + "`, which is neither an output nor a state");
        }
      }
      for (const auto& name : outputNames) {
        roots.push_back(outputNode(outputField(result, name, line, wrongType), wrongType, line));
      }
      for (const auto& [name, initial] : states) {
        if (keys.count(name) == 0) {
          throw LuaUnsupported(line,
                               "the state `" + name + "` gets no next value: return it by name");
        }
        stateRoots.push_back(
            outputNode(outputField(result, name, line, wrongType), wrongType, line));
      }
    } else {
      if (results.size() < outputNames.size()) {
        throw LuaError(line,
                       "`evaluate` returns " + std::to_string(results.size()) +
                           " values instead of " + std::to_string(outputNames.size()));
      }
      for (std::size_t i = 0; i < outputNames.size(); ++i) {
        roots.push_back(outputNode(results[i], wrongType, line));
      }
      for (const auto& [name, initial] : states) {
        const auto found = std::find(outputNames.begin(), outputNames.end(), name);
        if (found == outputNames.end()) {
          throw LuaUnsupported(line,
                               "the state `" + name +
                                   "` gets no next value: return an output of its name, or return "
                                   "a table to keep it out of the outputs");
        }
        stateRoots.push_back(roots[static_cast<std::size_t>(found - outputNames.begin())]);
      }
    }

    program_.addMeta(expr::MetaEntry{"generator", {expr::MetaValue{"seissol lua frontend", true}}});
    for (const auto& name : inputs) {
      program_.addInput(name, DataType::F64);
    }
    for (std::size_t i = 0; i < states.size(); ++i) {
      program_.addState(states[i].first, states[i].second, stateRoots[i]);
    }
    for (std::size_t i = 0; i < outputNames.size(); ++i) {
      program_.addOutput(outputNames[i], DataType::F64, roots[i]);
      for (const int at : returnLines_) {
        program_.addSource(
            expr::SourceOf{expr::SourceOf::What::Output,
                           NoNode,
                           outputNames[i],
                           expr::SourcePosition{origin_.file, origin_.firstLine + at - 1, 1}});
      }
    }
    for (const auto& check : checks_) {
      program_.addCheck(check.holds, check.message);
    }
  }

  private:
  // ------------------------------------------------------------- truths ---

  NodeId konst(double value) { return arena_.konst(value); }
  NodeId pw(Fn fn, NodeId a) { return arena_.pw(fn, a); }
  NodeId pw(Fn fn, NodeId a, NodeId b) { return arena_.pw(fn, a, b); }
  NodeId select(NodeId c, NodeId a, NodeId b) {
    return a == b ? a : arena_.pw(Fn::Select, c, a, b);
  }

  Truth truthOf(NodeId node) {
    truthNodes_.insert(node);
    return Truth{node, true};
  }
  NodeId nodeOf(Truth t) { return t.known() ? konst(t.value ? 1.0 : 0.0) : t.node; }

  bool isZero(NodeId node) const {
    const auto& n = arena_[node];
    return n.kind == expr::Kind::Const && n.value == 0.0;
  }
  /// The truth `t` negates, if it is a negation.
  std::optional<NodeId> negated(Truth t) const {
    if (t.known()) {
      return std::nullopt;
    }
    const auto& n = arena_[t.node];
    if (n.kind == expr::Kind::PW && n.fn == Fn::Eq && truthNodes_.count(n.a) != 0 && isZero(n.b)) {
      return n.a;
    }
    return std::nullopt;
  }
  bool opposite(Truth a, Truth b) const {
    const auto na = negated(a);
    const auto nb = negated(b);
    return (na.has_value() && *na == b.node) || (nb.has_value() && *nb == a.node);
  }

  Truth tNot(Truth t) {
    if (t.known()) {
      return known(!t.value);
    }
    if (const auto inner = negated(t)) {
      return Truth{*inner, true};
    }
    return truthOf(pw(Fn::Eq, t.node, konst(0.0)));
  }
  Truth tAnd(Truth a, Truth b) {
    if (a.isFalse() || b.isFalse()) {
      return False;
    }
    if (a.isTrue()) {
      return b;
    }
    if (b.isTrue() || a == b) {
      return a;
    }
    if (opposite(a, b)) {
      return False;
    }
    return truthOf(pw(Fn::And, a.node, b.node));
  }
  Truth tOr(Truth a, Truth b) {
    if (a.isTrue() || b.isTrue()) {
      return True;
    }
    if (a.isFalse()) {
      return b;
    }
    if (b.isFalse() || a == b) {
      return a;
    }
    if (opposite(a, b)) {
      return True;
    }
    return truthOf(pw(Fn::Or, a.node, b.node));
  }
  Truth tSelect(Truth c, Truth a, Truth b) {
    if (c.known()) {
      return c.value ? a : b;
    }
    if (a == b) {
      return a;
    }
    if (a.isTrue() && b.isFalse()) {
      return c;
    }
    if (a.isFalse() && b.isTrue()) {
      return tNot(c);
    }
    if (a.isTrue() || a == c) {
      return tOr(c, b);
    }
    if (b.isFalse() || b == c) {
      return tAnd(c, a);
    }
    if (a.isFalse()) {
      return tAnd(tNot(c), b);
    }
    if (b.isTrue()) {
      return tOr(tNot(c), a);
    }
    return truthOf(arena_.pw(Fn::Select, c.node, a.node, b.node));
  }

  /// `t`, where the way here decides it.
  Truth decide(Truth t) const {
    if (t.known()) {
      return t;
    }
    const auto inner = negated(t);
    for (auto it = state_.facts.rbegin(); it != state_.facts.rend(); ++it) {
      if (it->first == t.node) {
        return known(it->second);
      }
      if (inner.has_value() && it->first == *inner) {
        return known(!it->second);
      }
    }
    return t;
  }

  /// Where the current statement runs, among all points.
  Truth reach() { return tAnd(baseReach_, tAnd(state_.domain, state_.flow.live)); }

  // ------------------------------------------------------------ numbers ---

  NodeId nodeOf(const Number& n) { return n.known ? konst(n.k.asFloat()) : n.node; }
  static Truth intnessOf(const Number& n) { return n.known ? known(n.k.isInt) : n.intness; }

  /// The value `r` of an operation on integers, where `ints` holds: the zero of an integer of Lua
  /// has no sign.
  Number integral(NodeId r, Truth ints) {
    if (ints.isFalse()) {
      return Number::at(r);
    }
    const NodeId normal = pw(Fn::Add, r, konst(0.0));
    return Number::at(ints.isTrue() ? normal : select(nodeOf(ints), normal, r), ints);
  }

  Number selectNumber(Truth c, const Number& a, const Number& b) {
    if (sameNumber(a, b)) {
      return a;
    }
    return Number::at(select(nodeOf(c), nodeOf(a), nodeOf(b)),
                      tSelect(c, intnessOf(a), intnessOf(b)));
  }

  Constant foldArithmetic(Op op, const Constant& x, const Constant& y, int line) {
    const bool ints = x.isInt && y.isInt;
    const double a = x.asFloat();
    const double b = y.asFloat();
    const auto toInteger = [&](const Constant& c) -> std::int64_t {
      const auto value = exactInteger(c);
      if (!value.has_value()) {
        luaError("number has no integer representation", True, line);
        return 0;
      }
      return *value;
    };
    switch (op) {
    case Op::Add:
      return ints ? Constant::integer(
                        wrapped(static_cast<std::uint64_t>(x.i) + static_cast<std::uint64_t>(y.i)))
                  : Constant::real(a + b);
    case Op::Sub:
      return ints ? Constant::integer(
                        wrapped(static_cast<std::uint64_t>(x.i) - static_cast<std::uint64_t>(y.i)))
                  : Constant::real(a - b);
    case Op::Mul:
      return ints ? Constant::integer(
                        wrapped(static_cast<std::uint64_t>(x.i) * static_cast<std::uint64_t>(y.i)))
                  : Constant::real(a * b);
    case Op::Div:
      return Constant::real(a / b);
    case Op::Pow:
      return Constant::real(b == 2.0 ? a * a : std::pow(a, b));
    case Op::IntDiv:
      if (ints) {
        const auto q = floorDivide(x.i, y.i);
        if (!q.has_value()) {
          luaError("attempt to perform 'n//0'", True, line);
          return Constant::integer(0);
        }
        return Constant::integer(*q);
      }
      return Constant::real(std::floor(a / b));
    case Op::Mod:
      if (ints) {
        const auto r = floorModulo(x.i, y.i);
        if (!r.has_value()) {
          luaError("attempt to perform 'n%%0'", True, line);
          return Constant::integer(0);
        }
        return Constant::integer(*r);
      }
      return Constant::real(apply(Fn::Mod, a, b));
    case Op::BitAnd:
      return Constant::integer(toInteger(x) & toInteger(y));
    case Op::BitOr:
      return Constant::integer(toInteger(x) | toInteger(y));
    case Op::BitXor:
      return Constant::integer(toInteger(x) ^ toInteger(y));
    case Op::Shl:
      return Constant::integer(shiftLeft(toInteger(x), toInteger(y)));
    case Op::Shr:
      return Constant::integer(
          shiftLeft(toInteger(x), wrapped(0U - static_cast<std::uint64_t>(toInteger(y)))));
    default:
      throw std::logic_error("not an arithmetic operator");
    }
  }

  /// A division of integers by zero fails in Lua; where the divisor is zero, a check.
  void divisionByZero(Truth ints, const Number& divisor, const char* op, int line) {
    if (ints.isFalse()) {
      return;
    }
    if (divisor.known) {
      if (divisor.k.asFloat() == 0.0) {
        luaError(std::string("attempt to perform 'n") + op + "0'", ints, line);
      }
      return;
    }
    const Truth zero = truthOf(pw(Fn::Eq, divisor.node, konst(0.0)));
    luaError(std::string("attempt to perform 'n") + op + "0'", tAnd(ints, zero), line);
  }

  Number arithmetic(Op op, const Number& x, const Number& y, int line) {
    if (x.known && y.known) {
      return Number::of(foldArithmetic(op, x.k, y.k, line));
    }
    const NodeId a = nodeOf(x);
    const NodeId b = nodeOf(y);
    const Truth ints = tAnd(intnessOf(x), intnessOf(y));
    switch (op) {
    case Op::Add:
      return Number::at(pw(Fn::Add, a, b), ints);
    case Op::Sub:
      return Number::at(pw(Fn::Sub, a, b), ints);
    case Op::Mul:
      return integral(pw(Fn::Mul, a, b), ints);
    case Op::Div:
      return Number::at(pw(Fn::Div, a, b));
    case Op::Pow: {
      // Lua's luai_numpow: x * x for the exponent 2, pow() else
      const NodeId square = pw(Fn::Mul, a, a);
      if (y.known) {
        return Number::at(y.k.asFloat() == 2.0 ? square : pw(Fn::Pow, a, b));
      }
      return Number::at(select(pw(Fn::Eq, b, konst(2.0)), square, pw(Fn::Pow, a, b)));
    }
    case Op::IntDiv:
      divisionByZero(ints, y, "//", line);
      return integral(pw(Fn::Floor, pw(Fn::Div, a, b)), ints);
    case Op::Mod:
      divisionByZero(ints, y, "%%", line);
      return integral(pw(Fn::Mod, a, b), ints);
    default:
      throw LuaUnsupported(line, "a bitwise operation on a value of the point is not translated");
    }
  }

  Number negate(const Number& x) {
    if (x.known) {
      return Number::of(x.k.isInt
                            ? Constant::integer(wrapped(0U - static_cast<std::uint64_t>(x.k.i)))
                            : Constant::real(-x.k.f));
    }
    return integral(pw(Fn::Neg, x.node), x.intness);
  }

  /// `x < y` (`x <= y`) of numbers.
  Truth lessNumbers(const Number& x, const Number& y, bool orEqual) {
    if (x.known && y.known) {
      const auto order = compareConstants(x.k, y.k);
      return known(order.has_value() && (*order < 0 || (orEqual && *order == 0)));
    }
    return truthOf(pw(orEqual ? Fn::Le : Fn::Lt, nodeOf(x), nodeOf(y)));
  }

  Truth equalNumbers(const Number& x, const Number& y) {
    if (x.known && y.known) {
      const auto order = compareConstants(x.k, y.k);
      return known(order.has_value() && *order == 0);
    }
    return truthOf(pw(Fn::Eq, nodeOf(x), nodeOf(y)));
  }

  // ------------------------------------------------------------- values ---

  Value merge(Truth c, const Value& a, const Value& b) {
    if (c.isTrue()) {
      return a;
    }
    if (c.isFalse()) {
      return b;
    }
    if (identical(a, b)) {
      return a;
    }
    Value out;
    std::vector<bool> paired(b.alternatives.size(), false);
    for (const auto& x : a.alternatives) {
      Alternative merged = x;
      bool found = false;
      for (std::size_t j = 0; j < b.alternatives.size(); ++j) {
        const auto& y = b.alternatives[j];
        if (!paired[j] && x.sameIdentity(y)) {
          paired[j] = true;
          found = true;
          merged.guard = tSelect(c, x.guard, y.guard);
          if (x.type == Type::Boolean) {
            merged.boolean = tSelect(c, x.boolean, y.boolean);
          } else if (x.type == Type::Number) {
            merged.number = selectNumber(c, x.number, y.number);
          }
          break;
        }
      }
      if (!found) {
        merged.guard = tAnd(c, x.guard);
      }
      if (!merged.guard.isFalse()) {
        out.alternatives.push_back(std::move(merged));
      }
    }
    const Truth otherwise = tNot(c);
    for (std::size_t j = 0; j < b.alternatives.size(); ++j) {
      if (!paired[j]) {
        Alternative merged = b.alternatives[j];
        merged.guard = tAnd(otherwise, merged.guard);
        if (!merged.guard.isFalse()) {
          out.alternatives.push_back(std::move(merged));
        }
      }
    }
    if (out.alternatives.empty()) {
      return nil();
    }
    if (out.single()) {
      out.alternatives.front().guard = True;
    }
    return out;
  }

  std::vector<Value> mergeLists(Truth c, const std::vector<Value>& a, const std::vector<Value>& b) {
    std::vector<Value> out;
    const std::size_t count = std::max(a.size(), b.size());
    for (std::size_t i = 0; i < count; ++i) {
      out.push_back(merge(c, i < a.size() ? a[i] : nil(), i < b.size() ? b[i] : nil()));
    }
    return out;
  }

  /// The parts `parts` (values where their guards hold, the guards parting the points), as one.
  Value combine(const std::vector<std::pair<Truth, Value>>& parts) {
    if (parts.empty()) {
      return nil();
    }
    Value out = parts.back().second;
    for (std::size_t i = parts.size() - 1; i-- > 0;) {
      out = merge(parts[i].first, parts[i].second, out);
    }
    return out;
  }

  [[noreturn]] void refuseStale(const Alternative& alternative, int line) const {
    throw LuaUnsupported(alternative.id > 0 ? alternative.id : line,
                         alternative.text +
                             " keeps its value from one point to the next: the function writes it, "
                             "and reads it before it writes it at some points");
  }

  /// Lua's truth of `v`: all but nil and false.
  Truth truthy(const Value& v, int line) {
    Truth t = False;
    for (const auto& alternative : v.alternatives) {
      Truth here = True;
      if (alternative.type == Type::Nil) {
        here = False;
      } else if (alternative.type == Type::Boolean) {
        here = alternative.boolean;
      } else if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      }
      t = tOr(t, tAnd(alternative.guard, here));
    }
    return t;
  }

  /// The number `v` is, with an error of Lua where it is none.
  Number numberOf(const Value& v, int line, const std::string& doing) {
    const Alternative* found = nullptr;
    Truth wrong = False;
    std::string wrongType;
    for (const auto& alternative : v.alternatives) {
      if (alternative.type == Type::Number) {
        found = &alternative;
      } else if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      } else if (alternative.type == Type::String) {
        throw LuaUnsupported(line, "a string as a number is not translated");
      } else {
        wrong = tOr(wrong, alternative.guard);
        wrongType = typeName(alternative.type);
      }
    }
    if (!wrong.isFalse()) {
      luaError("attempt to " + doing + " a " + wrongType + " value", wrong, line);
    }
    return found == nullptr ? Number::of(Constant::real(Nan)) : found->number;
  }

  /// Where `a == b`.
  Truth equal(const Value& a, const Value& b, int line) {
    Truth t = False;
    for (const auto& x : a.alternatives) {
      for (const auto& y : b.alternatives) {
        if (x.type == Type::Stale) {
          refuseStale(x, line);
        }
        if (y.type == Type::Stale) {
          refuseStale(y, line);
        }
        Truth same = False;
        if (x.type == Type::Number && y.type == Type::Number) {
          same = equalNumbers(x.number, y.number);
        } else if (x.type == Type::Boolean && y.type == Type::Boolean) {
          if (x.boolean.known() && y.boolean.known()) {
            same = known(x.boolean.value == y.boolean.value);
          } else {
            same = truthOf(pw(Fn::Eq, nodeOf(x.boolean), nodeOf(y.boolean)));
          }
        } else {
          same = known(x.sameIdentity(y));
        }
        t = tOr(t, tAnd(tAnd(x.guard, y.guard), same));
      }
    }
    return t;
  }

  /// Where `a < b` (`a <= b`), with an error of Lua where they do not compare.
  Truth less(const Value& a, const Value& b, bool orEqual, int line) {
    Truth t = False;
    for (const auto& x : a.alternatives) {
      for (const auto& y : b.alternatives) {
        const Truth guard = tAnd(x.guard, y.guard);
        if (guard.isFalse()) {
          continue;
        }
        if (x.type == Type::Stale) {
          refuseStale(x, line);
        }
        if (y.type == Type::Stale) {
          refuseStale(y, line);
        }
        Truth result = False;
        if (x.type == Type::Number && y.type == Type::Number) {
          result = lessNumbers(x.number, y.number, orEqual);
        } else if (x.type == Type::String && y.type == Type::String) {
          const int order = x.text.compare(y.text);
          result = known(order < 0 || (orEqual && order == 0));
        } else {
          const std::string tx = typeName(x.type);
          const std::string ty = typeName(y.type);
          luaError(tx == ty ? "attempt to compare two " + tx + " values"
                            : "attempt to compare " + tx + " with " + ty,
                   guard,
                   line);
          continue;
        }
        t = tOr(t, tAnd(guard, result));
      }
    }
    return t;
  }

  // -------------------------------------------------------------- store ---

  int newCell(Value value, const std::string& name) {
    const int id = nextCell_++;
    state_.store.cells[id] = std::make_shared<Value>(std::move(value));
    cellNames_[id] = name;
    return id;
  }

  Value readCell(int id) const {
    const auto found = state_.store.cells.find(id);
    return found == state_.store.cells.end() ? nil() : *found->second;
  }

  void writeCell(int id, Value value) {
    if (!state_.flow.live.isTrue()) {
      value = merge(state_.flow.live, value, readCell(id));
    }
    if (recording_ && id < cellBase_) {
      leaks_.insert(Location{id, -1, Key{}});
    }
    auto& slot = state_.store.cells[id];
    if (slot != nullptr && slot.use_count() == 1) {
      *slot = std::move(value);
    } else {
      slot = std::make_shared<Value>(std::move(value));
    }
  }

  int newTable() {
    const int id = nextTable_++;
    state_.store.tables[id] = std::make_shared<Table>();
    return id;
  }

  const Table& tableAt(int id) const {
    static const Table Empty;
    const auto found = state_.store.tables.find(id);
    return found == state_.store.tables.end() ? Empty : *found->second;
  }

  Value readField(int id, const Key& key) const {
    const auto& fields = tableAt(id).fields;
    const auto found = fields.find(key);
    return found == fields.end() ? nil() : found->second;
  }

  Table& mutableTable(int id) {
    auto& slot = state_.store.tables[id];
    if (slot == nullptr) {
      slot = std::make_shared<Table>();
    } else if (slot.use_count() > 1) {
      slot = std::make_shared<Table>(*slot);
    }
    return *slot;
  }

  /// t[key] = value, without a guard: for a table no point sees yet, and the libraries.
  void rawSet(int id, const Key& key, Value value) {
    auto& fields = mutableTable(id).fields;
    if (definitelyNil(value)) {
      fields.erase(key);
    } else {
      fields[key] = std::move(value);
    }
  }

  /// t[key] = value where `guard` holds, among the points that run.
  void writeField(int id, const Key& key, Value value, Truth guard) {
    const Truth where = tAnd(state_.flow.live, guard);
    if (where.isFalse()) {
      return;
    }
    if (!where.isTrue()) {
      value = merge(where, value, readField(id, key));
    }
    if (recording_ && id < tableBase_) {
      leaks_.insert(Location{-1, id, key});
    }
    rawSet(id, key, std::move(value));
  }

  /// The key `k` is, if the script fixes it; none for nil and NaN, which no table holds.
  std::optional<Key> keyOf(const Value& k, int line) const {
    if (!k.single()) {
      throw LuaUnsupported(line, "a table index that depends on the point is not translated");
    }
    const auto& a = k.only();
    Key key;
    switch (a.type) {
    case Type::Nil:
      return std::nullopt;
    case Type::Boolean:
      if (!a.boolean.known()) {
        throw LuaUnsupported(line, "a table index that depends on the point is not translated");
      }
      key.kind = Key::Kind::Boolean;
      key.id = a.boolean.value ? 1 : 0;
      return key;
    case Type::Number: {
      if (!a.number.known) {
        throw LuaUnsupported(line, "a table index that depends on the point is not translated");
      }
      if (!a.number.k.isInt && std::isnan(a.number.k.f)) {
        return std::nullopt;
      }
      const auto whole = exactInteger(a.number.k);
      key.number = whole.has_value() ? Constant::integer(*whole) : a.number.k;
      return key;
    }
    case Type::String:
      key.kind = Key::Kind::String;
      key.text = a.text;
      return key;
    case Type::Table:
      key.kind = Key::Kind::Table;
      key.id = a.id;
      return key;
    case Type::Function:
      key.kind = Key::Kind::Function;
      key.id = a.id;
      return key;
    case Type::Builtin:
      key.kind = Key::Kind::Builtin;
      key.text = a.text;
      key.id = a.id;
      return key;
    case Type::Stale:
      refuseStale(a, line);
    }
    return std::nullopt;
  }

  /// t[k], with an error of Lua where t is no table.
  Value index(const Value& t, const Value& k, int line, const std::string& what = "") {
    std::vector<std::pair<Truth, Value>> parts;
    std::optional<Key> key;
    bool keyRead = false;
    for (const auto& alternative : t.alternatives) {
      if (alternative.type == Type::Table) {
        if (!keyRead) {
          key = keyOf(k, line);
          keyRead = true;
        }
        parts.emplace_back(alternative.guard,
                           key.has_value() ? readField(alternative.id, *key) : nil());
      } else if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      } else if (alternative.type == Type::String) {
        throw LuaUnsupported(line, "the methods of strings are not translated");
      } else {
        luaError("attempt to index a " + std::string(typeName(alternative.type)) + " value" + what,
                 alternative.guard,
                 line);
      }
    }
    return combine(parts);
  }

  /// t[k] = v.
  void setIndex(const Value& t, const Value& k, const Value& v, int line, const std::string& what) {
    std::optional<Key> key;
    bool keyRead = false;
    for (const auto& alternative : t.alternatives) {
      if (alternative.type == Type::Table) {
        if (!keyRead) {
          key = keyOf(k, line);
          keyRead = true;
          if (!key.has_value()) {
            const bool isNil = definitelyNil(k);
            luaError(isNil ? "table index is nil" : "table index is NaN", True, line);
            return;
          }
        }
        writeField(alternative.id, *key, v, alternative.guard);
      } else if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      } else {
        luaError("attempt to index a " + std::string(typeName(alternative.type)) + " value" + what,
                 alternative.guard,
                 line);
      }
    }
  }

  /// The length of the table `id` (a border): from 1 on, the keys up to it hold values, at every
  /// point.
  std::int64_t border(int id, int line) {
    const auto& fields = tableAt(id).fields;
    std::int64_t n = 0;
    while (true) {
      const auto found = fields.find(Key::integer(n + 1));
      if (found == fields.end()) {
        break;
      }
      for (const auto& alternative : found->second.alternatives) {
        if (alternative.type == Type::Nil) {
          throw LuaUnsupported(line, "the length of a table depends on the point");
        }
      }
      ++n;
    }
    // a key beyond the border makes the length Lua's choice
    for (const auto& [key, value] : fields) {
      if (key.kind == Key::Kind::Number && key.number.isInt && key.number.i > n + 1) {
        throw LuaUnsupported(line, "the length of a table with holes is not translated");
      }
    }
    return n;
  }

  /// The list of strings a table holds at the keys 1, 2, ...
  std::vector<std::string> stringList(const Value& v) {
    std::vector<std::string> out;
    if (!v.is(Type::Table)) {
      return out;
    }
    const auto count = border(v.only().id, 0);
    for (std::int64_t i = 1; i <= count; ++i) {
      const Value entry = readField(v.only().id, Key::integer(i));
      if (entry.is(Type::String)) {
        out.push_back(entry.only().text);
      }
    }
    return out;
  }

  // -------------------------------------------------------------- flows ---

  /// An error of Lua where `where` holds, among the points that run: a check of the program, or
  /// an error of the translation where every point fails. Points that fail stop.
  void luaError(const std::string& message, Truth where, int line) {
    const Truth at = tAnd(reach(), where);
    if (at.isFalse()) {
      return;
    }
    if (at.isTrue()) {
      throw LuaError(line, message);
    }
    checks_.push_back(LuaCheck{nodeOf(tNot(at)), origin_.at(line) + ": " + message});
    if (where.isTrue()) {
      state_.flow.returned = tOr(state_.flow.returned, state_.flow.live);
      state_.flow.live = False;
    }
  }

  Store mergeStores(Truth c, const Store& a, const Store& b) {
    Store out;
    auto ia = a.cells.begin();
    auto ib = b.cells.begin();
    while (ia != a.cells.end() || ib != b.cells.end()) {
      if (ib == b.cells.end() || (ia != a.cells.end() && ia->first < ib->first)) {
        out.cells.emplace_hint(out.cells.end(), *ia++);
      } else if (ia == a.cells.end() || ib->first < ia->first) {
        out.cells.emplace_hint(out.cells.end(), *ib++);
      } else {
        if (ia->second == ib->second || identical(*ia->second, *ib->second)) {
          out.cells.emplace_hint(out.cells.end(), *ia);
        } else {
          out.cells.emplace_hint(out.cells.end(),
                                 ia->first,
                                 std::make_shared<Value>(merge(c, *ia->second, *ib->second)));
        }
        ++ia;
        ++ib;
      }
    }
    auto ta = a.tables.begin();
    auto tb = b.tables.begin();
    while (ta != a.tables.end() || tb != b.tables.end()) {
      if (tb == b.tables.end() || (ta != a.tables.end() && ta->first < tb->first)) {
        out.tables.emplace_hint(out.tables.end(), *ta++);
      } else if (ta == a.tables.end() || tb->first < ta->first) {
        out.tables.emplace_hint(out.tables.end(), *tb++);
      } else {
        if (ta->second == tb->second) {
          out.tables.emplace_hint(out.tables.end(), *ta);
        } else {
          auto merged = std::make_shared<Table>();
          const auto& fa = ta->second->fields;
          const auto& fb = tb->second->fields;
          auto xa = fa.begin();
          auto xb = fb.begin();
          const auto put = [&](const Key& key, Value value) {
            if (!definitelyNil(value)) {
              merged->fields.emplace_hint(merged->fields.end(), key, std::move(value));
            }
          };
          while (xa != fa.end() || xb != fb.end()) {
            if (xb == fb.end() || (xa != fa.end() && xa->first < xb->first)) {
              put(xa->first, merge(c, xa->second, nil()));
              ++xa;
            } else if (xa == fa.end() || xb->first < xa->first) {
              put(xb->first, merge(c, nil(), xb->second));
              ++xb;
            } else {
              put(xa->first, merge(c, xa->second, xb->second));
              ++xa;
              ++xb;
            }
          }
          out.tables.emplace_hint(out.tables.end(), ta->first, std::move(merged));
        }
        ++ta;
        ++tb;
      }
    }
    return out;
  }

  /// The flow of `a` where `c` holds and of `b` elsewhere; `version` that of the results where the
  /// two parted.
  Flow mergeFlows(Truth c, const Flow& a, const Flow& b, int version) {
    Flow out;
    out.live = tSelect(c, a.live, b.live);
    out.returned = tSelect(c, a.returned, b.returned);
    out.broken = tSelect(c, a.broken, b.broken);
    if (a.version == b.version || !b.hasResults || b.version == version) {
      out.hasResults = a.hasResults;
      out.results = a.results;
      out.version = a.version;
    } else if (!a.hasResults || a.version == version) {
      out.hasResults = b.hasResults;
      out.results = b.results;
      out.version = b.version;
    } else {
      out.hasResults = true;
      out.results = mergeLists(c, a.results, b.results);
      out.version = ++nextVersion_;
    }
    return out;
  }

  /// The points of the pendings of `state` after the first `depth` back into it.
  void resolvePendings(State& state, std::size_t depth) {
    while (state.pendings.size() > depth) {
      const auto pending = state.pendings.back();
      state.pendings.pop_back();
      const Truth c = pending->condition;
      if (pending->stoppedOnTrue) {
        state.store = mergeStores(c, pending->store, state.store);
        state.flow = mergeFlows(c, pending->flow, state.flow, pending->version);
      } else {
        state.store = mergeStores(c, state.store, pending->store);
        state.flow = mergeFlows(c, state.flow, pending->flow, pending->version);
      }
      state.domain = pending->domain;
      state.facts.resize(pending->facts);
    }
  }
  void resolvePendings(std::size_t depth) { resolvePendings(state_, depth); }

  /// The pendings of `side` after the first `depth`, for the points of both sides of a branch:
  /// their stopped points are those where the side runs (`where`) as well. A pending further in
  /// may name points of one further out, which that one, resolved later, takes back.
  std::vector<std::shared_ptr<const Pending>>
      lifted(State& side, std::size_t depth, Truth where, Truth domain, std::size_t facts) {
    std::vector<std::shared_ptr<const Pending>> out;
    for (std::size_t i = depth; i < side.pendings.size(); ++i) {
      const auto& inner = *side.pendings[i];
      auto pending = std::make_shared<Pending>(inner);
      pending->condition =
          tAnd(where, inner.stoppedOnTrue ? inner.condition : tNot(inner.condition));
      pending->stoppedOnTrue = true;
      pending->domain = domain;
      pending->facts = facts;
      pending->version = -1;
      out.push_back(std::move(pending));
    }
    side.pendings.resize(depth);
    return out;
  }

  /// Runs `onTrue` where `c` holds and `onFalse` where not, each on a state of its own. A side
  /// that stops every point it runs (by a return, a break or an error) waits in a pending, and
  /// the evaluation goes on with the other; else the two are merged. What stops within a side
  /// that goes on waits as well, so that the state is of the points that run.
  template <typename OnTrue, typename OnFalse>
  void branch(Truth c, OnTrue&& onTrue, OnFalse&& onFalse) {
    const std::size_t depth = state_.pendings.size();
    const std::size_t facts = state_.facts.size();
    const Truth domain = state_.domain;
    const int version = state_.flow.version;
    State before = state_;

    state_.domain = tAnd(domain, c);
    state_.facts.emplace_back(c.node, true);
    onTrue();
    State whenTrue = std::move(state_);

    state_ = std::move(before);
    state_.domain = tAnd(domain, tNot(c));
    state_.facts.emplace_back(c.node, false);
    onFalse();

    const bool trueStopped = whenTrue.flow.live.isFalse();
    const bool falseStopped = state_.flow.live.isFalse();
    if (trueStopped) {
      resolvePendings(whenTrue, depth);
    }
    if (falseStopped) {
      resolvePendings(state_, depth);
    }
    if (trueStopped != falseStopped) {
      // the stopped side waits; the other goes on, in its own domain
      State& going = trueStopped ? state_ : whenTrue;
      State& stopped = trueStopped ? whenTrue : state_;
      auto pending = std::make_shared<Pending>();
      pending->condition = c;
      pending->stoppedOnTrue = trueStopped;
      pending->store = std::move(stopped.store);
      pending->flow = std::move(stopped.flow);
      pending->domain = domain;
      pending->facts = facts;
      pending->version = version;
      auto inner = lifted(going, depth, trueStopped ? tNot(c) : c, domain, facts);
      going.pendings.push_back(std::move(pending));
      for (auto& p : inner) {
        going.pendings.push_back(std::move(p));
      }
      if (!trueStopped) {
        state_ = std::move(whenTrue);
      }
      return;
    }
    auto innerTrue = trueStopped ? std::vector<std::shared_ptr<const Pending>>{}
                                 : lifted(whenTrue, depth, c, domain, facts);
    auto innerFalse = falseStopped ? std::vector<std::shared_ptr<const Pending>>{}
                                   : lifted(state_, depth, tNot(c), domain, facts);
    state_.store = mergeStores(c, whenTrue.store, state_.store);
    state_.flow = mergeFlows(c, whenTrue.flow, state_.flow, version);
    state_.domain = domain;
    state_.facts.resize(facts);
    for (auto& p : innerTrue) {
      state_.pendings.push_back(std::move(p));
    }
    for (auto& p : innerFalse) {
      state_.pendings.push_back(std::move(p));
    }
    separateStopped();
  }

  /// Where some points of the state stopped and others run on (after a loop, say, that some
  /// left by a return): the stopped ones wait in a pending, so that the state is of the points
  /// that run, and an assignment needs no select on the flow.
  void separateStopped() {
    const Truth live = state_.flow.live;
    if (live.known()) {
      return;
    }
    auto pending = std::make_shared<Pending>();
    pending->condition = tNot(live);
    pending->stoppedOnTrue = true;
    pending->store = state_.store;
    pending->flow = state_.flow;
    pending->flow.live = False;
    pending->domain = state_.domain;
    pending->facts = state_.facts.size();
    pending->version = -1;
    state_.pendings.push_back(std::move(pending));
    state_.domain = tAnd(state_.domain, live);
    state_.facts.emplace_back(live.node, true);
    state_.flow = Flow{};
    state_.flow.version = ++nextVersion_;
  }

  void step(int line) {
    if (++steps_ > options_.stepBudget) {
      throw LuaUnsupported(line,
                           "the evaluation takes more than " + std::to_string(options_.stepBudget) +
                               " steps; a loop that long is not unrolled");
    }
    if (arena_.size() > options_.nodeBudget) {
      throw LuaUnsupported(
          line, "the program grows beyond " + std::to_string(options_.nodeBudget) + " nodes");
    }
  }

  // ------------------------------------------------------------- scopes ---

  void declare(Scope& scope, const std::string& name, Value value, bool constant) {
    scope.names.emplace_back(name, newCell(std::move(value), name));
    scope.constants.push_back(constant);
    scope.compileTime.emplace_back();
  }

  /// The value Lua's compiler gives the local `name`, a <const> one, if it gives it one.
  static std::optional<Constant> compileTimeLocal(const Scope* scope, const std::string& name) {
    std::size_t limit = scope != nullptr ? scope->names.size() : 0;
    for (; scope != nullptr; scope = scope->parent.get()) {
      for (std::size_t i = std::min(limit, scope->names.size()); i-- > 0;) {
        if (scope->names[i].first == name) {
          return scope->compileTime[i];
        }
      }
      limit = scope->parentVisible;
    }
    return std::nullopt;
  }

  /// The number `e` is where Lua's compiler computes it (constfolding of lcode.c): a numeral, and
  /// an operation on such that gives an integer, or a float other than zero and NaN; through
  /// parentheses and <const> locals.
  std::optional<Constant> compileTimeConstant(const Expression& e, const ScopePtr& scope) {
    const auto folded = [](const Constant& c) -> std::optional<Constant> {
      if (!c.isInt && (std::isnan(c.f) || c.f == 0.0)) {
        return std::nullopt;
      }
      return c;
    };
    switch (e.kind) {
    case ExprKind::Number:
      return e.isInteger ? Constant::integer(e.integer) : Constant::real(e.number);
    case ExprKind::Paren:
      return compileTimeConstant(*e.a, scope);
    case ExprKind::Name:
      return compileTimeLocal(scope.get(), e.text);
    case ExprKind::Unary: {
      if (e.op != Op::Negate && e.op != Op::BitNot) {
        return std::nullopt;
      }
      const auto a = compileTimeConstant(*e.a, scope);
      if (!a.has_value()) {
        return std::nullopt;
      }
      if (e.op == Op::Negate) {
        return folded(a->isInt ? Constant::integer(wrapped(0U - static_cast<std::uint64_t>(a->i)))
                               : Constant::real(-a->f));
      }
      const auto whole = exactInteger(*a);
      return whole.has_value() ? std::optional<Constant>(Constant::integer(~*whole)) : std::nullopt;
    }
    case ExprKind::Binary: {
      const bool bitwise = e.op == Op::BitAnd || e.op == Op::BitOr || e.op == Op::BitXor ||
                           e.op == Op::Shl || e.op == Op::Shr;
      const bool arithmetic = e.op == Op::Add || e.op == Op::Sub || e.op == Op::Mul ||
                              e.op == Op::Div || e.op == Op::IntDiv || e.op == Op::Mod ||
                              e.op == Op::Pow;
      if (!bitwise && !arithmetic) {
        return std::nullopt;
      }
      const auto a = compileTimeConstant(*e.a, scope);
      const auto b = compileTimeConstant(*e.b, scope);
      if (!a.has_value() || !b.has_value()) {
        return std::nullopt;
      }
      // validop: no conversion errors, no division by zero
      if (bitwise && (!exactInteger(*a).has_value() || !exactInteger(*b).has_value())) {
        return std::nullopt;
      }
      if ((e.op == Op::Div || e.op == Op::IntDiv || e.op == Op::Mod) && b->asFloat() == 0.0) {
        return std::nullopt;
      }
      return folded(foldArithmetic(e.op, *a, *b, e.line));
    }
    default:
      return std::nullopt;
    }
  }

  /// The cell of the local `name` visible in `scope`, and whether it is <const>.
  static std::optional<std::pair<int, bool>> lookup(const Scope* scope, const std::string& name) {
    std::size_t limit = scope != nullptr ? scope->names.size() : 0;
    for (; scope != nullptr; scope = scope->parent.get()) {
      for (std::size_t i = std::min(limit, scope->names.size()); i-- > 0;) {
        if (scope->names[i].first == name) {
          return std::make_pair(scope->names[i].second, static_cast<bool>(scope->constants[i]));
        }
      }
      limit = scope->parentVisible;
    }
    return std::nullopt;
  }

  /// The cells of a scope no closure holds, given up at its end.
  void release(const ScopePtr& scope) {
    if (scope->captured) {
      return;
    }
    for (const auto& [name, cell] : scope->names) {
      state_.store.cells.erase(cell);
    }
  }

  Value readName(const std::string& name, const ScopePtr& scope, int line) {
    if (const auto local = lookup(scope.get(), name)) {
      return readAt(readCell(local->first), line);
    }
    return readAt(
        index(readName("_ENV", scope, line), string(name), line, " (global '" + name + "')"), line);
  }

  /// `v`, what of it is of an earlier point marked as read at `line`, where it was not read
  /// before: a refusal names the line that reads it, rather than the one that uses it.
  static Value readAt(Value v, int line) {
    for (auto& alternative : v.alternatives) {
      if (alternative.type == Type::Stale && alternative.id < 0) {
        alternative.id = line;
      }
    }
    return v;
  }

  void assignName(const std::string& name, const Value& value, const ScopePtr& scope, int line) {
    if (const auto local = lookup(scope.get(), name)) {
      if (local->second) {
        throw LuaError(line, "attempt to assign to const variable '" + name + "'");
      }
      writeCell(local->first, value);
      return;
    }
    setIndex(readName("_ENV", scope, line), string(name), value, line, " (global '" + name + "')");
  }

  // --------------------------------------------------------- statements ---

  struct Activation {
    const FunctionBody* body{nullptr};
    std::vector<Value> varargs;
    int loops{0};
  };

  void execBlock(const Block& block, const ScopePtr& parent) {
    auto scope = childOf(parent);
    execStatements(block, scope);
    release(scope);
  }

  void execStatements(const Block& block, const ScopePtr& scope) {
    for (const auto& statement : block.statements) {
      if (state_.flow.live.isFalse()) {
        return;
      }
      exec(statement, scope);
    }
  }

  void exec(const Statement& s, const ScopePtr& scope) {
    step(s.line);
    switch (s.kind) {
    case StatKind::Empty:
    case StatKind::Label:
      return;
    case StatKind::Local: {
      auto values = evalList(s.values, scope);
      // localstat of lparser.c: the last variable, <const> and of as many as there are
      // expressions, is a constant of the compiler where its expression is one
      const std::optional<Constant> compileTime =
          s.attributes.back() == "const" && s.names.size() == s.values.size()
              ? compileTimeConstant(*s.values.back(), scope)
              : std::nullopt;
      for (std::size_t i = 0; i < s.names.size(); ++i) {
        const auto& attribute = s.attributes[i];
        if (!attribute.empty() && attribute != "const") {
          if (attribute == "close") {
            throw LuaUnsupported(s.line, "a to-be-closed variable is not translated");
          }
          throw LuaError(s.line, "unknown attribute '" + attribute + "'");
        }
        declare(*scope, s.names[i], i < values.size() ? values[i] : nil(), attribute == "const");
      }
      scope->compileTime.back() = compileTime;
      return;
    }
    case StatKind::LocalFunction: {
      // `local function f` is `local f; f = function`: the function sees itself
      declare(*scope, s.names.front(), nil(), false);
      const int cell = scope->names.back().second;
      state_.store.cells[cell] = std::make_shared<Value>(closureOf(s.function, scope));
      return;
    }
    case StatKind::Function: {
      // `function a.b.c:m()`: a is a variable, the rest fields
      Value value = closureOf(s.function, scope);
      if (s.path.size() == 1 && s.method.empty()) {
        assignName(s.path.front(), value, scope, s.line);
        return;
      }
      Value target = readName(s.path.front(), scope, s.line);
      const std::size_t last = s.method.empty() ? s.path.size() - 1 : s.path.size();
      for (std::size_t i = 1; i < last; ++i) {
        target = index(target, string(s.path[i]), s.line, " (field '" + s.path[i] + "')");
      }
      const std::string field = s.method.empty() ? s.path.back() : s.method;
      setIndex(target, string(field), value, s.line, "");
      return;
    }
    case StatKind::Assign: {
      // the tables and keys of the targets, then the values, then the assignments
      struct Target {
        bool name{false};
        Value table;
        Value key;
      };
      std::vector<Target> targets;
      for (const auto& target : s.targets) {
        if (target->kind == ExprKind::Name) {
          targets.push_back(Target{true, nil(), string(target->text)});
        } else {
          Value t = eval(*target->a, scope);
          Value k = eval(*target->b, scope);
          targets.push_back(Target{false, std::move(t), std::move(k)});
        }
      }
      auto values = evalList(s.values, scope);
      for (std::size_t i = targets.size(); i-- > 0;) {
        const Value value = i < values.size() ? values[i] : nil();
        if (targets[i].name) {
          assignName(targets[i].key.only().text, value, scope, s.line);
        } else {
          setIndex(
              targets[i].table, targets[i].key, value, s.line, variableOf(*s.targets[i]->a, scope));
        }
      }
      return;
    }
    case StatKind::Call:
      static_cast<void>(callExpression(*s.call, scope));
      return;
    case StatKind::Do:
      execBlock(*s.body, scope);
      return;
    case StatKind::If:
      ifChain(s, 0, scope);
      return;
    case StatKind::While:
      whileLoop(s, scope);
      return;
    case StatKind::Repeat:
      repeatLoop(s, scope);
      return;
    case StatKind::NumericFor:
      numericFor(s, scope);
      return;
    case StatKind::GenericFor:
      genericFor(s, scope);
      return;
    case StatKind::Return: {
      auto values = evalList(s.values, scope);
      if (state_.flow.live.isFalse()) {
        return;
      }
      if (calls_.size() == entryDepth_) {
        returnLines_.push_back(s.line);
      }
      auto& flow = state_.flow;
      if (!flow.hasResults || flow.returned.isFalse()) {
        flow.results = std::move(values);
      } else {
        flow.results = mergeLists(flow.returned, flow.results, values);
      }
      flow.hasResults = true;
      flow.version = ++nextVersion_;
      flow.returned = tOr(flow.returned, flow.live);
      flow.live = False;
      return;
    }
    case StatKind::Break:
      if (calls_.back().loops == 0) {
        throw LuaError(s.line, "break outside a loop");
      }
      state_.flow.broken = tOr(state_.flow.broken, state_.flow.live);
      state_.flow.live = False;
      return;
    case StatKind::Goto:
      throw LuaUnsupported(s.line, "goto is not translated");
    }
  }

  void ifChain(const Statement& s, std::size_t k, const ScopePtr& scope) {
    if (k == s.conditions.size()) {
      if (s.otherwise != nullptr) {
        execBlock(*s.otherwise, scope);
      }
      return;
    }
    const Value condition = eval(*s.conditions[k], scope);
    const Truth t = decide(truthy(condition, s.conditions[k]->line));
    if (t.isTrue()) {
      execBlock(*s.blocks[k], scope);
      return;
    }
    if (t.isFalse()) {
      ifChain(s, k + 1, scope);
      return;
    }
    branch(t, [&] { execBlock(*s.blocks[k], scope); }, [&] { ifChain(s, k + 1, scope); });
  }

  /// A loop: `iteration` runs one, and says whether there is one more. Breaks are Lua's; the
  /// points that break wait until the end of the loop.
  template <typename Iteration>
  void loop(int line, Iteration&& iteration) {
    const std::size_t depth = state_.pendings.size();
    const Truth outerBroken = state_.flow.broken;
    const Truth liveBefore = state_.flow.live;
    state_.flow.broken = False;
    ++calls_.back().loops;
    while (!state_.flow.live.isFalse()) {
      step(line);
      if (!iteration()) {
        break;
      }
    }
    --calls_.back().loops;
    resolvePendings(depth);
    // the points that broke run on; those that returned or failed do not
    state_.flow.live = tAnd(liveBefore, tNot(state_.flow.returned));
    state_.flow.broken = outerBroken;
    separateStopped();
  }

  /// The truth of a loop's condition, which the script has to fix.
  bool loopCondition(const Expression& condition, const ScopePtr& scope) {
    const Truth t = decide(truthy(eval(condition, scope), condition.line));
    if (!t.known()) {
      throw LuaUnsupported(condition.line,
                           "a loop whose end depends on the point is not translated");
    }
    return t.value;
  }

  void whileLoop(const Statement& s, const ScopePtr& scope) {
    loop(s.line, [&] {
      if (!loopCondition(*s.conditions.front(), scope)) {
        return false;
      }
      execBlock(*s.body, scope);
      return true;
    });
  }

  void repeatLoop(const Statement& s, const ScopePtr& scope) {
    loop(s.line, [&] {
      // the condition sees the locals of the body
      auto inner = childOf(scope);
      execStatements(*s.body, inner);
      bool more = false;
      if (!state_.flow.live.isFalse()) {
        more = !loopCondition(*s.conditions.front(), inner);
      }
      release(inner);
      return more;
    });
  }

  Constant loopBound(const Expression& e, const ScopePtr& scope, const char* what) {
    const Value v = eval(e, scope);
    if (!v.is(Type::Number)) {
      if (v.single() && v.only().type != Type::Number) {
        throw LuaError(e.line, std::string("'for' ") + what + " must be a number");
      }
      throw LuaUnsupported(e.line, "a loop whose bounds depend on the point is not translated");
    }
    if (!v.only().number.known) {
      throw LuaUnsupported(e.line, "a loop whose bounds depend on the point is not translated");
    }
    return v.only().number.k;
  }

  void numericFor(const Statement& s, const ScopePtr& scope) {
    const Constant start = loopBound(*s.start, scope, "initial value");
    const Constant stop = loopBound(*s.stop, scope, "limit");
    const Constant stepValue =
        s.step != nullptr ? loopBound(*s.step, scope, "step") : Constant::integer(1);
    const auto body = [&](const Constant& value) {
      auto inner = childOf(scope);
      declare(*inner, s.names.front(), number(value), false);
      execStatements(*s.body, inner);
      release(inner);
    };
    if (start.isInt && stepValue.isInt) {
      // forprep of lvm.c: the limit as an integer, the count of iterations up front
      const std::int64_t init = start.i;
      const std::int64_t increment = stepValue.i;
      if (increment == 0) {
        throw LuaError(s.line, "'for' step is zero");
      }
      std::int64_t limit = 0;
      if (stop.isInt) {
        limit = stop.i;
      } else {
        const double f = increment < 0 ? std::ceil(stop.f) : std::floor(stop.f);
        const auto whole = exactInteger(f);
        if (whole.has_value()) {
          limit = *whole;
        } else if (std::isnan(stop.f)) {
          return;
        } else if (stop.f > 0) {
          if (increment < 0) {
            return;
          }
          limit = std::numeric_limits<std::int64_t>::max();
        } else {
          if (increment > 0) {
            return;
          }
          limit = std::numeric_limits<std::int64_t>::min();
        }
      }
      if (increment > 0 ? init > limit : init < limit) {
        return;
      }
      const std::uint64_t count =
          increment > 0 ? (static_cast<std::uint64_t>(limit) - static_cast<std::uint64_t>(init)) /
                              static_cast<std::uint64_t>(increment)
                        : (static_cast<std::uint64_t>(init) - static_cast<std::uint64_t>(limit)) /
                              (static_cast<std::uint64_t>(-(increment + 1)) + 1U);
      std::uint64_t done = 0;
      std::int64_t value = init;
      loop(s.line, [&] {
        body(Constant::integer(value));
        if (done++ == count) {
          return false;
        }
        value = wrapped(static_cast<std::uint64_t>(value) + static_cast<std::uint64_t>(increment));
        return true;
      });
      return;
    }
    const double init = start.asFloat();
    const double limit = stop.asFloat();
    const double increment = stepValue.asFloat();
    if (increment == 0.0) {
      throw LuaError(s.line, "'for' step is zero");
    }
    double value = init;
    if (!(increment > 0 ? value <= limit : value >= limit)) {
      return;
    }
    loop(s.line, [&] {
      body(Constant::real(value));
      value += increment;
      return increment > 0 ? value <= limit : value >= limit;
    });
  }

  void genericFor(const Statement& s, const ScopePtr& scope) {
    auto values = evalList(s.values, scope);
    values.resize(std::max<std::size_t>(values.size(), 3), nil());
    const Value iterator = values[0];
    const Value invariant = values[1];
    Value control = values[2];
    loop(s.line, [&] {
      auto results = call(iterator, {invariant, control}, s.line, " (for iterator)");
      if (state_.flow.live.isFalse()) {
        return false;
      }
      results.resize(std::max(results.size(), s.names.size()), nil());
      const Value& first = results.front();
      if (definitelyNil(first)) {
        return false;
      }
      for (const auto& alternative : first.alternatives) {
        if (alternative.type == Type::Nil) {
          throw LuaUnsupported(s.line, "a loop whose end depends on the point is not translated");
        }
      }
      control = first;
      auto inner = childOf(scope);
      for (std::size_t i = 0; i < s.names.size(); ++i) {
        declare(*inner, s.names[i], results[i], false);
      }
      execStatements(*s.body, inner);
      release(inner);
      return true;
    });
  }

  // -------------------------------------------------------- expressions ---

  static bool isMulti(const Expression& e) {
    return e.kind == ExprKind::Call || e.kind == ExprKind::Method || e.kind == ExprKind::Vararg;
  }

  /// Whether `e` cannot write: no calls, no functions, no tables.
  bool pure(const Expression& e) {
    const auto found = purity_.find(&e);
    if (found != purity_.end()) {
      return found->second;
    }
    bool result = true;
    switch (e.kind) {
    case ExprKind::Call:
    case ExprKind::Method:
    case ExprKind::Function:
    case ExprKind::Table:
      result = false;
      break;
    default:
      result = (e.a == nullptr || pure(*e.a)) && (e.b == nullptr || pure(*e.b));
    }
    purity_[&e] = result;
    return result;
  }

  std::string variableOf(const Expression& e, const ScopePtr& scope) const {
    if (e.kind == ExprKind::Name) {
      return lookup(scope.get(), e.text).has_value() ? " (local '" + e.text + "')"
                                                     : " (global '" + e.text + "')";
    }
    if (e.kind == ExprKind::Index && e.b->kind == ExprKind::String) {
      return " (field '" + e.b->text + "')";
    }
    if (e.kind == ExprKind::Method) {
      return " (method '" + e.text + "')";
    }
    return "";
  }

  Value closureOf(const std::shared_ptr<FunctionBody>& body, const ScopePtr& scope) {
    for (Scope* s = scope.get(); s != nullptr; s = s->parent.get()) {
      s->captured = true;
    }
    auto closure = std::make_shared<Closure>();
    closure->body = body;
    closure->scope = scope;
    closure->visible = scope->names.size();
    closure->id = nextClosure_++;
    return function(std::move(closure));
  }

  const std::vector<Value>& varargs(int line) const {
    const auto& activation = calls_.back();
    if (activation.body == nullptr || !activation.body->vararg) {
      throw LuaError(line, "cannot use '...' outside a vararg function");
    }
    return activation.varargs;
  }

  Value eval(const Expression& e, const ScopePtr& scope) {
    switch (e.kind) {
    case ExprKind::Nil:
      return nil();
    case ExprKind::True:
      return boolean(True);
    case ExprKind::False:
      return boolean(False);
    case ExprKind::Number:
      return number(e.isInteger ? Constant::integer(e.integer) : Constant::real(e.number));
    case ExprKind::String:
      return string(e.text);
    case ExprKind::Vararg: {
      const auto& values = varargs(e.line);
      return values.empty() ? nil() : values.front();
    }
    case ExprKind::Function:
      return closureOf(e.function, scope);
    case ExprKind::Table:
      return construct(e, scope);
    case ExprKind::Binary:
      return binary(e, scope);
    case ExprKind::Unary:
      return unary(e, scope);
    case ExprKind::Name:
      return readName(e.text, scope, e.line);
    case ExprKind::Index: {
      const Value t = eval(*e.a, scope);
      const Value k = eval(*e.b, scope);
      return readAt(index(t, k, e.line, variableOf(*e.a, scope)), e.line);
    }
    case ExprKind::Call:
    case ExprKind::Method: {
      auto values = callExpression(e, scope);
      return values.empty() ? nil() : values.front();
    }
    case ExprKind::Paren:
      return eval(*e.a, scope);
    }
    return nil();
  }

  std::vector<Value> evalMulti(const Expression& e, const ScopePtr& scope) {
    if (e.kind == ExprKind::Call || e.kind == ExprKind::Method) {
      return callExpression(e, scope);
    }
    if (e.kind == ExprKind::Vararg) {
      return varargs(e.line);
    }
    return {eval(e, scope)};
  }

  std::vector<Value> evalList(const std::vector<ExpressionPtr>& list, const ScopePtr& scope) {
    std::vector<Value> out;
    for (std::size_t i = 0; i < list.size(); ++i) {
      if (i + 1 == list.size()) {
        auto rest = evalMulti(*list[i], scope);
        for (auto& value : rest) {
          out.push_back(std::move(value));
        }
      } else {
        out.push_back(eval(*list[i], scope));
      }
    }
    return out;
  }

  Value construct(const Expression& e, const ScopePtr& scope) {
    const int id = newTable();
    std::vector<Value> positional;
    for (std::size_t i = 0; i < e.fields.size(); ++i) {
      const auto& field = e.fields[i];
      if (field.key != nullptr) {
        const Value k = eval(*field.key, scope);
        const Value v = eval(*field.value, scope);
        const auto key = keyOf(k, e.line);
        if (!key.has_value()) {
          luaError(definitelyNil(k) ? "table index is nil" : "table index is NaN", True, e.line);
          continue;
        }
        rawSet(id, *key, v);
      } else if (i + 1 == e.fields.size() && isMulti(*field.value)) {
        auto values = evalMulti(*field.value, scope);
        for (auto& value : values) {
          positional.push_back(std::move(value));
        }
      } else {
        positional.push_back(eval(*field.value, scope));
      }
    }
    for (std::size_t i = 0; i < positional.size(); ++i) {
      rawSet(id, Key::integer(static_cast<std::int64_t>(i + 1)), positional[i]);
    }
    return table(id);
  }

  /// `a and b`, `a or b`: b only where a does not decide.
  Value logical(const Expression& e, const ScopePtr& scope, bool isAnd) {
    const Value a = eval(*e.a, scope);
    const Truth t = decide(truthy(a, e.line));
    if (t.known()) {
      return t.value == isAnd ? eval(*e.b, scope) : a;
    }
    const Truth when = isAnd ? t : tNot(t);
    Value b;
    if (pure(*e.b)) {
      // no writes: the store stays; failures are checks where b runs
      const Truth domain = state_.domain;
      const Flow flow = state_.flow;
      state_.domain = tAnd(domain, when);
      state_.facts.emplace_back(when.node, true);
      b = eval(*e.b, scope);
      state_.facts.pop_back();
      state_.domain = domain;
      state_.flow = flow;
    } else {
      branch(when, [&] { b = eval(*e.b, scope); }, [] {});
    }
    return merge(when, b, a);
  }

  Value binary(const Expression& e, const ScopePtr& scope) {
    if (e.op == Op::And || e.op == Op::Or) {
      return logical(e, scope, e.op == Op::And);
    }
    const Value a = eval(*e.a, scope);
    const Value b = eval(*e.b, scope);
    if (e.op == Op::Sub) {
      // Lua's code generator computes `a - k`, for an integer constant k from -127 to 127, as
      // `a + -k` (OP_ADDI): the same, but at k = 0, where -0.0 - 0 is -0.0 and -0.0 + 0 is 0.0
      const auto k = compileTimeConstant(*e.b, scope);
      if (k.has_value() && k->isInt && k->i >= -127 && k->i <= 127) {
        return number(arithmetic(Op::Add,
                                 numberOf(a, e.line, "perform arithmetic on"),
                                 Number::of(Constant::integer(-k->i)),
                                 e.line));
      }
    }
    switch (e.op) {
    case Op::Eq:
      return boolean(equal(a, b, e.line));
    case Op::Ne:
      return boolean(tNot(equal(a, b, e.line)));
    case Op::Lt:
      return boolean(less(a, b, false, e.line));
    case Op::Le:
      return boolean(less(a, b, true, e.line));
    case Op::Gt:
      return boolean(less(b, a, false, e.line));
    case Op::Ge:
      return boolean(less(b, a, true, e.line));
    case Op::Concat:
      return string(text(a, e.line) + text(b, e.line));
    case Op::BitAnd:
    case Op::BitOr:
    case Op::BitXor:
    case Op::Shl:
    case Op::Shr:
      return number(arithmetic(e.op,
                               numberOf(a, e.line, "perform bitwise operation on"),
                               numberOf(b, e.line, "perform bitwise operation on"),
                               e.line));
    default:
      return number(arithmetic(e.op,
                               numberOf(a, e.line, "perform arithmetic on"),
                               numberOf(b, e.line, "perform arithmetic on"),
                               e.line));
    }
  }

  /// The text of a string or a number the script fixes, for `..`.
  std::string text(const Value& v, int line) {
    if (v.is(Type::String)) {
      return v.only().text;
    }
    if (v.is(Type::Number) && v.only().number.known) {
      return numberText(v.only().number.k);
    }
    for (const auto& alternative : v.alternatives) {
      if (alternative.type == Type::Number || alternative.type == Type::String) {
        throw LuaUnsupported(line, "a string that depends on the point is not translated");
      }
      if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      }
    }
    luaError("attempt to concatenate a " + std::string(typeName(v.alternatives.front().type)) +
                 " value",
             True,
             line);
    return "";
  }

  Value unary(const Expression& e, const ScopePtr& scope) {
    const Value a = eval(*e.a, scope);
    switch (e.op) {
    case Op::Not:
      return boolean(tNot(truthy(a, e.line)));
    case Op::Negate:
      return number(negate(numberOf(a, e.line, "perform arithmetic on")));
    case Op::BitNot: {
      const Number n = numberOf(a, e.line, "perform bitwise operation on");
      if (!n.known) {
        throw LuaUnsupported(e.line,
                             "a bitwise operation on a value of the point is not translated");
      }
      const auto whole = exactInteger(n.k);
      if (!whole.has_value()) {
        luaError("number has no integer representation", True, e.line);
        return number(Constant::integer(0));
      }
      return number(Constant::integer(~*whole));
    }
    case Op::Length:
      return length(a, e.line);
    default:
      throw std::logic_error("not a unary operator");
    }
  }

  Value length(const Value& v, int line) {
    if (v.is(Type::String)) {
      return number(Constant::integer(static_cast<std::int64_t>(v.only().text.size())));
    }
    if (v.is(Type::Table)) {
      return number(Constant::integer(border(v.only().id, line)));
    }
    for (const auto& alternative : v.alternatives) {
      if (alternative.type == Type::Table || alternative.type == Type::String) {
        throw LuaUnsupported(line, "the length of a value that depends on the point");
      }
      if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      }
    }
    luaError("attempt to get length of a " + std::string(typeName(v.alternatives.front().type)) +
                 " value",
             True,
             line);
    return number(Constant::integer(0));
  }

  // -------------------------------------------------------------- calls ---

  std::vector<Value> callExpression(const Expression& e, const ScopePtr& scope) {
    Value callee;
    std::vector<Value> args;
    std::string what;
    if (e.kind == ExprKind::Method) {
      const Value object = eval(*e.a, scope);
      callee = index(object, string(e.text), e.line, " (method '" + e.text + "')");
      args.push_back(object);
      what = " (method '" + e.text + "')";
    } else {
      callee = eval(*e.a, scope);
      what = variableOf(*e.a, scope);
    }
    // print and io.write do nothing, with their arguments
    if (callee.is(Type::Builtin) &&
        (callee.only().text == "print" || callee.only().text == "io.write")) {
      return {};
    }
    auto rest = evalList(e.args, scope);
    for (auto& value : rest) {
      args.push_back(std::move(value));
    }
    return call(callee, std::move(args), e.line, what);
  }

  std::vector<Value>
      call(const Value& callee, std::vector<Value> args, int line, const std::string& what) {
    const Alternative* callable = nullptr;
    Truth wrong = False;
    std::string wrongType;
    for (const auto& alternative : callee.alternatives) {
      if (alternative.type == Type::Function || alternative.type == Type::Builtin) {
        if (callable != nullptr) {
          throw LuaUnsupported(line, "a function that depends on the point is not translated");
        }
        callable = &alternative;
      } else if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      } else {
        wrong = tOr(wrong, alternative.guard);
        wrongType = typeName(alternative.type);
      }
    }
    if (!wrong.isFalse()) {
      luaError("attempt to call a " + wrongType + " value" + what, wrong, line);
    }
    if (callable == nullptr || state_.flow.live.isFalse()) {
      return {};
    }
    if (callable->type == Type::Function) {
      return callClosure(*callable->closure, std::move(args), line);
    }
    return callBuiltin(callable->text, callable->id, args, line);
  }

  std::vector<Value> callClosure(const Closure& closure, std::vector<Value> args, int line) {
    const FunctionBody* body = closure.body.get();
    for (const auto& activation : calls_) {
      if (activation.body == body) {
        throw LuaUnsupported(line,
                             "recursion is not translated: " +
                                 (body->name.empty() ? std::string("a function") : body->name) +
                                 " calls itself");
      }
    }
    if (calls_.size() >= MaxCallDepth) {
      throw LuaUnsupported(line, "the calls nest too deep");
    }
    auto scope = std::make_shared<Scope>();
    scope->parent = closure.scope;
    scope->parentVisible = closure.visible;
    for (std::size_t i = 0; i < body->parameters.size(); ++i) {
      declare(*scope, body->parameters[i], i < args.size() ? args[i] : nil(), false);
    }
    Activation activation;
    activation.body = body;
    if (body->vararg && args.size() > body->parameters.size()) {
      activation.varargs.assign(args.begin() + static_cast<std::ptrdiff_t>(body->parameters.size()),
                                args.end());
    }
    calls_.push_back(std::move(activation));
    Flow caller = std::move(state_.flow);
    state_.flow = Flow{};
    state_.flow.live = caller.live;
    state_.flow.version = ++nextVersion_;
    const std::size_t depth = state_.pendings.size();

    execStatements(body->block, scope);
    resolvePendings(depth);

    std::vector<Value> results;
    const Flow& flow = state_.flow;
    if (flow.hasResults && !flow.returned.isFalse()) {
      results = flow.live.isFalse() ? flow.results : mergeLists(flow.returned, flow.results, {});
    }
    state_.flow = std::move(caller);
    calls_.pop_back();
    release(scope);
    return results;
  }

  void beginEntry() {
    cellBase_ = nextCell_;
    tableBase_ = nextTable_;
    recording_ = true;
    entryDepth_ = calls_.size() + 1;
    if (stale_ == nullptr) {
      return;
    }
    for (const auto& location : *stale_) {
      if (location.cell >= 0) {
        state_.store.cells[location.cell] =
            std::make_shared<Value>(stale("the variable `" + cellNames_[location.cell] + "`"));
      } else {
        const std::string what = location.table == globals_
                                     ? "the global `" + keyText(location.key) + "`"
                                     : "the field `" + keyText(location.key) + "`";
        mutableTable(location.table).fields[location.key] = stale(what);
      }
    }
  }

  void endEntry() {
    recording_ = false;
    entryDepth_ = 0;
  }

  /// A field of a result, where the result is a table; an error of Lua where not.
  Value outputField(const Value& result,
                    const std::string& name,
                    int line,
                    const std::string& message) {
    std::vector<std::pair<Truth, Value>> parts;
    Truth wrong = False;
    for (const auto& alternative : result.alternatives) {
      if (alternative.type == Type::Table) {
        parts.emplace_back(alternative.guard, readField(alternative.id, Key::string(name)));
      } else if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      } else {
        wrong = tOr(wrong, alternative.guard);
      }
    }
    if (!wrong.isFalse()) {
      luaError(message, wrong, line);
    }
    return combine(parts);
  }

  /// The node of an output, a number at every point; a check with `message` where not.
  NodeId outputNode(const Value& v, const std::string& message, int line) {
    const Alternative* found = nullptr;
    Truth wrong = False;
    for (const auto& alternative : v.alternatives) {
      if (alternative.type == Type::Number) {
        found = &alternative;
      } else if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      } else if (alternative.type == Type::String) {
        throw LuaUnsupported(line, "a string as an output is not translated");
      } else {
        wrong = tOr(wrong, alternative.guard);
      }
    }
    if (!wrong.isFalse()) {
      luaError(message, wrong, line);
    }
    return found == nullptr ? konst(Nan) : nodeOf(found->number);
  }

  // ------------------------------------------------------------ library ---

  void makeGlobals() {
    globals_ = newTable();
    const auto set = [&](int id, const std::string& name, const Value& value) {
      rawSet(id, Key::string(name), value);
    };
    const auto library = [&](const std::string& name, const std::vector<std::string>& functions) {
      const int id = newTable();
      for (const auto& function : functions) {
        set(id, function, builtin(name + "." + function));
      }
      set(globals_, name, table(id));
      return id;
    };
    set(globals_, "_G", table(globals_));
    set(globals_, "_VERSION", string("Lua 5.4"));
    for (const char* name :
         {"assert",   "error",  "ipairs",         "next",         "pairs",        "print",
          "rawequal", "rawget", "rawlen",         "rawset",       "select",       "tonumber",
          "tostring", "type",   "collectgarbage", "dofile",       "getmetatable", "load",
          "loadfile", "pcall",  "require",        "setmetatable", "warn",         "xpcall"}) {
      set(globals_, name, builtin(name));
    }
    const int math =
        library("math", {"abs",  "ceil", "cos",       "deg",  "exp", "floor",  "fmod",      "log",
                         "max",  "min",  "modf",      "rad",  "sin", "sqrt",   "tan",       "acos",
                         "asin", "atan", "tointeger", "type", "ult", "random", "randomseed"});
    set(math, "pi", number(Constant::real(Pi)));
    set(math, "huge", number(Constant::real(Inf)));
    set(math, "maxinteger", number(Constant::integer(std::numeric_limits<std::int64_t>::max())));
    set(math, "mininteger", number(Constant::integer(std::numeric_limits<std::int64_t>::min())));
    if (options_.lmathx) {
      // `for fname,fval in pairs(mathx) do math[fname] = fval end`, as easi and SeisSol load it
      const std::vector<std::string> functions = {
          "fabs",      "acos",      "acosh",    "asin",  "asinh",     "atan",  "atan2",  "atanh",
          "cbrt",      "ceil",      "copysign", "cos",   "cosh",      "deg",   "erf",    "erfc",
          "exp",       "exp2",      "expm1",    "fdim",  "floor",     "fma",   "fmax",   "fmin",
          "fmod",      "frexp",     "gamma",    "hypot", "isfinite",  "isinf", "isnan",  "isnormal",
          "ldexp",     "lgamma",    "log",      "log10", "log1p",     "log2",  "logb",   "modf",
          "nearbyint", "nextafter", "pow",      "rad",   "remainder", "round", "scalbn", "sin",
          "sinh",      "sqrt",      "tan",      "tanh",  "trunc"};
      const int mathx = library("mathx", functions);
      set(mathx, "version", string("mathx library for Lua 5.3 / Jun 2015"));
      set(mathx, "inf", number(Constant::real(Inf)));
      set(mathx, "nan", number(Constant::real(Nan)));
      set(mathx, "pi", number(Constant::real(Pi)));
      for (const auto& [key, value] : tableAt(mathx).fields) {
        rawSet(math, key, value);
      }
    }
    library("string",
            {"byte",
             "char",
             "dump",
             "find",
             "format",
             "gmatch",
             "gsub",
             "len",
             "lower",
             "match",
             "rep",
             "reverse",
             "sub",
             "upper",
             "pack",
             "packsize",
             "unpack"});
    library("table", {"concat", "insert", "move", "pack", "remove", "sort", "unpack"});
    library("os",
            {"clock",
             "date",
             "difftime",
             "execute",
             "exit",
             "getenv",
             "remove",
             "rename",
             "setlocale",
             "time",
             "tmpname"});
    library("io",
            {"close",
             "flush",
             "input",
             "lines",
             "open",
             "output",
             "popen",
             "read",
             "tmpfile",
             "type",
             "write"});
    library("coroutine",
            {"close", "create", "isyieldable", "resume", "running", "status", "wrap", "yield"});
    library("utf8", {"char", "codes", "codepoint", "len", "offset"});
    library("debug", {"traceback", "getinfo", "sethook", "gethook", "getlocal"});
    library("package", {"searchpath"});
    if (module_) {
      const int ssol = library("ssol",
                               {"lt",
                                "le",
                                "gt",
                                "ge",
                                "eq",
                                "land",
                                "lor",
                                "lnot",
                                "select",
                                "sign",
                                "mod",
                                "atan2",
                                "pow"});
      set(ssol, "pi", number(Constant::real(M_PI)));
      set(ssol, "g", number(Constant::real(9.80665)));
    }
  }

  /// The grid of a field spec, and its object in `fields`.
  void gridOf(const Value& spec, int fields) {
    if (!spec.is(Type::Table)) {
      return;
    }
    const int id = spec.only().id;
    const auto text = [&](const char* name) {
      const Value v = readField(id, Key::string(name));
      return v.is(Type::String) ? v.only().text : std::string();
    };
    const std::string name = text("name");
    const std::string kindName = text("kind");
    const auto parameters = stringList(readField(id, Key::string("parameters")));
    if (parameters.empty()) {
      throw LuaError(
          0, "field `" + name + "` declares no parameters; they give the number of its components");
    }
    const auto kind = datafield::parseGridKind(kindName);
    if (!kind.has_value()) {
      throw LuaError(0, "field `" + name + "` declares the unknown kind \"" + kindName + "\"");
    }
    const std::string interpolationName =
        text("interpolation").empty() ? std::string("linear") : text("interpolation");
    if (!datafield::parseInterpolation(interpolationName).has_value()) {
      throw LuaError(0,
                     "field `" + name + "` declares the unknown interpolation \"" +
                         interpolationName + "\"");
    }
    // the grid as sderiv declares it, with the defaults of a field spec: the variable of an ASAGI
    // file, the interpolation, and the boundary
    datafield::GridSpec declared;
    declared.kind = kindName;
    declared.path = text("file");
    declared.variable = text("variable");
    if (*kind == datafield::GridKind::AsagiLite && declared.variable.empty()) {
      declared.variable = datafield::DefaultGridVariable;
    }
    declared.components = parameters;
    declared.interpolation = {interpolationName};
    declared.boundary = {"clamp"};
    datafield::GridDesc desc;
    try {
      desc = datafield::gridDescOf(declared);
    } catch (const std::invalid_argument& error) {
      throw LuaError(0, "field `" + name + "`: " + error.what());
    }
    const int slot = static_cast<int>(grids_.size());
    grids_.push_back(Grid{program_.internGrid(desc), static_cast<int>(parameters.size()), name});
    const int object = newTable();
    rawSet(object, Key::string("sample"), builtin("sample", slot));
    rawSet(fields, Key::string(name), table(object));
  }

  /// The number argument `i` (from 0) of a builtin, as luaL_checknumber reads it.
  Number
      argument(const std::vector<Value>& args, std::size_t i, const std::string& name, int line) {
    if (i >= args.size()) {
      luaError("bad argument #" + std::to_string(i + 1) + " to '" + name +
                   "' (number expected, got no value)",
               True,
               line);
      return Number::of(Constant::real(Nan));
    }
    return numberOf(args[i], line, "call '" + name + "' with");
  }

  /// The argument `i` as an integer the script fixes (luaL_checkinteger).
  std::int64_t integerArgument(const std::vector<Value>& args,
                               std::size_t i,
                               const std::string& name,
                               int line) {
    const Number n = argument(args, i, name, line);
    if (!n.known) {
      throw LuaUnsupported(line, "`" + name + "` with an argument that depends on the point");
    }
    const auto whole = exactInteger(n.k);
    if (!whole.has_value()) {
      luaError("bad argument #" + std::to_string(i + 1) + " to '" + name +
                   "' (number has no integer representation)",
               True,
               line);
      return 0;
    }
    return *whole;
  }

  /// A function of floats: Lua's math of the C library, or lmathx's, one argument.
  Value unaryFloat(const std::vector<Value>& args, const std::string& name, int line, Fn fn) {
    const Number x = argument(args, 0, name, line);
    if (x.known) {
      return number(Constant::real(apply(fn, x.k.asFloat())));
    }
    return number(Number::at(pw(fn, x.node)));
  }

  /// A function of the C library the program has no operation for: computed where the script fixes
  /// its arguments, else refused.
  template <typename F>
  Value constantOnly(
      const std::vector<Value>& args, std::size_t count, const std::string& name, int line, F&& f) {
    std::vector<double> values;
    for (std::size_t i = 0; i < count; ++i) {
      const Number n = argument(args, i, name, line);
      if (!n.known) {
        throw LuaUnsupported(line,
                             name + " of a value of the point is not translated: the program has "
                                    "no operation for it");
      }
      values.push_back(n.k.asFloat());
    }
    return number(Constant::real(f(values)));
  }

  std::vector<Value>
      callBuiltin(const std::string& name, int datum, std::vector<Value>& args, int line) {
    const auto arg = [&](std::size_t i) { return i < args.size() ? args[i] : nil(); };
    if (name.rfind("math.", 0) == 0 || name.rfind("mathx.", 0) == 0) {
      return {mathFunction(name, args, line)};
    }
    if (name.rfind("ssol.", 0) == 0) {
      return {ssolFunction(name.substr(5), args, line)};
    }
    if (name == "sample") {
      return sample(datum, args, line);
    }
    if (name == "print" || name == "io.write") {
      return {};
    }
    if (name == "assert") {
      const Truth holds = truthy(arg(0), line);
      if (!holds.isTrue()) {
        const Value message = arg(1);
        luaError(message.is(Type::String) ? message.only().text : "assertion failed!",
                 tNot(holds),
                 line);
      }
      return args;
    }
    if (name == "error") {
      const Value message = arg(0);
      std::string text = "error";
      if (message.is(Type::String)) {
        text = message.only().text;
      } else if (message.is(Type::Number) && message.only().number.known) {
        text = numberText(message.only().number.k);
      }
      luaError(text, True, line);
      return {};
    }
    if (name == "type") {
      const Value v = arg(0);
      if (args.empty()) {
        luaError("bad argument #1 to 'type' (value expected)", True, line);
        return {nil()};
      }
      const std::string first = typeName(v.alternatives.front().type);
      for (const auto& alternative : v.alternatives) {
        if (alternative.type == Type::Stale) {
          refuseStale(alternative, line);
        }
        if (typeName(alternative.type) != first) {
          throw LuaUnsupported(line, "the type of a value depends on the point");
        }
      }
      return {string(first)};
    }
    if (name == "tostring") {
      const Value v = arg(0);
      if (v.is(Type::String)) {
        return {v};
      }
      if (v.is(Type::Number) && v.only().number.known) {
        return {string(numberText(v.only().number.k))};
      }
      if (v.is(Type::Nil)) {
        return {string("nil")};
      }
      if (v.is(Type::Boolean) && v.only().boolean.known()) {
        return {string(v.only().boolean.value ? "true" : "false")};
      }
      throw LuaUnsupported(line, "tostring of a value that depends on the point or of a reference");
    }
    if (name == "tonumber") {
      if (args.size() > 1 && !definitelyNil(args[1])) {
        throw LuaUnsupported(line, "tonumber with a base is not translated");
      }
      const Value v = arg(0);
      if (v.is(Type::String)) {
        const auto parsed = stringToNumber(v.only().text);
        return {parsed.has_value() ? number(*parsed) : nil()};
      }
      std::vector<std::pair<Truth, Value>> parts;
      for (const auto& alternative : v.alternatives) {
        if (alternative.type == Type::Stale) {
          refuseStale(alternative, line);
        }
        if (alternative.type == Type::String) {
          throw LuaUnsupported(line, "tonumber of a string that depends on the point");
        }
        Alternative copy = alternative;
        copy.guard = True;
        parts.emplace_back(alternative.guard,
                           alternative.type == Type::Number ? atom(copy) : nil());
      }
      return {combine(parts)};
    }
    if (name == "ipairs") {
      return {builtin("ipairs.iterator"), arg(0), number(Constant::integer(0))};
    }
    if (name == "ipairs.iterator") {
      const std::int64_t i = integerArgument(args, 1, "ipairs", line) + 1;
      const Value v = index(arg(0), number(Constant::integer(i)), line);
      if (definitelyNil(v)) {
        return {nil()};
      }
      for (const auto& alternative : v.alternatives) {
        if (alternative.type == Type::Nil) {
          throw LuaUnsupported(line, "a loop whose end depends on the point is not translated");
        }
      }
      return {number(Constant::integer(i)), v};
    }
    if (name == "pairs") {
      if (!arg(0).is(Type::Table)) {
        luaError("bad argument #1 to 'for iterator' (table expected)", True, line);
        return {nil()};
      }
      return {builtin("next"), arg(0), nil()};
    }
    if (name == "next") {
      return next(arg(0), arg(1), line);
    }
    if (name == "select") {
      const Value n = arg(0);
      if (n.is(Type::String) && n.only().text == "#") {
        return {number(Constant::integer(static_cast<std::int64_t>(args.size()) - 1))};
      }
      std::int64_t i = integerArgument(args, 0, "select", line);
      const auto count = static_cast<std::int64_t>(args.size()) - 1;
      if (i < 0) {
        i = count + i;
        if (i < 0) {
          luaError("bad argument #1 to 'select' (index out of range)", True, line);
          return {};
        }
      } else if (i == 0) {
        luaError("bad argument #1 to 'select' (index out of range)", True, line);
        return {};
      } else {
        i -= 1;
      }
      if (i >= count) {
        return {};
      }
      return std::vector<Value>(args.begin() + 1 + i, args.end());
    }
    if (name == "rawget") {
      return {index(arg(0), arg(1), line)};
    }
    if (name == "rawset") {
      setIndex(arg(0), arg(1), arg(2), line, "");
      return {arg(0)};
    }
    if (name == "rawequal") {
      return {boolean(equal(arg(0), arg(1), line))};
    }
    if (name == "rawlen") {
      return {length(arg(0), line)};
    }
    if (name == "table.insert") {
      return tableInsert(args, line);
    }
    if (name == "table.remove") {
      return tableRemove(args, line);
    }
    if (name == "table.unpack") {
      if (!arg(0).is(Type::Table)) {
        throw LuaUnsupported(line, "table.unpack of a value that is no table at every point");
      }
      const int id = arg(0).only().id;
      const std::int64_t first =
          args.size() > 1 && !definitelyNil(args[1]) ? integerArgument(args, 1, "unpack", line) : 1;
      const std::int64_t last = args.size() > 2 && !definitelyNil(args[2])
                                    ? integerArgument(args, 2, "unpack", line)
                                    : border(id, line);
      std::vector<Value> out;
      for (std::int64_t i = first; i <= last; ++i) {
        out.push_back(readField(id, Key::integer(i)));
      }
      return out;
    }
    if (name == "table.pack") {
      const int id = newTable();
      for (std::size_t i = 0; i < args.size(); ++i) {
        rawSet(id, Key::integer(static_cast<std::int64_t>(i + 1)), args[i]);
      }
      rawSet(
          id, Key::string("n"), number(Constant::integer(static_cast<std::int64_t>(args.size()))));
      return {table(id)};
    }
    if (name == "collectgarbage") {
      return {number(Constant::integer(0))};
    }
    throw LuaUnsupported(line, name + " is not translated");
  }

  std::vector<Value> next(const Value& t, const Value& k, int line) {
    if (!t.is(Type::Table)) {
      throw LuaUnsupported(line, "next over a value that is no table at every point");
    }
    const auto& fields = tableAt(t.only().id).fields;
    auto it = fields.begin();
    if (!definitelyNil(k)) {
      const auto key = keyOf(k, line);
      it = key.has_value() ? fields.find(*key) : fields.end();
      if (it == fields.end()) {
        luaError("invalid key to 'next'", True, line);
        return {nil()};
      }
      ++it;
    }
    if (it == fields.end()) {
      return {nil()};
    }
    for (const auto& alternative : it->second.alternatives) {
      if (alternative.type == Type::Nil) {
        throw LuaUnsupported(line, "the keys of a table depend on the point");
      }
    }
    Value key;
    switch (it->first.kind) {
    case Key::Kind::Number:
      key = number(it->first.number);
      break;
    case Key::Kind::String:
      key = string(it->first.text);
      break;
    case Key::Kind::Boolean:
      key = boolean(known(it->first.id != 0));
      break;
    case Key::Kind::Table:
      key = table(it->first.id);
      break;
    default:
      throw LuaUnsupported(line, "a table with functions as keys is not translated");
    }
    return {key, it->second};
  }

  std::vector<Value> tableInsert(const std::vector<Value>& args, int line) {
    if (args.empty() || !args[0].is(Type::Table)) {
      throw LuaUnsupported(line, "table.insert into a value that is no table at every point");
    }
    const int id = args[0].only().id;
    const std::int64_t size = border(id, line);
    if (args.size() == 2) {
      writeField(id, Key::integer(size + 1), args[1], True);
      return {};
    }
    if (args.size() != 3) {
      luaError("wrong number of arguments to 'insert'", True, line);
      return {};
    }
    const std::int64_t position = integerArgument(args, 1, "insert", line);
    if (position < 1 || position > size + 1) {
      luaError("bad argument #2 to 'insert' (position out of bounds)", True, line);
      return {};
    }
    for (std::int64_t i = size; i >= position; --i) {
      writeField(id, Key::integer(i + 1), readField(id, Key::integer(i)), True);
    }
    writeField(id, Key::integer(position), args[2], True);
    return {};
  }

  std::vector<Value> tableRemove(const std::vector<Value>& args, int line) {
    if (args.empty() || !args[0].is(Type::Table)) {
      throw LuaUnsupported(line, "table.remove from a value that is no table at every point");
    }
    const int id = args[0].only().id;
    const std::int64_t size = border(id, line);
    std::int64_t position = size;
    if (args.size() > 1 && !definitelyNil(args[1])) {
      position = integerArgument(args, 1, "remove", line);
      if (position != size &&
          static_cast<std::uint64_t>(position) - 1U > static_cast<std::uint64_t>(size)) {
        luaError("bad argument #2 to 'remove' (position out of bounds)", True, line);
        return {};
      }
    }
    const Value removed = readField(id, Key::integer(position));
    for (; position < size; ++position) {
      writeField(id, Key::integer(position), readField(id, Key::integer(position + 1)), True);
    }
    writeField(id, Key::integer(position), nil(), True);
    return {removed};
  }

  std::vector<Value> sample(int slot, const std::vector<Value>& args, int line) {
    const auto& grid = grids_[static_cast<std::size_t>(slot)];
    const auto coordinates = args.size() > 0 ? args.size() - 1 : 0;
    if (coordinates < 1 || coordinates > static_cast<std::size_t>(expr::MaxLookupDimension)) {
      throw LuaUnsupported(line,
                           "field `" + grid.name + "` sampled with " + std::to_string(coordinates) +
                               " coordinates; a grid has 1 to " +
                               std::to_string(expr::MaxLookupDimension));
    }
    std::vector<NodeId> nodes;
    for (std::size_t i = 1; i < args.size(); ++i) {
      nodes.push_back(nodeOf(argument(args, i, "sample", line)));
    }
    std::vector<Value> out;
    for (int component = 0; component < grid.components; ++component) {
      out.push_back(number(Number::at(arena_.lookup(grid.id, component, nodes))));
    }
    return out;
  }

  /// The truth of an argument of ssol: a boolean, a number other than zero (lua_tonumber).
  Truth ssolTruth(const Value& v, int line) {
    Truth t = False;
    for (const auto& alternative : v.alternatives) {
      Truth here = False;
      if (alternative.type == Type::Boolean) {
        here = alternative.boolean;
      } else if (alternative.type == Type::Number) {
        here = tNot(equalNumbers(alternative.number, Number::of(Constant::real(0.0))));
      } else if (alternative.type == Type::Stale) {
        refuseStale(alternative, line);
      } else if (alternative.type == Type::String) {
        throw LuaUnsupported(line, "a string as a condition of ssol is not translated");
      }
      t = tOr(t, tAnd(alternative.guard, here));
    }
    return t;
  }

  Value ssolFunction(const std::string& name, const std::vector<Value>& args, int line) {
    const auto arg = [&](std::size_t i) { return i < args.size() ? args[i] : nil(); };
    const auto two = [&](Fn fn) -> Value {
      const Number x = argument(args, 0, "ssol." + name, line);
      const Number y = argument(args, 1, "ssol." + name, line);
      if (x.known && y.known) {
        return number(Constant::real(apply(fn, x.k.asFloat(), y.k.asFloat())));
      }
      return number(Number::at(pw(fn, nodeOf(x), nodeOf(y))));
    };
    if (name == "lt" || name == "le" || name == "gt" || name == "ge" || name == "eq") {
      const Number x = argument(args, 0, "ssol." + name, line);
      const Number y = argument(args, 1, "ssol." + name, line);
      const Number fx = x.known ? Number::of(Constant::real(x.k.asFloat())) : x;
      const Number fy = y.known ? Number::of(Constant::real(y.k.asFloat())) : y;
      if (name == "eq") {
        return boolean(equalNumbers(fx, fy));
      }
      const bool swap = name == "gt" || name == "ge";
      return boolean(lessNumbers(swap ? fy : fx, swap ? fx : fy, name == "le" || name == "ge"));
    }
    if (name == "land") {
      return boolean(tAnd(ssolTruth(arg(0), line), ssolTruth(arg(1), line)));
    }
    if (name == "lor") {
      return boolean(tOr(ssolTruth(arg(0), line), ssolTruth(arg(1), line)));
    }
    if (name == "lnot") {
      return boolean(tNot(ssolTruth(arg(0), line)));
    }
    if (name == "select") {
      return merge(ssolTruth(arg(0), line), arg(1), arg(2));
    }
    if (name == "sign") {
      return unaryFloat(args, "ssol.sign", line, Fn::Sign);
    }
    if (name == "mod") {
      return two(Fn::Mod);
    }
    if (name == "atan2") {
      return two(Fn::Atan2);
    }
    if (name == "pow") {
      return two(Fn::Pow);
    }
    throw LuaUnsupported(line, "ssol." + name + " is not translated");
  }

  Value mathFunction(const std::string& qualified, const std::vector<Value>& args, int line) {
    const bool lmathx = qualified.rfind("mathx.", 0) == 0;
    const std::string name = qualified.substr(qualified.find('.') + 1);
    const auto x = [&]() { return argument(args, 0, name, line); };
    const auto floatOf = [&](const Number& n) { return nodeOf(n); };
    static const std::map<std::string, Fn> Unary = {{"sqrt", Fn::Sqrt},
                                                    {"exp", Fn::Exp},
                                                    {"sin", Fn::Sin},
                                                    {"cos", Fn::Cos},
                                                    {"tan", Fn::Tan},
                                                    {"asin", Fn::Asin},
                                                    {"acos", Fn::Acos},
                                                    {"log2", Fn::Log2},
                                                    {"log10", Fn::Log10},
                                                    {"sinh", Fn::Sinh},
                                                    {"cosh", Fn::Cosh},
                                                    {"tanh", Fn::Tanh},
                                                    {"asinh", Fn::Asinh},
                                                    {"acosh", Fn::Acosh},
                                                    {"atanh", Fn::Atanh},
                                                    {"erf", Fn::Erf},
                                                    {"fabs", Fn::Abs},
                                                    {"round", Fn::Round}};
    if (const auto found = Unary.find(name); found != Unary.end()) {
      return unaryFloat(args, name, line, found->second);
    }
    if (name == "floor" || name == "ceil") {
      const Fn fn = name == "floor" ? Fn::Floor : Fn::Ceil;
      const Number n = x();
      if (lmathx) {
        // lmathx: C's floor, a float
        if (n.known) {
          return number(Constant::real(apply(fn, n.k.asFloat())));
        }
        return number(Number::at(pw(fn, n.node)));
      }
      // Lua's: an integer where it fits one (pushnumint)
      if (n.known) {
        if (n.k.isInt) {
          return number(n.k);
        }
        const double value = apply(fn, n.k.f);
        const auto whole = exactInteger(value);
        return number(whole.has_value() ? Constant::integer(*whole) : Constant::real(value));
      }
      const NodeId rounded = pw(fn, n.node);
      const Truth fits = tAnd(truthOf(pw(Fn::Le, konst(-TwoTo63), rounded)),
                              truthOf(pw(Fn::Lt, rounded, konst(TwoTo63))));
      return number(Number::at(pw(Fn::Add, rounded, konst(0.0)), tOr(n.intness, fits)));
    }
    if (name == "abs") {
      const Number n = x();
      if (n.known) {
        return number(n.k.isInt
                          ? Constant::integer(
                                n.k.i < 0 ? wrapped(0U - static_cast<std::uint64_t>(n.k.i)) : n.k.i)
                          : Constant::real(std::fabs(n.k.f)));
      }
      return number(Number::at(pw(Fn::Abs, n.node), n.intness));
    }
    if (name == "max" || name == "min") {
      // Lua's own: the first argument no later one is above (below)
      Number m = x();
      for (std::size_t i = 1; i < args.size(); ++i) {
        const Number next = argument(args, i, name, line);
        const Truth better =
            name == "max" ? lessNumbers(m, next, false) : lessNumbers(next, m, false);
        m = better.known() ? (better.value ? next : m) : selectNumber(better, next, m);
      }
      return number(m);
    }
    if (name == "fmax" || name == "fmin") {
      const Fn fn = name == "fmax" ? Fn::Max : Fn::Min;
      Number m = x();
      if (m.known) {
        m = Number::of(Constant::real(m.k.asFloat()));
      } else {
        m.intness = False;
      }
      for (std::size_t i = 1; i < args.size(); ++i) {
        const Number next = argument(args, i, name, line);
        if (m.known && next.known) {
          m = Number::of(Constant::real(apply(fn, m.k.asFloat(), next.k.asFloat())));
        } else {
          m = Number::at(pw(fn, floatOf(m), floatOf(next)));
        }
      }
      return number(m);
    }
    if (name == "log") {
      const Number n = x();
      const bool withBase =
          lmathx ? args.size() >= 2 : (args.size() >= 2 && !definitelyNil(args[1]));
      if (!withBase) {
        return unaryFloat(args, name, line, Fn::Log);
      }
      const Number base = argument(args, 1, name, line);
      if (base.known) {
        const double b = base.k.asFloat();
        const Fn fn = b == 10.0 ? Fn::Log10 : (b == 2.0 ? Fn::Log2 : Fn::Log);
        if (n.known) {
          const double value = apply(fn, n.k.asFloat());
          return number(Constant::real(fn == Fn::Log ? value / std::log(b) : value));
        }
        const NodeId value = pw(fn, n.node);
        return number(Number::at(fn == Fn::Log ? pw(Fn::Div, value, konst(std::log(b))) : value));
      }
      // the order of the tests differs between lmathx and Lua; the values do not
      const NodeId v = floatOf(n);
      const NodeId b = base.node;
      const NodeId other = pw(Fn::Div, pw(Fn::Log, v), pw(Fn::Log, b));
      return number(Number::at(select(pw(Fn::Eq, b, konst(10.0)),
                                      pw(Fn::Log10, v),
                                      select(pw(Fn::Eq, b, konst(2.0)), pw(Fn::Log2, v), other))));
    }
    if (name == "atan") {
      if (lmathx) {
        if (args.size() < 2) {
          return unaryFloat(args, name, line, Fn::Atan);
        }
      } else if (args.size() < 2 || definitelyNil(args[1])) {
        // Lua's own: atan2(y, 1)
        const Number y = x();
        if (y.known) {
          return number(Constant::real(std::atan2(y.k.asFloat(), 1.0)));
        }
        return number(Number::at(pw(Fn::Atan2, y.node, konst(1.0))));
      }
    }
    if (name == "atan" || name == "atan2" || name == "pow") {
      const Fn fn = name == "pow" ? Fn::Pow : Fn::Atan2;
      const Number a = argument(args, 0, name, line);
      const Number b = argument(args, 1, name, line);
      if (a.known && b.known) {
        return number(Constant::real(apply(fn, a.k.asFloat(), b.k.asFloat())));
      }
      return number(Number::at(pw(fn, floatOf(a), floatOf(b))));
    }
    if (name == "deg" || name == "rad") {
      const Number n = x();
      const double factor = name == "deg" ? 180.0 / Pi : Pi / 180.0;
      if (n.known) {
        return number(Constant::real(n.k.asFloat() * factor));
      }
      return number(Number::at(pw(Fn::Mul, n.node, konst(factor))));
    }
    if (name == "trunc") {
      const Number n = x();
      if (n.known) {
        return number(Constant::real(std::trunc(n.k.asFloat())));
      }
      return number(Number::at(
          select(pw(Fn::Lt, n.node, konst(0.0)), pw(Fn::Ceil, n.node), pw(Fn::Floor, n.node))));
    }
    if (name == "fdim") {
      const Number a = argument(args, 0, name, line);
      const Number b = argument(args, 1, name, line);
      if (a.known && b.known) {
        return number(Constant::real(std::fdim(a.k.asFloat(), b.k.asFloat())));
      }
      const NodeId u = floatOf(a);
      const NodeId v = floatOf(b);
      return number(Number::at(select(pw(Fn::Le, u, v), konst(0.0), pw(Fn::Sub, u, v))));
    }
    if (name == "isnan" || name == "isinf" || name == "isfinite" || name == "isnormal") {
      const Number n = x();
      if (n.known) {
        const double v = n.k.asFloat();
        const bool result = name == "isnan"      ? std::isnan(v)
                            : name == "isinf"    ? std::isinf(v)
                            : name == "isfinite" ? std::isfinite(v)
                                                 : std::isnormal(v);
        return boolean(known(result));
      }
      const NodeId magnitude = pw(Fn::Abs, n.node);
      if (name == "isnan") {
        return boolean(tNot(truthOf(pw(Fn::Eq, n.node, n.node))));
      }
      if (name == "isinf") {
        return boolean(truthOf(pw(Fn::Eq, magnitude, konst(Inf))));
      }
      const Truth finite = truthOf(pw(Fn::Lt, magnitude, konst(Inf)));
      if (name == "isfinite") {
        return boolean(finite);
      }
      return boolean(tAnd(finite, truthOf(pw(Fn::Le, konst(DBL_MIN), magnitude))));
    }
    if (name == "ldexp" || name == "scalbn") {
      const Number n = x();
      std::int64_t e = 0;
      if (name == "ldexp") {
        e = integerArgument(args, 1, name, line);
      } else {
        const Number exponent = argument(args, 1, name, line);
        if (!exponent.known) {
          throw LuaUnsupported(line, "scalbn with an exponent that depends on the point");
        }
        // lmathx passes the float on to an int parameter: it is truncated
        e = static_cast<std::int64_t>(std::trunc(exponent.k.asFloat()));
      }
      const int exponent = static_cast<int>(std::clamp<std::int64_t>(e, -100000, 100000));
      if (n.known) {
        return number(Constant::real(std::ldexp(n.k.asFloat(), exponent)));
      }
      if (exponent < -1074 || exponent > 1023) {
        throw LuaUnsupported(line,
                             name + " of a value of the point by 2^" + std::to_string(exponent));
      }
      // x * 2^e is rounded once, as ldexp rounds
      return number(Number::at(pw(Fn::Mul, n.node, konst(std::ldexp(1.0, exponent)))));
    }
    if (name == "tointeger") {
      const Value v = args.empty() ? nil() : args[0];
      std::vector<std::pair<Truth, Value>> parts;
      for (const auto& alternative : v.alternatives) {
        if (alternative.type == Type::Stale) {
          refuseStale(alternative, line);
        }
        if (alternative.type == Type::String) {
          throw LuaUnsupported(line, "math.tointeger of a string is not translated");
        }
        if (alternative.type != Type::Number) {
          parts.emplace_back(alternative.guard, nil());
          continue;
        }
        const Number& n = alternative.number;
        if (n.known) {
          const auto whole = exactInteger(n.k);
          parts.emplace_back(alternative.guard,
                             whole.has_value() ? number(Constant::integer(*whole)) : nil());
          continue;
        }
        // an integer where the float is one
        const NodeId rounded = pw(Fn::Floor, n.node);
        const Truth whole = tOr(n.intness,
                                tAnd(truthOf(pw(Fn::Eq, rounded, n.node)),
                                     tAnd(truthOf(pw(Fn::Le, konst(-TwoTo63), n.node)),
                                          truthOf(pw(Fn::Lt, n.node, konst(TwoTo63))))));
        parts.emplace_back(
            alternative.guard,
            merge(whole, number(Number::at(pw(Fn::Add, n.node, konst(0.0)), True)), nil()));
      }
      return combine(parts);
    }
    if (name == "type") {
      const Value v = args.empty() ? nil() : args[0];
      if (v.is(Type::Number)) {
        const Truth intness = intnessOf(v.only().number);
        if (!intness.known()) {
          throw LuaUnsupported(line, "math.type of a number whose subtype depends on the point");
        }
        return string(intness.value ? "integer" : "float");
      }
      for (const auto& alternative : v.alternatives) {
        if (alternative.type == Type::Number) {
          throw LuaUnsupported(line, "math.type of a value that depends on the point");
        }
      }
      return nil();
    }
    if (name == "ult") {
      const std::int64_t a = integerArgument(args, 0, name, line);
      const std::int64_t b = integerArgument(args, 1, name, line);
      return boolean(known(static_cast<std::uint64_t>(a) < static_cast<std::uint64_t>(b)));
    }
    if (name == "random" || name == "randomseed") {
      throw LuaUnsupported(
          line, "math." + name + " is not translated: the program would depend on the run");
    }
    // the C library, where the script fixes the arguments
    const auto one = [&](double (*f)(double)) {
      return constantOnly(
          args, 1, "math." + name, line, [&](const std::vector<double>& v) { return f(v[0]); });
    };
    const auto pair = [&](double (*f)(double, double)) {
      return constantOnly(args, 2, "math." + name, line, [&](const std::vector<double>& v) {
        return f(v[0], v[1]);
      });
    };
    if (name == "fmod") {
      if (!lmathx) {
        const Number a = argument(args, 0, name, line);
        const Number b = argument(args, 1, name, line);
        if (a.known && b.known && a.k.isInt && b.k.isInt) {
          if (b.k.i == 0) {
            luaError("bad argument #2 to 'fmod' (zero)", True, line);
            return number(Constant::integer(0));
          }
          return number(Constant::integer(b.k.i == -1 ? 0 : a.k.i % b.k.i));
        }
      }
      return pair([](double a, double b) { return std::fmod(a, b); });
    }
    if (name == "modf") {
      throw LuaUnsupported(line, "math.modf is not translated");
    }
    if (name == "cbrt") {
      return one([](double v) { return std::cbrt(v); });
    }
    if (name == "erfc") {
      return one([](double v) { return std::erfc(v); });
    }
    if (name == "exp2") {
      return one([](double v) { return std::exp2(v); });
    }
    if (name == "expm1") {
      return one([](double v) { return std::expm1(v); });
    }
    if (name == "log1p") {
      return one([](double v) { return std::log1p(v); });
    }
    if (name == "logb") {
      return one([](double v) { return std::logb(v); });
    }
    if (name == "gamma") {
      return one([](double v) { return std::tgamma(v); });
    }
    if (name == "lgamma") {
      return one([](double v) { return std::lgamma(v); });
    }
    if (name == "nearbyint") {
      return one([](double v) { return std::nearbyint(v); });
    }
    if (name == "copysign") {
      return pair([](double a, double b) { return std::copysign(a, b); });
    }
    if (name == "hypot") {
      return pair([](double a, double b) { return std::hypot(a, b); });
    }
    if (name == "nextafter") {
      return pair([](double a, double b) { return std::nextafter(a, b); });
    }
    if (name == "remainder") {
      return pair([](double a, double b) { return std::remainder(a, b); });
    }
    if (name == "fma") {
      return constantOnly(args, 3, "math.fma", line, [](const std::vector<double>& v) {
        return std::fma(v[0], v[1], v[2]);
      });
    }
    throw LuaUnsupported(line, "math." + name + " is not translated");
  }

  struct Grid {
    expr::GridId id;
    int components;
    std::string name;
  };

  expr::Program& program_;
  expr::Arena& arena_;
  const LuaOptions& options_;
  const LuaOrigin& origin_;
  const std::set<Location>* stale_;
  bool module_;

  Truth baseReach_{True};
  State state_;
  std::vector<LuaCheck> checks_;
  std::unordered_set<NodeId> truthNodes_;
  std::unordered_map<const Expression*, bool> purity_;
  std::map<int, std::string> cellNames_;
  std::vector<Activation> calls_;
  std::vector<Grid> grids_;
  int globals_{-1};
  int nextCell_{0};
  int nextTable_{0};
  int nextClosure_{0};
  int nextVersion_{0};
  std::size_t steps_{0};

  // the call of the entry function (f, evaluate): what it finds, writes and returns
  bool recording_{false};
  int cellBase_{0};
  int tableBase_{0};
  std::size_t entryDepth_{0};
  std::set<Location> leaks_;
  std::vector<int> returnLines_;
};

/// Rejects what Lua rejects when it compiles a chunk, before anything runs: `goto` (which the
/// frontend refuses) and `break` outside a loop.
void checkStatically(const Block& block, int loops) {
  for (const auto& s : block.statements) {
    switch (s.kind) {
    case StatKind::Goto:
      throw LuaUnsupported(s.line, "goto is not translated");
    case StatKind::Break:
      if (loops == 0) {
        throw LuaError(s.line, "break outside a loop");
      }
      break;
    default:
      break;
    }
    if (s.body != nullptr) {
      const bool isLoop = s.kind == StatKind::While || s.kind == StatKind::Repeat ||
                          s.kind == StatKind::NumericFor || s.kind == StatKind::GenericFor;
      checkStatically(*s.body, isLoop ? loops + 1 : loops);
    }
    for (const auto& b : s.blocks) {
      checkStatically(*b, loops);
    }
    if (s.otherwise != nullptr) {
      checkStatically(*s.otherwise, loops);
    }
    if (s.function != nullptr) {
      checkStatically(s.function->block, 0);
    }
  }
}

std::shared_ptr<FunctionBody> parse(const std::string& code) {
  auto chunk = parseChunk(code);
  checkStatically(chunk->block, 0);
  return chunk;
}

} // namespace

LuaMapTranslation translateLuaMap(expr::Program& program,
                                  const std::string& code,
                                  const std::vector<std::pair<std::string, expr::NodeId>>& in,
                                  const std::set<std::string>& returns,
                                  expr::NodeId reach,
                                  const LuaOrigin& origin,
                                  const LuaOptions& options) {
  const auto chunk = parse(code);
  // a value f keeps from one point to the next: easi runs the chunk once per thread, and f at
  // every point; a second run marks what f writes of that as stale, which f may not read first
  std::set<Location> stale;
  for (int run = 0;; ++run) {
    Evaluator evaluator(program, options, origin, reach, run == 0 ? nullptr : &stale, false);
    auto translation = evaluator.luaMap(chunk, in, returns);
    if (run == 0 && !evaluator.leaks().empty()) {
      stale = evaluator.leaks();
      continue;
    }
    return translation;
  }
}

expr::Program translateLuaModule(const std::string& code,
                                 const LuaOrigin& origin,
                                 const LuaOptions& options) {
  const auto chunk = parse(code);
  std::set<Location> stale;
  for (int run = 0;; ++run) {
    expr::Program program;
    Evaluator evaluator(program, options, origin, NoNode, run == 0 ? nullptr : &stale, true);
    evaluator.module(chunk);
    if (run == 0 && !evaluator.leaks().empty()) {
      stale = evaluator.leaks();
      continue;
    }
    try {
      expr::validate(program);
    } catch (const std::invalid_argument& error) {
      throw LuaError(0, std::string("the translated program is invalid: ") + error.what());
    }
    return program;
  }
}

} // namespace seissol::reader::scripting::lua
