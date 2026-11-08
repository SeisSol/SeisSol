// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/SderivFrontend.h"

#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

// ============================================================ tokenizer =====
enum class TokenKind : std::uint8_t {
  Num,
  Str,
  Pow,
  Op,
  Name,
  Let,
  In,
  Def,
  Out,
  Grid,
  State,
  Eof
};

struct Token {
  TokenKind kind{};
  std::string value;
  int position{0};
};

const char* describe(TokenKind kind) {
  switch (kind) {
  case TokenKind::Num:
    return "a number";
  case TokenKind::Str:
    return "a quoted string";
  case TokenKind::Name:
    return "a name";
  case TokenKind::Eof:
    return "end of input";
  default:
    return "a token";
  }
}

std::vector<Token> tokenize(const std::string& source) {
  std::vector<Token> tokens;
  const std::size_t n = source.size();
  auto nameStart = [](unsigned char c) { return std::isalpha(c) != 0 || c == '_'; };
  auto nameCont = [](unsigned char c) { return std::isalnum(c) != 0 || c == '_'; };

  std::size_t i = 0;
  while (i < n) {
    const auto c = static_cast<unsigned char>(source[i]);
    if (std::isspace(c) != 0) {
      ++i;
      continue;
    }
    if (c == '#') {
      while (i < n && source[i] != '\n') {
        ++i;
      }
      continue;
    }
    const std::size_t start = i;

    // Only \" and \\ are escapes; everything else is literal. That is enough to
    // keep a Windows path intact without growing an escape language, and an
    // unterminated string is reported at its OPENING quote, because reporting
    // it at end of input points at the wrong line.
    if (c == '"') {
      ++i;
      std::string value;
      bool closed = false;
      while (i < n) {
        if (source[i] == '\\' && i + 1 < n && (source[i + 1] == '"' || source[i + 1] == '\\')) {
          value.push_back(source[i + 1]);
          i += 2;
          continue;
        }
        if (source[i] == '"') {
          ++i;
          closed = true;
          break;
        }
        if (source[i] == '\n') {
          break;
        }
        value.push_back(source[i]);
        ++i;
      }
      if (!closed) {
        throw SderivError("lex", "unterminated string", static_cast<int>(start));
      }
      tokens.push_back({TokenKind::Str, std::move(value), static_cast<int>(start)});
      continue;
    }

    if (std::isdigit(c) != 0) {
      while (i < n && std::isdigit(static_cast<unsigned char>(source[i])) != 0) {
        ++i;
      }
      if (i + 1 < n && source[i] == '.' &&
          std::isdigit(static_cast<unsigned char>(source[i + 1])) != 0) {
        ++i;
        while (i < n && std::isdigit(static_cast<unsigned char>(source[i])) != 0) {
          ++i;
        }
      }
      if (i < n && (source[i] == 'e' || source[i] == 'E')) {
        std::size_t j = i + 1;
        if (j < n && (source[j] == '+' || source[j] == '-')) {
          ++j;
        }
        if (j < n && std::isdigit(static_cast<unsigned char>(source[j])) != 0) {
          i = j + 1;
          while (i < n && std::isdigit(static_cast<unsigned char>(source[i])) != 0) {
            ++i;
          }
        }
      }
      tokens.push_back({TokenKind::Num, source.substr(start, i - start), static_cast<int>(start)});
      continue;
    }
    if (c == '*' && i + 1 < n && source[i + 1] == '*') {
      tokens.push_back({TokenKind::Pow, "**", static_cast<int>(start)});
      i += 2;
      continue;
    }
    if (c == '+' || c == '-' || c == '*' || c == '/' || c == '(' || c == ')' || c == ',' ||
        c == '=') {
      tokens.push_back(
          {TokenKind::Op, std::string(1, static_cast<char>(c)), static_cast<int>(start)});
      ++i;
      continue;
    }
    if (nameStart(c)) {
      ++i;
      while (i < n && nameCont(static_cast<unsigned char>(source[i]))) {
        ++i;
      }
      std::string value = source.substr(start, i - start);
      TokenKind kind = TokenKind::Name;
      if (value == "let") {
        kind = TokenKind::Let;
      } else if (value == "in") {
        kind = TokenKind::In;
      } else if (value == "def") {
        kind = TokenKind::Def;
      } else if (value == "out") {
        kind = TokenKind::Out;
      } else if (value == "grid") {
        kind = TokenKind::Grid;
      } else if (value == "state") {
        kind = TokenKind::State;
      }
      tokens.push_back({kind, std::move(value), static_cast<int>(start)});
      continue;
    }
    throw SderivError("lex",
                      "unexpected character '" + std::string(1, static_cast<char>(c)) + "'",
                      static_cast<int>(i));
  }
  tokens.push_back({TokenKind::Eof, "", static_cast<int>(n)});
  return tokens;
}

// ========================================================== surface AST =====
using SurfaceId = int;
enum class SurfaceKind : std::uint8_t { Num, Name, Call, Binary, Unary, Let, Def };

struct SurfaceNode {
  SurfaceKind kind{};
  double number{0};
  std::string text; // Name.id / Call.fn / Binary.op / Unary.op / Let.name / Def.name
  std::vector<SurfaceId> args;
  std::vector<std::string> params;
  SurfaceId a{-1};
  SurfaceId b{-1};
  int position{0};
};

class SurfaceArena {
  public:
  SurfaceId add(SurfaceNode&& node) {
    nodes_.push_back(std::move(node));
    return static_cast<SurfaceId>(nodes_.size()) - 1;
  }
  const SurfaceNode& operator[](SurfaceId id) const { return nodes_[id]; }

  private:
  std::vector<SurfaceNode> nodes_;
};

// Grid declarations are a table rather than AST nodes, so no AST walk needs a
// new case for them.
struct GridDeclaration {
  std::string name;
  std::string kind;
  std::string file;
  std::string variable;
  std::string interpolation;
  std::vector<std::string> components;
  /// what a sample beyond the volume does (GridDesc::boundary), and the axis that is time, if any
  /// (GridDesc::timeAxis): only the form with attributes states them
  std::string boundary{"clamp"};
  std::optional<int> timeAxis;
  int position{0};
};

// A value kept per point across calls. Its next value is the definition of the same name, and
// every reference to the name reads the value from the previous call.
struct StateDeclaration {
  std::string name;
  double initial{0.0};
  int position{0};
  // `state f64 name = ...`: the type of the state (see ValueType)
  ValueType type{ValueType::Real};
};

struct ParsedProgram {
  std::vector<StateDeclaration> states;
  std::vector<SurfaceId> defs;
  // The subset of `defs` marked `out`, in declaration order. A subset rather
  // than a separate list: an exported definition IS an ordinary definition and
  // resolves through the same symbol table, so a later output can read an
  // earlier one by name. Keeping them in two containers would have made that a
  // coincidence rather than a property.
  std::vector<SurfaceId> exports;
  std::vector<GridDeclaration> grids;
  SurfaceId expression{-1}; // -1 when the module declares its outputs itself
};

// =============================================== catalogs (spec-defined) ====
const std::unordered_map<std::string, double> Constants = {{"pi", M_PI}, {"g", 9.80665}};

// Spelled exactly as the Lua-side ssol.* names, so a model transliterates.
const std::unordered_map<std::string, Fn> Builtins = {
    {"sqrt", Fn::Sqrt},   {"abs", Fn::Abs},      {"exp", Fn::Exp},     {"log", Fn::Log},
    {"log2", Fn::Log2},   {"log10", Fn::Log10},  {"sign", Fn::Sign},   {"floor", Fn::Floor},
    {"ceil", Fn::Ceil},   {"round", Fn::Round},  {"sin", Fn::Sin},     {"cos", Fn::Cos},
    {"tan", Fn::Tan},     {"asin", Fn::Asin},    {"acos", Fn::Acos},   {"sinh", Fn::Sinh},
    {"cosh", Fn::Cosh},   {"tanh", Fn::Tanh},    {"asinh", Fn::Asinh}, {"acosh", Fn::Acosh},
    {"atanh", Fn::Atanh}, {"erf", Fn::Erf},      {"min", Fn::Min},     {"max", Fn::Max},
    {"pow", Fn::Pow},     {"atan2", Fn::Atan2},  {"mod", Fn::Mod},     {"lt", Fn::Lt},
    {"le", Fn::Le},       {"eq", Fn::Eq},        {"land", Fn::And},    {"lor", Fn::Or},
    {"atan", Fn::Atan},   {"select", Fn::Select}};

// The conversions between the types of values (see ValueType): `f64(x)` and `f32(x)` give x in
// binary64 and binary32, `real(x)` in the compute type of the program. The same names give the
// type of a state: `state f64 sum = 0.0`.
const std::unordered_map<std::string, ValueType> Conversions = {
    {"f64", ValueType::F64}, {"f32", ValueType::F32}, {"real", ValueType::Real}};

/// Whether `name` is a function: a definition, a state or an input of that name would be read as
/// the function instead.
bool isFunction(const std::string& name) {
  return Builtins.count(name) != 0 || Conversions.count(name) != 0;
}

// The names are spelled exactly as the Lua-side ssol.* ones, so a model
// transliterates. `gt`, `ge` and `lnot` are rewrites, as in the Lua tracer: the
// IR has Lt and Le and no Gt or Ge, because a > b IS b < a and a second node
// kind would be a second thing for every consumer to handle.
const std::unordered_map<std::string, Fn> SwappedComparisons = {{"gt", Fn::Lt}, {"ge", Fn::Le}};

// The grid vocabularies are closed, so a swapped file/interpolation pair fails
// at parse time rather than at grid-load time on one rank of a large job. They
// are NOT restated here: datafield::parseGridKind / parseInterpolation are the
// vocabulary, and a scheme added to Grid.h becomes available to this frontend
// without a second edit that someone can forget.

// The dataset inside the file. Optional, and it is the closed interpolation
// vocabulary above that makes the slot unambiguous: in
//   "asagi", FILE, X, ...
// X is the interpolation if it names one and the variable otherwise. The
// default matches AsagiLiteGrid::open's own default rather than restating a
// convention in a second place.
//
// The alternative -- a required slot -- was rejected because the Lua field_spec
// table has to gain the same field, and there it MUST be optional: field_specs
// is a pre-existing ABI with scripts already written against it. Two spellings
// of the same option in the two frontends is exactly the divergence the shared
// Program is supposed to remove.
// Variable default and the two vocabulary parsers now live in Grid.h, so the
// Lua path cannot disagree with this one about what "linear" means.
using reader::datafield::DefaultGridVariable;

// ============================================================== parser ======
class Parser {
  public:
  Parser(const std::vector<Token>& tokens, SurfaceArena& arena) : tokens_(tokens), arena_(arena) {}

  ParsedProgram program() {
    ParsedProgram parsed;
    // def and grid interleave freely and need no separator, for the same reason
    // def already needed none: each starts with its own token kind, so the
    // greedy expression in a def body stops cleanly at the next declaration.
    while (peek().kind == TokenKind::Def || peek().kind == TokenKind::Grid ||
           peek().kind == TokenKind::Out || peek().kind == TokenKind::State) {
      if (peek().kind == TokenKind::Grid) {
        parsed.grids.push_back(gridDeclaration());
        continue;
      }
      if (peek().kind == TokenKind::State) {
        parsed.states.push_back(stateDeclaration());
        continue;
      }
      const bool exported = peek().kind == TokenKind::Out;
      if (exported) {
        eat(TokenKind::Out);
      }
      const SurfaceId id = definition(exported);
      parsed.defs.push_back(id);
      if (exported) {
        parsed.exports.push_back(id);
      }
    }

    // A module is EITHER a trailing expression, whose name the caller supplies,
    // OR a set of `out def`s that name themselves. Allowing both would leave it
    // open which one is "the" output, and a file that grew an `out def` while
    // keeping its trailing expression would silently gain a second one.
    if (parsed.exports.empty()) {
      parsed.expression = expression();
    }
    eat(TokenKind::Eof);
    return parsed;
  }

  private:
  const Token& peek() const { return tokens_[index_]; }
  const Token& next() { return tokens_[index_++]; }
  bool valueIs(const char* v) const { return peek().value == v; }

  const Token& eat(TokenKind kind, const char* value = nullptr) {
    const Token& token = peek();
    if (token.kind != kind || (value != nullptr && token.value != value)) {
      const std::string expected =
          value != nullptr ? std::string("'") + value + "'" : std::string(describe(kind));
      throw SderivError(
          "parse", "expected " + expected + ", got '" + token.value + "'", token.position);
    }
    return next();
  }

  // 'grid' NAME '=' KIND ',' FILE ',' [VARIABLE ','] INTERPOLATION {',' COMPONENT}, which clamps
  // and has no time axis; or, as in the core format (the sderiv spec), 'grid' NAME '(' ATTRIBUTE
  // {',' ATTRIBUTE} ')', with ATTRIBUTE = NAME '=' (STRING | INTEGER | '(' STRING {',' STRING}
  // ')').
  GridDeclaration gridDeclaration() {
    GridDeclaration grid;
    grid.position = peek().position;
    eat(TokenKind::Grid);
    grid.name = eat(TokenKind::Name).value;
    if (peek().kind == TokenKind::Op && valueIs("(")) {
      gridAttributes(grid);
      return grid;
    }
    eat(TokenKind::Op, "=");
    grid.kind = eat(TokenKind::Str).value;
    eat(TokenKind::Op, ",");
    grid.file = eat(TokenKind::Str).value;
    eat(TokenKind::Op, ",");
    const Token& slot = eat(TokenKind::Str);
    if (reader::datafield::parseInterpolation(slot.value).has_value()) {
      grid.variable = DefaultGridVariable;
      grid.interpolation = slot.value;
    } else {
      grid.variable = slot.value;
      eat(TokenKind::Op, ",");
      grid.interpolation = eat(TokenKind::Str).value;
    }
    while (valueIs(",")) {
      eat(TokenKind::Op, ",");
      grid.components.push_back(eat(TokenKind::Str).value);
    }
    return grid;
  }

  // The attributes of a grid, of which the core format requires all but `time`: it has no
  // defaults, so a declaration says what it means; the short form above is the one with them.
  void gridAttributes(GridDeclaration& grid) {
    const std::vector<std::string> names = {
        "kind", "path", "variable", "components", "interpolation", "boundary", "time"};
    std::set<std::string> given;
    eat(TokenKind::Op, "(");
    while (true) {
      const Token& key = eat(TokenKind::Name);
      if (std::find(names.begin(), names.end(), key.value) == names.end()) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: unknown attribute `" + key.value +
                              "`; a grid has kind, path, variable, components, interpolation, "
                              "boundary and time",
                          key.position);
      }
      if (!given.insert(key.value).second) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: the attribute `" + key.value +
                              "` is given twice",
                          key.position);
      }
      eat(TokenKind::Op, "=");
      if (key.value == "components") {
        if (peek().kind == TokenKind::Op && valueIs("(")) {
          eat(TokenKind::Op, "(");
          grid.components.push_back(eat(TokenKind::Str).value);
          while (valueIs(",")) {
            eat(TokenKind::Op, ",");
            grid.components.push_back(eat(TokenKind::Str).value);
          }
          eat(TokenKind::Op, ")");
        } else {
          grid.components.push_back(eat(TokenKind::Str).value);
        }
      } else if (key.value == "time") {
        const Token& axis = eat(TokenKind::Num);
        if (axis.value.find_first_not_of("0123456789") != std::string::npos ||
            axis.value.size() > 2 || std::stoi(axis.value) >= MaxLookupDimension) {
          throw SderivError("grid",
                            "grid `" + grid.name + "`: the time is an axis, 0 to " +
                                std::to_string(MaxLookupDimension - 1) + ", not " + axis.value,
                            axis.position);
        }
        grid.timeAxis = std::stoi(axis.value);
      } else {
        const std::string value = eat(TokenKind::Str).value;
        if (key.value == "kind") {
          grid.kind = value;
        } else if (key.value == "path") {
          grid.file = value;
        } else if (key.value == "variable") {
          grid.variable = value;
        } else if (key.value == "interpolation") {
          grid.interpolation = value;
        } else {
          grid.boundary = value;
        }
      }
      if (!valueIs(",")) {
        break;
      }
      eat(TokenKind::Op, ",");
    }
    eat(TokenKind::Op, ")");
    for (const auto& name : names) {
      if (name != "time" && given.count(name) == 0) {
        throw SderivError("grid",
                          "grid `" + grid.name + "` gives no " + name +
                              "; a grid with attributes gives kind, path, variable, components, "
                              "interpolation and boundary, and may give time",
                          grid.position);
      }
    }
  }

  // 'state' [TYPE] NAME '=' ['-'] NUM -- the initial value is a literal: it is set before any
  // point has been seen, so there is nothing for an expression to read. TYPE is f64, f32 or real.
  StateDeclaration stateDeclaration() {
    StateDeclaration state;
    state.position = peek().position;
    eat(TokenKind::State);
    state.name = eat(TokenKind::Name).value;
    if (peek().kind == TokenKind::Name) {
      const auto type = Conversions.find(state.name);
      if (type == Conversions.end()) {
        throw SderivError("parse",
                          "`" + state.name +
                              "` is no type of a state; the types are f64, f32 "
                              "and real",
                          state.position);
      }
      state.type = type->second;
      state.name = eat(TokenKind::Name).value;
    }
    eat(TokenKind::Op, "=");
    const bool negative = valueIs("-");
    if (negative) {
      eat(TokenKind::Op, "-");
    }
    const double magnitude = std::stod(eat(TokenKind::Num).value);
    state.initial = negative ? -magnitude : magnitude;
    return state;
  }

  SurfaceId definition(bool exported) {
    const int position = peek().position;
    eat(TokenKind::Def);
    SurfaceNode node;
    node.kind = SurfaceKind::Def;
    node.position = position;
    node.text = eat(TokenKind::Name).value;
    if (valueIs("(")) {
      eat(TokenKind::Op, "(");
      if (!valueIs(")")) {
        node.params.push_back(eat(TokenKind::Name).value);
        while (valueIs(",")) {
          eat(TokenKind::Op, ",");
          node.params.push_back(eat(TokenKind::Name).value);
        }
      }
      eat(TokenKind::Op, ")");
    }
    // An output is a value per point, not a function -- there is no call site
    // for it to take arguments from. Rejected here rather than at lowering so
    // the diagnostic carries the declaration's position.
    if (exported && !node.params.empty()) {
      throw SderivError("parse", "an output definition takes no parameters", position);
    }
    eat(TokenKind::Op, "=");
    node.a = expression();
    return arena_.add(std::move(node));
  }

  /// How deep expressions may nest, in parentheses, arguments, lets and signs: deep enough for any
  /// source written or generated (the specification asks for at least 256 levels), and shallow
  /// enough that the parser, which recurses into them, cannot overflow its stack.
  static constexpr int MaxNesting = 2048;

  /// One more level of nesting while it lives.
  class Nesting {
public:
    explicit Nesting(Parser& parser) : parser_(parser) {
      if (++parser_.nesting_ > MaxNesting) {
        throw SderivError("parse",
                          "the expression is nested deeper than " + std::to_string(MaxNesting) +
                              " levels",
                          parser_.peek().position);
      }
    }
    ~Nesting() { --parser_.nesting_; }
    Nesting(const Nesting&) = delete;
    Nesting& operator=(const Nesting&) = delete;
    Nesting(Nesting&&) = delete;
    Nesting& operator=(Nesting&&) = delete;

private:
    Parser& parser_;
  };

  SurfaceId expression() {
    const Nesting nesting(*this);
    if (peek().kind == TokenKind::Let) {
      SurfaceNode node;
      node.kind = SurfaceKind::Let;
      node.position = peek().position;
      eat(TokenKind::Let);
      node.text = eat(TokenKind::Name).value;
      eat(TokenKind::Op, "=");
      node.a = expression();
      eat(TokenKind::In);
      node.b = expression();
      return arena_.add(std::move(node));
    }
    return additive();
  }

  SurfaceId additive() {
    SurfaceId lhs = multiplicative();
    while (valueIs("+") || valueIs("-")) {
      SurfaceNode node;
      node.kind = SurfaceKind::Binary;
      node.position = peek().position;
      node.text = next().value;
      node.a = lhs;
      node.b = multiplicative();
      lhs = arena_.add(std::move(node));
    }
    return lhs;
  }

  SurfaceId multiplicative() {
    SurfaceId lhs = unary();
    while (valueIs("*") || valueIs("/")) {
      SurfaceNode node;
      node.kind = SurfaceKind::Binary;
      node.position = peek().position;
      node.text = next().value;
      node.a = lhs;
      node.b = unary();
      lhs = arena_.add(std::move(node));
    }
    return lhs;
  }

  SurfaceId unary() {
    const Nesting nesting(*this);
    if (valueIs("-")) {
      SurfaceNode node;
      node.kind = SurfaceKind::Unary;
      node.position = peek().position;
      node.text = "neg";
      next();
      node.a = unary();
      return arena_.add(std::move(node));
    }
    return power();
  }

  SurfaceId power() {
    SurfaceId base = primary();
    if (peek().kind == TokenKind::Pow) {
      SurfaceNode node;
      node.kind = SurfaceKind::Binary;
      node.position = peek().position;
      node.text = "pow";
      next();
      node.a = base;
      node.b = unary(); // right-associative
      return arena_.add(std::move(node));
    }
    return base;
  }

  SurfaceId primary() {
    const Token& token = peek();
    if (token.kind == TokenKind::Num) {
      next();
      SurfaceNode node;
      node.kind = SurfaceKind::Num;
      node.position = token.position;
      node.number = std::stod(token.value);
      return arena_.add(std::move(node));
    }
    // The whole point of confining strings: no case here, and the message says
    // where they ARE allowed rather than "unexpected token".
    if (token.kind == TokenKind::Str) {
      throw SderivError(
          "parse", "string literals are only allowed in a `grid` declaration", token.position);
    }
    if (token.kind == TokenKind::Name) {
      const int position = token.position;
      std::string name = next().value;
      if (valueIs("(")) {
        eat(TokenKind::Op, "(");
        SurfaceNode node;
        node.kind = SurfaceKind::Call;
        node.position = position;
        node.text = std::move(name);
        if (!valueIs(")")) {
          node.args.push_back(expression());
          while (valueIs(",")) {
            eat(TokenKind::Op, ",");
            node.args.push_back(expression());
          }
        }
        eat(TokenKind::Op, ")");
        return arena_.add(std::move(node));
      }
      SurfaceNode node;
      node.kind = SurfaceKind::Name;
      node.position = position;
      node.text = std::move(name);
      return arena_.add(std::move(node));
    }
    if (token.value == "(") {
      eat(TokenKind::Op, "(");
      SurfaceId inner = expression();
      eat(TokenKind::Op, ")");
      return inner;
    }
    throw SderivError("parse", "unexpected '" + token.value + "'", token.position);
  }

  const std::vector<Token>& tokens_;
  SurfaceArena& arena_;
  std::size_t index_{0};
  int nesting_{0};
};

// ================================================= grid component functions =
struct ComponentFunction {
  std::size_t grid{0};
  std::int32_t component{0};
};

// std::map, not unordered: every diagnostic below iterates it, and a message
// whose ordering depends on a hash is a message that differs between builds.
using ComponentTable = std::map<std::string, ComponentFunction>;

ComponentTable checkGrids(const std::vector<GridDeclaration>& grids,
                          const std::set<std::string>& defNames) {
  ComponentTable table;
  std::set<std::string> gridNames;

  for (std::size_t g = 0; g < grids.size(); ++g) {
    const GridDeclaration& grid = grids[g];
    if (!gridNames.insert(grid.name).second) {
      throw SderivError("grid", "duplicate grid `" + grid.name + "`", grid.position);
    }
    if (!reader::datafield::parseGridKind(grid.kind).has_value()) {
      throw SderivError(
          "grid", "grid `" + grid.name + "`: unknown kind \"" + grid.kind + "\"", grid.position);
    }
    if (!reader::datafield::parseInterpolation(grid.interpolation).has_value()) {
      throw SderivError("grid",
                        "grid `" + grid.name + "`: unknown interpolation \"" + grid.interpolation +
                            "\"; did the file and the interpolation get swapped?",
                        grid.position);
    }
    if (!reader::datafield::parseBoundaryMode(grid.boundary).has_value()) {
      throw SderivError("grid",
                        "grid `" + grid.name + "`: unknown boundary \"" + grid.boundary +
                            "\"; it is \"clamp\" or \"fail\"",
                        grid.position);
    }
    if (grid.components.empty()) {
      throw SderivError("grid", "grid `" + grid.name + "` declares no components", grid.position);
    }
    std::set<std::string> seen;
    for (std::size_t c = 0; c < grid.components.size(); ++c) {
      const std::string& component = grid.components[c];
      if (!seen.insert(component).second) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: duplicate component \"" + component + "\"",
                          grid.position);
      }
      const std::string generated = grid.name + "_" + component;
      // Collisions are errors rather than shadowing: silently overriding a
      // builtin or a def is the exact failure mode this frontend exists to
      // remove from the reader path.
      if (Constants.count(generated) != 0 || isFunction(generated) ||
          defNames.count(generated) != 0 || table.count(generated) != 0) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: the generated name `" + generated +
                              "` collides with an existing name",
                          grid.position);
      }
      table[generated] = {g, static_cast<std::int32_t>(c)};
    }
  }

  // A grid name is never itself an expression, so this is not formally
  // ambiguous -- but nobody can read `m_rho(...)` and tell which declaration it
  // came from. This check is what underscore flattening costs, and the reason
  // it is still cheaper than adding a `.` token to the language.
  for (const auto& grid : grids) {
    if (table.count(grid.name) != 0) {
      throw SderivError("grid",
                        "grid `" + grid.name +
                            "` has the same name as a component function generated by another grid",
                        grid.position);
    }
  }
  return table;
}

// ============================================================== lowering ====
class Lowering {
  public:
  Lowering(const SurfaceArena& surface,
           const ParsedProgram& parsed,
           const ComponentTable& components,
           const std::vector<GridId>& gridIds,
           const std::set<std::string>& inputs,
           Program& program)
      : surface_(surface), parsed_(parsed), components_(components), gridIds_(gridIds),
        inputs_(inputs), program_(program) {
    for (const SurfaceId id : parsed_.defs) {
      defs_[surface_[id].text] = id;
    }
    for (const auto& state : parsed_.states) {
      states_.insert(state.name);
    }
  }

  NodeId lower(SurfaceId id) {
    Environment empty;
    return lower(id, empty);
  }

  /// Lowers the definitions without parameters that the expressions `roots` read, directly or
  /// through others, each once, and in an order in which a definition comes after those it reads;
  /// the expressions themselves are lowered after, with lower(). So no definition is lowered again
  /// where it is read once more, which was exponential in the depth of a chain of shared
  /// definitions, and no definition is lowered inside another: the length of a chain of
  /// definitions takes no stack. A definition that reads itself is an error, as before. It also
  /// fixes the order of the channels, that of their first read in the source, from the first
  /// expression on (channels()).
  void prepare(const std::vector<SurfaceId>& roots) {
    // the definitions the roots reach, and which each reads
    std::map<std::string, std::vector<std::string>> reads;
    std::vector<std::string> pending;
    for (const SurfaceId root : roots) {
      for (auto& name : definitionsRead(root, {})) {
        pending.push_back(std::move(name));
      }
    }
    while (!pending.empty()) {
      const std::string name = std::move(pending.back());
      pending.pop_back();
      if (reads.count(name) != 0) {
        continue;
      }
      const SurfaceNode& definition = surface_[defs_.at(name)];
      auto read = definitionsRead(
          definition.a, std::set<std::string>(definition.params.begin(), definition.params.end()));
      pending.insert(pending.end(), read.begin(), read.end());
      reads.emplace(name, std::move(read));
    }

    // a definition after those it reads: a depth-first walk, without recursion
    enum class Mark : std::uint8_t { New, Open, Done };
    std::map<std::string, Mark> marks;
    std::vector<std::string> order;
    for (const auto& [start, _] : reads) {
      if (marks[start] != Mark::New) {
        continue;
      }
      std::vector<std::pair<std::string, std::size_t>> stack{{start, 0}};
      marks[start] = Mark::Open;
      while (!stack.empty()) {
        auto& [name, next] = stack.back();
        const auto& read = reads.at(name);
        if (next < read.size()) {
          const std::string child = read[next++];
          auto& mark = marks[child];
          if (mark == Mark::Open) {
            fail(surface_[defs_.at(child)], "`" + child + "` is defined in terms of itself");
          }
          if (mark == Mark::New) {
            mark = Mark::Open;
            stack.emplace_back(child, 0);
          }
          continue;
        }
        marks[name] = Mark::Done;
        order.push_back(name);
        stack.pop_back();
      }
    }
    for (const auto& name : order) {
      const SurfaceNode& definition = surface_[defs_.at(name)];
      if (definition.params.empty() && lowered_.count(name) == 0) {
        Environment empty;
        lowered_.emplace(name, lower(definition.a, empty));
      }
    }

    channels_ = channelsInOrder(roots);
  }

  // Channels the roots read, in the order of their first read in the source, from the first root
  // on, as a walk that reads a definition where it is named meets them (see prepare()).
  // Deterministic because the walk is, which matters: the signature order feeds the binding and,
  // through it, the tile layout.
  [[nodiscard]] const std::vector<std::string>& channels() const { return channels_; }

  private:
  using Environment = std::map<std::string, NodeId>;

  /// What a name resolves to where it is read, besides a binding of a let or a parameter (which
  /// comes first): an input the consumer names, a state, a constant, a definition, or else (where
  /// it is not called) a channel. The order is that of lower().
  enum class Resolution : std::uint8_t { Input, State, Constant, Definition, Other };

  [[nodiscard]] Resolution resolve(const std::string& name) const {
    if (inputs_.count(name) != 0) {
      return Resolution::Input;
    }
    if (states_.count(name) != 0) {
      return Resolution::State;
    }
    if (Constants.count(name) != 0) {
      return Resolution::Constant;
    }
    if (defs_.count(name) != 0) {
      return Resolution::Definition;
    }
    return Resolution::Other;
  }

  /// The names of the definitions the expression `root` reads directly, where the names `bound`
  /// are bound; without recursion.
  [[nodiscard]] std::vector<std::string> definitionsRead(SurfaceId root,
                                                         const std::set<std::string>& bound) const {
    std::vector<std::string> found;
    using Scope = std::shared_ptr<const std::set<std::string>>;
    std::vector<std::pair<SurfaceId, Scope>> stack{
        {root, std::make_shared<const std::set<std::string>>(bound)}};
    while (!stack.empty()) {
      const auto [id, scope] = stack.back();
      stack.pop_back();
      const SurfaceNode& node = surface_[id];
      switch (node.kind) {
      case SurfaceKind::Unary:
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Binary:
        stack.emplace_back(node.a, scope);
        stack.emplace_back(node.b, scope);
        break;
      case SurfaceKind::Let: {
        stack.emplace_back(node.a, scope);
        auto inner = std::make_shared<std::set<std::string>>(*scope);
        inner->insert(node.text);
        stack.emplace_back(node.b, std::move(inner));
        break;
      }
      case SurfaceKind::Name:
        if (scope->count(node.text) == 0 && resolve(node.text) == Resolution::Definition) {
          found.push_back(node.text);
        }
        break;
      case SurfaceKind::Call:
        for (const SurfaceId arg : node.args) {
          stack.emplace_back(arg, scope);
        }
        if (scope->count(node.text) == 0 && components_.count(node.text) == 0 &&
            defs_.count(node.text) != 0) {
          found.push_back(node.text);
        }
        break;
      case SurfaceKind::Num:
      case SurfaceKind::Def:
        break;
      }
    }
    return found;
  }

  /// The channels the roots read, as described at channels(): a walk in the order of the source,
  /// which reads a definition where it is first named; without recursion.
  [[nodiscard]] std::vector<std::string>
      channelsInOrder(const std::vector<SurfaceId>& roots) const {
    std::vector<std::string> channels;
    std::set<std::string> seen;
    std::set<std::string> expanded;
    std::set<std::pair<std::string, std::vector<SurfaceId>>> calls;
    using Scope = std::shared_ptr<const std::set<std::string>>;
    const auto empty = std::make_shared<const std::set<std::string>>();
    std::vector<std::pair<SurfaceId, Scope>> stack;
    // the last root first, so that the first one is walked first
    for (auto root = roots.rbegin(); root != roots.rend(); ++root) {
      stack.emplace_back(*root, empty);
    }
    const auto add = [&](const std::string& name) {
      if (seen.insert(name).second) {
        channels.push_back(name);
      }
    };
    while (!stack.empty()) {
      const auto [id, scope] = stack.back();
      stack.pop_back();
      const SurfaceNode& node = surface_[id];
      switch (node.kind) {
      case SurfaceKind::Unary:
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Binary:
        // the left operand first, so it is pushed last
        stack.emplace_back(node.b, scope);
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Let: {
        auto inner = std::make_shared<std::set<std::string>>(*scope);
        inner->insert(node.text);
        stack.emplace_back(node.b, std::move(inner));
        stack.emplace_back(node.a, scope);
        break;
      }
      case SurfaceKind::Name: {
        if (scope->count(node.text) != 0) {
          break;
        }
        const auto resolution = resolve(node.text);
        if (resolution == Resolution::Input || resolution == Resolution::Other) {
          add(node.text);
        } else if (resolution == Resolution::Definition && expanded.insert(node.text).second) {
          stack.emplace_back(surface_[defs_.at(node.text)].a, empty);
        }
        break;
      }
      case SurfaceKind::Call: {
        const auto definition = defs_.find(node.text);
        const bool callsDefinition = scope->count(node.text) == 0 &&
                                     components_.count(node.text) == 0 && definition != defs_.end();
        // the body of a definition after its arguments, which lower() binds first; the body of one
        // with parameters once per call of other arguments
        if (callsDefinition && calls.emplace(node.text, node.args).second) {
          const SurfaceNode& callee = surface_[definition->second];
          stack.emplace_back(callee.a,
                             std::make_shared<const std::set<std::string>>(callee.params.begin(),
                                                                           callee.params.end()));
        }
        for (auto arg = node.args.rbegin(); arg != node.args.rend(); ++arg) {
          stack.emplace_back(*arg, scope);
        }
        break;
      }
      case SurfaceKind::Num:
      case SurfaceKind::Def:
        break;
      }
    }
    return channels;
  }

  [[noreturn]] void fail(const SurfaceNode& node, const std::string& message) const {
    throw SderivError("resolve", message, node.position);
  }

  /// The definition of `node` while it is being lowered: one that is reached again reads itself.
  class Expansion {
public:
    Expansion(Lowering& lowering, const SurfaceNode& node) : lowering_(lowering) {
      auto& expanding = lowering_.expanding_;
      if (std::find(expanding.begin(), expanding.end(), node.text) != expanding.end()) {
        lowering_.fail(node, "`" + node.text + "` is defined in terms of itself");
      }
      expanding.push_back(node.text);
    }
    ~Expansion() { lowering_.expanding_.pop_back(); }
    Expansion(const Expansion&) = delete;
    Expansion& operator=(const Expansion&) = delete;
    Expansion(Expansion&&) = delete;
    Expansion& operator=(Expansion&&) = delete;

private:
    Lowering& lowering_;
  };

  NodeId channel(const SurfaceNode& node) { return program_.arena().field(node.text); }

  static Fn binaryFn(const std::string& op) {
    if (op == "+") {
      return Fn::Add;
    }
    if (op == "-") {
      return Fn::Sub;
    }
    if (op == "*") {
      return Fn::Mul;
    }
    if (op == "/") {
      return Fn::Div;
    }
    return Fn::Pow;
  }

  NodeId lowerLookup(const SurfaceNode& node, const ComponentFunction& fn, Environment& env) {
    if (node.args.empty() || node.args.size() > 6) {
      fail(node,
           "`" + node.text + "`: sampled with " + std::to_string(node.args.size()) +
               " coordinates, the IR allows 1..6");
    }
    const auto& timeAxis = parsed_.grids[fn.grid].timeAxis;
    if (timeAxis.has_value() && static_cast<std::size_t>(*timeAxis) >= node.args.size()) {
      fail(node,
           "`" + node.text + "`: sampled with " + std::to_string(node.args.size()) +
               " coordinates, but axis " + std::to_string(*timeAxis) + " is the time");
    }
    std::vector<NodeId> coords;
    coords.reserve(node.args.size());
    for (const SurfaceId arg : node.args) {
      coords.push_back(lower(arg, env));
    }
    return program_.arena().lookup(gridIds_[fn.grid], fn.component, coords);
  }

  NodeId lower(SurfaceId id, Environment& env) {
    const SurfaceNode& node = surface_[id];
    switch (node.kind) {
    case SurfaceKind::Num:
      return program_.arena().konst(node.number);

    case SurfaceKind::Unary:
      return program_.arena().pw(Fn::Neg, lower(node.a, env));

    case SurfaceKind::Binary: {
      // A sum of many terms is a long chain of binary nodes along their left operands, as the
      // parser builds it; it is walked without recursion, from its first term on.
      std::vector<SurfaceId> spine{id};
      while (surface_[surface_[spine.back()].a].kind == SurfaceKind::Binary) {
        spine.push_back(surface_[spine.back()].a);
      }
      NodeId value = lower(surface_[spine.back()].a, env);
      for (auto link = spine.rbegin(); link != spine.rend(); ++link) {
        const SurfaceNode& binary = surface_[*link];
        const NodeId b = lower(binary.b, env);
        value = program_.arena().pw(binaryFn(binary.text), value, b);
      }
      return value;
    }

    case SurfaceKind::Let: {
      const NodeId value = lower(node.a, env);
      Environment inner = env;
      inner[node.text] = value;
      // Sharing rather than duplication falls out of interning: the bound id is
      // substituted, so both uses reach the same node without a CSE pass.
      return lower(node.b, inner);
    }

    case SurfaceKind::Name: {
      const auto bound = env.find(node.text);
      if (bound != env.end()) {
        return bound->second;
      }
      // An input the consumer names reads the input, as a state reads its previous value.
      if (inputs_.count(node.text) != 0) {
        return channel(node);
      }
      // A state reads its value from the previous call -- also where the definition of the same
      // name, its next value, is in scope. That is what makes the update a parallel assignment.
      if (states_.count(node.text) != 0) {
        return program_.arena().field(node.text);
      }
      const auto constant = Constants.find(node.text);
      if (constant != Constants.end()) {
        return program_.arena().konst(constant->second);
      }
      const auto def = defs_.find(node.text);
      if (def != defs_.end()) {
        const SurfaceNode& definition = surface_[def->second];
        if (!definition.params.empty()) {
          fail(node,
               "`" + node.text + "` expects " + std::to_string(definition.params.size()) +
                   " arguments");
        }
        // it reads no environment, so it has one value wherever it is read
        if (const auto known = lowered_.find(node.text); known != lowered_.end()) {
          return known->second;
        }
        const Expansion expansion(*this, node);
        Environment inner;
        const NodeId value = lower(definition.a, inner);
        lowered_.emplace(node.text, value);
        return value;
      }
      if (components_.count(node.text) != 0) {
        fail(node, "`" + node.text + "` is a grid component and must be called with coordinates");
      }
      if (isFunction(node.text)) {
        fail(node, "`" + node.text + "` is a function and must be called");
      }
      // Everything left is an input channel. No declaration, no catalogue.
      return channel(node);
    }

    case SurfaceKind::Call: {
      if (env.count(node.text) != 0) {
        fail(node, "`" + node.text + "` is not callable");
      }

      const auto component = components_.find(node.text);
      if (component != components_.end()) {
        return lowerLookup(node, component->second, env);
      }

      const auto def = defs_.find(node.text);
      if (def != defs_.end()) {
        const SurfaceNode& definition = surface_[def->second];
        if (node.args.size() != definition.params.size()) {
          fail(node,
               "`" + node.text + "` expects " + std::to_string(definition.params.size()) +
                   " arguments, got " + std::to_string(node.args.size()));
        }
        Environment inner;
        std::vector<NodeId> arguments;
        arguments.reserve(definition.params.size());
        for (std::size_t i = 0; i < definition.params.size(); ++i) {
          arguments.push_back(lower(node.args[i], env));
          inner[definition.params[i]] = arguments.back();
        }
        // the value of a call depends on its arguments alone
        auto key = std::make_pair(node.text, std::move(arguments));
        if (const auto known = calls_.find(key); known != calls_.end()) {
          return known->second;
        }
        const Expansion expansion(*this, node);
        const NodeId value = lower(definition.a, inner);
        calls_.emplace(std::move(key), value);
        return value;
      }

      // `lnot` has no IR op on purpose: negating a 0/1 value is 1-c, and an op
      // for it would be an op every backend has to implement. The rewrite lives
      // here so that promise is kept in exactly one place.
      if (node.text == "lnot") {
        if (node.args.size() != 1) {
          fail(node, "`lnot` takes 1 argument, got " + std::to_string(node.args.size()));
        }
        const NodeId one = program_.arena().konst(1.0);
        return program_.arena().pw(Fn::Sub, one, lower(node.args[0], env));
      }

      // a > b is b < a. Rewritten here rather than given its own Fn, so every
      // consumer of the IR keeps one comparison shape to handle.
      if (const auto swapped = SwappedComparisons.find(node.text);
          swapped != SwappedComparisons.end()) {
        if (node.args.size() != 2) {
          fail(node,
               "`" + node.text + "` takes 2 arguments, got " + std::to_string(node.args.size()));
        }
        // in the order of the source, which the order of the evaluation of the arguments of a
        // call would leave to the compiler
        const NodeId left = lower(node.args[0], env);
        const NodeId right = lower(node.args[1], env);
        return program_.arena().pw(swapped->second, right, left);
      }
      // A condition is 0.0 or 1.0, so negation is 1 - c. Same rewrite the Lua
      // tracer performs, and it keeps conditions in the arithmetic domain
      // rather than introducing a boolean one.
      if (node.text == "lnot") {
        if (node.args.size() != 1) {
          fail(node, "`lnot` takes 1 argument, got " + std::to_string(node.args.size()));
        }
        return program_.arena().pw(Fn::Sub, program_.arena().konst(1.0), lower(node.args[0], env));
      }

      // a conversion to another type: f64(x), f32(x), real(x)
      const auto conversion = Conversions.find(node.text);
      if (conversion != Conversions.end()) {
        if (node.args.size() != 1) {
          fail(node,
               "`" + node.text + "` takes 1 argument, got " + std::to_string(node.args.size()));
        }
        return program_.arena().cast(conversion->second, lower(node.args[0], env));
      }

      const auto builtin = Builtins.find(node.text);
      if (builtin == Builtins.end()) {
        fail(node, "unknown function `" + node.text + "`");
      }
      const Fn fn = builtin->second;
      const int wanted = arity(fn);

      // min and max are variadic in every language a model author is likely to
      // come from, so fold rather than reject; the IR op stays binary.
      if ((fn == Fn::Min || fn == Fn::Max) && node.args.size() > 2) {
        NodeId acc = lower(node.args[0], env);
        for (std::size_t i = 1; i < node.args.size(); ++i) {
          acc = program_.arena().pw(fn, acc, lower(node.args[i], env));
        }
        return acc;
      }
      if (static_cast<int>(node.args.size()) != wanted) {
        fail(node,
             "`" + node.text + "` takes " + std::to_string(wanted) + " arguments, got " +
                 std::to_string(node.args.size()));
      }
      if (wanted == 1) {
        return program_.arena().pw(fn, lower(node.args[0], env));
      }
      if (wanted == 2) {
        const NodeId a = lower(node.args[0], env);
        const NodeId b = lower(node.args[1], env);
        return program_.arena().pw(fn, a, b);
      }
      const NodeId a = lower(node.args[0], env);
      const NodeId b = lower(node.args[1], env);
      const NodeId c = lower(node.args[2], env);
      return program_.arena().pw(fn, a, b, c);
    }

    case SurfaceKind::Def:
      break;
    }
    fail(node, "cannot lower this node");
  }

  const SurfaceArena& surface_;
  const ParsedProgram& parsed_;
  const ComponentTable& components_;
  const std::vector<GridId>& gridIds_;
  const std::set<std::string>& inputs_;
  Program& program_;
  std::map<std::string, SurfaceId> defs_;
  std::set<std::string> states_;
  std::vector<std::string> channels_;
  /// the value of every definition without parameters lowered so far
  std::map<std::string, NodeId> lowered_;
  /// the value of every call of a definition with parameters so far, by its arguments
  std::map<std::pair<std::string, std::vector<NodeId>>, NodeId> calls_;
  /// the definitions being lowered, innermost last
  std::vector<std::string> expanding_;
};

} // namespace

namespace {

// Compile one source into `program`. `externalName` is set for the trailing-
// expression form, where the caller names the single output; it is null for a
// module that declares its own outputs with `out def`.
//
// Shared by both entry points rather than duplicated, because everything from
// the grid interning down is identical and the only difference is where the
// output names come from.
void compileSource(const std::string& source,
                   const std::string* externalName,
                   reader::scripting::DataType type,
                   const SderivOptions& options,
                   Program& program,
                   std::vector<std::string>& channelOrder) {
  SurfaceArena surface;
  const auto tokens = tokenize(source);
  Parser parser(tokens, surface);
  const ParsedProgram parsed = parser.program();

  if (externalName != nullptr && parsed.expression < 0) {
    throw SderivError("parse",
                      "the caller named the output `" + *externalName +
                          "`, but this module declares its own outputs with `out def`",
                      0);
  }
  if (externalName == nullptr && parsed.exports.empty()) {
    throw SderivError(
        "parse", "the module declares no outputs; mark at least one definition `out def`", 0);
  }

  std::set<std::string> defNames;
  for (const SurfaceId id : parsed.defs) {
    // validate() catches a duplicate OUTPUT name, but a duplicate plain `def`
    // would silently take the last one, because the symbol table is a map. That
    // is a shadowing bug the author cannot see, so it is caught here.
    if (!defNames.insert(surface[id].text).second) {
      throw SderivError(
          "parse", "duplicate definition `" + surface[id].text + "`", surface[id].position);
    }
    // Neither may a definition take the name of a constant or a function: a read of the name
    // would resolve to that instead, and the definition be ignored without a word. (That of a
    // state is reported for the state, below.)
    const bool ofState =
        std::any_of(parsed.states.begin(), parsed.states.end(), [&](const auto& state) {
          return state.name == surface[id].text;
        });
    if (!ofState && (Constants.count(surface[id].text) != 0 || isFunction(surface[id].text))) {
      throw SderivError("parse",
                        "the definition `" + surface[id].text +
                            "` has the name of a constant or function",
                        surface[id].position);
    }
  }
  const ComponentTable components = checkGrids(parsed.grids, defNames);

  // Every state needs exactly one definition of the same name, its next value; and a state name
  // must not be anything else a name can resolve to.
  std::set<std::string> stateNames;
  for (const auto& state : parsed.states) {
    if (!stateNames.insert(state.name).second) {
      throw SderivError("state", "duplicate state `" + state.name + "`", state.position);
    }
    if (Constants.count(state.name) != 0 || isFunction(state.name) ||
        components.count(state.name) != 0) {
      throw SderivError("state",
                        "the state `" + state.name +
                            "` has the name of a constant, function or "
                            "grid component",
                        state.position);
    }
    if (defNames.count(state.name) == 0) {
      throw SderivError("state",
                        "the state `" + state.name +
                            "` is never updated; define its next value "
                            "with `def " +
                            state.name + " = ...` or `out def " + state.name + " = ...`",
                        state.position);
    }
  }
  for (const SurfaceId id : parsed.defs) {
    if (stateNames.count(surface[id].text) != 0 && !surface[id].params.empty()) {
      throw SderivError("state",
                        "the next value of the state `" + surface[id].text +
                            "` cannot take parameters",
                        surface[id].position);
    }
  }

  std::vector<GridId> gridIds;
  gridIds.reserve(parsed.grids.size());
  for (const auto& grid : parsed.grids) {
    reader::datafield::GridDesc desc;
    // Grid.h is the contract; the surface vocabularies map onto it here.
    // `components` deliberately does NOT reach the desc: Grid.h declines to
    // carry a component-name table, because Kind::Lookup already holds a
    // resolved integer and a second name-resolution path is the defect class
    // this frontend exists to remove. The names stay in ComponentTable, which
    // is frontend-local and already built above.
    //
    // GridDesc::boundary and GridDesc::timeAxis come from the form with attributes; the short
    // form clamps and has no time axis.
    desc.kind = *reader::datafield::parseGridKind(grid.kind);
    desc.path = grid.file;
    desc.variable = grid.variable;
    desc.interpolation = *reader::datafield::parseInterpolation(grid.interpolation);
    desc.boundary = *reader::datafield::parseBoundaryMode(grid.boundary);
    desc.timeAxis = grid.timeAxis;
    // internGrid dedupes, so two outputs naming the same grid share a GridId
    // and the backend loads the file once.
    gridIds.push_back(program.internGrid(desc));
  }

  Lowering lowering(surface, parsed, components, gridIds, options.inputs, program);

  // the expressions of the outputs and of the next values of the states, whose definitions are
  // lowered first
  std::vector<SurfaceId> rootExpressions;
  if (externalName != nullptr) {
    rootExpressions.push_back(parsed.expression);
  } else {
    for (const SurfaceId id : parsed.exports) {
      rootExpressions.push_back(surface[id].a);
    }
  }
  for (const auto& state : parsed.states) {
    for (const SurfaceId id : parsed.defs) {
      if (surface[id].text == state.name) {
        rootExpressions.push_back(surface[id].a);
      }
    }
  }
  lowering.prepare(rootExpressions);

  // Roots first, names after, so `channels()` has seen every output before the
  // channel order is folded in.
  std::vector<std::pair<std::string, NodeId>> roots;
  if (externalName != nullptr) {
    roots.emplace_back(*externalName, lowering.lower(parsed.expression));
  } else {
    for (const SurfaceId id : parsed.exports) {
      // Lowered through the ordinary Name path, so an `out def` reads exactly
      // like any other definition -- including from a later output that names
      // it. Interning makes the shared subexpression one node either way.
      roots.emplace_back(surface[id].text, lowering.lower(surface[id].a));
    }
  }

  // The next values of the states, lowered before the channels are collected: a definition that
  // is not an output may read channels nothing else does.
  std::vector<NodeId> stateRoots;
  for (const auto& state : parsed.states) {
    for (const SurfaceId id : parsed.defs) {
      if (surface[id].text == state.name) {
        stateRoots.push_back(lowering.lower(surface[id].a));
      }
    }
  }

  for (const auto& name : lowering.channels()) {
    if (std::find(channelOrder.begin(), channelOrder.end(), name) == channelOrder.end()) {
      channelOrder.push_back(name);
    }
  }
  for (auto& [name, root] : roots) {
    program.addOutput(name, type, root);
  }
  for (std::size_t i = 0; i < parsed.states.size(); ++i) {
    program.addState(
        parsed.states[i].name, parsed.states[i].initial, stateRoots[i], parsed.states[i].type);
  }
}

} // namespace

Program compileSderiv(const std::vector<SderivOutput>& outputs) {
  Program program;
  std::vector<std::string> channelOrder;
  for (const auto& output : outputs) {
    compileSource(output.source, &output.name, output.type, {}, program, channelOrder);
  }
  for (const auto& name : channelOrder) {
    program.addInput(name, reader::scripting::DataType::F64);
  }
  validate(program);
  return program;
}

Program compileSderiv(const std::string& source, const std::string& outputName) {
  return compileSderiv({SderivOutput{outputName, source, reader::scripting::DataType::F64}});
}

Program compileSderivModule(const std::string& source) {
  return compileSderivModule(source, SderivOptions{});
}

Program compileSderivModule(const std::string& source, const SderivOptions& options) {
  Program program;
  std::vector<std::string> channelOrder;
  compileSource(source, nullptr, reader::scripting::DataType::F64, options, program, channelOrder);
  for (const auto& name : channelOrder) {
    program.addInput(name, reader::scripting::DataType::F64);
  }
  validate(program);
  return program;
}

} // namespace seissol::expr
