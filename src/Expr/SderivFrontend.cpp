// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/SderivFrontend.h"

#include "Expr/Core.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
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
  Check,
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

    const auto digit = [&](std::size_t at) {
      return at < n && std::isdigit(static_cast<unsigned char>(source[at])) != 0;
    };
    if (c == '0' && i + 1 < n && (source[i + 1] == 'x' || source[i + 1] == 'X')) {
      // hexadecimal floating point, as in C: 0x1.8p+1
      i += 2;
      while (i < n &&
             (std::isxdigit(static_cast<unsigned char>(source[i])) != 0 || source[i] == '.')) {
        ++i;
      }
      if (i < n && (source[i] == 'p' || source[i] == 'P')) {
        ++i;
        if (i < n && (source[i] == '+' || source[i] == '-')) {
          ++i;
        }
        while (digit(i)) {
          ++i;
        }
      }
      tokens.push_back({TokenKind::Num, source.substr(start, i - start), static_cast<int>(start)});
      continue;
    }
    if (std::isdigit(c) != 0 || (c == '.' && digit(i + 1))) {
      while (digit(i)) {
        ++i;
      }
      // `2.5`, `.5` and `2.`, but not the point of `2.x`
      if (i < n && source[i] == '.' &&
          (digit(i + 1) || i + 1 >= n ||
           std::isalpha(static_cast<unsigned char>(source[i + 1])) == 0 || source[i + 1] == 'e' ||
           source[i + 1] == 'E')) {
        ++i;
        while (digit(i)) {
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
        c == '=' || c == '.' || c == '{' || c == '}' || c == '[' || c == ']') {
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
      if (value == "inf" || value == "nan") {
        kind = TokenKind::Num;
      } else if (value == "let") {
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
      } else if (value == "check") {
        kind = TokenKind::Check;
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
// Record: the fields `params` of the values `args`; Array: the elements `args`; Field: the field
// `text` of `a`; Index: the element `b` of `a`.
enum class SurfaceKind : std::uint8_t {
  Num,
  Name,
  Call,
  Binary,
  Unary,
  Let,
  Def,
  Record,
  Array,
  Field,
  Index
};

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
  /// the columns of a table that are its axes
  std::vector<std::string> coordinates;
  /// what a sample beyond the volume does (GridDesc::boundary), and the axis that is time, if any
  /// (GridDesc::timeAxis)
  std::string boundary{"clamp"};
  std::optional<int> timeAxis;
  /// an interpolation and a boundary per axis (where there is more than one), the kinds of the axes
  /// of a file, and the knots and values of an inline grid (flat, as GridData has them)
  std::vector<std::string> interpolations;
  std::vector<std::string> boundaries;
  std::vector<std::string> axes;
  std::vector<std::vector<double>> knots;
  std::vector<double> values;
  int position{0};

  /// What GridDesc reads, before it is checked (datafield::gridDescOf()).
  [[nodiscard]] reader::datafield::GridSpec spec() const {
    reader::datafield::GridSpec spec;
    spec.kind = kind;
    spec.path = file;
    spec.variable = variable;
    spec.components = components;
    spec.coordinates = coordinates;
    spec.interpolation =
        interpolations.empty() ? std::vector<std::string>{interpolation} : interpolations;
    spec.boundary = boundaries.empty() ? std::vector<std::string>{boundary} : boundaries;
    spec.axes = axes;
    spec.knots = knots;
    spec.values = values;
    spec.timeAxis = timeAxis;
    return spec;
  }
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
  // `out e`: expressions whose records give an output per field
  std::vector<SurfaceId> spreads;
  // the outputs in the order of the source: an `out def` (an index of `exports`) or an `out e`
  // (one of `spreads`)
  std::vector<std::pair<bool, std::size_t>> outputOrder;
  std::vector<GridDeclaration> grids;
  // `check c, "message"`: the conditions every point has to meet (Program::checks())
  std::vector<std::pair<SurfaceId, std::string>> checks;
  SurfaceId expression{-1}; // -1 when the module declares its outputs itself
};

// =============================================== catalogs (spec-defined) ====
const std::unordered_map<std::string, double> Constants = {{"pi", M_PI}, {"g", 9.80665}};

// Spelled exactly as the Lua-side ssol.* names, so a model transliterates.
const std::unordered_map<std::string, Fn> Builtins = {
    {"sqrt", Fn::Sqrt},   {"abs", Fn::Abs},       {"exp", Fn::Exp},     {"log", Fn::Log},
    {"log2", Fn::Log2},   {"log10", Fn::Log10},   {"sign", Fn::Sign},   {"floor", Fn::Floor},
    {"ceil", Fn::Ceil},   {"round", Fn::Round},   {"sin", Fn::Sin},     {"cos", Fn::Cos},
    {"tan", Fn::Tan},     {"asin", Fn::Asin},     {"acos", Fn::Acos},   {"sinh", Fn::Sinh},
    {"cosh", Fn::Cosh},   {"tanh", Fn::Tanh},     {"asinh", Fn::Asinh}, {"acosh", Fn::Acosh},
    {"atanh", Fn::Atanh}, {"erf", Fn::Erf},       {"min", Fn::Min},     {"max", Fn::Max},
    {"pow", Fn::Pow},     {"atan2", Fn::Atan2},   {"mod", Fn::Mod},     {"lt", Fn::Lt},
    {"le", Fn::Le},       {"eq", Fn::Eq},         {"land", Fn::And},    {"lor", Fn::Or},
    {"atan", Fn::Atan},   {"select", Fn::Select}, {"rcp", Fn::Rcp}};

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

// The dataset inside an ASAGI file where a declaration names none, and the boundary where it gives
// none: the defaults of the full language, which the core format does not have. The variable is
// AsagiLiteGrid::open's own default (and that of a Lua field_spec), stated once in Grid.h.
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
           peek().kind == TokenKind::Out || peek().kind == TokenKind::State ||
           peek().kind == TokenKind::Check) {
      if (peek().kind == TokenKind::Grid) {
        parsed.grids.push_back(gridDeclaration());
        continue;
      }
      if (peek().kind == TokenKind::Check) {
        eat(TokenKind::Check);
        const SurfaceId condition = expression();
        eat(TokenKind::Op, ",");
        const Token& message = eat(TokenKind::Str);
        parsed.checks.emplace_back(condition, message.value);
        continue;
      }
      if (peek().kind == TokenKind::State) {
        parsed.states.push_back(stateDeclaration());
        continue;
      }
      const bool exported = peek().kind == TokenKind::Out;
      if (exported) {
        eat(TokenKind::Out);
        // `out e`: the fields of a record, each an output of its name
        if (peek().kind != TokenKind::Def) {
          parsed.outputOrder.emplace_back(true, parsed.spreads.size());
          parsed.spreads.push_back(expression());
          continue;
        }
      }
      const SurfaceId id = definition(exported);
      parsed.defs.push_back(id);
      if (exported) {
        parsed.outputOrder.emplace_back(false, parsed.exports.size());
        parsed.exports.push_back(id);
      }
    }

    // A module is EITHER a trailing expression, whose name the caller supplies,
    // OR a set of `out def`s that name themselves. Allowing both would leave it
    // open which one is "the" output, and a file that grew an `out def` while
    // keeping its trailing expression would silently gain a second one.
    if (parsed.exports.empty() && parsed.spreads.empty()) {
      parsed.expression = expression();
    }
    eat(TokenKind::Eof);
    return parsed;
  }

  private:
  const Token& peek() const { return tokens_[index_]; }
  const Token& next() { return tokens_[index_++]; }

  /// The value of a number, with correct rounding (Core.h); an infinity or zero beyond the range of
  /// binary64.
  static double numberAt(const Token& token) {
    if (token.value == "inf") {
      return std::numeric_limits<double>::infinity();
    }
    if (token.value == "nan") {
      return std::numeric_limits<double>::quiet_NaN();
    }
    try {
      return numberOf(token.value);
    } catch (const std::invalid_argument&) {
      throw SderivError("lex", "`" + token.value + "` is no number", token.position);
    }
  }
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

  // 'grid' NAME '(' ATTRIBUTE {',' ATTRIBUTE} ')', with ATTRIBUTE = NAME '=' (STRING | NUMBER |
  // TUPLE): the attributes of the core format (the sderiv spec), of which the variable of an ASAGI
  // file and the boundary have defaults here.
  GridDeclaration gridDeclaration() {
    GridDeclaration grid;
    grid.position = peek().position;
    eat(TokenKind::Grid);
    grid.name = eat(TokenKind::Name).value;
    if (!valueIs("(")) {
      throw SderivError("grid",
                        "grid `" + grid.name + "`: its attributes follow in parentheses, as grid " +
                            grid.name +
                            " (kind = \"asagi\", path = \"model.nc\", components = (\"rho\"), "
                            "interpolation = \"linear\")",
                        peek().position);
    }
    gridAttributes(grid);
    return grid;
  }

  /// A string, or strings in parentheses.
  std::vector<std::string> stringsHere() {
    std::vector<std::string> values;
    if (peek().kind == TokenKind::Op && valueIs("(")) {
      eat(TokenKind::Op, "(");
      values.push_back(eat(TokenKind::Str).value);
      while (valueIs(",")) {
        eat(TokenKind::Op, ",");
        values.push_back(eat(TokenKind::Str).value);
      }
      eat(TokenKind::Op, ")");
    } else {
      values.push_back(eat(TokenKind::Str).value);
    }
    return values;
  }

  /// Numbers in parentheses, each with its sign.
  std::vector<double> numbersHere() {
    std::vector<double> values;
    eat(TokenKind::Op, "(");
    while (true) {
      const bool negative = peek().kind == TokenKind::Op && valueIs("-");
      if (negative) {
        eat(TokenKind::Op, "-");
      }
      const double value = numberAt(eat(TokenKind::Num));
      values.push_back(negative ? -value : value);
      if (!valueIs(",")) {
        break;
      }
      eat(TokenKind::Op, ",");
    }
    eat(TokenKind::Op, ")");
    return values;
  }

  /// A number with its sign, or a tuple of values: `-2.5`, `(1, (2, 3))`. Its mark is its
  /// position in the source.
  reader::datafield::SpelledValue spelledValue(int depth = 0) {
    if (depth > MaxNesting) {
      throw SderivError("parse",
                        "the values are nested deeper than " + std::to_string(MaxNesting) +
                            " levels",
                        peek().position);
    }
    reader::datafield::SpelledValue value;
    value.mark = static_cast<std::size_t>(peek().position);
    if (valueIs("(")) {
      eat(TokenKind::Op, "(");
      value.tuple = true;
      while (true) {
        value.items.push_back(spelledValue(depth + 1));
        if (!valueIs(",")) {
          break;
        }
        eat(TokenKind::Op, ",");
      }
      eat(TokenKind::Op, ")");
      return value;
    }
    const bool negative = peek().kind == TokenKind::Op && valueIs("-");
    if (negative) {
      eat(TokenKind::Op, "-");
    }
    const double magnitude = numberAt(eat(TokenKind::Num));
    value.number = negative ? -magnitude : magnitude;
    return value;
  }

  void gridAttributes(GridDeclaration& grid) {
    const std::vector<std::string> names = {"kind",
                                            "path",
                                            "variable",
                                            "coordinates",
                                            "components",
                                            "axes",
                                            "knots",
                                            "values",
                                            "interpolation",
                                            "boundary",
                                            "time"};
    // the attributes given, where they are
    std::map<std::string, int> given;
    std::optional<reader::datafield::SpelledValue> values;
    eat(TokenKind::Op, "(");
    while (true) {
      const Token& key = eat(TokenKind::Name);
      if (std::find(names.begin(), names.end(), key.value) == names.end()) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: unknown attribute `" + key.value +
                              "`; a grid has kind, path, variable, coordinates, components, axes, "
                              "knots, values, interpolation, boundary and time",
                          key.position);
      }
      if (!given.emplace(key.value, key.position).second) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: the attribute `" + key.value +
                              "` is given twice",
                          key.position);
      }
      eat(TokenKind::Op, "=");
      if (key.value == "components") {
        grid.components = stringsHere();
      } else if (key.value == "coordinates") {
        grid.coordinates = stringsHere();
      } else if (key.value == "interpolation") {
        grid.interpolations = stringsHere();
        grid.interpolation = grid.interpolations.front();
      } else if (key.value == "boundary") {
        grid.boundaries = stringsHere();
        grid.boundary = grid.boundaries.front();
      } else if (key.value == "axes") {
        grid.axes = stringsHere();
      } else if (key.value == "knots") {
        eat(TokenKind::Op, "(");
        while (true) {
          grid.knots.push_back(numbersHere());
          if (!valueIs(",")) {
            break;
          }
          eat(TokenKind::Op, ",");
        }
        eat(TokenKind::Op, ")");
      } else if (key.value == "values") {
        // nested as the knots and the components say, which may come after
        if (!valueIs("(")) {
          throw SderivError("grid",
                            "grid `" + grid.name + "`: expected `(` and the values at the nodes",
                            peek().position);
        }
        values = spelledValue();
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
        } else {
          grid.variable = value;
        }
      }
      if (!valueIs(",")) {
        break;
      }
      eat(TokenKind::Op, ",");
    }
    eat(TokenKind::Op, ")");

    // every grid gives its kind, components and interpolation, and what its kind gives; the
    // variable of an ASAGI file and the boundary have defaults
    if (given.count("kind") == 0) {
      throw SderivError("grid", "grid `" + grid.name + "` gives no kind", grid.position);
    }
    const auto kind = reader::datafield::parseGridKind(grid.kind);
    if (!kind.has_value()) {
      throw SderivError("grid",
                        "grid `" + grid.name + "`: unknown kind \"" + grid.kind +
                            "\"; a grid is \"asagi\", \"scec\", \"table\" or \"inline\"",
                        given.at("kind"));
    }
    const auto& attributes = reader::datafield::gridAttributesOf(*kind);
    const auto& common = reader::datafield::commonGridAttributes();
    for (const auto& [attribute, position] : given) {
      if (std::find(common.begin(), common.end(), attribute) == common.end() &&
          !attributes.has(attribute)) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: a grid of the kind \"" + grid.kind +
                              "\" has no attribute `" + attribute + "`",
                          position);
      }
    }
    if (*kind == reader::datafield::GridKind::AsagiLite && given.count("variable") == 0) {
      grid.variable = DefaultGridVariable;
      given.emplace("variable", grid.position);
    }
    given.emplace("boundary", grid.position);
    std::vector<std::string> required = common;
    required.insert(required.end(), attributes.required.begin(), attributes.required.end());
    for (const auto& name : required) {
      if (given.count(name) == 0) {
        throw SderivError("grid",
                          "grid `" + grid.name + "` gives no " + name + "; a grid of the kind \"" +
                              grid.kind + "\" gives " + requiredList(*kind),
                          grid.position);
      }
    }
    if (values.has_value()) {
      if (const auto problem = reader::datafield::flattenGridValues(
              *values, grid.knots, grid.components, grid.values)) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: the values: " + problem->what,
                          static_cast<int>(problem->mark));
      }
    }
  }

  /// What a grid of the kind `kind` gives, but the defaults, as a list for a message.
  static std::string requiredList(reader::datafield::GridKind kind) {
    std::vector<std::string> names = {"kind", "components", "interpolation"};
    for (const auto& name : reader::datafield::gridAttributesOf(kind).required) {
      if (!(kind == reader::datafield::GridKind::AsagiLite && name == "variable")) {
        names.push_back(name);
      }
    }
    std::string list;
    for (std::size_t i = 0; i < names.size(); ++i) {
      list += (i == 0 ? "" : i + 1 == names.size() ? " and " : ", ") + names[i];
    }
    return list;
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
    const double magnitude = numberAt(eat(TokenKind::Num));
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
    SurfaceId base = postfix();
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

  /// A primary and what selects from it: `p.rho`, `a[0]`, `f(x).rho`, `m[1][2]`.
  SurfaceId postfix() {
    SurfaceId value = primary();
    while (peek().kind == TokenKind::Op && (valueIs(".") || valueIs("["))) {
      const Nesting nesting(*this);
      SurfaceNode node;
      node.position = peek().position;
      node.a = value;
      if (valueIs(".")) {
        eat(TokenKind::Op, ".");
        node.kind = SurfaceKind::Field;
        node.text = eat(TokenKind::Name).value;
      } else {
        eat(TokenKind::Op, "[");
        node.kind = SurfaceKind::Index;
        node.b = expression();
        eat(TokenKind::Op, "]");
      }
      value = arena_.add(std::move(node));
    }
    return value;
  }

  /// The token `ahead` tokens after the next one, or the end.
  [[nodiscard]] const Token& peekAt(std::size_t ahead) const {
    return tokens_[std::min(index_ + ahead, tokens_.size() - 1)];
  }

  SurfaceId primary() {
    const Token& token = peek();
    if (token.kind == TokenKind::Num) {
      next();
      SurfaceNode node;
      node.kind = SurfaceKind::Num;
      node.position = token.position;
      node.number = numberAt(token);
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
      // `m.rho(x, y, z)`, a sample of the component rho of the grid m: a call of `m.rho`, a name no
      // declaration can have; without the call, `p.rho` is the field of a record (postfix())
      if (peek().kind == TokenKind::Op && valueIs(".") && peekAt(1).kind == TokenKind::Name &&
          peekAt(2).kind == TokenKind::Op && peekAt(2).value == "(") {
        eat(TokenKind::Op, ".");
        name += "." + eat(TokenKind::Name).value;
      }
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
    if (token.kind == TokenKind::Op && token.value == "(") {
      eat(TokenKind::Op, "(");
      SurfaceId inner = expression();
      eat(TokenKind::Op, ")");
      return inner;
    }
    // a record, {rho = 2600, mu = 1e10}, its fields in their order, each once
    if (token.kind == TokenKind::Op && token.value == "{") {
      const Nesting nesting(*this);
      SurfaceNode node;
      node.kind = SurfaceKind::Record;
      node.position = token.position;
      eat(TokenKind::Op, "{");
      while (true) {
        const Token& field = eat(TokenKind::Name);
        if (std::find(node.params.begin(), node.params.end(), field.value) != node.params.end()) {
          throw SderivError(
              "parse", "the record has the field `" + field.value + "` twice", field.position);
        }
        node.params.push_back(field.value);
        eat(TokenKind::Op, "=");
        node.args.push_back(expression());
        if (!valueIs(",")) {
          break;
        }
        eat(TokenKind::Op, ",");
      }
      eat(TokenKind::Op, "}");
      return arena_.add(std::move(node));
    }
    // an array, [1, 2, 3]; a matrix is an array of arrays
    if (token.kind == TokenKind::Op && token.value == "[") {
      const Nesting nesting(*this);
      SurfaceNode node;
      node.kind = SurfaceKind::Array;
      node.position = token.position;
      eat(TokenKind::Op, "[");
      while (true) {
        node.args.push_back(expression());
        if (!valueIs(",")) {
          break;
        }
        eat(TokenKind::Op, ",");
      }
      eat(TokenKind::Op, "]");
      return arena_.add(std::move(node));
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
    // one namespace for grids and definitions, as in the core format; a bare name of a grid is
    // none of an input either (Lowering)
    if (defNames.count(grid.name) != 0) {
      throw SderivError(
          "grid", "grid `" + grid.name + "` has the name of a definition", grid.position);
    }
    if (!reader::datafield::parseInterpolation(grid.interpolation).has_value()) {
      throw SderivError("grid",
                        "grid `" + grid.name + "`: unknown interpolation \"" + grid.interpolation +
                            "\"",
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
    for (const auto& component : grid.components) {
      if (!seen.insert(component).second) {
        throw SderivError("grid",
                          "grid `" + grid.name + "`: duplicate component \"" + component + "\"",
                          grid.position);
      }
    }
    // and all a grid has to be, the data of an inline one included
    try {
      static_cast<void>(reader::datafield::gridDescOf(grid.spec()));
    } catch (const std::invalid_argument& error) {
      throw SderivError("grid", "grid `" + grid.name + "`: " + error.what(), grid.position);
    }
    // `m.rho`, which no other declaration can be called
    for (std::size_t c = 0; c < grid.components.size(); ++c) {
      table[grid.name + "." + grid.components[c]] = {g, static_cast<std::int32_t>(c)};
    }
  }
  return table;
}

// ============================================================== lowering ====

// A value of the full language: a number, which is a node of the program, a record of named values,
// or an array of values (a matrix is an array of arrays). Records and arrays exist only here: the
// program has the numbers they hold, each its own, and an operation on them is one per number.
struct Value {
  enum class Shape : std::uint8_t { Number, Record, Array };
  Shape shape{Shape::Number};
  NodeId node{NoNode};
  /// the names of the fields of a record, in their order
  std::vector<std::string> names;
  /// the fields of a record, or the elements of an array
  std::vector<Value> items;

  static Value number(NodeId node) {
    Value value;
    value.node = node;
    return value;
  }
  static Value array(std::vector<Value> items) {
    Value value;
    value.shape = Shape::Array;
    value.items = std::move(items);
    return value;
  }
  [[nodiscard]] bool isNumber() const { return shape == Shape::Number; }
  /// The field `name` of a record, or null.
  [[nodiscard]] const Value* field(const std::string& name) const {
    for (std::size_t i = 0; i < names.size(); ++i) {
      if (names[i] == name) {
        return &items[i];
      }
    }
    return nullptr;
  }
};

/// `value` as briefly as it reads back.
std::string shortest(double value) {
  std::array<char, 32> text{};
  const auto result = std::to_chars(text.data(), text.data() + text.size(), value);
  return {text.data(), result.ptr};
}

/// What `value` is, for a message: "a number", "a record of rho, mu", "an array of 3".
std::string described(const Value& value) {
  switch (value.shape) {
  case Value::Shape::Number:
    return "a number";
  case Value::Shape::Record: {
    std::string text = "a record of ";
    for (std::size_t i = 0; i < value.names.size(); ++i) {
      text += (i == 0 ? "" : ", ") + value.names[i];
    }
    return text;
  }
  case Value::Shape::Array:
    break;
  }
  return "an array of " + std::to_string(value.items.size());
}

// The intrinsics of arrays, of their numbers of arguments. Unlike a builtin, a definition may take
// the name of one, which it then hides: `sum` and `norm` are names a script may well have.
const std::map<std::string, std::size_t> Intrinsics = {{"sum", 1},
                                                       {"dot", 2},
                                                       {"norm", 1},
                                                       {"cross", 2},
                                                       {"matmul", 2},
                                                       {"transpose", 1},
                                                       {"trace", 1},
                                                       {"det", 1},
                                                       {"eye", 1},
                                                       {"len", 1},
                                                       {"range", 1}};

/// The most elements `range(n)` makes, and the largest matrix of `eye(n)`: arrays are unrolled.
constexpr std::size_t MostRange = std::size_t{1} << 16U;
constexpr std::size_t MostEye = 64;

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

  /// The number `id` is, `what` for a message: an output, a next value, a condition.
  NodeId lower(SurfaceId id, const std::string& what) {
    Environment empty;
    return numberOf(surface_[id], lower(id, empty), what);
  }

  /// The value `id` is: a number, a record or an array.
  Value lowerValue(SurfaceId id) {
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
  using Environment = std::map<std::string, Value>;

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
      case SurfaceKind::Record:
      case SurfaceKind::Array:
        for (const SurfaceId arg : node.args) {
          stack.emplace_back(arg, scope);
        }
        break;
      case SurfaceKind::Field:
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Index:
        stack.emplace_back(node.a, scope);
        stack.emplace_back(node.b, scope);
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
      case SurfaceKind::Record:
      case SurfaceKind::Array:
        // the first field first, so it is pushed last
        for (auto arg = node.args.rbegin(); arg != node.args.rend(); ++arg) {
          stack.emplace_back(*arg, scope);
        }
        break;
      case SurfaceKind::Field:
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Index:
        stack.emplace_back(node.b, scope);
        stack.emplace_back(node.a, scope);
        break;
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

  /// The coordinates of a sample of the grid `grid` (an index of the declarations): as many as the
  /// IR allows, and its time axis among them.
  std::vector<NodeId> coordinatesOf(const SurfaceNode& node, std::size_t grid, Environment& env) {
    if (node.args.empty() || node.args.size() > 6) {
      fail(node,
           "`" + node.text + "`: sampled with " + std::to_string(node.args.size()) +
               " coordinates, the IR allows 1..6");
    }
    const auto& timeAxis = parsed_.grids[grid].timeAxis;
    if (timeAxis.has_value() && static_cast<std::size_t>(*timeAxis) >= node.args.size()) {
      fail(node,
           "`" + node.text + "`: sampled with " + std::to_string(node.args.size()) +
               " coordinates, but axis " + std::to_string(*timeAxis) + " is the time");
    }
    std::vector<NodeId> coords;
    coords.reserve(node.args.size());
    for (const SurfaceId arg : node.args) {
      coords.push_back(number(arg, env, "a coordinate of a sample"));
    }
    return coords;
  }

  /// The number `value` is, of `at`; `what` is it, for a message.
  NodeId numberOf(const SurfaceNode& at, const Value& value, const std::string& what) const {
    if (!value.isNumber()) {
      fail(at, what + " is a number, not " + described(value));
    }
    return value.node;
  }

  NodeId number(SurfaceId id, Environment& env, const std::string& what) {
    return numberOf(surface_[id], lower(id, env), what);
  }

  /// `op` of the numbers of `args`, leaf by leaf: a number with every number of a record or an
  /// array, records field by field (those of the same names, in the order of the first), arrays
  /// element by element (those of the same length). `what` names the operation for a message.
  Value leafwise(const SurfaceNode& at,
                 const std::string& what,
                 const std::vector<Value>& args,
                 const std::function<NodeId(const std::vector<NodeId>&)>& op) {
    const Value* form = nullptr;
    for (const auto& arg : args) {
      if (!arg.isNumber()) {
        form = &arg;
        break;
      }
    }
    if (form == nullptr) {
      std::vector<NodeId> nodes;
      nodes.reserve(args.size());
      for (const auto& arg : args) {
        nodes.push_back(arg.node);
      }
      return Value::number(op(nodes));
    }
    Value out;
    out.shape = form->shape;
    out.names = form->names;
    out.items.reserve(form->items.size());
    for (std::size_t i = 0; i < form->items.size(); ++i) {
      std::vector<Value> parts;
      parts.reserve(args.size());
      for (const auto& arg : args) {
        if (arg.isNumber()) {
          parts.push_back(arg);
          continue;
        }
        const Value* part = nullptr;
        if (arg.shape == form->shape && arg.items.size() == form->items.size()) {
          part = arg.shape == Value::Shape::Array ? &arg.items[i] : arg.field(form->names[i]);
        }
        if (part == nullptr) {
          fail(at,
               "the operands of " + what + " do not fit: " + described(*form) + " and " +
                   described(arg));
        }
        parts.push_back(*part);
      }
      out.items.push_back(leafwise(at, what, parts, op));
    }
    return out;
  }

  /// The operation `fn` of `args`, leaf by leaf.
  Value pointwise(const SurfaceNode& at, const std::string& what, Fn fn, std::vector<Value> args) {
    return leafwise(at, what, args, [&](const std::vector<NodeId>& n) {
      auto& arena = program_.arena();
      return n.size() == 1   ? arena.pw(fn, n[0])
             : n.size() == 2 ? arena.pw(fn, n[0], n[1])
                             : arena.pw(fn, n[0], n[1], n[2]);
    });
  }

  /// The value of `node` where the program fixes it -- a number, or operations of such -- or
  /// nothing where it depends on the point.
  [[nodiscard]] std::optional<double> fixedValueOf(NodeId node) const {
    const Node& n = program_.arena()[node];
    if (n.kind == Kind::Const) {
      return n.value;
    }
    if (n.kind != Kind::PW || n.a == NoNode) {
      return std::nullopt;
    }
    const auto a = fixedValueOf(n.a);
    if (!a.has_value()) {
      return std::nullopt;
    }
    switch (n.fn) {
    case Fn::Neg:
      return -*a;
    case Fn::Floor:
      return std::floor(*a);
    case Fn::Ceil:
      return std::ceil(*a);
    default:
      break;
    }
    const auto b = n.b == NoNode ? std::nullopt : fixedValueOf(n.b);
    if (!b.has_value()) {
      return std::nullopt;
    }
    switch (n.fn) {
    case Fn::Add:
      return *a + *b;
    case Fn::Sub:
      return *a - *b;
    case Fn::Mul:
      return *a * *b;
    case Fn::Min:
      return std::min(*a, *b);
    case Fn::Max:
      return std::max(*a, *b);
    default:
      return std::nullopt;
    }
  }

  /// The whole number from `least` to `limit` - 1 that `value` of `at` is, which the program
  /// fixes; `what` is it, for a message.
  std::size_t wholeNumberOf(const SurfaceNode& at,
                            const Value& value,
                            std::size_t least,
                            std::size_t limit,
                            const std::string& what) const {
    const auto fixed = fixedValueOf(numberOf(at, value, what));
    if (!fixed.has_value()) {
      fail(at, what + " is a whole number the program fixes, not one of the point");
    }
    if (!(*fixed >= static_cast<double>(least) && *fixed < static_cast<double>(limit)) ||
        std::floor(*fixed) != *fixed) {
      fail(at,
           what + " is a whole number from " + std::to_string(least) + " to " +
               std::to_string(limit - 1) + ", not " + shortest(*fixed));
    }
    return static_cast<std::size_t>(*fixed);
  }

  /// The numbers of the vector `value`, an array of numbers; `what` takes it.
  std::vector<NodeId>
      vectorOf(const SurfaceNode& at, const std::string& what, const Value& value) const {
    std::vector<NodeId> numbers;
    if (value.shape == Value::Shape::Array) {
      for (const auto& item : value.items) {
        if (!item.isNumber()) {
          numbers.clear();
          break;
        }
        numbers.push_back(item.node);
      }
    }
    if (numbers.empty()) {
      fail(at, what + " takes a vector, an array of numbers, not " + described(value));
    }
    return numbers;
  }

  /// The rows of the matrix `value`, an array of vectors of one length; `what` takes it.
  std::vector<std::vector<NodeId>>
      matrixOf(const SurfaceNode& at, const std::string& what, const Value& value) const {
    const auto refused = [&]() {
      fail(at,
           what + " takes a matrix, an array of arrays of numbers of one length, not " +
               described(value));
    };
    if (value.shape != Value::Shape::Array) {
      refused();
    }
    std::vector<std::vector<NodeId>> rows;
    for (const auto& item : value.items) {
      if (item.shape != Value::Shape::Array) {
        refused();
      }
      rows.push_back(vectorOf(at, what, item));
      if (rows.back().size() != rows.front().size()) {
        refused();
      }
    }
    return rows;
  }

  static Value vectorValue(const std::vector<NodeId>& numbers) {
    std::vector<Value> items;
    items.reserve(numbers.size());
    for (const NodeId number : numbers) {
      items.push_back(Value::number(number));
    }
    return Value::array(std::move(items));
  }

  static Value matrixValue(const std::vector<std::vector<NodeId>>& rows) {
    std::vector<Value> items;
    items.reserve(rows.size());
    for (const auto& row : rows) {
      items.push_back(vectorValue(row));
    }
    return Value::array(std::move(items));
  }

  /// a[0] b[0] + a[1] b[1] + ..., summed in that order.
  NodeId dotOf(const std::vector<NodeId>& a, const std::vector<NodeId>& b) {
    auto& arena = program_.arena();
    NodeId sum = arena.pw(Fn::Mul, a[0], b[0]);
    for (std::size_t i = 1; i < a.size(); ++i) {
      sum = arena.pw(Fn::Add, sum, arena.pw(Fn::Mul, a[i], b[i]));
    }
    return sum;
  }

  /// The elements of the array `array` folded by `fn` from the first on, leaf by leaf.
  Value folded(const SurfaceNode& at, const std::string& what, const Value& array, Fn fn) {
    if (array.shape != Value::Shape::Array) {
      fail(at, what + " takes an array, not " + described(array));
    }
    Value acc = array.items.front();
    for (std::size_t i = 1; i < array.items.size(); ++i) {
      acc = pointwise(at, what, fn, {acc, array.items[i]});
    }
    return acc;
  }

  /// The intrinsic of arrays `node` calls (Intrinsics), of the values `args`.
  Value intrinsic(const SurfaceNode& node, const std::vector<Value>& args) {
    auto& arena = program_.arena();
    const std::string& name = node.text;
    const std::string what = "`" + name + "`";
    if (name == "sum") {
      return folded(node, what, args[0], Fn::Add);
    }
    if (name == "len") {
      if (args[0].shape != Value::Shape::Array) {
        fail(node, what + " takes an array, not " + described(args[0]));
      }
      return Value::number(arena.konst(static_cast<double>(args[0].items.size())));
    }
    if (name == "range" || name == "eye") {
      const std::size_t n = wholeNumberOf(
          node, args[0], 1, (name == "range" ? MostRange : MostEye) + 1, "the argument of " + what);
      if (name == "range") {
        std::vector<Value> items;
        for (std::size_t i = 0; i < n; ++i) {
          items.push_back(Value::number(arena.konst(static_cast<double>(i))));
        }
        return Value::array(std::move(items));
      }
      std::vector<std::vector<NodeId>> rows(n, std::vector<NodeId>(n, arena.konst(0.0)));
      for (std::size_t i = 0; i < n; ++i) {
        rows[i][i] = arena.konst(1.0);
      }
      return matrixValue(rows);
    }
    if (name == "dot" || name == "cross") {
      const auto a = vectorOf(node, what, args[0]);
      const auto b = vectorOf(node, what, args[1]);
      if (a.size() != b.size() || (name == "cross" && a.size() != 3)) {
        fail(node,
             what + " takes two vectors of " + (name == "cross" ? "3" : "one length") + ", not " +
                 described(args[0]) + " and " + described(args[1]));
      }
      if (name == "dot") {
        return Value::number(dotOf(a, b));
      }
      const auto term = [&](std::size_t i, std::size_t j) {
        return arena.pw(Fn::Sub, arena.pw(Fn::Mul, a[i], b[j]), arena.pw(Fn::Mul, a[j], b[i]));
      };
      return vectorValue({term(1, 2), term(2, 0), term(0, 1)});
    }
    if (name == "norm") {
      const auto a = vectorOf(node, what, args[0]);
      return Value::number(arena.pw(Fn::Sqrt, dotOf(a, a)));
    }
    if (name == "matmul") {
      // a matrix or a vector on either side: a vector on the left is a row, on the right a column
      const auto isMatrix = [](const Value& value) {
        return value.shape == Value::Shape::Array && !value.items.empty() &&
               value.items.front().shape == Value::Shape::Array;
      };
      const bool leftMatrix = isMatrix(args[0]);
      const bool rightMatrix = isMatrix(args[1]);
      if (!leftMatrix && !rightMatrix) {
        fail(node, what + " of two vectors: their product is dot(a, b)");
      }
      const auto left = leftMatrix
                            ? matrixOf(node, what, args[0])
                            : std::vector<std::vector<NodeId>>{vectorOf(node, what, args[0])};
      std::vector<std::vector<NodeId>> right;
      if (rightMatrix) {
        right = matrixOf(node, what, args[1]);
      } else {
        for (const NodeId element : vectorOf(node, what, args[1])) {
          right.push_back({element});
        }
      }
      if (left.front().size() != right.size()) {
        fail(node,
             what + " of a " + std::to_string(left.size()) + " by " +
                 std::to_string(left.front().size()) + " and a " + std::to_string(right.size()) +
                 " by " + std::to_string(right.front().size()) + " operand");
      }
      std::vector<std::vector<NodeId>> product(left.size());
      for (std::size_t i = 0; i < left.size(); ++i) {
        for (std::size_t j = 0; j < right.front().size(); ++j) {
          std::vector<NodeId> column;
          for (const auto& row : right) {
            column.push_back(row[j]);
          }
          product[i].push_back(dotOf(left[i], column));
        }
      }
      if (!rightMatrix) {
        std::vector<NodeId> vector;
        for (const auto& row : product) {
          vector.push_back(row.front());
        }
        return vectorValue(vector);
      }
      return leftMatrix ? matrixValue(product) : vectorValue(product.front());
    }
    // the intrinsics of a matrix
    const auto rows = matrixOf(node, what, args[0]);
    if (name == "transpose") {
      std::vector<std::vector<NodeId>> transposed(rows.front().size());
      for (const auto& row : rows) {
        for (std::size_t j = 0; j < row.size(); ++j) {
          transposed[j].push_back(row[j]);
        }
      }
      return matrixValue(transposed);
    }
    if (rows.size() != rows.front().size()) {
      fail(node,
           what + " takes a square matrix, not " + std::to_string(rows.size()) + " by " +
               std::to_string(rows.front().size()));
    }
    if (name == "trace") {
      NodeId sum = rows[0][0];
      for (std::size_t i = 1; i < rows.size(); ++i) {
        sum = arena.pw(Fn::Add, sum, rows[i][i]);
      }
      return Value::number(sum);
    }
    // det: of 1 by 1 to 3 by 3, by the first row
    const auto mul = [&](NodeId a, NodeId b) { return arena.pw(Fn::Mul, a, b); };
    const auto sub = [&](NodeId a, NodeId b) { return arena.pw(Fn::Sub, a, b); };
    const auto& m = rows;
    switch (rows.size()) {
    case 1:
      return Value::number(m[0][0]);
    case 2:
      return Value::number(sub(mul(m[0][0], m[1][1]), mul(m[0][1], m[1][0])));
    case 3: {
      const NodeId minor0 = sub(mul(m[1][1], m[2][2]), mul(m[1][2], m[2][1]));
      const NodeId minor1 = sub(mul(m[1][0], m[2][2]), mul(m[1][2], m[2][0]));
      const NodeId minor2 = sub(mul(m[1][0], m[2][1]), mul(m[1][1], m[2][0]));
      return Value::number(
          arena.pw(Fn::Add, sub(mul(m[0][0], minor0), mul(m[0][1], minor1)), mul(m[0][2], minor2)));
    }
    default:
      fail(node,
           what + " takes a matrix of 1 by 1 to 3 by 3, not " + std::to_string(rows.size()) +
               " by " + std::to_string(rows.size()));
    }
  }

  /// The numbers of a value, in its order, and its form ("n", "{rho:n}", "[n,n]"): what a call
  /// of a definition depends on.
  static void keyOf(const Value& value, std::vector<NodeId>& leaves, std::string& form) {
    switch (value.shape) {
    case Value::Shape::Number:
      leaves.push_back(value.node);
      form += 'n';
      return;
    case Value::Shape::Record:
      form += '{';
      for (std::size_t i = 0; i < value.items.size(); ++i) {
        form += value.names[i] + ':';
        keyOf(value.items[i], leaves, form);
        form += ',';
      }
      form += '}';
      return;
    case Value::Shape::Array:
      form += '[';
      for (const auto& item : value.items) {
        keyOf(item, leaves, form);
        form += ',';
      }
      form += ']';
      return;
    }
  }

  Value lower(SurfaceId id, Environment& env) {
    const SurfaceNode& node = surface_[id];
    switch (node.kind) {
    case SurfaceKind::Num:
      return Value::number(program_.arena().konst(node.number));

    case SurfaceKind::Unary:
      return pointwise(node, "`-`", Fn::Neg, {lower(node.a, env)});

    case SurfaceKind::Binary: {
      // A sum of many terms is a long chain of binary nodes along their left operands, as the
      // parser builds it; it is walked without recursion, from its first term on.
      std::vector<SurfaceId> spine{id};
      while (surface_[surface_[spine.back()].a].kind == SurfaceKind::Binary) {
        spine.push_back(surface_[spine.back()].a);
      }
      Value value = lower(surface_[spine.back()].a, env);
      for (auto link = spine.rbegin(); link != spine.rend(); ++link) {
        const SurfaceNode& binary = surface_[*link];
        Value b = lower(binary.b, env);
        value = pointwise(binary,
                          "`" + std::string(binary.text == "pow" ? "**" : binary.text) + "`",
                          binaryFn(binary.text),
                          {std::move(value), std::move(b)});
      }
      return value;
    }

    case SurfaceKind::Let: {
      Value value = lower(node.a, env);
      Environment inner = env;
      inner[node.text] = std::move(value);
      // Sharing rather than duplication falls out of interning: the bound id is
      // substituted, so both uses reach the same node without a CSE pass.
      return lower(node.b, inner);
    }

    case SurfaceKind::Record: {
      Value record;
      record.shape = Value::Shape::Record;
      record.names = node.params;
      for (const SurfaceId arg : node.args) {
        record.items.push_back(lower(arg, env));
      }
      return record;
    }

    case SurfaceKind::Array: {
      std::vector<Value> items;
      items.reserve(node.args.size());
      for (const SurfaceId arg : node.args) {
        items.push_back(lower(arg, env));
      }
      return Value::array(std::move(items));
    }

    case SurfaceKind::Field: {
      // `m.rho` of a grid m, without its coordinates
      if (const SurfaceNode& base = surface_[node.a]; base.kind == SurfaceKind::Name &&
                                                      env.count(base.text) == 0 &&
                                                      gridNamed(base.text) != nullptr) {
        const std::string sample = base.text + "." + node.text;
        fail(node,
             "`" + sample + "` is a component of a grid, sampled with coordinates: " + sample +
                 "(x, y, z)");
      }
      const Value record = lower(node.a, env);
      if (record.shape != Value::Shape::Record) {
        fail(node, "`." + node.text + "` is a field of a record, not of " + described(record));
      }
      const Value* field = record.field(node.text);
      if (field == nullptr) {
        fail(node, described(record) + " has no field `" + node.text + "`");
      }
      return *field;
    }

    case SurfaceKind::Index: {
      const Value array = lower(node.a, env);
      if (array.shape != Value::Shape::Array) {
        fail(node, "an index selects an element of an array, not of " + described(array));
      }
      const std::size_t index =
          wholeNumberOf(node, lower(node.b, env), 0, array.items.size(), "the index");
      return array.items[index];
    }

    case SurfaceKind::Name: {
      const auto bound = env.find(node.text);
      if (bound != env.end()) {
        return bound->second;
      }
      // An input the consumer names reads the input, as a state reads its previous value.
      if (inputs_.count(node.text) != 0) {
        return Value::number(channel(node));
      }
      // A state reads its value from the previous call -- also where the definition of the same
      // name, its next value, is in scope. That is what makes the update a parallel assignment.
      if (states_.count(node.text) != 0) {
        return Value::number(program_.arena().field(node.text));
      }
      const auto constant = Constants.find(node.text);
      if (constant != Constants.end()) {
        return Value::number(program_.arena().konst(constant->second));
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
        Value value = lower(definition.a, inner);
        lowered_.emplace(node.text, value);
        return value;
      }
      if (gridNamed(node.text) != nullptr) {
        fail(node,
             "`" + node.text + "` is a grid; a component of it is sampled as " + node.text +
                 ".c(x, y, z), all of them as " + node.text + "(x, y, z)");
      }
      // (the name of an intrinsic of arrays is a name like any other where it is not called)
      if (isFunction(node.text)) {
        fail(node, "`" + node.text + "` is a function and must be called");
      }
      // Everything left is an input channel. No declaration, no catalogue.
      return Value::number(channel(node));
    }

    case SurfaceKind::Call: {
      if (env.count(node.text) != 0) {
        fail(node, "`" + node.text + "` is not callable");
      }

      const auto component = components_.find(node.text);
      if (component != components_.end()) {
        const auto coords = coordinatesOf(node, component->second.grid, env);
        return Value::number(program_.arena().lookup(
            gridIds_[component->second.grid], component->second.component, coords));
      }
      if (const auto dot = node.text.find('.'); dot != std::string::npos) {
        const std::string gridName = node.text.substr(0, dot);
        const GridDeclaration* grid = gridNamed(gridName);
        if (grid == nullptr) {
          fail(node, "`" + gridName + "` is no grid, so `" + node.text + "` samples nothing");
        }
        std::string list;
        for (const auto& each : grid->components) {
          list += (list.empty() ? "" : ", ") + each;
        }
        fail(node,
             "the grid `" + gridName + "` has no component `" + node.text.substr(dot + 1) +
                 "`; it has " + list);
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
        std::vector<NodeId> leaves;
        std::string form;
        for (std::size_t i = 0; i < definition.params.size(); ++i) {
          Value argument = lower(node.args[i], env);
          keyOf(argument, leaves, form);
          form += ';';
          inner[definition.params[i]] = std::move(argument);
        }
        // the value of a call depends on its arguments alone
        auto key = std::make_tuple(node.text, std::move(form), std::move(leaves));
        if (const auto known = calls_.find(key); known != calls_.end()) {
          return known->second;
        }
        const Expansion expansion(*this, node);
        Value value = lower(definition.a, inner);
        calls_.emplace(std::move(key), value);
        return value;
      }

      // m(x, y, z): every component of the grid m, a record of their names
      if (const GridDeclaration* grid = gridNamed(node.text); grid != nullptr) {
        const auto g = static_cast<std::size_t>(grid - parsed_.grids.data());
        const auto coords = coordinatesOf(node, g, env);
        Value record;
        record.shape = Value::Shape::Record;
        for (std::size_t c = 0; c < grid->components.size(); ++c) {
          record.names.push_back(grid->components[c]);
          record.items.push_back(Value::number(
              program_.arena().lookup(gridIds_[g], static_cast<std::int32_t>(c), coords)));
        }
        return record;
      }

      // the arguments, in the order of the source, which the order of the evaluation of the
      // arguments of a call would leave to the compiler
      std::vector<Value> args;
      args.reserve(node.args.size());
      for (const SurfaceId arg : node.args) {
        args.push_back(lower(arg, env));
      }
      const auto arguments = [&](std::size_t wanted) {
        if (args.size() != wanted) {
          fail(node,
               "`" + node.text + "` takes " + std::to_string(wanted) + " argument" +
                   (wanted == 1 ? "" : "s") + ", got " + std::to_string(args.size()));
        }
      };
      const std::string what = "`" + node.text + "`";

      // `lnot` has no operation of its own: it is eq(x, 0), true where x is false, as the
      // specification has it -- also for values other than 0 and 1, and for a NaN, which is true
      // and so not false.
      if (node.text == "lnot") {
        arguments(1);
        return pointwise(node, what, Fn::Eq, {args[0], Value::number(program_.arena().konst(0.0))});
      }

      // a > b is b < a. Rewritten here rather than given its own Fn, so every
      // consumer of the IR keeps one comparison shape to handle.
      if (const auto swapped = SwappedComparisons.find(node.text);
          swapped != SwappedComparisons.end()) {
        arguments(2);
        return pointwise(node, what, swapped->second, {args[1], args[0]});
      }

      // a conversion to another type: f64(x), f32(x), real(x)
      const auto conversion = Conversions.find(node.text);
      if (conversion != Conversions.end()) {
        arguments(1);
        return leafwise(node, what, args, [&](const std::vector<NodeId>& n) {
          return program_.arena().cast(conversion->second, n[0]);
        });
      }

      if (const auto intrinsicArity = Intrinsics.find(node.text);
          intrinsicArity != Intrinsics.end()) {
        arguments(intrinsicArity->second);
        return intrinsic(node, args);
      }

      const auto builtin = Builtins.find(node.text);
      if (builtin == Builtins.end()) {
        fail(node, "unknown function `" + node.text + "`");
      }
      const Fn fn = builtin->second;

      // min and max are variadic in every language a model author is likely to
      // come from, so fold rather than reject; the IR op stays binary. Of one array, they are
      // its least and its largest element.
      if (fn == Fn::Min || fn == Fn::Max) {
        if (args.size() == 1) {
          return folded(node, what, args[0], fn);
        }
        if (args.size() > 2) {
          Value acc = args[0];
          for (std::size_t i = 1; i < args.size(); ++i) {
            acc = pointwise(node, what, fn, {acc, args[i]});
          }
          return acc;
        }
      }
      arguments(static_cast<std::size_t>(arity(fn)));
      return pointwise(node, what, fn, args);
    }

    case SurfaceKind::Def:
      break;
    }
    fail(node, "cannot lower this node");
  }

  /// The grid of the name `name`, or null.
  [[nodiscard]] const GridDeclaration* gridNamed(const std::string& name) const {
    for (const auto& grid : parsed_.grids) {
      if (grid.name == name) {
        return &grid;
      }
    }
    return nullptr;
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
  std::map<std::string, Value> lowered_;

  public:
  /// The node of the definition `name` without parameters, where it was lowered.
  [[nodiscard]] std::optional<NodeId> definition(const std::string& name) const {
    const auto found = lowered_.find(name);
    if (found == lowered_.end() || !found->second.isNumber()) {
      return std::nullopt;
    }
    return found->second.node;
  }

  private:
  /// the value of every call of a definition with parameters so far, by its arguments: their form
  /// and their numbers
  std::map<std::tuple<std::string, std::string, std::vector<NodeId>>, Value> calls_;
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
  if (externalName == nullptr && parsed.exports.empty() && parsed.spreads.empty()) {
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
    auto spec = grid.spec();
    // relative to the directory of the file of the module, as easi and the parameter file have it
    if (spec.kind != "inline") {
      std::filesystem::path path(grid.file);
      if (path.is_relative() && !options.file.empty()) {
        path = (std::filesystem::path(options.file).parent_path() / path).lexically_normal();
      }
      spec.path = path.string();
    }
    desc = reader::datafield::gridDescOf(spec);
    // internGrid dedupes, so two outputs naming the same grid share a GridId
    // and the backend loads the file once. The names are what the core format prints.
    gridIds.push_back(program.internGrid(desc));
    program.nameGrid(gridIds.back(), GridNames{grid.name, grid.components});
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
    for (const SurfaceId id : parsed.spreads) {
      rootExpressions.push_back(id);
    }
  }
  for (const auto& state : parsed.states) {
    for (const SurfaceId id : parsed.defs) {
      if (surface[id].text == state.name) {
        rootExpressions.push_back(surface[id].a);
      }
    }
  }
  for (const auto& check : parsed.checks) {
    rootExpressions.push_back(check.first);
  }
  lowering.prepare(rootExpressions);

  // Roots first, names after, so `channels()` has seen every output before the
  // channel order is folded in.
  std::vector<std::pair<std::string, NodeId>> roots;
  // and where they are, for their `meta line`
  std::vector<int> rootPositions;
  if (externalName != nullptr) {
    roots.emplace_back(*externalName,
                       lowering.lower(parsed.expression, "the expression of the module"));
  } else {
    for (const auto& [spread, index] : parsed.outputOrder) {
      if (!spread) {
        // Lowered through the ordinary Name path, so an `out def` reads exactly
        // like any other definition -- including from a later output that names
        // it. Interning makes the shared subexpression one node either way.
        const SurfaceId id = parsed.exports[index];
        roots.emplace_back(surface[id].text,
                           lowering.lower(surface[id].a, "the output `" + surface[id].text + "`"));
        rootPositions.push_back(surface[id].position);
        continue;
      }
      // `out e`: a field of the record e is an output of its name
      const SurfaceId id = parsed.spreads[index];
      const Value record = lowering.lowerValue(id);
      if (record.shape != Value::Shape::Record) {
        throw SderivError("resolve",
                          "`out` of an expression gives the fields of a record as outputs, and "
                          "this is " +
                              described(record) + "; `out def v = e` gives one output",
                          surface[id].position);
      }
      for (std::size_t i = 0; i < record.items.size(); ++i) {
        if (!record.items[i].isNumber()) {
          throw SderivError("resolve",
                            "the output `" + record.names[i] + "` is a number, not " +
                                described(record.items[i]),
                            surface[id].position);
        }
        roots.emplace_back(record.names[i], record.items[i].node);
        rootPositions.push_back(surface[id].position);
      }
    }
  }

  // The next values of the states, lowered before the channels are collected: a definition that
  // is not an output may read channels nothing else does.
  std::vector<NodeId> stateRoots;
  for (const auto& state : parsed.states) {
    for (const SurfaceId id : parsed.defs) {
      if (surface[id].text == state.name) {
        stateRoots.push_back(
            lowering.lower(surface[id].a, "the next value of the state `" + state.name + "`"));
      }
    }
  }
  // and those of the checks
  std::vector<NodeId> checkRoots;
  for (const auto& check : parsed.checks) {
    checkRoots.push_back(lowering.lower(check.first, "the condition of a check"));
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
  for (std::size_t i = 0; i < parsed.checks.size(); ++i) {
    program.addCheck(checkRoots[i], parsed.checks[i].second);
  }

  // where they come from (`meta line` of the core format): the definitions, the states, the
  // outputs and the conditions of the checks
  const auto at = [&](int offset) {
    const CorePosition position =
        positionOf(source, static_cast<std::size_t>(std::max(offset, 0)), options.file);
    return SourcePosition{position.file, position.line, position.column};
  };
  for (const SurfaceId id : parsed.defs) {
    if (!surface[id].params.empty()) {
      continue;
    }
    if (const auto node = lowering.definition(surface[id].text); node.has_value()) {
      program.addSource(
          SourceOf{SourceOf::What::Node, *node, std::string(), at(surface[id].position)});
    }
  }
  if (externalName == nullptr) {
    for (std::size_t i = 0; i < roots.size(); ++i) {
      program.addSource(
          SourceOf{SourceOf::What::Output, NoNode, roots[i].first, at(rootPositions[i])});
    }
  }
  for (const auto& state : parsed.states) {
    program.addSource(SourceOf{SourceOf::What::State, NoNode, state.name, at(state.position)});
  }
  for (std::size_t i = 0; i < parsed.checks.size(); ++i) {
    program.addSource(SourceOf{SourceOf::What::Node,
                               checkRoots[i],
                               std::string(),
                               at(surface[parsed.checks[i].first].position)});
  }
}

} // namespace

/// `compile`, with the byte offset of an error turned into its line and column in `source`.
template <typename F>
void positioned(const std::string& source, const std::string& file, F&& compile) {
  try {
    compile();
  } catch (const CoreError&) {
    throw;
  } catch (const SderivError& error) {
    if (error.position() < 0) {
      throw;
    }
    throw CoreError(positionOf(source, static_cast<std::size_t>(error.position()), file),
                    error.what());
  }
}

Program compileSderiv(const std::vector<SderivOutput>& outputs) {
  Program program;
  std::vector<std::string> channelOrder;
  for (const auto& output : outputs) {
    positioned(output.source, "", [&] {
      compileSource(output.source, &output.name, output.type, {}, program, channelOrder);
    });
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
  if (isCoreText(source)) {
    return loadCore(source, options.file);
  }
  Program program;
  std::vector<std::string> channelOrder;
  positioned(source, options.file, [&] {
    compileSource(
        source, nullptr, reader::scripting::DataType::F64, options, program, channelOrder);
  });
  for (const auto& name : channelOrder) {
    program.addInput(name, reader::scripting::DataType::F64);
  }
  validate(program);
  return program;
}

} // namespace seissol::expr
