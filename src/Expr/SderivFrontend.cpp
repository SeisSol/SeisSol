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
#include "Expr/SderivLibrary.h"
#include "Reader/Datafield/Grid.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
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
  Dataset,
  State,
  Check,
  Type,
  For,
  Import,
  If,
  Then,
  Else,
  And,
  Or,
  Not,
  Meta,
  Next,
  Eof
};

struct Token {
  TokenKind kind{};
  std::string value;
  int position{0};
  /// a name in backquotes, which is never a keyword
  bool escaped{false};
};

const char* describe(TokenKind kind) {
  switch (kind) {
  case TokenKind::Num:
    return "a number";
  case TokenKind::Str:
    return "a quoted string";
  case TokenKind::Name:
    return "a name";
  case TokenKind::In:
    return "`in`";
  case TokenKind::Then:
    return "`then`";
  case TokenKind::Else:
    return "`else`";
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

    // A name in backquotes, as the core format escapes it: any character but a backquote and a line
    // end, never a keyword -- `in` is a name there.
    if (c == '`') {
      ++i;
      const std::size_t begin = i;
      while (i < n && source[i] != '`' && source[i] != '\n' && source[i] != '\r') {
        ++i;
      }
      if (i >= n || source[i] != '`') {
        throw SderivError(
            "lex", "an escaped name ends with ` on its line", static_cast<int>(start));
      }
      if (i == begin) {
        throw SderivError("lex", "an escaped name is not empty", static_cast<int>(start));
      }
      tokens.push_back(
          {TokenKind::Name, source.substr(begin, i - begin), static_cast<int>(start), true});
      ++i;
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
    // the comparisons: `<`, `<=`, `>`, `>=`, `==`, `!=`
    if (c == '<' || c == '>' || ((c == '=' || c == '!') && i + 1 < n && source[i + 1] == '=')) {
      const bool equal = i + 1 < n && source[i + 1] == '=';
      tokens.push_back(
          {TokenKind::Op, source.substr(start, equal ? 2 : 1), static_cast<int>(start)});
      i += equal ? 2 : 1;
      continue;
    }
    if (c == '!') {
      throw SderivError(
          "lex", "unexpected character '!'; the negation is `not`", static_cast<int>(i));
    }
    if (c == '+' || c == '-' || c == '*' || c == '/' || c == '(' || c == ')' || c == ',' ||
        c == '=' || c == '.' || c == '{' || c == '}' || c == '[' || c == ']' || c == ':' ||
        c == '@') {
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
      } else if (value == "dataset") {
        kind = TokenKind::Dataset;
      } else if (value == "state") {
        kind = TokenKind::State;
      } else if (value == "check") {
        kind = TokenKind::Check;
      } else if (value == "type") {
        kind = TokenKind::Type;
      } else if (value == "for") {
        kind = TokenKind::For;
      } else if (value == "import") {
        kind = TokenKind::Import;
      } else if (value == "if") {
        kind = TokenKind::If;
      } else if (value == "then") {
        kind = TokenKind::Then;
      } else if (value == "else") {
        kind = TokenKind::Else;
      } else if (value == "and") {
        kind = TokenKind::And;
      } else if (value == "or") {
        kind = TokenKind::Or;
      } else if (value == "not") {
        kind = TokenKind::Not;
      } else if (value == "meta") {
        kind = TokenKind::Meta;
      } else if (value == "next") {
        kind = TokenKind::Next;
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
// `text` of `a`; Index: the element `b` of `a`; Ascription: `a` as the type `type`; Comprehension:
// the array of `a` for every element `text` of the array `b`; Conditional: `if args[0] then
// args[1] else args[2]`.
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
  Index,
  Ascription,
  Comprehension,
  Conditional
};

using TypeId = int;
constexpr TypeId NoType = -1;

struct SurfaceNode {
  SurfaceKind kind{};
  double number{0};
  std::string text; // Name.id / Call.fn / Binary.op / Unary.op / Let.name / Def.name
  std::vector<SurfaceId> args;
  std::vector<std::string> params;
  SurfaceId a{-1};
  SurfaceId b{-1};
  int position{0};
  /// of an Ascription its type; of a Def the type of its value, and of its parameters
  TypeId type{NoType};
  std::vector<TypeId> paramTypes;
  /// the module the node is of: 0 the module compiled, else one it imports (Lowering)
  int module{0};
};

// A type as a module spells it (`type t = ...`, `e : t`): a number, the type of a name, a record of
// typed fields, or an array of a length of a type.
struct TypeExpression {
  enum class Kind : std::uint8_t { Number, Name, Record, Array };
  Kind kind{Kind::Number};
  std::string name;
  /// of a record its fields, with `items` their types; of an array, `items` its element type
  std::vector<std::string> fields;
  std::vector<TypeId> items;
  std::size_t length{0};
  int position{0};
  int module{0};
};

class SurfaceArena {
  public:
  SurfaceId add(SurfaceNode&& node) {
    node.module = module_;
    nodes_.push_back(std::move(node));
    return static_cast<SurfaceId>(nodes_.size()) - 1;
  }
  /// The module the nodes and types added from now on are of.
  void setModule(int module) { module_ = module; }
  const SurfaceNode& operator[](SurfaceId id) const { return nodes_[id]; }
  TypeId addType(TypeExpression&& type) {
    type.module = module_;
    types_.push_back(std::move(type));
    return static_cast<TypeId>(types_.size()) - 1;
  }
  [[nodiscard]] const TypeExpression& type(TypeId id) const {
    return types_[static_cast<std::size_t>(id)];
  }

  private:
  std::vector<SurfaceNode> nodes_;
  std::vector<TypeExpression> types_;
  int module_{0};
};

// `import stress` (a module of the library of SeisSol) or `import "file.sderiv" as name` (a file,
// relative to the importing one): the definitions and types of the module, as name.f, name.t.
struct ImportDeclaration {
  /// the name it is imported as
  std::string name;
  /// the file, or the module of the library
  std::string path;
  bool file{false};
  int position{0};
};

// `type NAME = ...`: a name for a type, which stands for it wherever it is named.
struct TypeDeclaration {
  std::string name;
  TypeId type{NoType};
  int position{0};
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
    if (kind == "cells") {
      // data on the mesh, which has none: they are an error of gridDescOf() where they are given
      spec.interpolation = interpolations;
      spec.boundary = boundaries;
    } else {
      spec.interpolation =
          interpolations.empty() ? std::vector<std::string>{interpolation} : interpolations;
      spec.boundary = boundaries.empty() ? std::vector<std::string>{boundary} : boundaries;
    }
    spec.axes = axes;
    spec.knots = knots;
    spec.values = values;
    spec.timeAxis = timeAxis;
    return spec;
  }
};

// A value kept per point across calls. Its next value is the definition of the same name, and
// every reference to the name reads the value from the previous call.
// The initial value of a state: a number, or a record or an array of them, which make a state of
// the program per number (`s[0]`, `s.v`).
struct InitialValue {
  enum class Shape : std::uint8_t { Number, Record, Array };
  Shape shape{Shape::Number};
  double number{0.0};
  std::vector<std::string> names;
  std::vector<InitialValue> items;
};

struct StateDeclaration {
  std::string name;
  /// its initial value, an expression of numbers the program fixes (Lowering::initializeStates())
  SurfaceId initial{-1};
  int position{0};
  // `state f64 name = ...`: the type of the state (see ValueType)
  ValueType type{ValueType::Real};
};

// `in x, y`: an input, declared as the core format declares it.
struct InputDeclaration {
  std::string name;
  int position{0};
};

// `out v = e`: an output of the name v, which is no definition, as the core format has it.
struct NamedOutput {
  std::string name;
  SurfaceId value{-1};
  int position{0};
};

// `meta key = "text", 12`: a line of provenance, as the core format has it.
struct MetaDeclaration {
  MetaEntry entry;
  int position{0};
};

// What an output of the source is: an `out def`, an `out e` of a record, or an `out v = e`.
enum class OutputKind : std::uint8_t { Definition, Spread, Named };

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
  // `out v = e`
  std::vector<NamedOutput> named;
  // the outputs in the order of the source: an `out def` (an index of `exports`), an `out e` (one
  // of `spreads`) or an `out v = e` (one of `named`)
  std::vector<std::pair<OutputKind, std::size_t>> outputOrder;
  // `in x, y`: the inputs the module declares, in their order; a module that declares any has no
  // others
  std::vector<InputDeclaration> inputs;
  // `next s = e`, which is `def s = e` of a state s: the names, and where
  std::vector<std::pair<std::string, int>> nexts;
  // `meta`, the lines of provenance before the declarations
  std::vector<MetaDeclaration> meta;
  std::vector<GridDeclaration> grids;
  std::vector<TypeDeclaration> types;
  std::vector<ImportDeclaration> imports;
  // `check c, "message"`: the conditions every point has to meet (Program::checks())
  std::vector<std::pair<SurfaceId, std::string>> checks;
  SurfaceId expression{-1}; // -1 when the module declares its outputs itself
};

// =============================================== catalogs (spec-defined) ====
// `g` is none: the gravitational acceleration is that of the run, an input
// (ScriptSettings::runInputs)
const std::unordered_map<std::string, double> Constants = {{"pi", M_PI}};

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

  /// A parser of a module imported as `qualifier`: a library, of definitions, types and imports
  /// only, whose definitions and types are named `qualifier.name`.
  Parser(const std::vector<Token>& tokens, SurfaceArena& arena, std::string qualifier)
      : tokens_(tokens), arena_(arena), qualifier_(std::move(qualifier)) {}

  ParsedProgram program() {
    ParsedProgram parsed;
    // def and grid interleave freely and need no separator, for the same reason
    // def already needed none: each starts with its own token kind, so the
    // greedy expression in a def body stops cleanly at the next declaration.
    bool metaAllowed = true;
    while (peek().kind == TokenKind::Def || peek().kind == TokenKind::Dataset ||
           peek().kind == TokenKind::Out || peek().kind == TokenKind::State ||
           peek().kind == TokenKind::Check || peek().kind == TokenKind::Type ||
           peek().kind == TokenKind::Import || peek().kind == TokenKind::In ||
           peek().kind == TokenKind::Next || peek().kind == TokenKind::Meta) {
      if (!qualifier_.empty() && peek().kind != TokenKind::Def && peek().kind != TokenKind::Type &&
          peek().kind != TokenKind::Import) {
        throw SderivError("parse",
                          "a library declares definitions, types and imports, not `" +
                              peek().value + "`",
                          peek().position);
      }
      // the declarations of the core format: `meta` before the others, `in`, `next`
      if (peek().kind == TokenKind::Meta) {
        if (!metaAllowed) {
          throw SderivError("parse", "meta lines come before the declarations", peek().position);
        }
        parsed.meta.push_back(metaDeclaration());
        continue;
      }
      metaAllowed = false;
      if (peek().kind == TokenKind::In) {
        // `in x, y`
        eat(TokenKind::In);
        while (true) {
          const Token& input = eat(TokenKind::Name);
          parsed.inputs.push_back(InputDeclaration{input.value, input.position});
          if (!(peek().kind == TokenKind::Op && valueIs(","))) {
            break;
          }
          eat(TokenKind::Op, ",");
        }
        continue;
      }
      if (peek().kind == TokenKind::Next) {
        // `next s = e`: the next value of the state s, as `def s = e`
        SurfaceNode node;
        node.kind = SurfaceKind::Def;
        node.position = peek().position;
        eat(TokenKind::Next);
        node.text = eat(TokenKind::Name).value;
        eat(TokenKind::Op, "=");
        node.a = expression();
        parsed.nexts.emplace_back(node.text, node.position);
        parsed.defs.push_back(arena_.add(std::move(node)));
        continue;
      }
      if (peek().kind == TokenKind::Import) {
        // `import stress`, `import "lib.sderiv" as lib`
        ImportDeclaration declaration;
        declaration.position = peek().position;
        eat(TokenKind::Import);
        if (peek().kind == TokenKind::Str) {
          declaration.file = true;
          declaration.path = eat(TokenKind::Str).value;
          const Token& as = eat(TokenKind::Name);
          if (as.value != "as" || as.escaped) {
            throw SderivError("parse",
                              "a file is imported as a name: import \"" + declaration.path +
                                  "\" as lib",
                              as.position);
          }
          declaration.name = eat(TokenKind::Name).value;
        } else {
          // `import stress`, or `import stress as s`
          declaration.path = eat(TokenKind::Name).value;
          declaration.name = declaration.path;
          if (peek().kind == TokenKind::Name && valueIs("as") && !peek().escaped) {
            next();
            declaration.name = eat(TokenKind::Name).value;
          }
        }
        parsed.imports.push_back(std::move(declaration));
        continue;
      }
      if (peek().kind == TokenKind::Dataset) {
        parsed.grids.push_back(gridDeclaration());
        continue;
      }
      if (peek().kind == TokenKind::Type) {
        // `type elastic = {rho, mu, lambda}`
        TypeDeclaration declaration;
        declaration.position = peek().position;
        eat(TokenKind::Type);
        declaration.name = qualified(eat(TokenKind::Name).value);
        eat(TokenKind::Op, "=");
        declaration.type = typeExpression();
        parsed.types.push_back(std::move(declaration));
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
        // `out v = e`: the output v, which is no definition (`=` is no operator of an expression)
        if (peek().kind == TokenKind::Name && peekAt(1).kind == TokenKind::Op &&
            peekAt(1).value == "=") {
          const Token& output = next();
          eat(TokenKind::Op, "=");
          parsed.outputOrder.emplace_back(OutputKind::Named, parsed.named.size());
          parsed.named.push_back(NamedOutput{output.value, expression(), output.position});
          continue;
        }
        // `out e`: the fields of a record, each an output of its name
        if (peek().kind != TokenKind::Def) {
          parsed.outputOrder.emplace_back(OutputKind::Spread, parsed.spreads.size());
          parsed.spreads.push_back(expression());
          continue;
        }
      }
      const SurfaceId id = definition(exported);
      parsed.defs.push_back(id);
      if (exported) {
        parsed.outputOrder.emplace_back(OutputKind::Definition, parsed.exports.size());
        parsed.exports.push_back(id);
      }
    }

    // A module is EITHER a trailing expression, whose name the caller supplies,
    // OR a set of `out def`s that name themselves. Allowing both would leave it
    // open which one is "the" output, and a file that grew an `out def` while
    // keeping its trailing expression would silently gain a second one.
    if (parsed.exports.empty() && parsed.spreads.empty() && parsed.named.empty() &&
        qualifier_.empty()) {
      parsed.expression = expression();
    }
    if (!qualifier_.empty() && peek().kind != TokenKind::Eof) {
      throw SderivError(
          "parse", "a library declares definitions, types and imports only", peek().position);
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
      const bool word = token.kind != TokenKind::Num && token.kind != TokenKind::Str &&
                        token.kind != TokenKind::Pow && token.kind != TokenKind::Op &&
                        token.kind != TokenKind::Name && token.kind != TokenKind::Eof;
      throw SderivError("parse",
                        "expected " + expected + ", got '" + token.value + "'" +
                            (kind == TokenKind::Name && word
                                 ? "; `" + token.value +
                                       "` is a reserved word, which an escaped name, `" +
                                       token.value + "`, avoids"
                                 : std::string()),
                        token.position);
    }
    return next();
  }

  /// `meta key = "text", 12`: a line of provenance, its values strings or whole numbers.
  MetaDeclaration metaDeclaration() {
    MetaDeclaration declaration;
    declaration.position = peek().position;
    eat(TokenKind::Meta);
    declaration.entry.key = eat(TokenKind::Name).value;
    eat(TokenKind::Op, "=");
    while (true) {
      const Token& value = peek();
      const bool whole = value.kind == TokenKind::Num && !value.value.empty() &&
                         std::all_of(value.value.begin(), value.value.end(), [](char c) {
                           return std::isdigit(static_cast<unsigned char>(c)) != 0;
                         });
      if (value.kind != TokenKind::Str && !whole) {
        throw SderivError(
            "parse", "the value of a meta line is a string or a whole number", value.position);
      }
      declaration.entry.values.push_back(MetaValue{next().value, value.kind == TokenKind::Str});
      if (!(peek().kind == TokenKind::Op && valueIs(","))) {
        break;
      }
      eat(TokenKind::Op, ",");
    }
    return declaration;
  }

  // 'grid' NAME '(' ATTRIBUTE {',' ATTRIBUTE} ')', with ATTRIBUTE = NAME '=' (STRING | NUMBER |
  // TUPLE): the attributes of the core format (the sderiv spec), of which the variable of an ASAGI
  // file and the boundary have defaults here.
  GridDeclaration gridDeclaration() {
    GridDeclaration grid;
    grid.position = peek().position;
    eat(TokenKind::Dataset);
    grid.name = eat(TokenKind::Name).value;
    if (!valueIs("(")) {
      throw SderivError("dataset",
                        "dataset `" + grid.name +
                            "`: its attributes follow in parentheses, as dataset " + grid.name +
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
        throw SderivError(
            "dataset",
            "dataset `" + grid.name + "`: unknown attribute `" + key.value +
                "`; a dataset has kind, path, variable, coordinates, components, axes, "
                "knots, values, interpolation, boundary and time",
            key.position);
      }
      if (!given.emplace(key.value, key.position).second) {
        throw SderivError("dataset",
                          "dataset `" + grid.name + "`: the attribute `" + key.value +
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
          throw SderivError("dataset",
                            "dataset `" + grid.name + "`: expected `(` and the values at the nodes",
                            peek().position);
        }
        values = spelledValue();
      } else if (key.value == "time") {
        const Token& axis = eat(TokenKind::Num);
        if (axis.value.find_first_not_of("0123456789") != std::string::npos ||
            axis.value.size() > 2 || std::stoi(axis.value) >= MaxLookupDimension) {
          throw SderivError("dataset",
                            "dataset `" + grid.name + "`: the time is an axis, 0 to " +
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
      throw SderivError("dataset", "dataset `" + grid.name + "` gives no kind", grid.position);
    }
    const auto kind = reader::datafield::parseGridKind(grid.kind);
    if (!kind.has_value()) {
      throw SderivError("dataset",
                        "dataset `" + grid.name + "`: unknown kind \"" + grid.kind +
                            "\"; a dataset is \"asagi\", \"scec\", \"table\", \"inline\" or "
                            "\"cells\"",
                        given.at("kind"));
    }
    const auto& attributes = reader::datafield::gridAttributesOf(*kind);
    const auto& common = reader::datafield::commonGridAttributes(*kind);
    for (const auto& [attribute, position] : given) {
      if (std::find(common.begin(), common.end(), attribute) == common.end() &&
          !attributes.has(attribute)) {
        throw SderivError("dataset",
                          "dataset `" + grid.name + "`: a dataset of the kind \"" + grid.kind +
                              "\" has no attribute `" + attribute + "`",
                          position);
      }
    }
    if (*kind == reader::datafield::GridKind::AsagiLite && given.count("variable") == 0) {
      grid.variable = DefaultGridVariable;
      given.emplace("variable", grid.position);
    }
    if (reader::datafield::sampledAtCoordinates(*kind)) {
      given.emplace("boundary", grid.position);
    } else {
      // data on the mesh: the value of the cell, with nothing to interpolate
      grid.boundary.clear();
    }
    std::vector<std::string> required = common;
    required.insert(required.end(), attributes.required.begin(), attributes.required.end());
    for (const auto& name : required) {
      if (given.count(name) == 0) {
        throw SderivError("dataset",
                          "dataset `" + grid.name + "` gives no " + name +
                              "; a dataset of the kind \"" + grid.kind + "\" gives " +
                              requiredList(*kind),
                          grid.position);
      }
    }
    if (values.has_value()) {
      if (const auto problem = reader::datafield::flattenGridValues(
              *values, grid.knots, grid.components, grid.values)) {
        throw SderivError("dataset",
                          "dataset `" + grid.name + "`: the values: " + problem->what,
                          static_cast<int>(problem->mark));
      }
    }
  }

  /// What a dataset of the kind `kind` gives, but the defaults, as a list for a message.
  static std::string requiredList(reader::datafield::GridKind kind) {
    std::vector<std::string> names = {"kind", "components"};
    if (reader::datafield::sampledAtCoordinates(kind)) {
      names.emplace_back("interpolation");
    }
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
    state.initial = expression();
    return state;
  }

  SurfaceId definition(bool exported) {
    const int position = peek().position;
    eat(TokenKind::Def);
    SurfaceNode node;
    node.kind = SurfaceKind::Def;
    node.position = position;
    node.text = qualified(eat(TokenKind::Name).value);
    if (valueIs("(")) {
      eat(TokenKind::Op, "(");
      // a parameter, and its type: `p : velocities`
      const auto parameter = [&]() {
        node.params.push_back(eat(TokenKind::Name).value);
        node.paramTypes.push_back(NoType);
        if (peek().kind == TokenKind::Op && valueIs(":")) {
          eat(TokenKind::Op, ":");
          node.paramTypes.back() = typeExpression();
        }
      };
      if (!valueIs(")")) {
        parameter();
        while (valueIs(",")) {
          eat(TokenKind::Op, ",");
          parameter();
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
    // the type of its value: `def lame(p : velocities) : elastic = ...`
    if (peek().kind == TokenKind::Op && valueIs(":")) {
      eat(TokenKind::Op, ":");
      node.type = typeExpression();
    }
    eat(TokenKind::Op, "=");
    node.a = expression();
    return arena_.add(std::move(node));
  }

  /// A type: `number`, the name of a type, a record of fields `{rho, mu : number, v : vector}`,
  /// each followed by the lengths of the arrays of it, the outermost first: `number[3][3]`.
  TypeId typeExpression() {
    const Nesting nesting(*this);
    TypeExpression type;
    type.position = peek().position;
    if (peek().kind == TokenKind::Op && valueIs("{")) {
      eat(TokenKind::Op, "{");
      type.kind = TypeExpression::Kind::Record;
      while (true) {
        const Token& field = eat(TokenKind::Name);
        if (std::find(type.fields.begin(), type.fields.end(), field.value) != type.fields.end()) {
          throw SderivError(
              "parse", "the type has the field `" + field.value + "` twice", field.position);
        }
        type.fields.push_back(field.value);
        TypeId fieldType = NoType;
        if (peek().kind == TokenKind::Op && valueIs(":")) {
          eat(TokenKind::Op, ":");
          fieldType = typeExpression();
        } else {
          // a field without a type is a number
          TypeExpression number;
          number.position = field.position;
          fieldType = arena_.addType(std::move(number));
        }
        type.items.push_back(fieldType);
        if (!valueIs(",")) {
          break;
        }
        eat(TokenKind::Op, ",");
      }
      eat(TokenKind::Op, "}");
    } else {
      const Token& name = peek();
      if (name.kind != TokenKind::Name) {
        throw SderivError("parse",
                          "expected a type: number, the name of a type, or a record {...}, not '" +
                              name.value + "'",
                          name.position);
      }
      const std::string spelled = next().value;
      if (spelled != "number") {
        type.kind = TypeExpression::Kind::Name;
        type.name = spelled;
        // `stress.normalized`, a type of a module imported as stress
        if (peek().kind == TokenKind::Op && valueIs(".") && peekAt(1).kind == TokenKind::Name) {
          eat(TokenKind::Op, ".");
          type.name += "." + eat(TokenKind::Name).value;
        }
      }
    }
    TypeId id = arena_.addType(std::move(type));
    // arrays of it, the first length the outermost
    std::vector<std::pair<std::size_t, int>> lengths;
    while (peek().kind == TokenKind::Op && valueIs("[")) {
      eat(TokenKind::Op, "[");
      const Token& length = eat(TokenKind::Num);
      if (length.value.find_first_not_of("0123456789") != std::string::npos ||
          length.value.size() > 6 || std::stoul(length.value) == 0) {
        throw SderivError("parse",
                          "the length of an array is a whole number from 1, not " + length.value,
                          length.position);
      }
      eat(TokenKind::Op, "]");
      lengths.emplace_back(std::stoul(length.value), length.position);
    }
    for (auto length = lengths.rbegin(); length != lengths.rend(); ++length) {
      TypeExpression array;
      array.kind = TypeExpression::Kind::Array;
      array.length = length->first;
      array.items = {id};
      array.position = length->second;
      id = arena_.addType(std::move(array));
    }
    return id;
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
    // `if c then a else b`, its last branch as far as an expression goes, as the body of a let
    if (peek().kind == TokenKind::If) {
      SurfaceNode node;
      node.kind = SurfaceKind::Conditional;
      node.position = peek().position;
      eat(TokenKind::If);
      node.args.push_back(expression());
      eat(TokenKind::Then);
      node.args.push_back(expression());
      eat(TokenKind::Else);
      node.args.push_back(expression());
      return arena_.add(std::move(node));
    }
    const SurfaceId value = disjunction();
    // `e : elastic`: e as a value of the type, its fields those of the type
    if (peek().kind == TokenKind::Op && valueIs(":")) {
      SurfaceNode node;
      node.kind = SurfaceKind::Ascription;
      node.position = peek().position;
      eat(TokenKind::Op, ":");
      node.a = value;
      node.type = typeExpression();
      return arena_.add(std::move(node));
    }
    return value;
  }

  /// `a or b`, `a and b` and `not a`, which bind less than the comparisons, as in Python.
  SurfaceId disjunction() {
    SurfaceId lhs = conjunction();
    while (peek().kind == TokenKind::Or) {
      lhs = logical("or", lhs, [this] { return conjunction(); });
    }
    return lhs;
  }

  SurfaceId conjunction() {
    SurfaceId lhs = negation();
    while (peek().kind == TokenKind::And) {
      lhs = logical("and", lhs, [this] { return negation(); });
    }
    return lhs;
  }

  template <typename F>
  SurfaceId logical(const char* op, SurfaceId lhs, const F& operand) {
    SurfaceNode node;
    node.kind = SurfaceKind::Binary;
    node.position = peek().position;
    next();
    node.text = op;
    node.a = lhs;
    node.b = operand();
    return arena_.add(std::move(node));
  }

  SurfaceId negation() {
    const Nesting nesting(*this);
    if (peek().kind == TokenKind::Not) {
      SurfaceNode node;
      node.kind = SurfaceKind::Unary;
      node.position = peek().position;
      node.text = "not";
      next();
      node.a = negation();
      return arena_.add(std::move(node));
    }
    return comparison();
  }

  [[nodiscard]] bool atComparison() const {
    return peek().kind == TokenKind::Op && (valueIs("<") || valueIs("<=") || valueIs(">") ||
                                            valueIs(">=") || valueIs("==") || valueIs("!="));
  }

  /// `a < b` and the other comparisons, of two sums; they do not chain.
  SurfaceId comparison() {
    const SurfaceId lhs = additive();
    if (!atComparison()) {
      return lhs;
    }
    SurfaceNode node;
    node.kind = SurfaceKind::Binary;
    node.position = peek().position;
    node.text = next().value;
    node.a = lhs;
    node.b = additive();
    if (atComparison()) {
      throw SderivError(
          "parse", "comparisons do not chain: `a < b < c` is `a < b and b < c`", peek().position);
    }
    return arena_.add(std::move(node));
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
    // `a @ b` is matmul(a, b), as in Python
    while (peek().kind == TokenKind::Op && (valueIs("*") || valueIs("/") || valueIs("@"))) {
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
          "parse", "string literals are only allowed in a `dataset` declaration", token.position);
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
      node.args.push_back(expression());
      // a comprehension, [f(i) for i in range(3)]: an element per element of an array
      if (peek().kind == TokenKind::For) {
        eat(TokenKind::For);
        node.kind = SurfaceKind::Comprehension;
        node.a = node.args.front();
        node.args.clear();
        node.text = eat(TokenKind::Name).value;
        eat(TokenKind::In);
        node.b = expression();
        eat(TokenKind::Op, "]");
        return arena_.add(std::move(node));
      }
      while (valueIs(",")) {
        eat(TokenKind::Op, ",");
        node.args.push_back(expression());
      }
      eat(TokenKind::Op, "]");
      return arena_.add(std::move(node));
    }
    throw SderivError("parse", "unexpected '" + token.value + "'", token.position);
  }

  /// The name a declaration of the module declares, qualified where the module is imported.
  [[nodiscard]] std::string qualified(const std::string& name) const {
    return qualifier_.empty() ? name : qualifier_ + "." + name;
  }

  const std::vector<Token>& tokens_;
  SurfaceArena& arena_;
  std::string qualifier_;
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
      throw SderivError("dataset", "duplicate dataset `" + grid.name + "`", grid.position);
    }
    // one namespace for grids and definitions, as in the core format; a bare name of a grid is
    // none of an input either (Lowering)
    if (defNames.count(grid.name) != 0) {
      throw SderivError(
          "dataset", "dataset `" + grid.name + "` has the name of a definition", grid.position);
    }
    const bool cells = grid.kind == "cells";
    if (!cells && !reader::datafield::parseInterpolation(grid.interpolation).has_value()) {
      throw SderivError("dataset",
                        "dataset `" + grid.name + "`: unknown interpolation \"" +
                            grid.interpolation + "\"",
                        grid.position);
    }
    if (!cells && !reader::datafield::parseBoundaryMode(grid.boundary).has_value()) {
      throw SderivError("dataset",
                        "dataset `" + grid.name + "`: unknown boundary \"" + grid.boundary +
                            "\"; it is \"clamp\" or \"fail\"",
                        grid.position);
    }
    if (grid.components.empty()) {
      throw SderivError(
          "dataset", "dataset `" + grid.name + "` declares no components", grid.position);
    }
    std::set<std::string> seen;
    for (const auto& component : grid.components) {
      if (!seen.insert(component).second) {
        throw SderivError("dataset",
                          "dataset `" + grid.name + "`: duplicate component \"" + component + "\"",
                          grid.position);
      }
    }
    // and all a grid has to be, the data of an inline one included
    try {
      static_cast<void>(reader::datafield::gridDescOf(grid.spec()));
    } catch (const std::invalid_argument& error) {
      throw SderivError("dataset", "dataset `" + grid.name + "`: " + error.what(), grid.position);
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

// A type, as a value is checked against it (`e : t`): a number, a record of typed fields, or an
// array of a length of a type; and how a message names it.
struct Type {
  enum class Kind : std::uint8_t { Number, Record, Array };
  Kind kind{Kind::Number};
  std::vector<std::string> names;
  /// the types of the fields of a record, or the element type of an array
  std::vector<Type> items;
  std::size_t length{0};
  std::string name{"number"};
};

/// The names of the leaves of a value, a state of the program each: `s`, `s[0]`, `s.v`, `s[1].v`.
void leafNames(const InitialValue& value, const std::string& name, std::vector<std::string>& out) {
  switch (value.shape) {
  case InitialValue::Shape::Number:
    out.push_back(name);
    return;
  case InitialValue::Shape::Record:
    for (std::size_t i = 0; i < value.items.size(); ++i) {
      leafNames(value.items[i], name + "." + value.names[i], out);
    }
    return;
  case InitialValue::Shape::Array:
    for (std::size_t i = 0; i < value.items.size(); ++i) {
      leafNames(value.items[i], name + "[" + std::to_string(i) + "]", out);
    }
    return;
  }
}

/// The numbers of an initial value, in the order of leafNames().
void leafNumbers(const InitialValue& value, std::vector<double>& out) {
  if (value.shape == InitialValue::Shape::Number) {
    out.push_back(value.number);
    return;
  }
  for (const auto& item : value.items) {
    leafNumbers(item, out);
  }
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

/// Whether a call of `name` is a function of the language: a builtin, a conversion, an intrinsic
/// of arrays, a comparison it swaps, `lnot`, `percentile` or `median`.
bool isCalledFunction(const std::string& name) {
  return isFunction(name) || Intrinsics.count(name) != 0 || SwappedComparisons.count(name) != 0 ||
         name == "lnot" || name == "percentile" || name == "median";
}

/// The material of the point: `material.rho` reads the input `rho`, and `material : elastic` the
/// inputs of the fields of the type, as the consumer supplies the parameters of the material.
constexpr const char* MaterialName = "material";

/// The most elements `range(n)` makes, and the largest matrix of `eye(n)`: arrays are unrolled.
constexpr std::size_t MostRange = std::size_t{1} << 16U;
constexpr std::size_t MostEye = 64;

// A module of a compilation: the qualifier of the names of its declarations (empty for the module
// compiled), and the modules it imports, by the names it imports them as.
struct ModuleInfo {
  std::string qualifier;
  std::map<std::string, int> imports;
};

class Lowering {
  public:
  Lowering(const SurfaceArena& surface,
           const ParsedProgram& parsed,
           const ComponentTable& components,
           const std::vector<GridId>& gridIds,
           const std::set<std::string>& inputs,
           const std::vector<ModuleInfo>& modules,
           Program& program)
      : surface_(surface), parsed_(parsed), components_(components), gridIds_(gridIds),
        inputs_(inputs), modules_(modules), program_(program) {
    for (const SurfaceId id : parsed_.defs) {
      defs_[surface_[id].text] = id;
    }
    for (const auto& state : parsed_.states) {
      states_.insert(state.name);
    }
    for (const auto& input : parsed_.inputs) {
      declared_.insert(input.name);
    }
    builtinTypes_ = recordTypes();
  }

  /// The number `id` is, `what` for a message: an output, a next value, a condition.
  NodeId lower(SurfaceId id, const std::string& what) {
    Environment empty;
    return numberOf(surface_[id], lower(id, empty), what);
  }

  /// Evaluates the initial values of the states, numbers or records or arrays of numbers the
  /// program fixes; they give the shapes of the states.
  void initializeStates() {
    for (const auto& state : parsed_.states) {
      Environment empty;
      const SurfaceNode& at = surface_[state.initial];
      stateShapes_.emplace(state.name,
                           initialOf(at,
                                     lower(state.initial, empty),
                                     "the initial value of the state `" + state.name + "`"));
    }
  }

  /// The initial value of the state `name` (initializeStates()).
  [[nodiscard]] const InitialValue& initialOfState(const std::string& name) const {
    return stateShapes_.at(name);
  }

  /// Resolves the types the module declares, so that one that names no type is an error also where
  /// nothing uses it.
  void resolveTypes() const {
    for (const auto& declaration : parsed_.types) {
      static_cast<void>(typeOf(declaration.type));
    }
  }

  /// The value `id` is: a number, a record or an array.
  Value lowerValue(SurfaceId id) {
    Environment empty;
    return lower(id, empty);
  }

  /// The number the definition `id` without parameters gives, as a value of its type where it
  /// gives one; `what` it is, for a message.
  NodeId lowerDefinition(SurfaceId id, const std::string& what) {
    Environment empty;
    return numberOf(surface_[id], definitionValue(surface_[id], empty), what);
  }

  /// The numbers of the value of the definition `id` without parameters, which has the shape of
  /// `shape`, in the order of leafNames(); `what` it is, for a message.
  std::vector<NodeId>
      leavesOfShape(SurfaceId id, const InitialValue& shape, const std::string& what) {
    Environment empty;
    const Value value = definitionValue(surface_[id], empty);
    std::vector<NodeId> leaves;
    collectLeaves(surface_[id], value, shape, what, leaves);
    return leaves;
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
        lowered_.emplace(name, definitionValue(definition, empty));
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
  /// comes first): an input the module declares or the consumer names, a state, a definition, a
  /// constant, or else (where it is not called) a channel. The order is that of lower(): what the
  /// module declares comes before what the language has, so that `in g` is an input, not the
  /// constant, as in the core format.
  enum class Resolution : std::uint8_t { Input, State, Constant, Definition, Other };

  /// The name the name `text` of a node of the module `module` reads, as definitions and types are
  /// registered: in a library its own names qualified (`stress.phi`), `b.f` of a module it imports
  /// as b the definition f of that module.
  [[nodiscard]] std::string qualifiedName(int module, const std::string& text) const {
    const ModuleInfo& info = modules_[static_cast<std::size_t>(module)];
    if (const auto dot = text.find('.'); dot != std::string::npos) {
      const auto imported = info.imports.find(text.substr(0, dot));
      if (imported != info.imports.end()) {
        return modules_[static_cast<std::size_t>(imported->second)].qualifier + "." +
               text.substr(dot + 1);
      }
      return text;
    }
    if (!info.qualifier.empty()) {
      std::string own = info.qualifier + "." + text;
      if (defs_.count(own) != 0 || declaresType(own)) {
        return own;
      }
    }
    return text;
  }

  /// Whether a module declares the type `name` (qualified).
  [[nodiscard]] bool declaresType(const std::string& name) const {
    return std::any_of(parsed_.types.begin(), parsed_.types.end(), [&](const auto& declaration) {
      return declaration.name == name;
    });
  }

  /// Whether `name` is the name a node of the module `module` imports a module as.
  [[nodiscard]] bool importsAs(int module, const std::string& name) const {
    return modules_[static_cast<std::size_t>(module)].imports.count(name) != 0;
  }

  /// The definition `alias.field` of a Field node of a Name `alias` a module is imported as, or
  /// empty.
  [[nodiscard]] std::string importedValue(const SurfaceNode& field,
                                          const std::set<std::string>& bound) const {
    const SurfaceNode& base = surface_[field.a];
    if (base.kind != SurfaceKind::Name || bound.count(base.text) != 0 ||
        !importsAs(base.module, base.text)) {
      return {};
    }
    return qualifiedName(base.module, base.text + "." + field.text);
  }

  [[nodiscard]] Resolution resolve(const std::string& name) const {
    if (declared_.count(name) != 0 || inputs_.count(name) != 0) {
      return Resolution::Input;
    }
    if (states_.count(name) != 0) {
      return Resolution::State;
    }
    if (defs_.count(name) != 0) {
      return Resolution::Definition;
    }
    if (Constants.count(name) != 0) {
      return Resolution::Constant;
    }
    return Resolution::Other;
  }

  /// What the name `name` of `node` resolves to: in a library, a constant before a definition of
  /// the module compiled, as lower() reads it there.
  [[nodiscard]] Resolution resolveAt(const SurfaceNode& node, const std::string& name) const {
    if (node.module != 0 && Constants.count(node.text) != 0) {
      return Resolution::Constant;
    }
    return resolve(name);
  }

  /// Whether a call of `name` reaches a definition of that name: one with parameters, or one
  /// without where the language has no function of the name. A definition without parameters is
  /// read where its name is, and a call of the name is the function, as in the core format, where a
  /// program may have a definition `max` and call `max(a, b)`.
  [[nodiscard]] bool callsDefinition(const std::string& name, const std::string& text) const {
    const auto definition = defs_.find(name);
    return definition != defs_.end() &&
           (!surface_[definition->second].params.empty() || !isCalledFunction(text));
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
      case SurfaceKind::Name: {
        const std::string name = qualifiedName(node.module, node.text);
        if (scope->count(node.text) == 0 && resolveAt(node, name) == Resolution::Definition) {
          found.push_back(name);
        }
        break;
      }
      case SurfaceKind::Call: {
        for (const SurfaceId arg : node.args) {
          stack.emplace_back(arg, scope);
        }
        const std::string name = qualifiedName(node.module, node.text);
        if (scope->count(node.text) == 0 && components_.count(node.text) == 0 &&
            callsDefinition(name, node.text)) {
          found.push_back(name);
        }
        break;
      }
      case SurfaceKind::Record:
      case SurfaceKind::Array:
        for (const SurfaceId arg : node.args) {
          stack.emplace_back(arg, scope);
        }
        break;
      case SurfaceKind::Field:
        // `lib.v`, a value of a module imported as lib
        if (const auto value = importedValue(node, *scope);
            !value.empty() && defs_.count(value) != 0) {
          found.push_back(value);
          break;
        }
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Ascription:
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Index:
        stack.emplace_back(node.a, scope);
        stack.emplace_back(node.b, scope);
        break;
      case SurfaceKind::Comprehension: {
        stack.emplace_back(node.b, scope);
        auto inner = std::make_shared<std::set<std::string>>(*scope);
        inner->insert(node.text);
        stack.emplace_back(node.a, std::move(inner));
        break;
      }
      case SurfaceKind::Conditional:
        for (const SurfaceId part : node.args) {
          stack.emplace_back(part, scope);
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
    const auto bound = [](const Scope& scope) {
      return [scope](const std::string& name) { return scope->count(name) != 0; };
    };
    // the inputs of the fields of a type of `material`
    const auto addFields = [&](TypeId id) {
      const Type type = typeOf(id);
      for (const auto& name : type.names) {
        add(name);
      }
    };
    // the types `material` is given as by the arguments of a call, met after its arguments
    std::vector<TypeId> pending;
    // the body of a definition: of a type and `material`, the inputs of the fields of the type
    const auto walkBody = [&](const SurfaceNode& definition, const Scope& inner) {
      if (definition.type != NoType && isMaterial(surface_[definition.a], bound(inner))) {
        addFields(definition.type);
        return;
      }
      stack.emplace_back(definition.a, inner);
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
        if (scope->count(node.text) != 0 || isMaterial(node, bound(scope))) {
          break;
        }
        const std::string name = qualifiedName(node.module, node.text);
        const auto resolution = resolveAt(node, name);
        if ((resolution == Resolution::Input || resolution == Resolution::Other) &&
            node.module == 0) {
          add(node.text);
        } else if (resolution == Resolution::Definition && expanded.insert(name).second) {
          walkBody(surface_[defs_.at(name)], empty);
        }
        break;
      }
      case SurfaceKind::Call: {
        const std::string name = qualifiedName(node.module, node.text);
        const auto definition = defs_.find(name);
        const bool callsDefinition = scope->count(node.text) == 0 &&
                                     components_.count(node.text) == 0 &&
                                     this->callsDefinition(name, node.text);
        // the body of a definition after its arguments, which lower() binds first; the body of one
        // with parameters once per call of other arguments
        if (callsDefinition && calls.emplace(name, node.args).second) {
          const SurfaceNode& callee = surface_[definition->second];
          walkBody(callee,
                   std::make_shared<const std::set<std::string>>(callee.params.begin(),
                                                                 callee.params.end()));
        }
        for (std::size_t i = node.args.size(); i-- > 0;) {
          // `material` as the argument of a parameter of a type: the inputs of its fields
          if (callsDefinition && i < surface_[definition->second].paramTypes.size() &&
              surface_[definition->second].paramTypes[i] != NoType &&
              isMaterial(surface_[node.args[i]], bound(scope))) {
            pending.emplace_back(surface_[definition->second].paramTypes[i]);
            continue;
          }
          stack.emplace_back(node.args[i], scope);
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
        // `material.rho`: the input `rho`
        if (isMaterial(surface_[node.a], bound(scope))) {
          add(node.text);
          break;
        }
        // `d.c` of data on the mesh: an input channel after all others (compileSderivModule())
        if (const SurfaceNode& base = surface_[node.a]; base.kind == SurfaceKind::Name &&
                                                        scope->count(base.text) == 0 &&
                                                        gridNamed(base.text) != nullptr) {
          break;
        }
        // `lib.v`, a value of a module imported as lib
        if (const auto value = importedValue(node, *scope);
            !value.empty() && defs_.count(value) != 0) {
          if (expanded.insert(value).second) {
            walkBody(surface_[defs_.at(value)], empty);
          }
          break;
        }
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Ascription:
        if (isMaterial(surface_[node.a], bound(scope))) {
          addFields(node.type);
          break;
        }
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Index:
        stack.emplace_back(node.b, scope);
        stack.emplace_back(node.a, scope);
        break;
      case SurfaceKind::Comprehension: {
        auto inner = std::make_shared<std::set<std::string>>(*scope);
        inner->insert(node.text);
        stack.emplace_back(node.a, std::move(inner));
        stack.emplace_back(node.b, scope);
        break;
      }
      case SurfaceKind::Conditional:
        // the condition first, so it is pushed last
        for (auto part = node.args.rbegin(); part != node.args.rend(); ++part) {
          stack.emplace_back(*part, scope);
        }
        break;
      case SurfaceKind::Num:
      case SurfaceKind::Def:
        break;
      }
      // the inputs of the fields of a type `material` was given as, in their order
      for (const TypeId type : pending) {
        addFields(type);
      }
      pending.clear();
    }
    return channels;
  }

  [[noreturn]] void fail(const SurfaceNode& node, const std::string& message) const {
    throw SderivError("resolve", message, node.position);
  }

  /// The definition of `node` while it is being lowered: one that is reached again reads itself.
  class Expansion {
public:
    Expansion(Lowering& lowering, const SurfaceNode& node, const std::string& name)
        : lowering_(lowering) {
      auto& expanding = lowering_.expanding_;
      if (std::find(expanding.begin(), expanding.end(), name) != expanding.end()) {
        lowering_.fail(node, "`" + node.text + "` is defined in terms of itself");
      }
      expanding.push_back(name);
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

  /// Fails for the input `node` reads where the module declares its inputs and not this one.
  void undeclared(const SurfaceNode& node) const {
    if (!declared_.empty()) {
      fail(node,
           "`" + node.text + "` is not declared; the module declares its inputs, as `in " +
               node.text + "`");
    }
  }

  [[nodiscard]] static bool isComparisonOrLogic(const std::string& op) {
    return op == "<" || op == "<=" || op == ">" || op == ">=" || op == "==" || op == "!=" ||
           op == "and" || op == "or";
  }

  /// Whether `node` is a truth value, 1 or 0: a comparison, a logical operation, or a choice
  /// between such.
  [[nodiscard]] bool isTruthValue(NodeId node) const {
    const Node& n = program_.arena()[node];
    if (n.kind == Kind::Const) {
      return n.value == 0.0 || n.value == 1.0;
    }
    if (n.kind != Kind::PW) {
      return false;
    }
    switch (n.fn) {
    case Fn::Lt:
    case Fn::Le:
    case Fn::Eq:
    case Fn::And:
    case Fn::Or:
      return true;
    case Fn::Select:
      return isTruthValue(n.b) && isTruthValue(n.c);
    default:
      return false;
    }
  }

  /// A comparison of `a` and `b`, as lt, le and eq: `a > b` is `b < a`, `a != b` is not `a == b`.
  /// `a and b` and `a or b` are what land and lor give, where `b` decides only: `b` as a truth
  /// value where `a` is true, else 0; 1 where `a` is true, else `b` as a truth value -- a select,
  /// so that a value `b` alone reads is used only where it decides.
  Value comparisonOrLogic(const SurfaceNode& at, Value a, Value b) {
    auto& arena = program_.arena();
    const std::string what = "`" + at.text + "`";
    const auto& op = at.text;
    if (op == "<" || op == "<=" || op == ">" || op == ">=") {
      const Fn fn = op == "<" || op == ">" ? Fn::Lt : Fn::Le;
      const bool swapped = op == ">" || op == ">=";
      return pointwise(at, what, fn, {swapped ? b : a, swapped ? a : b});
    }
    if (op == "==" || op == "!=") {
      return leafwise(at, what, {a, b}, [&](const std::vector<NodeId>& n) {
        const NodeId equal = arena.pw(Fn::Eq, n[0], n[1]);
        return op == "==" ? equal : arena.pw(Fn::Eq, equal, arena.konst(0.0));
      });
    }
    // b as a truth value: 1 where it is true (a NaN too), else 0
    const auto truth = [&](NodeId node) {
      return isTruthValue(node) ? node : arena.pw(Fn::Or, node, arena.konst(0.0));
    };
    return leafwise(at, what, {a, b}, [&](const std::vector<NodeId>& n) {
      return op == "and" ? arena.pw(Fn::Select, n[0], truth(n[1]), arena.konst(0.0))
                         : arena.pw(Fn::Select, n[0], arena.konst(1.0), truth(n[1]));
    });
  }

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
    case Fn::Div:
      return *a / *b;
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

  /// matmul(a, b): of matrices, or of a matrix and a vector, a vector on the left a row and on
  /// the right a column; each element a dot product.
  Value matmulOf(const SurfaceNode& at, const std::string& what, const Value& a, const Value& b) {
    // a matrix or a vector on either side: a vector on the left is a row, on the right a column
    const auto isMatrix = [](const Value& value) {
      return value.shape == Value::Shape::Array && !value.items.empty() &&
             value.items.front().shape == Value::Shape::Array;
    };
    const bool leftMatrix = isMatrix(a);
    const bool rightMatrix = isMatrix(b);
    if (!leftMatrix && !rightMatrix) {
      fail(at, what + " of two vectors: their product is dot(a, b)");
    }
    const auto left = leftMatrix ? matrixOf(at, what, a)
                                 : std::vector<std::vector<NodeId>>{vectorOf(at, what, a)};
    std::vector<std::vector<NodeId>> right;
    if (rightMatrix) {
      right = matrixOf(at, what, b);
    } else {
      for (const NodeId element : vectorOf(at, what, b)) {
        right.push_back({element});
      }
    }
    if (left.front().size() != right.size()) {
      fail(at,
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

  /// percentile(p, values...) and median(values...), the percentile 50: of numbers, or of the
  /// elements of one array of numbers (Kind::Percentile).
  Value percentileValue(const SurfaceNode& node,
                        const std::string& what,
                        const std::vector<Value>& args) {
    auto& arena = program_.arena();
    const bool median = node.text == "median";
    const std::size_t first = median ? 0 : 1;
    const std::string takes =
        what + " takes " + (median ? "" : "the percentile and ") + "numbers, or one array of them";
    if (args.size() <= first) {
      fail(node, takes);
    }
    if (!median && !args[0].isNumber()) {
      fail(node, "the percentile of " + what + " is a number, not " + described(args[0]));
    }
    const NodeId p = median ? arena.konst(50.0) : args[0].node;
    const bool array = args.size() == first + 1 && args[first].shape == Value::Shape::Array;
    const std::vector<Value> values(array ? args[first].items.begin() : args.begin() + first,
                                    array ? args[first].items.end() : args.end());
    std::vector<NodeId> nodes;
    for (const auto& value : values) {
      if (!value.isNumber()) {
        fail(node, takes + ", not " + (array ? described(args[first]) : described(value)));
      }
      nodes.push_back(value.node);
    }
    if (nodes.empty() || nodes.size() > static_cast<std::size_t>(MaxPercentileValues)) {
      fail(node,
           what + " takes 1 to " + std::to_string(MaxPercentileValues) + " values, not " +
               std::to_string(nodes.size()));
    }
    return Value::number(arena.percentile(p, nodes));
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
      return matmulOf(node, what, args[0], args[1]);
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

  /// The value of the definition `name` without parameters, which `node` reads: the same wherever
  /// it is read, since it reads no environment.
  Value definitionNamed(const SurfaceNode& node, const std::string& name) {
    const SurfaceNode& definition = surface_[defs_.at(name)];
    if (!definition.params.empty()) {
      fail(node,
           "`" + node.text + "` expects " + std::to_string(definition.params.size()) +
               " arguments");
    }
    if (const auto known = lowered_.find(name); known != lowered_.end()) {
      return known->second;
    }
    const Expansion expansion(*this, node, name);
    Environment inner;
    Value value = definitionValue(definition, inner);
    lowered_.emplace(name, value);
    return value;
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

  [[noreturn]] void failAt(int position, const std::string& message) const {
    throw SderivError("resolve", message, position);
  }

  /// The type `id`, the names of types replaced by what they stand for.
  [[nodiscard]] Type typeOf(TypeId id) const {
    std::vector<std::string> resolving;
    return resolved(id, resolving);
  }

  [[nodiscard]] Type resolved(TypeId id, std::vector<std::string>& resolving) const {
    const TypeExpression& expression = surface_.type(id);
    Type type;
    switch (expression.kind) {
    case TypeExpression::Kind::Number:
      return type;
    case TypeExpression::Kind::Record: {
      type.kind = Type::Kind::Record;
      type.names = expression.fields;
      type.name = "{";
      for (std::size_t i = 0; i < expression.fields.size(); ++i) {
        type.items.push_back(resolved(expression.items[i], resolving));
        type.name += (i == 0 ? "" : ", ") + expression.fields[i];
      }
      type.name += "}";
      return type;
    }
    case TypeExpression::Kind::Array:
      type.kind = Type::Kind::Array;
      type.length = expression.length;
      type.items.push_back(resolved(expression.items.front(), resolving));
      type.name = type.items.front().name + "[" + std::to_string(type.length) + "]";
      return type;
    case TypeExpression::Kind::Name:
      break;
    }
    // a type a module declares, else one of SeisSol
    const std::string qualified = qualifiedName(expression.module, expression.name);
    for (const auto& declaration : parsed_.types) {
      if (declaration.name != qualified) {
        continue;
      }
      if (std::find(resolving.begin(), resolving.end(), declaration.name) != resolving.end()) {
        failAt(declaration.position,
               "the type `" + declaration.name + "` is defined in terms of itself");
      }
      resolving.push_back(declaration.name);
      type = resolved(declaration.type, resolving);
      resolving.pop_back();
      type.name = declaration.name;
      return type;
    }
    const auto builtin = builtinTypes_.find(expression.name);
    if (builtin == builtinTypes_.end()) {
      std::string known;
      for (const auto& [name, fields] : builtinTypes_) {
        known += (known.empty() ? "" : ", ") + name;
      }
      failAt(expression.position,
             "unknown type `" + expression.name +
                 "`; a type is number, one a module declares with `type`, or one of SeisSol" +
                 (known.empty() ? std::string() : " (" + known + ")"));
    }
    type.kind = Type::Kind::Record;
    type.names = builtin->second;
    type.items.assign(type.names.size(), Type{});
    type.name = expression.name;
    return type;
  }

  /// `value`, of `at`, as a value of `type`: checked, and of the fields of the type, in its order;
  /// a record may have more, which the value leaves out. `what` it is, for a message.
  Value conformed(const SurfaceNode& at,
                  const Value& value,
                  const Type& type,
                  const std::string& what) const {
    const auto refuse = [&](const std::string& why) {
      fail(at, what + " is no `" + type.name + "`: " + why);
    };
    switch (type.kind) {
    case Type::Kind::Number:
      if (!value.isNumber()) {
        refuse("it is " + described(value));
      }
      return value;
    case Type::Kind::Array: {
      if (value.shape != Value::Shape::Array || value.items.size() != type.length) {
        refuse("it is " + described(value));
      }
      std::vector<Value> items;
      for (std::size_t i = 0; i < value.items.size(); ++i) {
        items.push_back(conformed(
            at, value.items[i], type.items.front(), what + ", element " + std::to_string(i) + ","));
      }
      return Value::array(std::move(items));
    }
    case Type::Kind::Record:
      break;
    }
    if (value.shape != Value::Shape::Record) {
      refuse("it is " + described(value));
    }
    Value record;
    record.shape = Value::Shape::Record;
    record.names = type.names;
    for (std::size_t i = 0; i < type.names.size(); ++i) {
      const Value* field = value.field(type.names[i]);
      if (field == nullptr) {
        refuse("it has no field `" + type.names[i] + "`");
      }
      record.items.push_back(
          conformed(at, *field, type.items[i], what + ", field " + type.names[i] + ","));
    }
    return record;
  }

  /// Whether `node` is `material`, the material of the point, where `bound` tells the names a let
  /// or a parameter binds.
  template <typename Bound>
  [[nodiscard]] bool isMaterial(const SurfaceNode& node, const Bound& bound) const {
    return node.kind == SurfaceKind::Name && node.module == 0 && node.text == MaterialName &&
           !bound(node.text) && resolve(node.text) == Resolution::Other &&
           gridNamed(node.text) == nullptr;
  }

  /// `material` as a value of `type`, a record of numbers: the inputs of the names of its fields.
  Value materialAs(const SurfaceNode& at, const Type& type) {
    if (type.kind != Type::Kind::Record) {
      fail(at,
           "`material` is a record of the parameters of the material, not a `" + type.name + "`");
    }
    Value record;
    record.shape = Value::Shape::Record;
    record.names = type.names;
    for (std::size_t i = 0; i < type.names.size(); ++i) {
      if (type.items[i].kind != Type::Kind::Number) {
        fail(at,
             "a parameter of the material is a number, and the field `" + type.names[i] + "` of `" +
                 type.name + "` is not");
      }
      record.items.push_back(Value::number(program_.arena().field(type.names[i])));
    }
    return record;
  }

  /// The value of `id` as a value of `type` (conformed()); of `material`, the record of its inputs.
  Value ascribed(const SurfaceNode& at,
                 SurfaceId id,
                 const Type& type,
                 Environment& env,
                 const std::string& what) {
    if (isMaterial(surface_[id], [&](const std::string& name) { return env.count(name) != 0; })) {
      return materialAs(at, type);
    }
    return conformed(at, lower(id, env), type, what);
  }

  /// The value of the definition `definition`, its parameters bound in `inner`, as a value of its
  /// type where it gives one.
  Value definitionValue(const SurfaceNode& definition, Environment& inner) {
    if (definition.type == NoType) {
      return lower(definition.a, inner);
    }
    return ascribed(definition,
                    definition.a,
                    typeOf(definition.type),
                    inner,
                    "the value of " + definition.text);
  }

  /// The numbers of `value`, which has the shape of `shape`, in the order of leafNames(); `what` it
  /// is, for a message.
  void collectLeaves(const SurfaceNode& at,
                     const Value& value,
                     const InitialValue& shape,
                     const std::string& what,
                     std::vector<NodeId>& leaves) const {
    const auto misfit = [&]() {
      fail(at, what + " has the shape of its initial value, and " + described(value) + " has not");
    };
    switch (shape.shape) {
    case InitialValue::Shape::Number:
      if (!value.isNumber()) {
        misfit();
      }
      leaves.push_back(value.node);
      return;
    case InitialValue::Shape::Record:
      if (value.shape != Value::Shape::Record || value.items.size() != shape.items.size()) {
        misfit();
      }
      for (std::size_t i = 0; i < shape.items.size(); ++i) {
        const Value* field = value.field(shape.names[i]);
        if (field == nullptr) {
          misfit();
        }
        collectLeaves(at, *field, shape.items[i], what, leaves);
      }
      return;
    case InitialValue::Shape::Array:
      if (value.shape != Value::Shape::Array || value.items.size() != shape.items.size()) {
        misfit();
      }
      for (std::size_t i = 0; i < shape.items.size(); ++i) {
        collectLeaves(at, value.items[i], shape.items[i], what, leaves);
      }
      return;
    }
  }

  /// `value`, of `at`, as an initial value: a number the program fixes, or a record or an array of
  /// such. `what` it is, for a message.
  InitialValue initialOf(const SurfaceNode& at, const Value& value, const std::string& what) const {
    InitialValue initial;
    switch (value.shape) {
    case Value::Shape::Number: {
      const auto fixed = fixedValueOf(value.node);
      if (!fixed.has_value()) {
        fail(at, what + " is a number the program fixes, not one of the point");
      }
      initial.number = *fixed;
      return initial;
    }
    case Value::Shape::Record:
      initial.shape = InitialValue::Shape::Record;
      initial.names = value.names;
      break;
    case Value::Shape::Array:
      initial.shape = InitialValue::Shape::Array;
      break;
    }
    for (const auto& item : value.items) {
      initial.items.push_back(initialOf(at, item, what));
    }
    return initial;
  }

  /// The value of the state `name` of the shape `shape`: its value of the previous call, a state of
  /// the program per number.
  Value stateValue(const InitialValue& shape, const std::string& name) {
    switch (shape.shape) {
    case InitialValue::Shape::Number:
      return Value::number(program_.arena().field(name));
    case InitialValue::Shape::Record: {
      Value record;
      record.shape = Value::Shape::Record;
      record.names = shape.names;
      for (std::size_t i = 0; i < shape.items.size(); ++i) {
        record.items.push_back(stateValue(shape.items[i], name + "." + shape.names[i]));
      }
      return record;
    }
    case InitialValue::Shape::Array:
      break;
    }
    std::vector<Value> items;
    for (std::size_t i = 0; i < shape.items.size(); ++i) {
      items.push_back(stateValue(shape.items[i], name + "[" + std::to_string(i) + "]"));
    }
    return Value::array(std::move(items));
  }

  Value lower(SurfaceId id, Environment& env) {
    const SurfaceNode& node = surface_[id];
    switch (node.kind) {
    case SurfaceKind::Num:
      return Value::number(program_.arena().konst(node.number));

    case SurfaceKind::Unary:
      // `not a` is eq(a, 0), as lnot
      if (node.text == "not") {
        return pointwise(node,
                         "`not`",
                         Fn::Eq,
                         {lower(node.a, env), Value::number(program_.arena().konst(0.0))});
      }
      return pointwise(node, "`-`", Fn::Neg, {lower(node.a, env)});

    case SurfaceKind::Conditional: {
      std::vector<Value> parts;
      parts.reserve(node.args.size());
      for (const SurfaceId part : node.args) {
        parts.push_back(lower(part, env));
      }
      return pointwise(node, "`if`", Fn::Select, std::move(parts));
    }

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
        if (binary.text == "@") {
          value = matmulOf(binary, "`@`", value, b);
          continue;
        }
        if (isComparisonOrLogic(binary.text)) {
          value = comparisonOrLogic(binary, std::move(value), std::move(b));
          continue;
        }
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

    case SurfaceKind::Ascription:
      return ascribed(node, node.a, typeOf(node.type), env, "the value");

    case SurfaceKind::Comprehension: {
      // [e for v in a]: e with v bound to every element of a
      const Value array = lower(node.b, env);
      if (array.shape != Value::Shape::Array) {
        fail(node, "a comprehension runs over an array, not " + described(array));
      }
      std::vector<Value> items;
      items.reserve(array.items.size());
      for (const auto& element : array.items) {
        Environment inner = env;
        inner[node.text] = element;
        items.push_back(lower(node.a, inner));
      }
      return Value::array(std::move(items));
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
      // `m.rho` of a dataset m, without coordinates
      if (const SurfaceNode& base = surface_[node.a]; base.kind == SurfaceKind::Name &&
                                                      env.count(base.text) == 0 &&
                                                      gridNamed(base.text) != nullptr) {
        const GridDeclaration& grid = *gridNamed(base.text);
        const std::string sample = base.text + "." + node.text;
        if (std::find(grid.components.begin(), grid.components.end(), node.text) ==
            grid.components.end()) {
          std::string list;
          for (const auto& each : grid.components) {
            list += (list.empty() ? "" : ", ") + each;
          }
          fail(node,
               "the dataset `" + base.text + "` has no component `" + node.text + "`; it has " +
                   list);
        }
        if (grid.kind == "cells") {
          // data on the mesh: the value of the cell of the point, an input channel the consumer
          // binds per cell (and the program has after all others, compileSderivModule())
          return Value::number(program_.arena().field(sample));
        }
        fail(node,
             "`" + sample + "` is a component of a dataset, sampled with coordinates: " + sample +
                 "(x, y, z)");
      }
      // `material.rho`: the input `rho`, the parameter of the material of the point
      if (isMaterial(surface_[node.a],
                     [&](const std::string& name) { return env.count(name) != 0; })) {
        return Value::number(program_.arena().field(node.text));
      }
      // `lib.v`: the value v of a module imported as lib
      {
        std::set<std::string> bound;
        for (const auto& [name, value] : env) {
          bound.insert(name);
        }
        if (const auto value = importedValue(node, bound); !value.empty()) {
          if (defs_.count(value) == 0) {
            fail(node,
                 "the module `" + surface_[node.a].text + "` has no definition `" + node.text +
                     "`");
          }
          return definitionNamed(node, value);
        }
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
      // a library reads its parameters, the definitions of its module and constants only
      if (node.module != 0) {
        const auto constant = Constants.find(node.text);
        if (constant != Constants.end()) {
          return Value::number(program_.arena().konst(constant->second));
        }
        const std::string name = qualifiedName(node.module, node.text);
        if (defs_.count(name) != 0) {
          return definitionNamed(node, name);
        }
        fail(node,
             "a definition of a library reads its parameters, the definitions of its module and "
             "constants, and `" +
                 node.text + "` is none of them");
      }
      // An input the module declares is the input, whatever else its name could be.
      if (declared_.count(node.text) != 0) {
        return Value::number(channel(node));
      }
      // An input the consumer names reads the input, as a state reads its previous value; a module
      // that declares its inputs declares it as well.
      if (inputs_.count(node.text) != 0) {
        undeclared(node);
        return Value::number(channel(node));
      }
      // A state reads its value from the previous call -- also where the definition of the same
      // name, its next value, is in scope. That is what makes the update a parallel assignment.
      if (states_.count(node.text) != 0) {
        const auto shape = stateShapes_.find(node.text);
        if (shape == stateShapes_.end()) {
          fail(node, "the initial value of a state reads no state");
        }
        return stateValue(shape->second, node.text);
      }
      // a definition of the module before a constant of the language: `def g = ...` hides g
      if (defs_.count(node.text) != 0) {
        return definitionNamed(node, node.text);
      }
      const auto constant = Constants.find(node.text);
      if (constant != Constants.end()) {
        return Value::number(program_.arena().konst(constant->second));
      }
      if (importsAs(node.module, node.text)) {
        fail(node,
             "`" + node.text + "` is a module; its definitions are read as " + node.text + ".name");
      }
      if (const GridDeclaration* grid = gridNamed(node.text); grid != nullptr) {
        if (grid->kind == "cells") {
          fail(node,
               "`" + node.text +
                   "` is a dataset of a value per cell of the mesh; a component of "
                   "it is read as " +
                   node.text + ".c");
        }
        fail(node,
             "`" + node.text + "` is a dataset; a component of it is sampled as " + node.text +
                 ".c(x, y, z), all of them as " + node.text + "(x, y, z)");
      }
      if (node.text == MaterialName) {
        fail(node,
             "`material` is the material of the point: a parameter of it is read as "
             "material.rho, those of a type as material : elastic");
      }
      // (the name of an intrinsic of arrays is a name like any other where it is not called)
      if (isFunction(node.text)) {
        fail(node, "`" + node.text + "` is a function and must be called");
      }
      // Everything left is an input channel. No declaration, no catalogue -- unless the module
      // declares its inputs.
      undeclared(node);
      return Value::number(channel(node));
    }

    case SurfaceKind::Call: {
      if (env.count(node.text) != 0) {
        fail(node, "`" + node.text + "` is not callable");
      }

      // the grids are those of the module compiled
      const auto component = node.module == 0 ? components_.find(node.text) : components_.end();
      if (component != components_.end()) {
        if (parsed_.grids[component->second.grid].kind == "cells") {
          fail(node,
               "the dataset `" + parsed_.grids[component->second.grid].name +
                   "` holds a value per cell of the mesh, read without coordinates, as " +
                   node.text);
        }
        const auto coords = coordinatesOf(node, component->second.grid, env);
        return Value::number(program_.arena().lookup(
            gridIds_[component->second.grid], component->second.component, coords));
      }
      const std::string name = qualifiedName(node.module, node.text);
      if (const auto dot = node.text.find('.');
          dot != std::string::npos && defs_.count(name) == 0) {
        const std::string gridName = node.text.substr(0, dot);
        if (importsAs(node.module, gridName)) {
          fail(node,
               "the module `" + gridName + "` has no definition `" + node.text.substr(dot + 1) +
                   "`");
        }
        const GridDeclaration* grid = gridNamed(gridName);
        if (grid == nullptr) {
          fail(node, "`" + gridName + "` is no dataset, so `" + node.text + "` samples nothing");
        }
        std::string list;
        for (const auto& each : grid->components) {
          list += (list.empty() ? "" : ", ") + each;
        }
        fail(node,
             "the dataset `" + gridName + "` has no component `" + node.text.substr(dot + 1) +
                 "`; it has " + list);
      }

      const auto def = defs_.find(name);
      if (def != defs_.end() && callsDefinition(name, node.text)) {
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
          // an argument as a value of the type of its parameter, where it gives one
          Value argument =
              definition.paramTypes[i] == NoType
                  ? lower(node.args[i], env)
                  : ascribed(node,
                             node.args[i],
                             typeOf(definition.paramTypes[i]),
                             env,
                             "the argument " + definition.params[i] + " of " + node.text);
          keyOf(argument, leaves, form);
          form += ';';
          inner[definition.params[i]] = std::move(argument);
        }
        // the value of a call depends on its arguments alone
        auto key = std::make_tuple(name, std::move(form), std::move(leaves));
        if (const auto known = calls_.find(key); known != calls_.end()) {
          return known->second;
        }
        const Expansion expansion(*this, node, name);
        Value value = definitionValue(definition, inner);
        calls_.emplace(std::move(key), value);
        return value;
      }

      // m(x, y, z): every component of the grid m, a record of their names
      if (const GridDeclaration* grid = node.module == 0 ? gridNamed(node.text) : nullptr;
          grid != nullptr) {
        if (grid->kind == "cells") {
          fail(node,
               "the dataset `" + grid->name +
                   "` holds a value per cell of the mesh, read "
                   "without coordinates, a component as " +
                   grid->name + ".c");
        }
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

      // the percentile p, or the median, of numbers or of the elements of one array
      if (node.text == "percentile" || node.text == "median") {
        return percentileValue(node, what, args);
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
  const std::vector<ModuleInfo>& modules_;
  Program& program_;
  std::map<std::string, SurfaceId> defs_;
  std::set<std::string> states_;
  /// the inputs the module declares (`in x`); where it declares any, it has no others
  std::set<std::string> declared_;
  /// the initial value of every state, which gives its shape (initializeStates())
  std::map<std::string, InitialValue> stateShapes_;
  /// the types SeisSol offers by name (recordTypes())
  std::map<std::string, std::vector<std::string>> builtinTypes_;
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
// The text of a module of a compilation: its file, and the offset of its first character among
// the positions of all modules -- the module compiled first, at 0, and then those it imports.
struct ModuleSource {
  std::string file;
  std::string text;
  int base{0};
};

/// The module of the position `position`.
const ModuleSource& moduleAt(const std::vector<ModuleSource>& sources, int position) {
  for (auto source = sources.rbegin(); source != sources.rend(); ++source) {
    if (position >= source->base) {
      return *source;
    }
  }
  return sources.front();
}

/// The tokens of the module `text`, their positions past `base`.
std::vector<Token> tokensAt(const std::string& text, int base) {
  std::vector<Token> tokens;
  try {
    tokens = tokenize(text);
  } catch (const SderivError& error) {
    throw SderivError("lex", std::string(error.what()).substr(5), error.position() + base);
  }
  for (auto& token : tokens) {
    token.position += base;
  }
  return tokens;
}

/// Loads the modules `parsed` imports, and those they import, into `surface`: every module once,
/// its definitions and types added to `parsed` under their qualified names, and its imports
/// recorded in `modules`. A module of the library is `std:name`; a file is relative to the
/// directory of the file importing it.
void loadImports(ParsedProgram& parsed,
                 SurfaceArena& surface,
                 std::vector<ModuleInfo>& modules,
                 std::vector<ModuleSource>& sources) {
  std::map<std::string, int> loaded;
  std::vector<std::string> loading;
  // the module `importer` imports by `declaration`, from the directory `directory`
  const std::function<int(int, const ImportDeclaration&, const std::string&)> load =
      [&](int importer, const ImportDeclaration& declaration, const std::string& directory) {
        if (modules[static_cast<std::size_t>(importer)].imports.count(declaration.name) != 0) {
          throw SderivError("import",
                            "the module is imported as `" + declaration.name + "` twice",
                            declaration.position);
        }
        std::string key;
        std::string file;
        std::string text;
        if (declaration.file) {
          std::filesystem::path path(declaration.path);
          if (path.is_relative() && !directory.empty()) {
            path = std::filesystem::path(directory) / path;
          }
          file = path.lexically_normal().string();
          key = file;
          std::ifstream stream(file, std::ios::binary);
          if (!stream) {
            throw SderivError(
                "import", "cannot read the module `" + file + "`", declaration.position);
          }
          text.assign(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
        } else {
          const auto& library = libraryModules();
          const auto found = library.find(declaration.path);
          if (found == library.end()) {
            std::string names;
            for (const auto& [name, module] : library) {
              names += (names.empty() ? "" : ", ") + name;
            }
            throw SderivError("import",
                              "the library has no module `" + declaration.path + "`; it has " +
                                  names,
                              declaration.position);
          }
          key = "std:" + declaration.path;
          file = key;
          text = found->second;
        }
        // a module that imports one that imports it, directly or not
        if (std::find(loading.begin(), loading.end(), key) != loading.end()) {
          throw SderivError(
              "import", "the module `" + file + "` imports itself", declaration.position);
        }
        if (const auto known = loaded.find(key); known != loaded.end()) {
          return known->second;
        }
        // the qualifier of its names: the name it is imported as, unique among the modules
        std::string qualifier = declaration.name;
        for (int n = 2;
             std::any_of(modules.begin(),
                         modules.end(),
                         [&](const auto& module) { return module.qualifier == qualifier; });
             ++n) {
          qualifier = declaration.name + "~" + std::to_string(n);
        }
        const int index = static_cast<int>(modules.size());
        modules.push_back(ModuleInfo{qualifier, {}});
        loaded.emplace(key, index);
        const int base = sources.back().base + static_cast<int>(sources.back().text.size()) + 1;
        sources.push_back(ModuleSource{file, text, base});
        loading.push_back(key);
        const auto tokens = tokensAt(text, base);
        surface.setModule(index);
        ParsedProgram module = Parser(tokens, surface, qualifier).program();
        surface.setModule(0);
        const std::string own =
            declaration.file ? std::filesystem::path(file).parent_path().string() : std::string();
        for (const auto& nested : module.imports) {
          modules[static_cast<std::size_t>(index)].imports.emplace(nested.name,
                                                                   load(index, nested, own));
        }
        loading.pop_back();
        parsed.defs.insert(parsed.defs.end(), module.defs.begin(), module.defs.end());
        parsed.types.insert(parsed.types.end(), module.types.begin(), module.types.end());
        return index;
      };
  const std::string directory =
      sources.front().file.empty()
          ? std::string()
          : std::filesystem::path(sources.front().file).parent_path().string();
  for (const auto& declaration : parsed.imports) {
    const int index = load(0, declaration, directory);
    modules.front().imports.emplace(declaration.name, index);
  }
}

void compileModules(const std::string& source,
                    const std::string* externalName,
                    reader::scripting::DataType type,
                    const SderivOptions& options,
                    Program& program,
                    std::vector<std::string>& channelOrder,
                    std::vector<ModuleSource>& sources) {
  SurfaceArena surface;
  const auto tokens = tokenize(source);
  Parser parser(tokens, surface);
  ParsedProgram parsed = parser.program();
  // the modules it imports, whose definitions and types join its own
  std::vector<ModuleInfo> modules(1);
  loadImports(parsed, surface, modules, sources);

  if (externalName != nullptr && parsed.expression < 0) {
    throw SderivError("parse",
                      "the caller named the output `" + *externalName +
                          "`, but this module declares its own outputs with `out def`",
                      0);
  }
  if (externalName == nullptr && parsed.exports.empty() && parsed.spreads.empty() &&
      parsed.named.empty()) {
    throw SderivError(
        "parse", "the module declares no outputs; mark at least one definition `out def`", 0);
  }

  // `meta`: the provenance first, as in the core format, so that a stale source is the error;
  // `meta line`, where a construct comes from, once the module is compiled. Of this module only,
  // where the program has those of other sources already.
  {
    Program provenance;
    std::vector<CorePosition> positions;
    for (const auto& declaration : parsed.meta) {
      if (declaration.entry.key != "line") {
        provenance.addMeta(declaration.entry);
        positions.push_back(positionOf(
            source, static_cast<std::size_t>(std::max(declaration.position, 0)), options.file));
      }
    }
    checkProvenance(provenance, positions, source, options.file);
    for (const auto& entry : provenance.meta()) {
      program.addMeta(entry);
    }
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
    // A definition without parameters may take the name of a constant or a function: where the
    // name is read, it is the definition, and where it is called, the function -- as in the core
    // format, whose names of operations are no reserved words. One with parameters may not, as a
    // call of the name would be the definition, and the function out of reach.
    if (!surface[id].params.empty() &&
        (Constants.count(surface[id].text) != 0 || isFunction(surface[id].text))) {
      throw SderivError("parse",
                        "the definition `" + surface[id].text +
                            "` takes parameters and has the name of a constant or function",
                        surface[id].position);
    }
  }
  const ComponentTable components = checkGrids(parsed.grids, defNames);

  // the types a module declares: each once, and none of the name of a type of SeisSol
  {
    const auto builtins = recordTypes();
    std::set<std::string> typeNames;
    for (const auto& declaration : parsed.types) {
      if (declaration.name == "number" || builtins.count(declaration.name) != 0) {
        throw SderivError("type",
                          "the type `" + declaration.name +
                              "` is one of SeisSol; a type of the module has another name",
                          declaration.position);
      }
      if (!typeNames.insert(declaration.name).second) {
        throw SderivError(
            "type", "duplicate type `" + declaration.name + "`", declaration.position);
      }
    }
  }

  // Every state needs exactly one definition of the same name, its next value; and a state name
  // must not be anything else a name can resolve to.
  std::set<std::string> stateNames;
  for (const auto& state : parsed.states) {
    if (!stateNames.insert(state.name).second) {
      throw SderivError("state", "duplicate state `" + state.name + "`", state.position);
    }
    if (components.count(state.name) != 0) {
      throw SderivError("state",
                        "the state `" + state.name + "` has the name of a dataset component",
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
  for (const auto& [name, position] : parsed.nexts) {
    if (stateNames.count(name) == 0) {
      throw SderivError(
          "state", "`next` gives the next value of a state, and `" + name + "` is none", position);
    }
  }

  // the inputs the module declares: each once, and of a name nothing else of it has (one namespace
  // of inputs, states, definitions and grids, as in the core format)
  {
    std::set<std::string> inputNames;
    for (const auto& input : parsed.inputs) {
      if (!inputNames.insert(input.name).second) {
        throw SderivError(
            "parse", "the input `" + input.name + "` is declared twice", input.position);
      }
      const bool taken = defNames.count(input.name) != 0 || stateNames.count(input.name) != 0 ||
                         std::any_of(parsed.grids.begin(),
                                     parsed.grids.end(),
                                     [&](const auto& grid) { return grid.name == input.name; });
      if (taken) {
        throw SderivError("parse",
                          "`" + input.name +
                              "` is declared as an input and as a definition, a state or a dataset",
                          input.position);
      }
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
    if (!reader::datafield::sampledAtCoordinates(desc.kind)) {
      // its input channels are named after it, `d.c`: two names of one file are two datasets
      desc.name = grid.name;
    }
    // internGrid dedupes, so two outputs naming the same grid share a GridId
    // and the backend loads the file once. The names are what the core format prints.
    gridIds.push_back(program.internGrid(desc));
    program.nameGrid(gridIds.back(), GridNames{grid.name, grid.components});
  }

  // the names modules are imported as are no other names of the module
  for (const auto& declaration : parsed.imports) {
    const bool taken =
        defNames.count(declaration.name) != 0 || stateNames.count(declaration.name) != 0 ||
        std::any_of(parsed.inputs.begin(),
                    parsed.inputs.end(),
                    [&](const auto& input) { return input.name == declaration.name; }) ||
        std::any_of(parsed.grids.begin(), parsed.grids.end(), [&](const auto& grid) {
          return grid.name == declaration.name;
        });
    if (taken) {
      throw SderivError("import",
                        "the module is imported as `" + declaration.name +
                            "`, the name of an input, a definition, a state or a dataset",
                        declaration.position);
    }
  }

  Lowering lowering(surface, parsed, components, gridIds, options.inputs, modules, program);
  lowering.resolveTypes();
  lowering.initializeStates();

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
    for (const auto& output : parsed.named) {
      rootExpressions.push_back(output.value);
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
    for (const auto& [kind, index] : parsed.outputOrder) {
      if (kind == OutputKind::Named) {
        // `out v = e`: an output, no definition
        const NamedOutput& output = parsed.named[index];
        roots.emplace_back(output.name,
                           lowering.lower(output.value, "the output `" + output.name + "`"));
        rootPositions.push_back(output.position);
        continue;
      }
      if (kind == OutputKind::Definition) {
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
  // A state of an array or a record is a state of the program per number, `s[0]`, `s.v`, whose
  // next value is the number of the next value of the same place.
  struct StateLeaf {
    std::string name;
    double initial{0.0};
    NodeId next{NoNode};
    ValueType type{ValueType::Real};
    int position{0};
  };
  std::vector<StateLeaf> stateLeaves;
  for (const auto& state : parsed.states) {
    for (const SurfaceId id : parsed.defs) {
      if (surface[id].text != state.name) {
        continue;
      }
      const std::string what = "the next value of the state `" + state.name + "`";
      const InitialValue& initial = lowering.initialOfState(state.name);
      std::vector<std::string> names;
      leafNames(initial, state.name, names);
      std::vector<double> initials;
      leafNumbers(initial, initials);
      std::vector<NodeId> nexts;
      if (initial.shape == InitialValue::Shape::Number) {
        nexts.push_back(lowering.lowerDefinition(id, what));
      } else {
        nexts = lowering.leavesOfShape(id, initial, what);
      }
      for (std::size_t i = 0; i < names.size(); ++i) {
        stateLeaves.push_back(
            StateLeaf{names[i], initials[i], nexts[i], state.type, state.position});
      }
    }
  }
  // and those of the checks
  std::vector<NodeId> checkRoots;
  for (const auto& check : parsed.checks) {
    checkRoots.push_back(lowering.lower(check.first, "the condition of a check"));
  }

  // the inputs the module declares, in their order, those it does not read too -- the program of a
  // core file, which reads the same here -- and then the others in the order of their first read
  for (const auto& input : parsed.inputs) {
    if (std::find(channelOrder.begin(), channelOrder.end(), input.name) == channelOrder.end()) {
      channelOrder.push_back(input.name);
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
  for (const auto& leaf : stateLeaves) {
    program.addState(leaf.name, leaf.initial, leaf.next, leaf.type);
  }
  for (std::size_t i = 0; i < parsed.checks.size(); ++i) {
    program.addCheck(checkRoots[i], parsed.checks[i].second);
  }

  // where they come from (`meta line` of the core format): the definitions, the states, the
  // outputs and the conditions of the checks
  const auto at = [&](int offset) {
    const ModuleSource& module = moduleAt(sources, std::max(offset, 0));
    const CorePosition position =
        positionOf(module.text,
                   static_cast<std::size_t>(std::max(offset, 0) - module.base),
                   module.base == 0 ? options.file : module.file);
    return SourcePosition{position.file, position.line, position.column};
  };
  // a `meta line` gives where a definition, a state or an output comes from in place of its
  // position in this file, which a generated file is not the source of: a definition of the name,
  // else a state, else an output, as the core format resolves it
  std::vector<SourceOf> lines;
  std::set<std::pair<SourceOf::What, std::string>> lined;
  for (const auto& declaration : parsed.meta) {
    const auto& values = declaration.entry.values;
    if (declaration.entry.key != "line") {
      continue;
    }
    if (values.size() != 4 || !values[0].quoted || !values[1].quoted || values[2].quoted ||
        values[3].quoted) {
      throw SderivError(
          "parse",
          "`meta line` names a definition, state or output, then the file, the line "
          "and the column it comes from: `meta line = \"t4\", \"model.sderiv\", 12, 9`",
          declaration.position);
    }
    const std::string& name = values[0].text;
    const SourcePosition position{
        values[1].text, std::stoi(values[2].text), std::stoi(values[3].text)};
    if (defNames.count(name) != 0) {
      // of a definition the program does not read, nothing
      if (const auto node = lowering.definition(name); node.has_value()) {
        lines.push_back(SourceOf{SourceOf::What::Node, *node, std::string(), position});
      }
      lined.emplace(SourceOf::What::Node, name);
    } else if (std::any_of(stateLeaves.begin(), stateLeaves.end(), [&](const auto& leaf) {
                 return leaf.name == name;
               })) {
      lines.push_back(SourceOf{SourceOf::What::State, NoNode, name, position});
      lined.emplace(SourceOf::What::State, name);
    } else if (std::any_of(roots.begin(), roots.end(), [&](const auto& root) {
                 return root.first == name;
               })) {
      lines.push_back(SourceOf{SourceOf::What::Output, NoNode, name, position});
      lined.emplace(SourceOf::What::Output, name);
    } else {
      throw SderivError("parse",
                        "`meta line` names `" + name +
                            "`, which is no definition, state or output of the module",
                        declaration.position);
    }
  }
  for (const SurfaceId id : parsed.defs) {
    if (!surface[id].params.empty() || lined.count({SourceOf::What::Node, surface[id].text}) != 0) {
      continue;
    }
    if (const auto node = lowering.definition(surface[id].text); node.has_value()) {
      program.addSource(
          SourceOf{SourceOf::What::Node, *node, std::string(), at(surface[id].position)});
    }
  }
  if (externalName == nullptr) {
    for (std::size_t i = 0; i < roots.size(); ++i) {
      if (lined.count({SourceOf::What::Output, roots[i].first}) == 0) {
        program.addSource(
            SourceOf{SourceOf::What::Output, NoNode, roots[i].first, at(rootPositions[i])});
      }
    }
  }
  for (const auto& leaf : stateLeaves) {
    if (lined.count({SourceOf::What::State, leaf.name}) == 0) {
      program.addSource(SourceOf{SourceOf::What::State, NoNode, leaf.name, at(leaf.position)});
    }
  }
  for (auto& line : lines) {
    program.addSource(std::move(line));
  }
  for (std::size_t i = 0; i < parsed.checks.size(); ++i) {
    program.addSource(SourceOf{SourceOf::What::Node,
                               checkRoots[i],
                               std::string(),
                               at(surface[parsed.checks[i].first].position)});
  }
}

} // namespace

namespace {

void compileSource(const std::string& source,
                   const std::string* externalName,
                   reader::scripting::DataType type,
                   const SderivOptions& options,
                   Program& program,
                   std::vector<std::string>& channelOrder) {
  std::vector<ModuleSource> sources{ModuleSource{options.file, source, 0}};
  try {
    compileModules(source, externalName, type, options, program, channelOrder, sources);
  } catch (const CoreError&) {
    throw;
  } catch (const SderivError& error) {
    // an error in a module it imports, at its line and column there
    const ModuleSource& module = moduleAt(sources, std::max(error.position(), 0));
    if (error.position() < 0 || module.base == 0) {
      throw;
    }
    throw CoreError(positionOf(module.text,
                               static_cast<std::size_t>(error.position() - module.base),
                               module.file),
                    error.what());
  }
}

/// The input channels of the data on the mesh of `program`, `d.c`, after all its other inputs, in
/// the order of the datasets and their components: the order a core file gives them, where its
/// declarations make them.
void addCellChannels(Program& program) {
  for (const auto& desc : program.grids()) {
    if (reader::datafield::sampledAtCoordinates(desc.kind)) {
      continue;
    }
    for (const auto& component : desc.componentNames) {
      const std::string channel = desc.name + "." + component;
      for (const auto& input : program.inputs()) {
        if (input.name == channel) {
          throw SderivError("dataset",
                            "the input `" + channel +
                                "` has the name of a component of the "
                                "data on the mesh `" +
                                desc.name + "`",
                            -1);
        }
      }
      program.addInput(channel, reader::scripting::DataType::F64);
    }
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
  addCellChannels(program);
  validate(program);
  return program;
}

Program compileSderiv(const std::string& source, const std::string& outputName) {
  return compileSderiv({SderivOutput{outputName, source, reader::scripting::DataType::F64}});
}

Program compileSderivModule(const std::string& source) {
  return compileSderivModule(source, SderivOptions{});
}

namespace {

/// The record types registered, and what guards them.
std::mutex& recordTypesMutex() {
  static std::mutex mutex;
  return mutex;
}
std::map<std::string, std::vector<std::string>>& recordTypeTable() {
  static std::map<std::string, std::vector<std::string>> table;
  return table;
}

} // namespace

void registerRecordType(const std::string& name, const std::vector<std::string>& fields) {
  const std::lock_guard lock(recordTypesMutex());
  recordTypeTable()[name] = fields;
}

std::map<std::string, std::vector<std::string>> recordTypes() {
  const std::lock_guard lock(recordTypesMutex());
  return recordTypeTable();
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
  addCellChannels(program);
  validate(program);
  return program;
}

} // namespace seissol::expr
