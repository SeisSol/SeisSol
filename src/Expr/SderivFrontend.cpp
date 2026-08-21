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
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

// ============================================================ tokenizer =====
enum class TokenKind : std::uint8_t { Num, Str, Pow, Op, Name, Let, In, Def, Grid, Eof };

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
      } else if (value == "grid") {
        kind = TokenKind::Grid;
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
  std::string interpolation;
  std::vector<std::string> components;
  int position{0};
};

struct ParsedProgram {
  std::vector<SurfaceId> defs;
  std::vector<GridDeclaration> grids;
  SurfaceId expression{-1};
};

// =============================================== catalogs (spec-defined) ====
const std::unordered_map<std::string, double> Constants = {{"pi", M_PI}, {"g", 9.80665}};

// Spelled exactly as the Lua-side ssol.* names, so a model transliterates.
const std::unordered_map<std::string, Fn> Builtins = {
    {"sqrt", Fn::Sqrt},    {"abs", Fn::Abs},     {"exp", Fn::Exp},     {"log", Fn::Log},
    {"log2", Fn::Log2},    {"log10", Fn::Log10}, {"sign", Fn::Sign},   {"floor", Fn::Floor},
    {"ceil", Fn::Ceil},    {"round", Fn::Round}, {"sin", Fn::Sin},     {"cos", Fn::Cos},
    {"tan", Fn::Tan},      {"asin", Fn::Asin},   {"acos", Fn::Acos},   {"sinh", Fn::Sinh},
    {"cosh", Fn::Cosh},    {"tanh", Fn::Tanh},   {"asinh", Fn::Asinh}, {"acosh", Fn::Acosh},
    {"atanh", Fn::Atanh},  {"erf", Fn::Erf},     {"min", Fn::Min},     {"max", Fn::Max},
    {"pow", Fn::Pow},      {"atan2", Fn::Atan2}, {"mod", Fn::Mod},     {"lt", Fn::Lt},
    {"le", Fn::Le},        {"eq", Fn::Eq},       {"land", Fn::And},    {"lor", Fn::Or},
    {"select", Fn::Select}};

// Closed vocabularies, so a swapped file/interpolation pair fails at parse time
// rather than at grid-load time on one rank of a large job.
const std::unordered_set<std::string> GridKinds = {"asagi", "scec"};
const std::unordered_set<std::string> Interpolations = {"linear", "nearest"};

// ============================================================== parser ======
class Parser {
  public:
  Parser(const std::vector<Token>& tokens, SurfaceArena& arena) : tokens_(tokens), arena_(arena) {}

  ParsedProgram program() {
    ParsedProgram parsed;
    // def and grid interleave freely and need no separator, for the same reason
    // def already needed none: each starts with its own token kind, so the
    // greedy expression in a def body stops cleanly at the next declaration.
    while (peek().kind == TokenKind::Def || peek().kind == TokenKind::Grid) {
      if (peek().kind == TokenKind::Def) {
        parsed.defs.push_back(definition());
      } else {
        parsed.grids.push_back(gridDeclaration());
      }
    }
    parsed.expression = expression();
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

  GridDeclaration gridDeclaration() {
    GridDeclaration grid;
    grid.position = peek().position;
    eat(TokenKind::Grid);
    grid.name = eat(TokenKind::Name).value;
    eat(TokenKind::Op, "=");
    grid.kind = eat(TokenKind::Str).value;
    eat(TokenKind::Op, ",");
    grid.file = eat(TokenKind::Str).value;
    eat(TokenKind::Op, ",");
    grid.interpolation = eat(TokenKind::Str).value;
    while (valueIs(",")) {
      eat(TokenKind::Op, ",");
      grid.components.push_back(eat(TokenKind::Str).value);
    }
    return grid;
  }

  SurfaceId definition() {
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
    eat(TokenKind::Op, "=");
    node.a = expression();
    return arena_.add(std::move(node));
  }

  SurfaceId expression() {
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
    if (GridKinds.count(grid.kind) == 0) {
      throw SderivError(
          "grid", "grid `" + grid.name + "`: unknown kind \"" + grid.kind + "\"", grid.position);
    }
    if (Interpolations.count(grid.interpolation) == 0) {
      throw SderivError("grid",
                        "grid `" + grid.name + "`: unknown interpolation \"" + grid.interpolation +
                            "\"; did the file and the interpolation get swapped?",
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
      if (Constants.count(generated) != 0 || Builtins.count(generated) != 0 ||
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
           Program& program)
      : surface_(surface), parsed_(parsed), components_(components), gridIds_(gridIds),
        program_(program) {
    for (const SurfaceId id : parsed_.defs) {
      defs_[surface_[id].text] = id;
    }
  }

  NodeId lower(SurfaceId id) {
    Environment empty;
    return lower(id, empty);
  }

  // Channels discovered along the way, in first-use order. Deterministic
  // because the walk is, which matters: the signature order feeds the binding
  // and, through it, the tile layout.
  [[nodiscard]] const std::vector<std::string>& channels() const { return channels_; }

  private:
  using Environment = std::map<std::string, NodeId>;

  [[noreturn]] void fail(const SurfaceNode& node, const std::string& message) const {
    throw SderivError("resolve", message, node.position);
  }

  NodeId channel(const SurfaceNode& node) {
    const std::string& name = node.text;
    if (std::find(channels_.begin(), channels_.end(), name) == channels_.end()) {
      channels_.push_back(name);
    }
    return program_.arena().field(name);
  }

  NodeId lowerLookup(const SurfaceNode& node, const ComponentFunction& fn, Environment& env) {
    if (node.args.empty() || node.args.size() > 6) {
      fail(node,
           "`" + node.text + "`: sampled with " + std::to_string(node.args.size()) +
               " coordinates, the IR allows 1..6");
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
      const NodeId a = lower(node.a, env);
      const NodeId b = lower(node.b, env);
      if (node.text == "+") {
        return program_.arena().pw(Fn::Add, a, b);
      }
      if (node.text == "-") {
        return program_.arena().pw(Fn::Sub, a, b);
      }
      if (node.text == "*") {
        return program_.arena().pw(Fn::Mul, a, b);
      }
      if (node.text == "/") {
        return program_.arena().pw(Fn::Div, a, b);
      }
      return program_.arena().pw(Fn::Pow, a, b);
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
        Environment inner;
        return lower(definition.a, inner);
      }
      if (components_.count(node.text) != 0) {
        fail(node, "`" + node.text + "` is a grid component and must be called with coordinates");
      }
      if (Builtins.count(node.text) != 0) {
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
        for (std::size_t i = 0; i < definition.params.size(); ++i) {
          inner[definition.params[i]] = lower(node.args[i], env);
        }
        return lower(definition.a, inner);
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
  Program& program_;
  std::map<std::string, SurfaceId> defs_;
  std::vector<std::string> channels_;
};

} // namespace

Program compileSderiv(const std::vector<SderivOutput>& outputs) {
  Program program;
  std::vector<std::string> channelOrder;

  for (const auto& output : outputs) {
    SurfaceArena surface;
    const auto tokens = tokenize(output.source);
    Parser parser(tokens, surface);
    const ParsedProgram parsed = parser.program();

    std::set<std::string> defNames;
    for (const SurfaceId id : parsed.defs) {
      defNames.insert(surface[id].text);
    }
    const ComponentTable components = checkGrids(parsed.grids, defNames);

    std::vector<GridId> gridIds;
    gridIds.reserve(parsed.grids.size());
    for (const auto& grid : parsed.grids) {
      datafield::GridDesc desc;
      // TODO(verify): field names taken from the Lua FieldSpec. Check against
      // Reader/Datafield/Grid.h — the one type here this frontend does not own.
      desc.kind = grid.kind;
      desc.file = grid.file;
      desc.interpolation = grid.interpolation;
      desc.components = grid.components;
      // internGrid dedupes, so two outputs naming the same grid share a GridId
      // and the backend loads the file once.
      gridIds.push_back(program.internGrid(desc));
    }

    Lowering lowering(surface, parsed, components, gridIds, program);
    const NodeId root = lowering.lower(parsed.expression);

    for (const auto& name : lowering.channels()) {
      if (std::find(channelOrder.begin(), channelOrder.end(), name) == channelOrder.end()) {
        channelOrder.push_back(name);
      }
    }
    program.addOutput(output.name, output.type, root);
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

} // namespace seissol::expr
