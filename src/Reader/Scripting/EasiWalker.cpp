// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "EasiWalker.h"

#include "Expr/Core.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Reader/Datafield/AsagiLite.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/CompiledReader.h"
#include "Reader/Scripting/DataReader.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/LuaFrontend.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::reader::scripting {

namespace {

namespace fs = std::filesystem;
using expr::Fn;
using expr::NodeId;
using expr::NoNode;

/// A value a component gives: its node, and where in the files the components are that compute it.
struct Value {
  NodeId node{NoNode};
  std::vector<expr::SourcePosition> from;
};

/// The values of the dimensions or parameters of a component, by name: in the order of easi's sets.
using Values = std::map<std::string, Value>;

std::string listOf(const std::set<std::string>& names) {
  std::string text;
  for (const auto& name : names) {
    text += (text.empty() ? "" : ", ") + name;
  }
  return text;
}

std::set<std::string> keysOf(const Values& values) {
  std::set<std::string> keys;
  for (const auto& [name, value] : values) {
    keys.insert(name);
  }
  return keys;
}

class Walker {
  public:
  explicit Walker(EasiQuery query) : query_(std::move(query)) {
    for (const auto& dimension : query_.dimensions) {
      rootDimensions_.insert(dimension);
    }
  }

  expr::Program translate(const std::string& path) {
    const std::string file = loadPath(path, "");
    const YAML::Node root = load(file);
    program_.addMeta(expr::MetaEntry{"generator", {expr::MetaValue{"seissol easi walker", true}}});
    // the dimensions, in the order of easi's set, then the group where a filter reads it
    Values in;
    for (const auto& dimension : rootDimensions_) {
      program_.addInput(dimension, DataType::F64);
      in[dimension] = Value{arena().field(dimension), {}};
    }
    const std::set<std::string> wanted(query_.parameters.begin(), query_.parameters.end());
    files_.push_back(file);
    const Values values = evaluate(root, in, wanted, NoNode);
    files_.pop_back();
    for (const auto& parameter : query_.parameters) {
      const auto found = values.find(parameter);
      const NodeId node = found != values.end() && found->second.node != NoNode
                              ? found->second.node
                              : arena().konst(std::numeric_limits<double>::quiet_NaN());
      program_.addOutput(parameter, DataType::F64, node);
      if (found != values.end()) {
        for (const auto& position : found->second.from) {
          program_.addSource(
              expr::SourceOf{expr::SourceOf::What::Output, NoNode, parameter, position});
        }
      }
    }
    for (auto& [condition, message] : checks_) {
      program_.addCheck(condition, message);
    }
    return std::move(program_);
  }

  std::set<std::string> parameters(const std::string& path) {
    const std::string file = loadPath(path, "");
    files_.push_back(file);
    const auto supplied = suppliedBy(load(file), rootDimensions_);
    files_.pop_back();
    return supplied;
  }

  private:
  expr::Arena& arena() { return program_.arena(); }
  NodeId konst(double value) { return arena().konst(value); }
  NodeId pw(Fn fn, NodeId a, NodeId b) { return arena().pw(fn, a, b); }
  NodeId pw(Fn fn, NodeId a, NodeId b, NodeId c) { return arena().pw(fn, a, b, c); }

  /// `a and b`, where NoNode is true.
  NodeId both(NodeId a, NodeId b) {
    if (a == NoNode) {
      return b;
    }
    if (b == NoNode) {
      return a;
    }
    return pw(Fn::And, a, b);
  }
  /// `a or b`, where NoNode is true.
  NodeId either(NodeId a, NodeId b) {
    if (a == NoNode || b == NoNode) {
      return NoNode;
    }
    return pw(Fn::Or, a, b);
  }
  NodeId negation(NodeId a) { return pw(Fn::Eq, a, konst(0.0)); }

  [[nodiscard]] const std::string& file() const { return files_.back(); }

  /// easi's reference to a node: the file and the line.
  [[nodiscard]] std::string where(const YAML::Node& node) const {
    return file() + "@" + std::to_string(node.Mark().line + 1);
  }

  [[nodiscard]] expr::SourcePosition positionOf(const YAML::Node& node) const {
    return expr::SourcePosition{file(), node.Mark().line + 1, node.Mark().column + 1};
  }

  [[noreturn]] void fail(const YAML::Node& node, const std::string& what) const {
    throw EasiModelError(where(node) + ": " + what);
  }

  /// That every point which reaches here (`reach`, NoNode for every point) meets `condition`:
  /// a check of the program, with easi's message where it does not.
  void require(NodeId reach, NodeId condition, const std::string& message) {
    const NodeId holds = reach == NoNode ? condition : pw(Fn::Or, negation(reach), condition);
    checks_.emplace_back(holds, message);
  }

  /// The file `name` an !Include names, as easi finds it: next to the file that includes it, else
  /// as it is given.
  [[nodiscard]] std::string loadPath(const std::string& name, const std::string& current) const {
    const fs::path next(name);
    const fs::path next0 =
        next.is_relative() && !current.empty() ? fs::path(current).parent_path() / next : next;
    std::error_code error;
    if (fs::exists(next0, error)) {
      return fs::canonical(next0, error).string();
    }
    if (fs::exists(next, error)) {
      return fs::canonical(next, error).string();
    }
    return next.string();
  }

  /// The data file `name` of a component (ASAGI), as easi finds it: next to the current file,
  /// where it is there, else as it is given.
  [[nodiscard]] std::string dataPath(const std::string& name) const {
    const fs::path entry(name);
    if (entry.is_relative()) {
      const fs::path next = fs::path(file()).parent_path() / entry;
      std::error_code error;
      if (fs::exists(next, error)) {
        return fs::canonical(next, error).string();
      }
    }
    return entry.string();
  }

  YAML::Node load(const std::string& path) const {
    try {
      return YAML::LoadFile(path);
    } catch (const YAML::Exception& error) {
      throw EasiModelError(path + ": " + error.what());
    }
  }

  template <typename T>
  T as(const YAML::Node& node, const char* key) const {
    const YAML::Node value = node[key];
    if (!value) {
      fail(node, std::string("the component has no `") + key + "`");
    }
    try {
      return value.as<T>();
    } catch (const YAML::Exception& error) {
      fail(value, std::string("`") + key + "` is not of the form easi reads: " + error.what());
    }
  }

  /// The components of a composite: its `components`, one or a list.
  static std::vector<YAML::Node> childrenOf(const YAML::Node& node) {
    std::vector<YAML::Node> children;
    if (node.IsMap() && node["components"]) {
      const YAML::Node components = node["components"];
      if (components.IsSequence()) {
        for (const auto& child : components) {
          children.push_back(child);
        }
      } else {
        children.push_back(components);
      }
    }
    return children;
  }

  /// The tags the walker does not translate yet.
  [[noreturn]] void unsupported(const YAML::Node& node) const {
    throw EasiUnsupported(where(node) + ": the walker does not translate " + node.Tag() +
                          " (yet); easi evaluates the file");
  }

  /// The root of the file an !Include names, and the file.
  std::pair<YAML::Node, std::string> included(const YAML::Node& node) {
    const std::string path = loadPath(node.as<std::string>(), file());
    return {load(path), path};
  }

  /// The values a map gives at `in`, its own outputs, before its components: those of a
  /// ConstantMap, an AffineMap, a PolynomialMap, an EvalModel and an ASAGI grid.
  Values mapped(const YAML::Node& node, const Values& in, NodeId reach) {
    const std::string tag = node.Tag();
    Values out;
    if (tag == "!ConstantMap") {
      for (const auto& [name, value] : as<std::map<std::string, double>>(node, "map")) {
        out[name] = Value{konst(value), {positionOf(node)}};
      }
      return out;
    }
    if (tag == "!AffineMap") {
      const auto matrix = as<std::map<std::string, std::vector<double>>>(node, "matrix");
      const auto translation = as<std::map<std::string, double>>(node, "translation");
      if (matrix.size() != translation.size()) {
        fail(node, "Matrix and translation must have the same size in an affine map.");
      }
      // easi's product x * M from zero, a dimension after the other, then the translation; a
      // translation is paired with the row of the same position
      auto shift = translation.begin();
      for (const auto& [name, row] : matrix) {
        if (row.size() != in.size()) {
          fail(node, "Number of matrix entries does not match number of input parameters.");
        }
        NodeId value = konst(0.0);
        std::size_t k = 0;
        for (const auto& [dimension, coordinate] : in) {
          value = pw(Fn::Add, value, pw(Fn::Mul, coordinate.node, konst(row[k++])));
        }
        out[name] = Value{pw(Fn::Add, value, konst((shift++)->second)), {positionOf(node)}};
      }
      return out;
    }
    if (tag == "!PolynomialMap") {
      if (in.size() != 1) {
        fail(node, "Polynomial map requires 1D input.");
      }
      const NodeId x = in.begin()->second.node;
      const auto coefficients = as<std::map<std::string, std::vector<double>>>(node, "map");
      std::size_t count = std::numeric_limits<std::size_t>::max();
      for (const auto& [name, values] : coefficients) {
        if (count != std::numeric_limits<std::size_t>::max() && values.size() != count) {
          fail(node,
               "All parameters in a polynomial map must have the same number of coefficients.");
        }
        count = values.size();
      }
      for (const auto& [name, values] : coefficients) {
        if (values.empty()) {
          fail(node, "a polynomial map has coefficients");
        }
        // Horner's scheme as easi has it: c_k + y * x
        NodeId value = konst(values[0]);
        for (std::size_t k = 1; k < values.size(); ++k) {
          value = pw(Fn::Add, konst(values[k]), pw(Fn::Mul, value, x));
        }
        out[name] = Value{value, {positionOf(node)}};
      }
      return out;
    }
    if (tag == "!EvalModel") {
      if (!node["model"]) {
        fail(node, "the component has no `model`");
      }
      const auto names = as<std::vector<std::string>>(node, "parameters");
      const std::set<std::string> wanted(names.begin(), names.end());
      Values out = evaluate(node["model"], in, wanted, reach);
      // a parameter the model does not give fails a check there; NaN here
      for (const auto& name : wanted) {
        if (out.count(name) == 0) {
          out[name] = Value{konst(std::numeric_limits<double>::quiet_NaN()), {}};
        }
      }
      return out;
    }
    if (tag == "!ASAGI") {
      return grid(node, in);
    }
    if (tag == "!FunctionMap") {
      Values out;
      for (const auto& [name, code] : as<std::map<std::string, std::string>>(node, "map")) {
        ImpalaFunction function(*this, node, in, name);
        out[name] = Value{function.translate("{" + code + "}", reach), {positionOf(node)}};
      }
      return out;
    }
    if (tag == "!LuaMap") {
      return luaMapped(node, in, reach);
    }
    if (tag == "!Include") {
      fail(node, "an !Include is no map");
    }
    unsupported(node);
  }

  /// The line of the file the first line of the text of `scalar` is on: for a block scalar (`|`,
  /// `>`), the line after its indicator.
  [[nodiscard]] int firstLineOf(const YAML::Node& scalar) const {
    const auto mark = scalar.Mark();
    std::ifstream stream(file());
    const std::string text((std::istreambuf_iterator<char>(stream)),
                           std::istreambuf_iterator<char>());
    const bool block = mark.pos >= 0 && static_cast<std::size_t>(mark.pos) < text.size() &&
                       (text[mark.pos] == '|' || text[mark.pos] == '>');
    return mark.line + (block ? 2 : 1);
  }

  /// The values of an !LuaMap at `in`: its function `f`, translated by the Lua frontend as easi
  /// runs it (Lua 5.4 with lmathx, the chunk once and `f` at every point).
  Values luaMapped(const YAML::Node& node, const Values& in, NodeId reach) {
    std::set<std::string> returns;
    const YAML::Node given = node["returns"];
    if (!given) {
      fail(node, "the component has no `returns`");
    }
    try {
      if (given.IsScalar()) {
        returns.insert(given.as<std::string>());
      } else {
        const auto names = given.as<std::vector<std::string>>();
        returns.insert(names.begin(), names.end());
      }
    } catch (const YAML::Exception& error) {
      fail(given, std::string("`returns` is not of the form easi reads: ") + error.what());
    }
    const auto code = as<std::string>(node, "function");
    std::vector<std::pair<std::string, NodeId>> inputs;
    for (const auto& [name, value] : in) {
      inputs.emplace_back(name, value.node);
    }
    const lua::LuaOrigin origin{file(), firstLineOf(node["function"]), '@'};
    try {
      auto translation = lua::translateLuaMap(program_, code, inputs, returns, reach, origin);
      for (auto& check : translation.checks) {
        checks_.emplace_back(check.holds, std::move(check.message));
      }
      std::vector<expr::SourcePosition> from{positionOf(node)};
      for (const int line : translation.lines) {
        from.push_back(expr::SourcePosition{file(), origin.firstLine + line - 1, 1});
      }
      Values out;
      for (const auto& [name, value] : translation.values) {
        out[name] = Value{value, from};
      }
      return out;
    } catch (const lua::LuaUnsupported& error) {
      throw EasiUnsupported(origin.at(error.line()) +
                            ": the walker does not translate the Lua of this !LuaMap: " +
                            error.message() + "; easi evaluates the file");
    } catch (const lua::LuaError& error) {
      throw EasiModelError(origin.at(error.line()) + ": " + error.message());
    }
  }

  /// The file and the variable of an !ASAGI, and the axes ASAGI reads there, as many as `in` has
  /// dimensions.
  std::vector<datafield::AsagiAxis> asagiAxes(const YAML::Node& node, const Values& in) {
    const std::string path = dataPath(as<std::string>(node, "file"));
    const std::string variable = node["var"] ? as<std::string>(node, "var") : std::string("data");
    const auto axes = datafield::asagiAxesOf(path, variable);
    if (!axes.has_value()) {
      fail(node, "Could not open " + path + " with ASAGI.");
    }
    if (axes->size() != in.size()) {
      fail(node,
           "ASAGI requires " + std::to_string(axes->size()) + "D input (got " + listOf(keysOf(in)) +
               ").");
    }
    for (const auto& axis : *axes) {
      if (axis.nodes > 1 && !(axis.first < axis.last)) {
        throw EasiUnsupported(where(node) + ": a grid of ASAGI whose coordinates decrease");
      }
    }
    return *axes;
  }

  /// The Impala of a !FunctionMap, translated as easi reads it: easi >= 1.5 converts it into Lua
  /// (convertImpalaToLua()) and runs that with Lua 5.4 and lmathx. Its grammar is easi's, the
  /// grouping of its chains included: easi reads `a - b + c` as `a - (b + c)` and `a / b * c` as
  /// `a / (b * c)`, which the walker computes as well, with a warning that names the chain. A
  /// number without a point is an integer of Lua, whose zero has no sign; a function is that of
  /// Lua's math (with lmathx) of its name. A variable is Lua's, global: read where no assignment
  /// before gives it a value, it has that of an earlier point, which the walker refuses; a point
  /// where the code returns nothing stops easi, which a check of the program does here.
  class ImpalaFunction {
public:
    ImpalaFunction(Walker& walker, const YAML::Node& node, const Values& in, std::string name)
        : walker_(walker), node_(node), in_(in), name_(std::move(name)) {}

    NodeId translate(const std::string& code, NodeId reach) {
      tokens_ = tokenize(code);
      position_ = 0;
      done_ = no();
      value_ = walker_.konst(std::numeric_limits<double>::quiet_NaN());
      body(yes());
      // a point that returns nothing: easi stops there, the field not being a number
      if (!isConst(done_, 1.0)) {
        walker_.require(reach,
                        done_,
                        walker_.where(node_) + ": Tried key " + name_ +
                            " but it did not return a number.");
      }
      return value_;
    }

private:
    enum class Tok : std::uint8_t {
      Variable,
      Constant,
      Comma,
      Assign,
      If,
      Else,
      Return,
      Add,
      Sub,
      Mul,
      Div,
      Or,
      And,
      Eq,
      Ne,
      Lt,
      Le,
      Gt,
      Ge,
      Open,
      Close,
      BracketOpen,
      BracketClose,
      Semicolon,
      End
    };
    struct Token {
      Tok kind;
      std::string text;
    };
    /// A value of the code: its node, whether it is an integer of Lua, and the level of its
    /// outermost operation (1 a sum or a difference, 2 a product or a quotient, 0 another)
    struct Term {
      NodeId node;
      bool integer;
      int level;
    };

    /// easi's tokenizer: the words and operators first, longest first, by prefix -- `iffy` is
    /// `if` and `fy` --, then numbers without an exponent, names, and space.
    std::vector<Token> tokenize(const std::string& code) const {
      static const std::vector<std::pair<std::string, Tok>> Words = {
          {"return", Tok::Return}, {"else", Tok::Else},      {"if", Tok::If},
          {"==", Tok::Eq},         {"!=", Tok::Ne},          {"<=", Tok::Le},
          {">=", Tok::Ge},         {"&&", Tok::And},         {"||", Tok::Or},
          {"<", Tok::Lt},          {">", Tok::Gt},           {"+", Tok::Add},
          {"-", Tok::Sub},         {"*", Tok::Mul},          {"/", Tok::Div},
          {",", Tok::Comma},       {"=", Tok::Assign},       {";", Tok::Semicolon},
          {"(", Tok::BracketOpen}, {")", Tok::BracketClose}, {"{", Tok::Open},
          {"}", Tok::Close}};
      std::vector<Token> tokens;
      std::size_t at = 0;
      const auto digit = [&](std::size_t i) {
        return i < code.size() && std::isdigit(static_cast<unsigned char>(code[i])) != 0;
      };
      while (at < code.size()) {
        bool found = false;
        for (const auto& [word, kind] : Words) {
          if (code.compare(at, word.size(), word) == 0) {
            tokens.push_back(Token{kind, word});
            at += word.size();
            found = true;
            break;
          }
        }
        if (found) {
          continue;
        }
        const auto c = static_cast<unsigned char>(code[at]);
        if (digit(at)) {
          std::size_t end = at;
          while (digit(end)) {
            ++end;
          }
          if (end < code.size() && code[end] == '.') {
            ++end;
            while (digit(end)) {
              ++end;
            }
          }
          tokens.push_back(Token{Tok::Constant, code.substr(at, end - at)});
          at = end;
        } else if (std::isalpha(c) != 0) {
          std::size_t end = at + 1;
          while (end < code.size() &&
                 (std::isalnum(static_cast<unsigned char>(code[end])) != 0 || code[end] == '_')) {
            ++end;
          }
          tokens.push_back(Token{Tok::Variable, code.substr(at, end - at)});
          at = end;
        } else if (std::isspace(c) != 0) {
          ++at;
        } else {
          walker_.fail(node_, "Tokenization error: " + code.substr(at));
        }
      }
      tokens.push_back(Token{Tok::End, ""});
      return tokens;
    }

    [[nodiscard]] const Token& peek(std::size_t ahead = 0) const {
      return tokens_[std::min(position_ + ahead, tokens_.size() - 1)];
    }
    bool accept(Tok kind) {
      if (peek().kind == kind) {
        ++position_;
        return true;
      }
      return false;
    }
    std::string expect(Tok kind) {
      if (peek().kind != kind) {
        walker_.fail(node_,
                     "the Impala of `" + name_ + "` has `" + peek().text +
                         "` where easi's conversion expects another token");
      }
      return tokens_[position_++].text;
    }

    // truth values, folded where known
    NodeId yes() { return walker_.konst(1.0); }
    NodeId no() { return walker_.konst(0.0); }
    [[nodiscard]] bool isConst(NodeId id, double value) const {
      const auto& node = walker_.program_.arena()[id];
      return node.kind == expr::Kind::Const && node.value == value;
    }
    NodeId conjunction(NodeId a, NodeId b) {
      if (isConst(a, 0.0) || isConst(b, 0.0)) {
        return no();
      }
      if (isConst(a, 1.0)) {
        return b;
      }
      if (isConst(b, 1.0)) {
        return a;
      }
      return walker_.pw(Fn::And, a, b);
    }
    NodeId disjunction(NodeId a, NodeId b) {
      if (isConst(a, 1.0) || isConst(b, 1.0)) {
        return yes();
      }
      if (isConst(a, 0.0)) {
        return b;
      }
      if (isConst(b, 0.0)) {
        return a;
      }
      return walker_.pw(Fn::Or, a, b);
    }
    NodeId complement(NodeId a) {
      if (isConst(a, 0.0)) {
        return yes();
      }
      if (isConst(a, 1.0)) {
        return no();
      }
      return walker_.negation(a);
    }
    NodeId choice(NodeId condition, NodeId a, NodeId b) {
      if (isConst(condition, 1.0) || a == b) {
        return a;
      }
      if (isConst(condition, 0.0)) {
        return b;
      }
      return walker_.pw(Fn::Select, condition, a, b);
    }
    /// An integer of Lua has no negative zero: where an operation of integers gives one, zero.
    NodeId integral(NodeId node) { return walker_.pw(Fn::Add, node, walker_.konst(0.0)); }

    /// The value of the node `id` where Lua's compiler computes it: a numeral, and a sum,
    /// difference or product of such.
    std::optional<double> constantOf(NodeId id) const {
      const auto& node = walker_.arena()[id];
      if (node.kind == expr::Kind::Const) {
        return node.value;
      }
      if (node.kind != expr::Kind::PW ||
          (node.fn != Fn::Add && node.fn != Fn::Sub && node.fn != Fn::Mul)) {
        return std::nullopt;
      }
      const auto a = constantOf(node.a);
      const auto b = constantOf(node.b);
      if (!a.has_value() || !b.has_value()) {
        return std::nullopt;
      }
      return node.fn == Fn::Add ? *a + *b : (node.fn == Fn::Sub ? *a - *b : *a * *b);
    }

    /// The value of `term` where it is an integer constant of Lua from -127 to 127.
    std::optional<double> smallInteger(const Term& term) const {
      if (!term.integer) {
        return std::nullopt;
      }
      const auto value = constantOf(term.node);
      if (!value.has_value() || !(std::fabs(*value) <= 127.0)) {
        return std::nullopt;
      }
      return *value + 0.0;
    }

    /// `{ statements }` under `guard`. Whether every path through it returns.
    bool body(NodeId guard) {
      expect(Tok::Open);
      bool returns = false;
      while (!accept(Tok::Close)) {
        if (peek().kind == Tok::End) {
          expect(Tok::Close);
        }
        // statements after a return are never reached
        const NodeId active = returns ? no() : guard;
        if (peek().kind == Tok::If) {
          returns = ifStatement(active) || returns;
        } else if (accept(Tok::Return)) {
          const Term term = expression();
          expect(Tok::Semicolon);
          const NodeId taking = conjunction(active, complement(done_));
          value_ = choice(taking, term.node, value_);
          done_ = disjunction(done_, active);
          returns = true;
        } else if (peek().kind == Tok::Variable && peek(1).kind == Tok::Assign) {
          const std::string variable = expect(Tok::Variable);
          expect(Tok::Assign);
          const Term term = expression();
          expect(Tok::Semicolon);
          declared_.insert(variable);
          const NodeId taking = conjunction(active, complement(done_));
          const auto old = variables_.find(variable);
          if (old == variables_.end()) {
            variables_[variable] = Term{
                choice(taking, term.node, walker_.konst(std::numeric_limits<double>::quiet_NaN())),
                term.integer,
                0};
          } else {
            old->second = Term{choice(taking, term.node, old->second.node),
                               term.integer && old->second.integer,
                               0};
          }
          if (!returns) {
            assigned_.insert(variable);
          }
        } else {
          walker_.fail(node_,
                       "the Impala of `" + name_ +
                           "` has an expression for a statement, which the Lua easi makes of it "
                           "does not read");
        }
      }
      return returns;
    }

    /// `if (condition) { ... } else { ... }`. Whether both branches return.
    bool ifStatement(NodeId guard) {
      expect(Tok::If);
      const NodeId condition = boolean();
      const auto before = assigned_;
      const bool thenReturns = body(conjunction(guard, condition));
      auto afterThen = assigned_;
      assigned_ = before;
      bool elseReturns = false;
      if (accept(Tok::Else)) {
        elseReturns = body(conjunction(guard, complement(condition)));
      }
      // assigned after it: in both branches, or in the one that does not return
      if (thenReturns && !elseReturns) {
        // assigned_ is that of the else branch
      } else if (elseReturns && !thenReturns) {
        assigned_ = afterThen;
      } else {
        std::set<std::string> both;
        std::set_intersection(afterThen.begin(),
                              afterThen.end(),
                              assigned_.begin(),
                              assigned_.end(),
                              std::inserter(both, both.end()));
        assigned_ = std::move(both);
      }
      return thenReturns && elseReturns;
    }

    Term expression() { return sum(); }

    /// easi's chains group to the right: `a - b + c` is `a - (b + c)`.
    Term sum() {
      const Term left = product();
      for (const auto [kind, fn] : {std::pair{Tok::Add, Fn::Add}, std::pair{Tok::Sub, Fn::Sub}}) {
        if (accept(kind)) {
          const Term right = sum();
          if (kind == Tok::Sub && right.level == 1) {
            warnGrouping("-", "+ or -");
          }
          // Lua's code generator computes `a - k`, for an integer constant k from -127 to 127, as
          // `a + -k` (OP_ADDI): the same, but at k = 0, where -0.0 - 0 is -0.0, -0.0 + 0 is 0.0
          if (const auto k = smallInteger(right); kind == Tok::Sub && k.has_value()) {
            return Term{walker_.pw(Fn::Add, left.node, walker_.konst(-*k + 0.0)), left.integer, 1};
          }
          return Term{walker_.pw(fn, left.node, right.node), left.integer && right.integer, 1};
        }
      }
      return left;
    }

    Term product() {
      const Term left = unary();
      if (accept(Tok::Mul)) {
        const Term right = product();
        const bool integer = left.integer && right.integer;
        const NodeId node = walker_.pw(Fn::Mul, left.node, right.node);
        return Term{integer ? integral(node) : node, integer, 2};
      }
      if (accept(Tok::Div)) {
        const Term right = product();
        if (right.level == 2) {
          warnGrouping("/", "* or /");
        }
        return Term{walker_.pw(Fn::Div, left.node, right.node), false, 2};
      }
      return left;
    }

    Term unary() {
      if (accept(Tok::Add)) {
        walker_.fail(node_,
                     "the Impala of `" + name_ +
                         "` has a unary +, which the Lua easi makes of it does not read");
      }
      const bool negative = accept(Tok::Sub);
      const Term operand = atom();
      if (!negative) {
        return operand;
      }
      const NodeId node = walker_.arena().pw(Fn::Neg, operand.node);
      return Term{operand.integer ? integral(node) : node, operand.integer, 0};
    }

    Term atom() {
      if (accept(Tok::BracketOpen)) {
        const Term inner = expression();
        expect(Tok::BracketClose);
        return Term{inner.node, inner.integer, 0};
      }
      if (peek().kind == Tok::Variable && peek(1).kind == Tok::BracketOpen) {
        return call();
      }
      if (peek().kind == Tok::Variable) {
        const std::string name = expect(Tok::Variable);
        // an input first, as easi's conversion reads a name, whatever is assigned to it
        const auto input = in_.find(name);
        if (input != in_.end()) {
          return Term{input->second.node, false, 0};
        }
        if (declared_.count(name) == 0) {
          walker_.fail(node_, "Unknown variable: " + name);
        }
        if (assigned_.count(name) == 0) {
          throw EasiUnsupported(walker_.where(node_) + ": the Impala of `" + name_ + "` reads `" +
                                name +
                                "` where no assignment before gives it a value, where easi "
                                "takes that of an earlier point");
        }
        const Term& value = variables_.at(name);
        return Term{value.node, value.integer, 0};
      }
      if (peek().kind == Tok::Constant) {
        const std::string text = expect(Tok::Constant);
        double value = 0.0;
        try {
          value = expr::numberOf(text);
        } catch (const std::invalid_argument& error) {
          walker_.fail(node_, error.what());
        }
        return Term{walker_.konst(value), text.find('.') == std::string::npos, 0};
      }
      walker_.fail(node_, "Unexpected token when parsing expression: " + peek().text);
    }

    /// A function of Lua's math, with lmathx.
    Term call() {
      const std::string function = expect(Tok::Variable);
      std::vector<Term> args;
      expect(Tok::BracketOpen);
      if (!accept(Tok::BracketClose)) {
        args.push_back(expression());
        while (!accept(Tok::BracketClose)) {
          expect(Tok::Comma);
          args.push_back(expression());
        }
      }
      const auto arity = [&](std::size_t low, std::size_t high) {
        if (args.size() < low || args.size() > high) {
          throw EasiUnsupported(walker_.where(node_) + ": the Impala of `" + name_ + "` calls " +
                                function + " with " + std::to_string(args.size()) + " arguments");
        }
      };
      const auto unaryOp = [&](Fn fn) {
        arity(1, 1);
        return Term{walker_.arena().pw(fn, args[0].node), false, 0};
      };
      static const std::map<std::string, Fn> Unary = {
          {"sqrt", Fn::Sqrt},   {"exp", Fn::Exp},     {"log10", Fn::Log10}, {"log2", Fn::Log2},
          {"sin", Fn::Sin},     {"cos", Fn::Cos},     {"tan", Fn::Tan},     {"asin", Fn::Asin},
          {"acos", Fn::Acos},   {"sinh", Fn::Sinh},   {"cosh", Fn::Cosh},   {"tanh", Fn::Tanh},
          {"asinh", Fn::Asinh}, {"acosh", Fn::Acosh}, {"atanh", Fn::Atanh}, {"erf", Fn::Erf},
          {"floor", Fn::Floor}, {"ceil", Fn::Ceil},   {"round", Fn::Round}, {"fabs", Fn::Abs}};
      if (const auto found = Unary.find(function); found != Unary.end()) {
        return unaryOp(found->second);
      }
      if (function == "abs") {
        // Lua's own: of an integer an integer
        arity(1, 1);
        return Term{walker_.arena().pw(Fn::Abs, args[0].node), args[0].integer, 0};
      }
      if (function == "log") {
        arity(1, 2);
        if (args.size() == 1) {
          return unaryOp(Fn::Log);
        }
        // lmathx: of the base 10 log10, of 2 log2, else the quotient of the logarithms
        const NodeId x = args[0].node;
        const NodeId base = args[1].node;
        const NodeId other =
            walker_.pw(Fn::Div, walker_.arena().pw(Fn::Log, x), walker_.arena().pw(Fn::Log, base));
        return Term{choice(walker_.pw(Fn::Eq, base, walker_.konst(10.0)),
                           walker_.arena().pw(Fn::Log10, x),
                           choice(walker_.pw(Fn::Eq, base, walker_.konst(2.0)),
                                  walker_.arena().pw(Fn::Log2, x),
                                  other)),
                    false,
                    0};
      }
      if (function == "atan") {
        arity(1, 2);
        if (args.size() == 1) {
          return unaryOp(Fn::Atan);
        }
        return Term{walker_.pw(Fn::Atan2, args[0].node, args[1].node), false, 0};
      }
      if (function == "atan2" || function == "pow") {
        arity(2, 2);
        return Term{walker_.pw(function == "pow" ? Fn::Pow : Fn::Atan2, args[0].node, args[1].node),
                    false,
                    0};
      }
      if (function == "trunc") {
        arity(1, 1);
        const NodeId x = args[0].node;
        return Term{choice(walker_.pw(Fn::Lt, x, walker_.konst(0.0)),
                           walker_.arena().pw(Fn::Ceil, x),
                           walker_.arena().pw(Fn::Floor, x)),
                    false,
                    0};
      }
      if (function == "fmin" || function == "fmax") {
        arity(1, std::numeric_limits<std::size_t>::max());
        NodeId m = args[0].node;
        for (std::size_t i = 1; i < args.size(); ++i) {
          m = walker_.pw(function == "fmin" ? Fn::Min : Fn::Max, m, args[i].node);
        }
        return Term{m, false, 0};
      }
      if (function == "min" || function == "max") {
        // Lua's own: the first argument that no later one is below (above)
        arity(1, std::numeric_limits<std::size_t>::max());
        Term m = args[0];
        for (std::size_t i = 1; i < args.size(); ++i) {
          const NodeId better = function == "min" ? walker_.pw(Fn::Lt, args[i].node, m.node)
                                                  : walker_.pw(Fn::Lt, m.node, args[i].node);
          m = Term{choice(better, args[i].node, m.node), m.integer && args[i].integer, 0};
        }
        return m;
      }
      if (function == "deg" || function == "rad") {
        arity(1, 1);
        constexpr double Pi = 3.141592653589793238462643383279502884;
        return Term{walker_.pw(Fn::Mul,
                               args[0].node,
                               walker_.konst(function == "deg" ? 180.0 / Pi : Pi / 180.0)),
                    false,
                    0};
      }
      throw EasiUnsupported(walker_.where(node_) + ": the Impala of `" + name_ + "` calls " +
                            function + ", which the walker does not translate");
    }

    /// `( boolean )`, a comparison, and their `&&` and `||`, grouped to the right.
    NodeId boolean() {
      const NodeId left = conjunctionOfAtoms();
      if (accept(Tok::Or)) {
        return disjunction(left, boolean());
      }
      return left;
    }
    NodeId conjunctionOfAtoms() {
      const NodeId left = comparison();
      if (accept(Tok::And)) {
        return conjunction(left, conjunctionOfAtoms());
      }
      return left;
    }
    NodeId comparison() {
      if (accept(Tok::BracketOpen)) {
        const NodeId inner = boolean();
        expect(Tok::BracketClose);
        return inner;
      }
      const Term left = expression();
      const Tok op = peek().kind;
      if (op != Tok::Eq && op != Tok::Ne && op != Tok::Lt && op != Tok::Le && op != Tok::Gt &&
          op != Tok::Ge) {
        walker_.fail(node_, "No comparison operator given");
      }
      ++position_;
      const Term right = expression();
      switch (op) {
      case Tok::Eq:
        return walker_.pw(Fn::Eq, left.node, right.node);
      case Tok::Ne:
        return complement(walker_.pw(Fn::Eq, left.node, right.node));
      case Tok::Lt:
        return walker_.pw(Fn::Lt, left.node, right.node);
      case Tok::Le:
        return walker_.pw(Fn::Le, left.node, right.node);
      case Tok::Gt:
        return walker_.pw(Fn::Lt, right.node, left.node);
      default:
        return walker_.pw(Fn::Le, right.node, left.node);
      }
    }

    /// A warning where easi's grouping changes what a chain means, once per chain and run.
    void warnGrouping(const char* op, const char* following) const {
      static std::mutex mutex;
      static std::set<std::string> warned;
      const std::lock_guard<std::mutex> lock(mutex);
      if (!warned.insert(walker_.where(node_) + "/" + name_ + "/" + op).second) {
        return;
      }
      logWarning() << "easi:" << walker_.where(node_) << ": the Impala of" << name_
                   << "has a chain where" << op << "is followed by" << following
                   << "; easi's conversion into Lua groups it to the right (a - b + c as a - (b + "
                      "c)), and the walker computes what easi computes.";
    }

    Walker& walker_;
    const YAML::Node& node_;
    const Values& in_;
    std::string name_;
    std::vector<Token> tokens_;
    std::size_t position_{0};
    /// whether a return happened, and its value where it did
    NodeId done_{NoNode};
    NodeId value_{NoNode};
    /// the variables, as easi's conversion knows them (by an assignment before, in the text), as
    /// assigned on every path to here, and their values
    std::set<std::string> declared_;
    std::set<std::string> assigned_;
    std::map<std::string, Term> variables_;
  };

  /// A grid of ASAGI, sampled at the dimensions of `in` in their order, as easi samples it.
  Values grid(const YAML::Node& node, const Values& in) {
    static_cast<void>(asagiAxes(node, in));
    const auto parameters = as<std::vector<std::string>>(node, "parameters");
    std::string interpolation =
        node["interpolation"] ? as<std::string>(node, "interpolation") : std::string("linear");
    std::transform(interpolation.begin(),
                   interpolation.end(),
                   interpolation.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    datafield::GridSpec spec;
    spec.kind = "asagi";
    spec.path = dataPath(as<std::string>(node, "file"));
    spec.variable = node["var"] ? as<std::string>(node, "var") : std::string("data");
    spec.components = parameters;
    spec.interpolation = {interpolation};
    spec.boundary = {"clamp"};
    datafield::GridDesc desc;
    try {
      desc = datafield::gridDescOf(spec);
    } catch (const std::invalid_argument& error) {
      fail(node, error.what());
    }
    const auto id = program_.internGrid(desc);
    program_.nameGrid(id, expr::GridNames{"asagi" + std::to_string(id), parameters});
    std::vector<NodeId> coordinates;
    for (const auto& [dimension, value] : in) {
      coordinates.push_back(value.node);
    }
    Values out;
    for (std::size_t p = 0; p < parameters.size(); ++p) {
      out[parameters[p]] =
          Value{arena().lookup(id, static_cast<std::int32_t>(p), coordinates), {positionOf(node)}};
    }
    return out;
  }

  /// A layered model: the nodes of its parameters along the one dimension its map gives, a grid
  /// of them (Grid.h) -- lower and upper the steps previous and next, linear linear.
  Values layers(const YAML::Node& node,
                const Values& in,
                const std::set<std::string>& wanted,
                NodeId reach) {
    if (!node["map"]) {
      fail(node, "the component has no `map`");
    }
    const YAML::Node map = node["map"];
    std::string interpolation = as<std::string>(node, "interpolation");
    std::transform(interpolation.begin(),
                   interpolation.end(),
                   interpolation.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const auto nodes = as<std::map<double, std::vector<double>>>(node, "nodes");
    const auto parameters = as<std::vector<std::string>>(node, "parameters");
    if (interpolation != "lower" && interpolation != "upper" && interpolation != "linear") {
      fail(node, "Invalid interpolation type " + interpolation + ".");
    }
    if (nodes.empty()) {
      fail(node, "a layered model has nodes");
    }
    std::vector<double> knots;
    std::vector<double> data;
    bool infinite = false;
    for (const auto& [depth, values] : nodes) {
      if (values.size() != parameters.size()) {
        fail(node, "Number of parameters must match number of node entries.");
      }
      if (std::isnan(depth)) {
        fail(node, "Not a number nodes are forbidden.");
      }
      if (std::isinf(depth)) {
        if (interpolation == "linear") {
          fail(node, "Infinite nodes are forbidden for interpolation = linear.");
        }
        infinite = true;
      }
      knots.push_back(depth);
      data.insert(data.end(), values.begin(), values.end());
    }

    // the one dimension of its map: its own, or what its components supply
    Value along;
    const auto mapChildren = childrenOf(map);
    if (mapChildren.empty()) {
      const Values out = mapped(map, in, reach);
      if (out.size() != 1) {
        fail(node, "Layered model requires 1D input.");
      }
      along = out.begin()->second;
    } else {
      const auto supplied = suppliedBy(map, keysOf(in));
      if (supplied.size() != 1) {
        fail(node, "Layered model requires 1D input.");
      }
      along = evaluate(map, in, supplied, reach).begin()->second;
    }
    // easi's layers accept every depth but NaN
    const YAML::Node at = map;
    require(reach,
            pw(Fn::Eq, along.node, along.node),
            "Could not find model for point (" + where(at) + ").");

    if (infinite) {
      return layerSteps(node, knots, data, parameters, interpolation, along, wanted, reach);
    }
    datafield::GridSpec spec;
    spec.kind = "inline";
    spec.components = parameters;
    spec.knots = {knots};
    spec.values = data;
    spec.interpolation = {interpolation == "lower"   ? "previous"
                          : interpolation == "upper" ? "next"
                                                     : "linear"};
    spec.boundary = {"clamp"};
    datafield::GridDesc desc;
    try {
      desc = datafield::gridDescOf(spec);
    } catch (const std::invalid_argument& error) {
      fail(node, error.what());
    }
    const auto id = program_.internGrid(desc);
    program_.nameGrid(id, expr::GridNames{"layers" + std::to_string(id), parameters});
    Values out;
    std::set<std::string> missing;
    for (const auto& parameter : wanted) {
      const auto found = std::find(parameters.begin(), parameters.end(), parameter);
      if (found == parameters.end()) {
        missing.insert(parameter);
        continue;
      }
      out[parameter] = Value{
          arena().lookup(id, static_cast<std::int32_t>(found - parameters.begin()), {along.node}),
          {positionOf(node)}};
    }
    if (!missing.empty()) {
      require(reach,
              konst(0.0),
              where(node) + ": Component does not supply all required parameters. Missing: {" +
                  listOf(missing) + "}.");
    }
    return out;
  }

  /// The steps of a layered model with a node at an infinite depth, which a grid has not: easi's
  /// layers one after another, each from a node to the next (the first from minus infinity, the
  /// last to infinity), a depth on a node in the first that contains it, as a select between them.
  Values layerSteps(const YAML::Node& node,
                    const std::vector<double>& knots,
                    const std::vector<double>& data,
                    const std::vector<std::string>& parameters,
                    const std::string& interpolation,
                    const Value& along,
                    const std::set<std::string>& wanted,
                    NodeId reach) {
    constexpr double Infinity = std::numeric_limits<double>::infinity();
    const std::size_t count = parameters.size();
    // per layer its limits and the node whose values it takes
    struct Layer {
      double lower;
      double upper;
      std::size_t node;
    };
    std::vector<Layer> layers;
    for (std::size_t n = 0; n <= knots.size(); ++n) {
      const double lower = n == 0 ? -Infinity : knots[n - 1];
      const double upper = n == knots.size() ? Infinity : knots[n];
      const std::size_t taken = n == 0                     ? 0
                                : n == knots.size()        ? knots.size() - 1
                                : interpolation == "lower" ? n - 1
                                                           : n;
      layers.push_back(Layer{lower, upper, taken});
    }
    Values out;
    std::set<std::string> missing;
    const NodeId d = along.node;
    for (const auto& parameter : wanted) {
      const auto found = std::find(parameters.begin(), parameters.end(), parameter);
      if (found == parameters.end()) {
        missing.insert(parameter);
        continue;
      }
      const auto p = static_cast<std::size_t>(found - parameters.begin());
      NodeId value = konst(std::numeric_limits<double>::quiet_NaN());
      for (std::size_t l = layers.size(); l-- > 0;) {
        const NodeId within = pw(
            Fn::And, pw(Fn::Le, konst(layers[l].lower), d), pw(Fn::Le, d, konst(layers[l].upper)));
        value = pw(Fn::Select, within, konst(data[layers[l].node * count + p]), value);
      }
      out[parameter] = Value{value, {positionOf(node)}};
    }
    if (!missing.empty()) {
      require(reach,
              konst(0.0),
              where(node) + ": Component does not supply all required parameters. Missing: {" +
                  listOf(missing) + "}.");
    }
    return out;
  }

  /// Whether a component, as a component of a composite, accepts a point at `in`: NoNode where it
  /// accepts every point.
  NodeId accepts(const YAML::Node& node, const Values& in) {
    const std::string tag = node.Tag();
    if (tag == "!Include") {
      const auto [root, path] = included(node);
      files_.push_back(path);
      const NodeId accepted = accepts(root, in);
      files_.pop_back();
      return accepted;
    }
    if (tag == "!GroupFilter") {
      const YAML::Node groups = node["groups"];
      if (!groups) {
        fail(node, "the component has no `groups`");
      }
      std::vector<int> listed;
      try {
        if (groups.IsScalar()) {
          listed.push_back(groups.as<int>());
        } else {
          listed = groups.as<std::vector<int>>();
        }
      } catch (const YAML::Exception& error) {
        fail(groups, std::string("`groups` is not of the form easi reads: ") + error.what());
      }
      NodeId accepted = konst(0.0);
      for (const int group : std::set<int>(listed.begin(), listed.end())) {
        accepted = pw(Fn::Or, accepted, pw(Fn::Eq, groupNode(), konst(group)));
      }
      return accepted;
    }
    if (tag == "!AxisAlignedCuboidalDomainFilter") {
      const auto limits = as<std::map<std::string, std::pair<double, double>>>(node, "limits");
      matchDimensions(node, limits, in);
      NodeId accepted = NoNode;
      for (const auto& [dimension, range] : limits) {
        const NodeId x = in.at(dimension).node;
        accepted = both(
            accepted,
            pw(Fn::And, pw(Fn::Le, konst(range.first), x), pw(Fn::Le, x, konst(range.second))));
      }
      return accepted;
    }
    if (tag == "!SphericalDomainFilter") {
      const auto centre = as<std::map<std::string, double>>(node, "centre");
      const auto radius = as<double>(node, "radius");
      matchDimensions(node, centre, in);
      NodeId distance = konst(0.0);
      for (const auto& [dimension, coordinate] : centre) {
        const NodeId difference = pw(Fn::Sub, in.at(dimension).node, konst(coordinate));
        distance = pw(Fn::Add, distance, pw(Fn::Mul, difference, difference));
      }
      return pw(Fn::Le, distance, konst(radius * radius));
    }
    if (tag == "!ASAGI") {
      // within the grid, as ASAGI bounds it: from the first coordinate by its scaling, a dimension
      // of one node unbounded (Level::open())
      const auto axes = asagiAxes(node, in);
      NodeId accepted = NoNode;
      std::size_t d = 0;
      for (const auto& [dimension, value] : in) {
        const auto& axis = axes[d++];
        double lo = -std::numeric_limits<double>::infinity();
        double hi = std::numeric_limits<double>::infinity();
        if (axis.nodes > 1) {
          const double scaling = (axis.last - axis.first) / static_cast<double>(axis.nodes - 1);
          const double far = axis.first + scaling * static_cast<double>(axis.nodes - 1);
          lo = std::min(axis.first, far);
          hi = std::max(axis.first, far);
        }
        accepted =
            both(accepted,
                 pw(Fn::And, pw(Fn::Le, konst(lo), value.node), pw(Fn::Le, value.node, konst(hi))));
      }
      return accepted;
    }
    // the other maps, !Switch, !Any and !IdentityMap
    return NoNode;
  }

  /// The group of the point: the input `group`, or 0 for points without one.
  NodeId groupNode() {
    if (!query_.group) {
      return konst(0.0);
    }
    if (!groupDeclared_) {
      program_.addInput("group", DataType::F64);
      groupDeclared_ = true;
    }
    return arena().field("group");
  }

  /// That the dimensions a filter names are those of its input.
  template <typename Map>
  void matchDimensions(const YAML::Node& node, const Map& named, const Values& in) const {
    std::set<std::string> dimensions;
    for (const auto& [dimension, value] : named) {
      dimensions.insert(dimension);
    }
    if (dimensions != keysOf(in)) {
      fail(node,
           "The input parameters of a component ({" + listOf(dimensions) +
               "}) added to a composite are incompatible (should be {" + listOf(keysOf(in)) +
               "}).");
    }
  }

  /// The parameters `wanted` that the component gives at `in`, for the points that reach it
  /// (`reach`, NoNode for every point).
  Values evaluate(const YAML::Node& node,
                  const Values& in,
                  const std::set<std::string>& wanted,
                  NodeId reach) {
    const std::string tag = node.Tag();
    if (tag == "!Include") {
      // easi reads an included file with the dimensions of the root
      const auto [root, path] = included(node);
      if (keysOf(in) != rootDimensions_) {
        fail(node,
             "The input parameters of a component ({" + listOf(rootDimensions_) +
                 "}) added to a composite are incompatible (should be {" + listOf(keysOf(in)) +
                 "}).");
      }
      files_.push_back(path);
      Values values = evaluate(root, in, wanted, reach);
      files_.pop_back();
      return values;
    }
    if (tag == "!Switch") {
      return switched(node, in, wanted, reach);
    }
    if (tag == "!LayeredModel") {
      return layers(node, in, wanted, reach);
    }
    if (tag == "!Any" || tag == "!IdentityMap" || tag == "!GroupFilter" ||
        tag == "!AxisAlignedCuboidalDomainFilter" || tag == "!SphericalDomainFilter") {
      return composite(node, in, wanted, reach);
    }
    if (tag == "!ConstantMap" || tag == "!AffineMap" || tag == "!PolynomialMap" ||
        tag == "!EvalModel" || tag == "!ASAGI" || tag == "!FunctionMap" || tag == "!LuaMap") {
      return composite(node, mapped(node, in, reach), wanted, reach);
    }
    if (tag == "!SCECFile" || tag == "!STRESS_STR_DIP_SLIP_AM" || tag == "!OptimalStress" ||
        tag == "!AndersonianStress") {
      unsupported(node);
    }
    fail(node, "Unknown tag " + tag);
  }

  /// A !Switch: each parameter from the component that names it.
  Values switched(const YAML::Node& node,
                  const Values& in,
                  const std::set<std::string>& wanted,
                  NodeId reach) {
    std::set<std::string> all;
    std::vector<std::pair<std::set<std::string>, YAML::Node>> branches;
    for (const auto& entry : node) {
      std::set<std::string> restrictions;
      try {
        if (entry.first.IsSequence()) {
          const auto names = entry.first.as<std::vector<std::string>>();
          restrictions.insert(names.begin(), names.end());
        } else {
          restrictions.insert(entry.first.as<std::string>());
        }
      } catch (const YAML::Exception& error) {
        fail(node, std::string("a branch of a switch names parameters: ") + error.what());
      }
      std::set<std::string> overlap;
      std::set_intersection(restrictions.begin(),
                            restrictions.end(),
                            all.begin(),
                            all.end(),
                            std::inserter(overlap, overlap.end()));
      if (!overlap.empty()) {
        fail(node,
             "Switch requires a partition of parameters. The parameters {" + listOf(overlap) +
                 "} were specified in multiple branches.");
      }
      all.insert(restrictions.begin(), restrictions.end());
      branches.emplace_back(std::move(restrictions), entry.second);
    }
    Values out;
    std::set<std::string> missing;
    for (const auto& parameter : wanted) {
      if (all.count(parameter) == 0) {
        missing.insert(parameter);
      }
    }
    if (!missing.empty()) {
      require(reach,
              konst(0.0),
              where(node) + ": Switch is not complete with respect to request. Missing: {" +
                  listOf(missing) + "}.");
    }
    for (const auto& [restrictions, child] : branches) {
      std::set<std::string> sub;
      std::set_intersection(restrictions.begin(),
                            restrictions.end(),
                            wanted.begin(),
                            wanted.end(),
                            std::inserter(sub, sub.end()));
      if (sub.empty()) {
        continue;
      }
      for (auto& [name, value] : evaluate(child, in, sub, reach)) {
        out[name] = std::move(value);
      }
    }
    return out;
  }

  /// A composite of the values `values` (its input mapped): its own where it has no components,
  /// else those of the first of its components that accepts the point.
  Values composite(const YAML::Node& node,
                   const Values& values,
                   const std::set<std::string>& wanted,
                   NodeId reach) {
    const auto children = childrenOf(node);
    if (children.empty()) {
      Values out;
      std::set<std::string> missing;
      for (const auto& parameter : wanted) {
        const auto found = values.find(parameter);
        if (found == values.end()) {
          missing.insert(parameter);
        } else {
          out[parameter] = found->second;
        }
      }
      if (!missing.empty()) {
        require(reach,
                konst(0.0),
                where(node) + ": Component does not supply all required parameters. Missing: {" +
                    listOf(missing) + "}.");
      }
      return out;
    }
    const std::string message = "Could not find model for point (" + where(node) + ").";
    // the first component that accepts a point takes it; the later ones see the rest
    std::vector<NodeId> accepted;
    std::vector<Values> results;
    NodeId before = NoNode; // that no component before accepts
    NodeId any = konst(0.0);
    bool always = false;
    for (const auto& child : children) {
      const NodeId accepts = this->accepts(child, values);
      const NodeId reaches = both(reach, both(before, accepts));
      accepted.push_back(accepts);
      results.push_back(evaluate(child, values, wanted, reaches));
      if (accepts == NoNode) {
        always = true;
        break;
      }
      any = pw(Fn::Or, any, accepts);
      before = both(before, negation(accepts));
    }
    if (!always) {
      require(reach, any, message);
    }
    // from the last component back: a select of each that accepts, before the rest, a value it
    // does not give NaN (where a check fails)
    const NodeId nan = konst(std::numeric_limits<double>::quiet_NaN());
    Values out;
    for (const auto& parameter : wanted) {
      NodeId node = nan;
      std::vector<expr::SourcePosition> from;
      bool given = false;
      for (std::size_t c = results.size(); c-- > 0;) {
        const auto found = results[c].find(parameter);
        const bool gives = found != results[c].end() && found->second.node != NoNode;
        const NodeId here = gives ? found->second.node : nan;
        if (gives) {
          given = true;
          for (const auto& position : found->second.from) {
            if (std::find(from.begin(), from.end(), position) == from.end()) {
              from.push_back(position);
            }
          }
        }
        node = accepted[c] == NoNode ? here : pw(Fn::Select, accepted[c], here, node);
      }
      if (given) {
        out[parameter] = Value{node, std::move(from)};
      }
    }
    return out;
  }

  /// The parameters a component supplies at `in` (easi's suppliedParameters()): those of its
  /// components, or its own.
  std::set<std::string> suppliedBy(const YAML::Node& node, const std::set<std::string>& in) {
    const std::string tag = node.Tag();
    if (tag == "!Include") {
      const auto [root, path] = included(node);
      files_.push_back(path);
      auto supplied = suppliedBy(root, rootDimensions_);
      files_.pop_back();
      return supplied;
    }
    if (tag == "!Switch") {
      std::set<std::string> supplied;
      for (const auto& entry : node) {
        const auto sub = suppliedBy(entry.second, in);
        supplied.insert(sub.begin(), sub.end());
      }
      return supplied;
    }
    if (tag == "!LayeredModel") {
      const auto parameters = as<std::vector<std::string>>(node, "parameters");
      return {parameters.begin(), parameters.end()};
    }
    // the outputs of the component itself
    std::set<std::string> own;
    if (tag == "!Any" || tag == "!IdentityMap" || tag == "!GroupFilter") {
      own = in;
    } else if (tag == "!AxisAlignedCuboidalDomainFilter" || tag == "!SphericalDomainFilter") {
      own = in;
    } else if (tag == "!ConstantMap" || tag == "!PolynomialMap" || tag == "!FunctionMap") {
      for (const auto& [name, value] : as<std::map<std::string, YAML::Node>>(node, "map")) {
        own.insert(name);
      }
    } else if (tag == "!AffineMap") {
      for (const auto& [name, value] : as<std::map<std::string, YAML::Node>>(node, "matrix")) {
        own.insert(name);
      }
    } else if (tag == "!EvalModel" || tag == "!ASAGI") {
      const auto names = as<std::vector<std::string>>(node, "parameters");
      own.insert(names.begin(), names.end());
    } else if (tag == "!LuaMap") {
      const YAML::Node returns = node["returns"];
      if (returns && returns.IsScalar()) {
        own.insert(returns.as<std::string>());
      } else if (returns) {
        const auto names = returns.as<std::vector<std::string>>();
        own.insert(names.begin(), names.end());
      }
    } else {
      unsupported(node);
    }
    const auto children = childrenOf(node);
    if (children.empty()) {
      return own;
    }
    std::set<std::string> supplied;
    for (const auto& child : children) {
      const auto sub = suppliedBy(child, own);
      supplied.insert(sub.begin(), sub.end());
    }
    return supplied;
  }

  EasiQuery query_;
  std::set<std::string> rootDimensions_;
  expr::Program program_;
  std::vector<std::string> files_;
  std::vector<std::pair<NodeId, std::string>> checks_;
  bool groupDeclared_{false};
};

} // namespace

expr::Program easiProgram(const std::string& path, const EasiQuery& query) {
  return Walker(query).translate(path);
}

std::vector<std::string> easiParameters(const std::string& path,
                                        const std::vector<std::string>& dimensions) {
  EasiQuery query;
  query.dimensions = dimensions;
  const auto parameters = Walker(query).parameters(path);
  return {parameters.begin(), parameters.end()};
}

EasiProgramReader::EasiProgramReader(std::string path,
                                     const std::vector<std::string>& dimensions,
                                     ReferenceMaker makeReference)
    : path_(std::move(path)), makeReference_(makeReference) {
  // easi orders the dimensions by name
  const std::set<std::string> sorted(dimensions.begin(), dimensions.end());
  dimensions_.assign(sorted.begin(), sorted.end());
  parameters_ = easiParameters(path_, dimensions_);
}

EasiProgramReader::~EasiProgramReader() = default;

void EasiProgramReader::prepare(const DataTable& table) {
  std::vector<std::string> wanted;
  bool group = false;
  for (const auto& entry : table.dataEntries()) {
    if (entry.direction != Direction::In) {
      wanted.push_back(entry.name);
    } else if (entry.name == "group") {
      group = true;
    }
  }
  if (compiled_ == nullptr || wanted != madeFor_ || group != madeWithGroup_) {
    EasiQuery query;
    query.dimensions = dimensions_;
    query.parameters = wanted;
    query.group = group;
    compiled_ = std::make_unique<CompiledReader>(
        easiProgram(path_, query),
        makeReference_ != nullptr ? makeReference_(path_, dimensions_) : nullptr);
    madeFor_ = std::move(wanted);
    madeWithGroup_ = group;
  }
  compiled_->prepare(table);
}

void EasiProgramReader::call(const DataTable& table) {
  if (compiled_ == nullptr) {
    prepare(table);
  }
  compiled_->call(table);
}

} // namespace seissol::reader::scripting
