// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/Core.h"

#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

std::string located(const CorePosition& position) {
  return (position.file.empty() ? std::string("<core>") : position.file) + ":" +
         std::to_string(position.line) + ":" + std::to_string(position.column);
}

/// How deep parentheses, calls and negations may nest: the loader rejects deeper nesting with an
/// error rather than running out of stack.
constexpr int MaxNesting = 1024;

/// The reserved words of the core format; `inf` and `nan` are numbers.
const std::set<std::string>& coreWords() {
  static const std::set<std::string> Words = {
      "meta", "in", "grid", "def", "state", "next", "out", "check", "inf", "nan"};
  return Words;
}

/// The words the full language reserves besides: a printer escapes them as well, so that a core
/// file reads there too.
const std::set<std::string>& fullWords() {
  static const std::set<std::string> Words = {
      "let", "if", "then", "else", "and", "or", "not", "import"};
  return Words;
}

/// The operations of the program that are called by name, as the table of the specification
/// spells them; +, -, * and / and the negation are written as operators.
const std::vector<std::pair<std::string, Fn>>& operationSpellings() {
  static const std::vector<std::pair<std::string, Fn>> Spellings = {
      {"rcp", Fn::Rcp},     {"sqrt", Fn::Sqrt},   {"abs", Fn::Abs},      {"floor", Fn::Floor},
      {"ceil", Fn::Ceil},   {"round", Fn::Round}, {"sign", Fn::Sign},    {"exp", Fn::Exp},
      {"log", Fn::Log},     {"log2", Fn::Log2},   {"log10", Fn::Log10},  {"sin", Fn::Sin},
      {"cos", Fn::Cos},     {"tan", Fn::Tan},     {"asin", Fn::Asin},    {"acos", Fn::Acos},
      {"atan", Fn::Atan},   {"sinh", Fn::Sinh},   {"cosh", Fn::Cosh},    {"tanh", Fn::Tanh},
      {"asinh", Fn::Asinh}, {"acosh", Fn::Acosh}, {"atanh", Fn::Atanh},  {"erf", Fn::Erf},
      {"pow", Fn::Pow},     {"atan2", Fn::Atan2}, {"min", Fn::Min},      {"max", Fn::Max},
      {"mod", Fn::Mod},     {"lt", Fn::Lt},       {"le", Fn::Le},        {"eq", Fn::Eq},
      {"land", Fn::And},    {"lor", Fn::Or},      {"select", Fn::Select}};
  return Spellings;
}

/// The conversions between the types of values, spelled as in the full language.
const std::vector<std::pair<std::string, ValueType>>& conversionSpellings() {
  static const std::vector<std::pair<std::string, ValueType>> Spellings = {
      {"f64", ValueType::F64}, {"f32", ValueType::F32}, {"real", ValueType::Real}};
  return Spellings;
}

bool asciiLetter(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }
bool asciiDigit(char c) { return c >= '0' && c <= '9'; }
bool hexDigit(char c) { return asciiDigit(c) || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'); }
bool nameStart(char c) { return asciiLetter(c) || c == '_'; }
bool nameContinues(char c) { return asciiLetter(c) || asciiDigit(c) || c == '_'; }

/// Whether a number that from_chars finds beyond the range of binary64 overflows rather than
/// underflows: its magnitude is the base to the power of the position of its first digit that is
/// not zero, against the point, plus its exponent -- far from 0 either way.
bool overflows(const std::string& body, bool hex) {
  const auto exponentAt = body.find_first_of(hex ? "pP" : "eE");
  long long exponent = 0;
  if (exponentAt != std::string::npos) {
    std::size_t i = exponentAt + 1;
    bool negative = false;
    if (i < body.size() && (body[i] == '+' || body[i] == '-')) {
      negative = body[i] == '-';
      ++i;
    }
    for (; i < body.size(); ++i) {
      exponent = std::min(exponent * 10 + (body[i] - '0'), 1LL << 40);
    }
    exponent = negative ? -exponent : exponent;
  }
  const std::string mantissa = body.substr(0, exponentAt);
  const auto point = mantissa.find('.');
  const auto integerDigits =
      static_cast<long long>(point == std::string::npos ? mantissa.size() : point);
  long long digit = 0;
  for (const char c : mantissa) {
    if (c == '.') {
      continue;
    }
    if (c != '0') {
      break;
    }
    ++digit;
  }
  return (integerDigits - 1 - digit) * (hex ? 4 : 1) + exponent > 0;
}

/// The value of the spelling of a number, with correct rounding; beyond the range of binary64, an
/// infinity or zero, as the specification has it.
double numberOf(const std::string& body, bool hex) {
  double value = 0.0;
  const auto result = std::from_chars(body.data(),
                                      body.data() + body.size(),
                                      value,
                                      hex ? std::chars_format::hex : std::chars_format::general);
  if (result.ec == std::errc::result_out_of_range) {
    return overflows(body, hex) ? std::numeric_limits<double>::infinity() : 0.0;
  }
  if (result.ec != std::errc() || result.ptr != body.data() + body.size()) {
    throw std::logic_error("expr: the core lexer let through the number " + body);
  }
  return value;
}

/// The shortest decimal spelling that reads back as the same value (as std::to_chars writes it),
/// and inf, -inf and nan. The sign and the payload of a NaN are not kept.
std::string spellNumber(double value) {
  if (std::isnan(value)) {
    return "nan";
  }
  if (std::isinf(value)) {
    return value < 0 ? "-inf" : "inf";
  }
  std::array<char, 64> buffer{};
  const auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value);
  return std::string(buffer.data(), result.ptr);
}

std::string spellString(const std::string& text) {
  std::string out = "\"";
  for (const char c : text) {
    if (c == '"' || c == '\\') {
      out.push_back('\\');
    }
    if (c == '\n' || c == '\r') {
      throw std::invalid_argument("expr: the core format has no strings across lines: " + text);
    }
    out.push_back(c);
  }
  out.push_back('"');
  return out;
}

/// A name as the core format writes it: plain where it is a plain name and no reserved word,
/// escaped in backquotes otherwise.
std::string spellName(const std::string& name) {
  const bool plain = !name.empty() && nameStart(name.front()) &&
                     std::all_of(name.begin(), name.end(), nameContinues) &&
                     coreWords().count(name) == 0 && fullWords().count(name) == 0;
  if (plain) {
    return name;
  }
  if (name.empty() || name.find_first_of("`\n\r") != std::string::npos) {
    throw std::invalid_argument("expr: the core format cannot write the name `" + name + "`");
  }
  return "`" + name + "`";
}

const char* spellKind(reader::datafield::GridKind kind) {
  switch (kind) {
  case reader::datafield::GridKind::AsagiLite:
    return "asagi";
  case reader::datafield::GridKind::Scec:
    return "scec";
  case reader::datafield::GridKind::Distributed:
    break;
  }
  throw std::invalid_argument("expr: the core format has no kind of grid for a distributed one");
}

const char* spellInterpolation(reader::datafield::Interpolation interpolation) {
  switch (interpolation) {
  case reader::datafield::Interpolation::Nearest:
    return "nearest";
  case reader::datafield::Interpolation::Linear:
    return "linear";
  case reader::datafield::Interpolation::Cubic:
    return "cubic";
  }
  return "linear";
}

const char* spellBoundary(reader::datafield::BoundaryMode boundary) {
  return boundary == reader::datafield::BoundaryMode::Fail ? "fail" : "clamp";
}

// ================================================================ lexer =====

enum class Tok : std::uint8_t {
  Name,
  Word, // a reserved word
  Number,
  String,
  Plus,
  Minus,
  Star,
  Slash,
  Open,
  Close,
  Comma,
  Equals,
  Dot,
  End
};

struct Token {
  Tok kind{Tok::End};
  std::string text; // a name, a word, the content of a string, the spelling of a number
  double number{0.0};
  bool escaped{false}; // a name in backquotes
  bool whole{false};   // a number of decimal digits only
  int line{1};
  int column{1};
};

class Lexer {
  public:
  Lexer(const std::string& text, std::string file) : text_(text), file_(std::move(file)) {}

  /// Reads the header and gives the version it names.
  int header() {
    if (text_.compare(0, 3, "\xEF\xBB\xBF") == 0) {
      fail(1, 1, "a core file has no byte order mark");
    }
    const std::string magic = "#!sderiv";
    if (text_.compare(0, magic.size(), magic) != 0) {
      fail(1, 1, "a core file starts with the header `#!sderiv 1`");
    }
    for (std::size_t i = 0; i < magic.size(); ++i) {
      advance();
    }
    if (peek() != ' ' && peek() != '\t') {
      fail(line_, column_, "the header names the version of the core format: `#!sderiv 1`");
    }
    while (peek() == ' ' || peek() == '\t') {
      advance();
    }
    const int line = line_;
    const int column = column_;
    std::string digits;
    while (asciiDigit(peek()) && digits.size() < 10) {
      digits.push_back(peek());
      advance();
    }
    while (peek() == ' ' || peek() == '\t' || peek() == '\r') {
      advance();
    }
    if (digits.empty() || (pos_ < text_.size() && peek() != '\n')) {
      fail(line, column, "the header names the version as a whole number: `#!sderiv 1`");
    }
    const long long version = std::stoll(digits);
    if (version < 1) {
      fail(line, column, "there is no version 0 of the core format");
    }
    if (version > CoreVersion) {
      fail(line,
           column,
           "the file needs version " + digits + " of the core format; this SeisSol reads up to " +
               std::to_string(CoreVersion));
    }
    return static_cast<int>(version);
  }

  std::vector<Token> tokens() {
    std::vector<Token> out;
    while (true) {
      skipSpace();
      Token token;
      token.line = line_;
      token.column = column_;
      if (pos_ >= text_.size()) {
        out.push_back(token);
        return out;
      }
      const char c = peek();
      if (nameStart(c)) {
        const std::size_t start = pos_;
        while (nameContinues(peek())) {
          advance();
        }
        token.text = text_.substr(start, pos_ - start);
        token.kind = coreWords().count(token.text) != 0 ? Tok::Word : Tok::Name;
        if (token.text == "inf" || token.text == "nan") {
          token.kind = Tok::Number;
          token.number = token.text == "inf" ? std::numeric_limits<double>::infinity()
                                             : std::numeric_limits<double>::quiet_NaN();
        }
      } else if (c == '`') {
        advance();
        const std::size_t start = pos_;
        while (pos_ < text_.size() && peek() != '`' && peek() != '\n' && peek() != '\r') {
          advance();
        }
        if (peek() != '`') {
          fail(token.line, token.column, "an escaped name ends with ` on its line");
        }
        token.text = text_.substr(start, pos_ - start);
        advance();
        if (token.text.empty()) {
          fail(token.line, token.column, "an escaped name is not empty");
        }
        token.kind = Tok::Name;
        token.escaped = true;
      } else if (asciiDigit(c) || (c == '.' && asciiDigit(peek(1)))) {
        number(token);
      } else if (c == '"') {
        string(token);
      } else {
        advance();
        switch (c) {
        case '+':
          token.kind = Tok::Plus;
          break;
        case '-':
          token.kind = Tok::Minus;
          break;
        case '*':
          token.kind = Tok::Star;
          break;
        case '/':
          token.kind = Tok::Slash;
          break;
        case '(':
          token.kind = Tok::Open;
          break;
        case ')':
          token.kind = Tok::Close;
          break;
        case ',':
          token.kind = Tok::Comma;
          break;
        case '=':
          token.kind = Tok::Equals;
          break;
        case '.':
          token.kind = Tok::Dot;
          break;
        default:
          fail(token.line,
               token.column,
               std::string("unexpected character `") + c +
                   "`; the operators are + - * / ( ) , = .");
        }
      }
      out.push_back(std::move(token));
    }
  }

  [[noreturn]] void fail(int line, int column, const std::string& what) const {
    throw CoreError(CorePosition{file_, line, column}, what);
  }

  private:
  [[nodiscard]] char peek(std::size_t ahead = 0) const {
    return pos_ + ahead < text_.size() ? text_[pos_ + ahead] : '\0';
  }

  /// The column counts characters: a byte that continues one in UTF-8 does not move it.
  void advance() {
    const char c = text_[pos_++];
    if (c == '\n') {
      ++line_;
      column_ = 1;
    } else if ((static_cast<unsigned char>(c) & 0xC0U) != 0x80U) {
      ++column_;
    }
  }

  void skipSpace() {
    while (pos_ < text_.size()) {
      const char c = peek();
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
        advance();
      } else if (c == '#') {
        while (pos_ < text_.size() && peek() != '\n') {
          advance();
        }
      } else {
        return;
      }
    }
  }

  void digits(bool hex) {
    while (hex ? hexDigit(peek()) : asciiDigit(peek())) {
      advance();
    }
  }

  void exponent(const Token& token) {
    advance();
    if (peek() == '+' || peek() == '-') {
      advance();
    }
    if (!asciiDigit(peek())) {
      fail(token.line, token.column, "the exponent of a number has digits");
    }
    digits(false);
  }

  void number(Token& token) {
    const std::size_t start = pos_;
    token.kind = Tok::Number;
    if (peek() == '0' && (peek(1) == 'x' || peek(1) == 'X')) {
      advance();
      advance();
      const std::size_t body = pos_;
      digits(true);
      bool any = pos_ > body;
      if (peek() == '.') {
        advance();
        const std::size_t fraction = pos_;
        digits(true);
        any = any || pos_ > fraction;
      }
      if (!any) {
        fail(token.line, token.column, "a hexadecimal number has digits");
      }
      if (peek() != 'p' && peek() != 'P') {
        fail(
            token.line, token.column, "a hexadecimal number has a binary exponent, as in 0x1.8p+1");
      }
      exponent(token);
      token.number = numberOf(text_.substr(body, pos_ - body), true);
    } else {
      token.whole = true;
      digits(false);
      if (peek() == '.') {
        token.whole = false;
        advance();
        digits(false);
      }
      if (peek() == 'e' || peek() == 'E') {
        token.whole = false;
        exponent(token);
      }
      token.number = numberOf(text_.substr(start, pos_ - start), false);
    }
    token.text = text_.substr(start, pos_ - start);
    if (nameStart(peek()) || peek() == '`' || peek() == '.') {
      fail(line_, column_, "a number is followed by `" + std::string(1, peek()) + "`");
    }
  }

  void string(Token& token) {
    advance();
    token.kind = Tok::String;
    while (true) {
      if (pos_ >= text_.size() || peek() == '\n' || peek() == '\r') {
        fail(token.line, token.column, "a string ends with \" on its line");
      }
      const char c = peek();
      if (c == '"') {
        advance();
        return;
      }
      if (c == '\\' && (peek(1) == '"' || peek(1) == '\\')) {
        advance();
      }
      token.text.push_back(peek());
      advance();
    }
  }

  const std::string& text_;
  std::string file_;
  std::size_t pos_{0};
  int line_{1};
  int column_{1};
};

// =============================================================== loader =====

class Loader {
  public:
  Loader(std::vector<Token> tokens, std::string file)
      : tokens_(std::move(tokens)), file_(std::move(file)) {}

  Program run() {
    bool metaAllowed = true;
    while (peek().kind != Tok::End) {
      const Token& at = peek();
      if (at.kind != Tok::Word || at.text == "inf" || at.text == "nan") {
        fail(at, "a declaration starts with meta, in, grid, def, state, next, out or check");
      }
      if (at.text == "meta") {
        if (!metaAllowed) {
          fail(at, "meta lines come directly after the header");
        }
        meta();
        continue;
      }
      metaAllowed = false;
      if (at.text == "in") {
        input();
      } else if (at.text == "grid") {
        grid();
      } else if (at.text == "def") {
        definition();
      } else if (at.text == "state") {
        state();
      } else if (at.text == "next") {
        nextValue();
      } else if (at.text == "out") {
        output();
      } else {
        check();
      }
    }
    for (const auto& state : states_) {
      if (state.next == NoNode) {
        fail(state.at,
             "the state `" + state.name + "` has no next value; give it with `next " +
                 spellName(state.name) + " = ...`");
      }
      program_.addState(state.name, state.initial, state.next, state.type);
    }
    for (const auto& [name, root] : outputs_) {
      program_.addOutput(name, reader::scripting::DataType::F64, root);
    }
    for (auto& [root, message] : checks_) {
      program_.addCheck(root, std::move(message));
    }
    try {
      validate(program_);
    } catch (const std::invalid_argument& error) {
      throw CoreError(CorePosition{file_, 1, 1}, error.what());
    }
    return std::move(program_);
  }

  private:
  enum class What : std::uint8_t { Input, State, Definition, Grid };
  struct Entry {
    What what{What::Input};
    NodeId node{NoNode};
    std::size_t index{0}; // of a grid or a state
  };
  struct StateDeclaration {
    std::string name;
    double initial{0.0};
    ValueType type{ValueType::Real};
    NodeId next{NoNode};
    Token at;
  };
  struct GridDeclaration {
    GridId id{NoGrid};
    std::vector<std::string> components;
  };

  [[nodiscard]] Arena& arena() { return program_.arena(); }
  [[nodiscard]] const Token& peek() const { return tokens_[pos_]; }
  const Token& next() { return tokens_[pos_++]; }
  bool accept(Tok kind) {
    if (peek().kind == kind) {
      ++pos_;
      return true;
    }
    return false;
  }

  [[noreturn]] void fail(const Token& at, const std::string& what) const {
    throw CoreError(CorePosition{file_, at.line, at.column}, what);
  }

  const Token& expect(Tok kind, const std::string& what) {
    if (peek().kind != kind) {
      fail(peek(), "expected " + what);
    }
    return next();
  }

  /// A name, plain or escaped.
  std::string name(const std::string& what) {
    if (peek().kind != Tok::Name) {
      fail(peek(),
           "expected " + what +
               (peek().kind == Tok::Word
                    ? "; `" + peek().text + "` is a reserved word, which an escaped name, `" +
                          peek().text + "`, avoids"
                    : ""));
    }
    return next().text;
  }

  /// Declares `name` in the namespace of inputs, states, definitions and grids.
  void declare(const Token& at, const std::string& name, Entry entry) {
    if (names_.count(name) != 0) {
      fail(at, "`" + name + "` is declared twice");
    }
    names_.emplace(name, entry);
  }

  void meta() {
    next();
    MetaEntry entry;
    entry.key = name("the key of the meta line");
    expect(Tok::Equals, "`=` and the values of the meta line");
    do {
      const Token& value = peek();
      if (value.kind == Tok::String) {
        entry.values.push_back(MetaValue{next().text, true});
      } else if (value.kind == Tok::Number && value.whole) {
        entry.values.push_back(MetaValue{next().text, false});
      } else {
        fail(value, "the value of a meta line is a string or a whole number");
      }
    } while (accept(Tok::Comma));
    program_.addMeta(std::move(entry));
  }

  void input() {
    next();
    do {
      const Token at = peek();
      const std::string input = name("the name of an input");
      program_.addInput(input, reader::scripting::DataType::F64);
      declare(at, input, Entry{What::Input, arena().field(input), 0});
    } while (accept(Tok::Comma));
  }

  void grid() {
    next();
    const Token at = peek();
    const std::string gridName = name("the name of the grid");
    expect(Tok::Open, "`(` and the attributes of the grid");
    std::map<std::string, std::string> strings;
    std::vector<std::string> components;
    std::optional<int> time;
    std::set<std::string> given;
    do {
      const Token key = peek();
      const std::string attribute = name("an attribute of the grid");
      if (!given.insert(attribute).second) {
        fail(key, "the attribute `" + attribute + "` is given twice");
      }
      expect(Tok::Equals, "`=` and the value of the attribute");
      const Token value = peek();
      if (attribute == "components") {
        if (accept(Tok::Open)) {
          do {
            components.push_back(expect(Tok::String, "the name of a component, a string").text);
          } while (accept(Tok::Comma));
          expect(Tok::Close, "`)` after the components");
        } else {
          components.push_back(
              expect(Tok::String, "the names of the components, as (\"a\", \"b\")").text);
        }
      } else if (attribute == "time") {
        if (value.kind != Tok::Number || !value.whole) {
          fail(value, "the time axis of a grid is a whole number");
        }
        next();
        if (value.number >= reader::datafield::MaxGridDimensions) {
          fail(value,
               "a grid has at most " + std::to_string(reader::datafield::MaxGridDimensions) +
                   " axes");
        }
        time = static_cast<int>(value.number);
      } else if (attribute == "kind" || attribute == "path" || attribute == "variable" ||
                 attribute == "interpolation" || attribute == "boundary") {
        strings[attribute] = expect(Tok::String, "a string").text;
      } else {
        fail(key,
             "a grid has no attribute `" + attribute +
                 "`; it has kind, path, variable, components, interpolation, boundary and time");
      }
    } while (accept(Tok::Comma));
    expect(Tok::Close, "`)` after the attributes of the grid");

    for (const char* required : {"kind", "path", "variable", "interpolation", "boundary"}) {
      if (strings.count(required) == 0) {
        fail(at,
             std::string("the grid `") + gridName + "` gives no " + required +
                 "; the core format has no defaults");
      }
    }
    if (components.empty()) {
      fail(at, "the grid `" + gridName + "` names no components");
    }
    std::set<std::string> distinct(components.begin(), components.end());
    if (distinct.size() != components.size()) {
      fail(at, "the grid `" + gridName + "` names a component twice");
    }
    reader::datafield::GridDesc desc;
    const auto kind = reader::datafield::parseGridKind(strings["kind"]);
    const auto interpolation = reader::datafield::parseInterpolation(strings["interpolation"]);
    const auto boundary = reader::datafield::parseBoundaryMode(strings["boundary"]);
    if (!kind.has_value()) {
      fail(at, "the kind of a grid is \"asagi\" or \"scec\", not \"" + strings["kind"] + "\"");
    }
    if (!interpolation.has_value()) {
      fail(at,
           "the interpolation of a grid is \"nearest\", \"linear\" or \"cubic\", not \"" +
               strings["interpolation"] + "\"");
    }
    if (!boundary.has_value()) {
      fail(at,
           "the boundary of a grid is \"clamp\" or \"fail\", not \"" + strings["boundary"] + "\"");
    }
    desc.kind = *kind;
    desc.interpolation = *interpolation;
    desc.boundary = *boundary;
    desc.variable = strings["variable"];
    desc.timeAxis = time;
    // relative to the directory of the core file
    std::filesystem::path path(strings["path"]);
    if (path.is_relative() && !file_.empty()) {
      path = (std::filesystem::path(file_).parent_path() / path).lexically_normal();
    }
    desc.path = path.string();
    const GridId id = program_.internGrid(desc);
    program_.nameGrid(id, GridNames{gridName, components});
    declare(at, gridName, Entry{What::Grid, NoNode, grids_.size()});
    grids_.push_back(GridDeclaration{id, components});
  }

  void definition() {
    next();
    const Token at = peek();
    const std::string definitionName = name("the name of the definition");
    expect(Tok::Equals, "`=` and the value of the definition");
    const NodeId value = expression(0);
    declare(at, definitionName, Entry{What::Definition, value, 0});
  }

  void state() {
    next();
    const Token at = peek();
    std::string stateName = name("the name of the state");
    ValueType type = ValueType::Real;
    if (peek().kind == Tok::Name) {
      // `state f64 s = 0`: the type of the state, then its name
      const auto& conversions = conversionSpellings();
      const auto found =
          std::find_if(conversions.begin(), conversions.end(), [&](const auto& conversion) {
            return conversion.first == stateName;
          });
      if (at.escaped || found == conversions.end()) {
        fail(at, "the type of a state is f64, f32 or real, not `" + stateName + "`");
      }
      type = found->second;
      stateName = name("the name of the state");
    }
    expect(Tok::Equals, "`=` and the initial value of the state, a number");
    const bool negative = accept(Tok::Minus);
    const Token& value = expect(Tok::Number, "the initial value of the state, a number");
    declare(at, stateName, Entry{What::State, arena().field(stateName), states_.size()});
    states_.push_back(
        StateDeclaration{stateName, negative ? -value.number : value.number, type, NoNode, at});
  }

  void nextValue() {
    next();
    const Token at = peek();
    const std::string stateName = name("the name of a state");
    const auto found = names_.find(stateName);
    if (found == names_.end() || found->second.what != What::State) {
      fail(at, "`next` gives the next value of a state, and `" + stateName + "` is none");
    }
    StateDeclaration& state = states_[found->second.index];
    if (state.next != NoNode) {
      fail(at, "the state `" + stateName + "` has its next value already");
    }
    expect(Tok::Equals, "`=` and the next value");
    state.next = expression(0);
  }

  void output() {
    next();
    const Token at = peek();
    const std::string outputName = name("the name of the output");
    if (std::any_of(outputs_.begin(), outputs_.end(), [&](const auto& output) {
          return output.first == outputName;
        })) {
      fail(at, "the output `" + outputName + "` is given twice");
    }
    expect(Tok::Equals, "`=` and the value of the output");
    outputs_.emplace_back(outputName, expression(0));
  }

  void check() {
    next();
    const NodeId condition = expression(0);
    expect(Tok::Comma, "`,` and the message of the check, a string");
    checks_.emplace_back(condition, expect(Tok::String, "the message of the check, a string").text);
  }

  // ------------------------------------------------------------ expressions

  void nest(const Token& at, int depth) const {
    if (depth > MaxNesting) {
      fail(at, "the expression nests deeper than " + std::to_string(MaxNesting) + " levels");
    }
  }

  NodeId expression(int depth) {
    nest(peek(), depth);
    NodeId left = term(depth);
    while (peek().kind == Tok::Plus || peek().kind == Tok::Minus) {
      const Fn fn = next().kind == Tok::Plus ? Fn::Add : Fn::Sub;
      const NodeId right = term(depth);
      left = arena().pw(fn, left, right);
    }
    return left;
  }

  NodeId term(int depth) {
    NodeId left = factor(depth);
    while (peek().kind == Tok::Star || peek().kind == Tok::Slash) {
      const Fn fn = next().kind == Tok::Star ? Fn::Mul : Fn::Div;
      const NodeId right = factor(depth);
      left = arena().pw(fn, left, right);
    }
    return left;
  }

  NodeId factor(int depth) {
    if (peek().kind == Tok::Minus) {
      nest(next(), depth + 1);
      const NodeId operand = factor(depth + 1);
      const Node& node = arena()[operand];
      // negation applied to a number gives a number, since negation is exact
      if (node.kind == Kind::Const) {
        return arena().konst(-node.value, node.type);
      }
      return arena().pw(Fn::Neg, operand);
    }
    return primary(depth);
  }

  NodeId primary(int depth) {
    const Token at = peek();
    switch (at.kind) {
    case Tok::Number:
      next();
      return arena().konst(at.number);
    case Tok::Open: {
      next();
      const NodeId inner = expression(depth + 1);
      expect(Tok::Close, "`)`");
      return inner;
    }
    case Tok::Name:
      next();
      if (peek().kind == Tok::Open) {
        if (at.escaped) {
          fail(at, "an escaped name is no operation");
        }
        return call(at, depth);
      }
      if (peek().kind == Tok::Dot) {
        return sample(at, depth);
      }
      return reference(at);
    case Tok::Word:
      fail(at, "`" + at.text + "` is a reserved word; an escaped name, `" + at.text + "`, is not");
    default:
      fail(at, "expected a value: a number, a name, an operation or `(`");
    }
  }

  std::vector<NodeId> arguments(int depth) {
    expect(Tok::Open, "`(`");
    std::vector<NodeId> args;
    if (peek().kind != Tok::Close) {
      do {
        args.push_back(expression(depth + 1));
      } while (accept(Tok::Comma));
    }
    expect(Tok::Close, "`)` after the arguments");
    return args;
  }

  NodeId call(const Token& at, int depth) {
    const auto& operations = operationSpellings();
    const auto operation =
        std::find_if(operations.begin(), operations.end(), [&](const auto& spelling) {
          return spelling.first == at.text;
        });
    const auto& conversions = conversionSpellings();
    const auto conversion =
        std::find_if(conversions.begin(), conversions.end(), [&](const auto& spelling) {
          return spelling.first == at.text;
        });
    if (operation == operations.end() && conversion == conversions.end()) {
      fail(at,
           "there is no operation `" + at.text +
               "`; a definition with parameters exists only in the full language");
    }
    const std::vector<NodeId> args = arguments(depth);
    const int expected = operation != operations.end() ? arity(operation->second) : 1;
    if (static_cast<int>(args.size()) != expected) {
      fail(at,
           "`" + at.text + "` takes " + std::to_string(expected) + " argument" +
               (expected == 1 ? "" : "s") + ", not " + std::to_string(args.size()));
    }
    if (conversion != conversions.end()) {
      return arena().cast(conversion->second, args[0]);
    }
    switch (expected) {
    case 1:
      return arena().pw(operation->second, args[0]);
    case 2:
      return arena().pw(operation->second, args[0], args[1]);
    default:
      return arena().pw(operation->second, args[0], args[1], args[2]);
    }
  }

  NodeId sample(const Token& at, int depth) {
    const auto found = names_.find(at.text);
    if (found == names_.end() || found->second.what != What::Grid) {
      fail(at, "`" + at.text + "` is no grid, and only a grid is sampled with `.`");
    }
    next(); // the dot
    const Token component = peek();
    const std::string componentName = name("the name of a component of the grid");
    const GridDeclaration& grid = grids_[found->second.index];
    const auto index = std::find(grid.components.begin(), grid.components.end(), componentName);
    if (index == grid.components.end()) {
      fail(component, "the grid `" + at.text + "` has no component `" + componentName + "`");
    }
    const std::vector<NodeId> coordinates = arguments(depth);
    if (coordinates.empty() || coordinates.size() > static_cast<std::size_t>(MaxLookupDimension)) {
      fail(at,
           "a sample takes 1 to " + std::to_string(MaxLookupDimension) + " coordinates, not " +
               std::to_string(coordinates.size()));
    }
    return arena().lookup(
        grid.id, static_cast<std::int32_t>(index - grid.components.begin()), coordinates);
  }

  NodeId reference(const Token& at) {
    const auto found = names_.find(at.text);
    if (found == names_.end()) {
      const bool output = std::any_of(outputs_.begin(), outputs_.end(), [&](const auto& output) {
        return output.first == at.text;
      });
      fail(at,
           "`" + at.text + "` is not declared" +
               (output ? "; an output cannot be read, but a definition of its value can"
                       : "; an input is declared with `in " + spellName(at.text) + "`"));
    }
    if (found->second.what == What::Grid) {
      fail(at, "the grid `" + at.text + "` is sampled with a component, as `" + at.text + ".c(x)`");
    }
    return found->second.node;
  }

  std::vector<Token> tokens_;
  std::size_t pos_{0};
  std::string file_;
  Program program_;
  std::unordered_map<std::string, Entry> names_;
  std::vector<StateDeclaration> states_;
  std::vector<GridDeclaration> grids_;
  std::vector<std::pair<std::string, NodeId>> outputs_;
  std::vector<std::pair<NodeId, std::string>> checks_;
};

// ============================================================== printer =====

class Printer {
  public:
  explicit Printer(const Program& program) : program_(program), arena_(program.arena()) {}

  std::string run() {
    for (const auto& input : program_.inputs()) {
      taken_.insert(input.name);
    }
    for (const auto& state : program_.state()) {
      taken_.insert(state.name);
    }
    nameGrids();

    std::string out = "#!sderiv " + std::to_string(CoreVersion) + "\n";
    printMeta(out);

    std::vector<std::string> inputs;
    for (const auto& input : program_.inputs()) {
      inputs.push_back(input.name);
    }
    std::sort(inputs.begin(), inputs.end());
    if (!inputs.empty()) {
      out += "in ";
      for (std::size_t i = 0; i < inputs.size(); ++i) {
        out += (i == 0 ? "" : ", ") + spellName(inputs[i]);
      }
      out += "\n";
    }

    std::vector<std::size_t> grids(program_.grids().size());
    for (std::size_t g = 0; g < grids.size(); ++g) {
      grids[g] = g;
    }
    std::sort(grids.begin(), grids.end(), [&](std::size_t a, std::size_t b) {
      return gridNames_[a].name < gridNames_[b].name;
    });
    for (const std::size_t g : grids) {
      printGrid(out, g);
    }

    std::vector<std::size_t> states(program_.state().size());
    for (std::size_t s = 0; s < states.size(); ++s) {
      states[s] = s;
    }
    std::sort(states.begin(), states.end(), [&](std::size_t a, std::size_t b) {
      return program_.state()[a].name < program_.state()[b].name;
    });
    for (const std::size_t s : states) {
      const StateSpec& state = program_.state()[s];
      out += "state ";
      if (state.type != ValueType::Real) {
        out += std::string(typeSpelling(state.type)) + " ";
      }
      out += spellName(state.name) + " = " + spellNumber(state.initial) + "\n";
    }

    std::vector<std::size_t> outputs(program_.outputs().size());
    for (std::size_t o = 0; o < outputs.size(); ++o) {
      outputs[o] = o;
    }
    std::sort(outputs.begin(), outputs.end(), [&](std::size_t a, std::size_t b) {
      return program_.outputs()[a].name < program_.outputs()[b].name;
    });

    // the definitions, in the order a depth-first walk needs them: from the outputs, then from
    // the next values
    names_.assign(arena_.size(), std::string());
    for (const std::size_t o : outputs) {
      define(out, program_.roots()[o]);
    }
    for (const std::size_t s : states) {
      define(out, program_.state()[s].root);
    }
    for (const auto& check : program_.checks()) {
      define(out, check.root);
    }
    for (const std::size_t s : states) {
      out += "next " + spellName(program_.state()[s].name) + " = " +
             operand(program_.state()[s].root) + "\n";
    }
    for (const auto& check : program_.checks()) {
      out += "check " + operand(check.root) + ", " + spellString(check.message) + "\n";
    }
    for (const std::size_t o : outputs) {
      out += "out " + spellName(program_.outputs()[o].name) + " = " + operand(program_.roots()[o]) +
             "\n";
    }
    return out;
  }

  private:
  static const char* typeSpelling(ValueType type) {
    for (const auto& [spelling, value] : conversionSpellings()) {
      if (value == type) {
        return spelling.c_str();
      }
    }
    throw std::invalid_argument(std::string("expr: the core format has no values of type ") +
                                name(type));
  }

  void printMeta(std::string& out) const {
    // generator, source, depends (by path), command, body, line (by name, then position), then
    // keys the specification does not know, by key
    const std::vector<std::string> order = {
        "generator", "source", "depends", "command", "body", "line"};
    std::vector<MetaEntry> entries = program_.meta();
    const auto rank = [&](const std::string& key) {
      const auto found = std::find(order.begin(), order.end(), key);
      return static_cast<std::size_t>(found - order.begin());
    };
    const auto values = [](const MetaEntry& entry) {
      std::vector<std::string> texts;
      for (const auto& value : entry.values) {
        texts.push_back(value.text);
      }
      return texts;
    };
    std::stable_sort(entries.begin(), entries.end(), [&](const MetaEntry& a, const MetaEntry& b) {
      if (rank(a.key) != rank(b.key)) {
        return rank(a.key) < rank(b.key);
      }
      if (rank(a.key) == order.size() && a.key != b.key) {
        return a.key < b.key;
      }
      if (a.key == "depends" || a.key == "line") {
        if (a.key == "line" && !a.values.empty() && !b.values.empty() &&
            a.values.front().text != b.values.front().text) {
          return a.values.front().text < b.values.front().text;
        }
        return values(a) < values(b);
      }
      return false;
    });
    for (const auto& entry : entries) {
      out += "meta " + spellName(entry.key) + " = ";
      for (std::size_t i = 0; i < entry.values.size(); ++i) {
        out += (i == 0 ? "" : ", ") +
               (entry.values[i].quoted ? spellString(entry.values[i].text) : entry.values[i].text);
      }
      out += "\n";
    }
  }

  /// The names of the grids and their components: those the program has, else g0, g1, ... and
  /// c0, c1, ... up to the last component it samples.
  void nameGrids() {
    const std::size_t count = program_.grids().size();
    gridNames_.assign(count, GridNames{});
    std::vector<std::int32_t> components(count, 0);
    for (std::size_t id = 0; id < arena_.size(); ++id) {
      const Node& node = arena_[static_cast<NodeId>(id)];
      if (node.kind == Kind::Lookup && node.grid >= 0 &&
          static_cast<std::size_t>(node.grid) < count) {
        auto& used = components[static_cast<std::size_t>(node.grid)];
        used = std::max(used, node.comp + 1);
      }
    }
    for (std::size_t g = 0; g < count; ++g) {
      if (g < program_.gridNames().size() && !program_.gridNames()[g].name.empty()) {
        gridNames_[g] = program_.gridNames()[g];
      } else {
        std::string name = "g" + std::to_string(g);
        while (taken_.count(name) != 0) {
          name += "_";
        }
        gridNames_[g].name = name;
      }
      for (auto c = static_cast<std::int32_t>(gridNames_[g].components.size()); c < components[g];
           ++c) {
        gridNames_[g].components.push_back("c" + std::to_string(c));
      }
      if (gridNames_[g].components.empty()) {
        gridNames_[g].components.emplace_back("c0");
      }
      if (!taken_.insert(gridNames_[g].name).second) {
        throw std::invalid_argument("expr: the grid `" + gridNames_[g].name +
                                    "` has the name of another grid, an input or a state");
      }
    }
  }

  void printGrid(std::string& out, std::size_t g) const {
    const auto& desc = program_.grids()[g];
    const GridNames& names = gridNames_[g];
    out += "grid " + spellName(names.name) + " (kind = " + spellString(spellKind(desc.kind)) +
           ", path = " + spellString(desc.path) + ", variable = " + spellString(desc.variable) +
           ", components = (";
    for (std::size_t c = 0; c < names.components.size(); ++c) {
      out += (c == 0 ? "" : ", ") + spellString(names.components[c]);
    }
    out += "), interpolation = " + spellString(spellInterpolation(desc.interpolation)) +
           ", boundary = " + spellString(spellBoundary(desc.boundary));
    if (desc.timeAxis.has_value()) {
      out += ", time = " + std::to_string(*desc.timeAxis);
    }
    out += ")\n";
  }

  /// The next free name of a definition: t and a counter, past the names of the program.
  std::string freshName() {
    while (true) {
      std::string candidate = "t" + std::to_string(counter_++);
      if (taken_.count(candidate) == 0) {
        taken_.insert(candidate);
        return candidate;
      }
    }
  }

  /// The number `id` is: a constant, or a constant negated, which the loader reads as a number
  /// as well (negation is exact).
  [[nodiscard]] std::optional<double> number(NodeId id) const {
    const Node& node = arena_[id];
    if (node.kind == Kind::Const) {
      return node.value;
    }
    if (node.kind == Kind::PW && node.fn == Fn::Neg && !node.strict) {
      const auto negated = number(node.a);
      if (negated.has_value()) {
        return -*negated;
      }
    }
    return std::nullopt;
  }

  [[nodiscard]] bool leaf(NodeId id) const {
    return arena_[id].kind == Kind::Field || number(id).has_value();
  }

  /// An operand as a definition reads it: a number, the name of an input or a state, or the name
  /// of the definition of an operation.
  [[nodiscard]] std::string operand(NodeId id) const {
    const Node& node = arena_[id];
    if (const auto value = number(id); value.has_value()) {
      return spellNumber(*value);
    }
    if (node.kind == Kind::Field) {
      return spellName(arena_.channelName(node.ch));
    }
    return names_[static_cast<std::size_t>(id)];
  }

  /// Defines every operation under `root` that is not defined yet, children first.
  void define(std::string& out, NodeId root) {
    if (root == NoNode || leaf(root) || !names_[static_cast<std::size_t>(root)].empty()) {
      return;
    }
    // depth first, children in order, without recursion: a chain of definitions is as deep as it
    // is long
    std::vector<std::pair<NodeId, bool>> stack{{root, false}};
    std::vector<NodeId> children;
    while (!stack.empty()) {
      const auto [id, expanded] = stack.back();
      stack.pop_back();
      if (!names_[static_cast<std::size_t>(id)].empty()) {
        continue;
      }
      if (expanded) {
        const std::string defined = freshName();
        const std::string value = operation(id);
        names_[static_cast<std::size_t>(id)] = defined;
        out += "def " + defined + " = " + value + "\n";
        continue;
      }
      stack.emplace_back(id, true);
      arena_.children(id, children);
      for (auto child = children.rbegin(); child != children.rend(); ++child) {
        if (!leaf(*child) && names_[static_cast<std::size_t>(*child)].empty()) {
          stack.emplace_back(*child, false);
        }
      }
    }
  }

  /// The operation of `id`, on its operands.
  [[nodiscard]] std::string operation(NodeId id) const {
    const Node& node = arena_[id];
    if (node.type == ValueType::I64 || node.type == ValueType::Bool) {
      throw std::invalid_argument(
          "expr: the core format has no whole numbers and truth values of the lowering");
    }
    switch (node.kind) {
    case Kind::PW: {
      if (node.strict) {
        throw std::invalid_argument(
            "expr: the core format has no operations computed as written, which only the "
            "lowering makes");
      }
      switch (node.fn) {
      case Fn::Add:
        return operand(node.a) + " + " + operand(node.b);
      case Fn::Sub:
        return operand(node.a) + " - " + operand(node.b);
      case Fn::Mul:
        return operand(node.a) + " * " + operand(node.b);
      case Fn::Div:
        return operand(node.a) + " / " + operand(node.b);
      case Fn::Neg:
        return "-" + operand(node.a);
      default:
        break;
      }
      const auto& operations = operationSpellings();
      const auto spelling =
          std::find_if(operations.begin(), operations.end(), [&](const auto& candidate) {
            return candidate.second == node.fn;
          });
      if (spelling == operations.end()) {
        throw std::invalid_argument(std::string("expr: the core format has no spelling of ") +
                                    name(node.fn));
      }
      std::string text = spelling->first + "(" + operand(node.a);
      if (arity(node.fn) >= 2) {
        text += ", " + operand(node.b);
      }
      if (arity(node.fn) >= 3) {
        text += ", " + operand(node.c);
      }
      return text + ")";
    }
    case Kind::Cast:
      return std::string(typeSpelling(node.type)) + "(" + operand(node.a) + ")";
    case Kind::Lookup: {
      const GridNames& names = gridNames_[static_cast<std::size_t>(node.grid)];
      std::string text = spellName(names.name) + "." +
                         spellName(names.components[static_cast<std::size_t>(node.comp)]) + "(";
      const NodeId* coordinates = arena_.args(node);
      for (std::int32_t i = 0; i < node.argCount; ++i) {
        text += (i == 0 ? "" : ", ") + operand(coordinates[i]);
      }
      return text + ")";
    }
    case Kind::Contract:
      throw std::invalid_argument(
          "expr: the core format has no contractions, which a consumer introduces when it binds a "
          "program");
    case Kind::GridParam:
    case Kind::GridLoad:
      throw std::invalid_argument(
          "expr: the core format has no written-out lookups, which only the lowering makes");
    case Kind::Const:
    case Kind::Field:
      break;
    }
    return operand(id);
  }

  const Program& program_;
  const Arena& arena_;
  std::set<std::string> taken_;
  std::vector<GridNames> gridNames_;
  std::vector<std::string> names_;
  std::size_t counter_{0};
};

} // namespace

CoreError::CoreError(const CorePosition& position, const std::string& what)
    : SderivError(located(position), what, -1), where_(position) {}

bool isCoreText(const std::string& text) { return text.compare(0, 8, "#!sderiv") == 0; }

Program loadCore(const std::string& text, const std::string& file) {
  Lexer lexer(text, file);
  static_cast<void>(lexer.header());
  Loader loader(lexer.tokens(), file);
  return loader.run();
}

std::string printCore(const Program& program) { return Printer(program).run(); }

} // namespace seissol::expr
