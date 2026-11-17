// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "ImpalaToLua.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::reader::scripting {

namespace {

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
  Semicolon
};

struct Token {
  Tok kind;
  std::string text;
  /// the line of the text it is on, from 0
  int line;
};

bool isLetter(char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); }
bool isDigit(char c) { return c >= '0' && c <= '9'; }
bool isSpace(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v';
}

/// ImpalaJIT's tokens, as easi's conversion reads them: names [A-Za-z][A-Za-z0-9_]*, of which if,
/// else and return are keywords as whole words only (`ifactor` is a name); numbers [0-9]+(.[0-9]*)?
/// without an exponent; the operators, the longer first.
std::vector<Token> tokenize(const std::string& code, int firstLine, bool literal) {
  static const std::vector<std::pair<std::string, Tok>> Symbols = {{"==", Tok::Eq},
                                                                   {"!=", Tok::Ne},
                                                                   {"<=", Tok::Le},
                                                                   {">=", Tok::Ge},
                                                                   {"&&", Tok::And},
                                                                   {"||", Tok::Or},
                                                                   {"<", Tok::Lt},
                                                                   {">", Tok::Gt},
                                                                   {"+", Tok::Add},
                                                                   {"-", Tok::Sub},
                                                                   {"*", Tok::Mul},
                                                                   {"/", Tok::Div},
                                                                   {",", Tok::Comma},
                                                                   {"=", Tok::Assign},
                                                                   {";", Tok::Semicolon},
                                                                   {"(", Tok::BracketOpen},
                                                                   {")", Tok::BracketClose},
                                                                   {"{", Tok::Open},
                                                                   {"}", Tok::Close}};
  static const std::unordered_map<std::string, Tok> Keywords = {
      {"if", Tok::If}, {"else", Tok::Else}, {"return", Tok::Return}};
  std::vector<Token> tokens;
  std::size_t at = 0;
  int line = 0;
  while (at < code.size()) {
    const char c = code[at];
    if (isSpace(c)) {
      line += c == '\n' ? 1 : 0;
      ++at;
    } else if (isLetter(c)) {
      std::size_t end = at + 1;
      while (end < code.size() && (isLetter(code[end]) || isDigit(code[end]) || code[end] == '_')) {
        ++end;
      }
      const std::string word = code.substr(at, end - at);
      const auto keyword = Keywords.find(word);
      tokens.push_back(
          Token{keyword != Keywords.end() ? keyword->second : Tok::Variable, word, line});
      at = end;
    } else if (isDigit(c)) {
      std::size_t end = at + 1;
      while (end < code.size() && isDigit(code[end])) {
        ++end;
      }
      if (end < code.size() && code[end] == '.') {
        ++end;
        while (end < code.size() && isDigit(code[end])) {
          ++end;
        }
      }
      tokens.push_back(Token{Tok::Constant, code.substr(at, end - at), line});
      at = end;
    } else {
      bool found = false;
      for (const auto& [symbol, kind] : Symbols) {
        if (code.compare(at, symbol.size(), symbol) == 0) {
          tokens.push_back(Token{kind, symbol, line});
          at += symbol.size();
          found = true;
          break;
        }
      }
      if (!found) {
        throw ImpalaError(literal ? firstLine + line : firstLine,
                          "Tokenization error: " + code.substr(at));
      }
    }
  }
  return tokens;
}

/// The Lua name of a variable of Impala: apart from the keywords of Lua and the names of the
/// conversion itself.
std::string luaVariable(const std::string& name) { return "v_" + name; }

/// `text` as a string of Lua: a byte other than a letter, a digit or _ as its decimal escape.
std::string luaString(const std::string& text) {
  std::string quoted = "\"";
  for (const char c : text) {
    if (isLetter(c) || isDigit(c) || c == '_') {
      quoted += c;
    } else {
      const auto byte = std::to_string(static_cast<unsigned char>(c));
      quoted += "\\" + std::string(3 - byte.size(), '0') + byte;
    }
  }
  return quoted + "\"";
}

/// Whether `name` is a name of Impala, [A-Za-z][A-Za-z0-9_]*: an input of another name no code
/// reads.
bool isImpalaName(const std::string& name) {
  return !name.empty() && isLetter(name.front()) &&
         std::all_of(name.begin(), name.end(), [](char c) {
           return isLetter(c) || isDigit(c) || c == '_';
         });
}

/// A constant of Impala as a number of Lua with a fraction: ImpalaJIT computes in binary64, and an
/// integer of Lua would overflow.
std::string luaConstant(const std::string& constant) {
  if (constant.find('.') == std::string::npos) {
    return constant + ".0";
  }
  if (constant.back() == '.') {
    return constant + "0";
  }
  return constant;
}

/// The Lua being written, and the line of the file it is at.
class Layout {
  public:
  explicit Layout(int line) : line_(line) {}

  /// Breaks lines up to the line `line` of the file; never back.
  void at(int line) {
    for (; line_ < line; ++line_) {
      text_ += '\n';
    }
  }
  Layout& operator<<(const std::string& text) {
    text_ += text;
    return *this;
  }
  /// Appends `other`, which starts at the line this one is at.
  void append(const Layout& other) {
    text_ += other.text_;
    line_ = other.line_;
  }
  [[nodiscard]] const std::string& text() const { return text_; }

  /// What has been written so far, to go back to (restore()).
  struct Mark {
    std::size_t length;
    int line;
  };
  [[nodiscard]] Mark mark() const { return Mark{text_.size(), line_}; }
  void restore(const Mark& mark) {
    text_.resize(mark.length);
    line_ = mark.line;
  }

  private:
  std::string text_;
  int line_;
};

/// The parser of a function of ImpalaJIT (as the parser of easi's conversion), which writes its
/// Lua: + - and * / left-associative; a unary - on an atom; a condition a comparison of two
/// expressions or a condition in brackets, with && binding more strongly than ||. The inputs are
/// local copies, which an assignment changes; every other variable is local to the function, known
/// from its first assignment in the text on (not on its right-hand side), and nil until assigned.
/// A function is that of C's math library of its name: min, max and abs are fmin, fmax and fabs,
/// floor and ceil give floats.
class Converter {
  public:
  Converter(const ImpalaEntry& entry,
            const std::set<std::string>& inputs,
            std::vector<ImpalaChain>& chains)
      : entry_(entry), inputs_(inputs), chains_(chains), out_(entry.line),
        tokens_(tokenize("{" + entry.code + "}", entry.line, entry.literal)) {}

  /// Writes the function `function` of the entry to `target`, on the lines of its code.
  void function(Layout& target, const std::string& function) {
    body();
    if (position_ != tokens_.size()) {
      fail("Unexpected code after the end of the function");
    }
    target.at(entry_.line);
    target << " function " + function + "(__in) ";
    for (const auto& input : inputs_) {
      if (isImpalaName(input)) {
        target << "local " + luaVariable(input) + " = __in[\"" + input + "\"] ";
      }
    }
    for (const auto& local : localsInOrder_) {
      target << "local " + luaVariable(local) + " ";
    }
    target.append(out_);
    target << " end ";
  }

  private:
  [[nodiscard]] int lineOf(const Token& token) const {
    return entry_.literal ? entry_.line + token.line : entry_.line;
  }

  [[noreturn]] void fail(const std::string& what) const {
    const int line =
        tokens_.empty() ? entry_.line : lineOf(tokens_[std::min(position_, tokens_.size() - 1)]);
    throw ImpalaError(line, what);
  }

  /// Consumes the next token: what is written next is on its line.
  const Token& consume() {
    const Token& token = tokens_[position_++];
    out_.at(lineOf(token));
    return token;
  }

  [[nodiscard]] bool probe(std::initializer_list<Tok> kinds) const {
    if (tokens_.size() < position_ + kinds.size()) {
      return false;
    }
    std::size_t i = position_;
    return std::all_of(
        kinds.begin(), kinds.end(), [&](Tok kind) { return tokens_[i++].kind == kind; });
  }
  bool probeConsume(Tok kind) {
    if (probe({kind})) {
      consume();
      return true;
    }
    return false;
  }
  std::string expect(Tok kind) {
    if (tokens_.size() <= position_) {
      fail("Unexpected end of the function");
    }
    if (tokens_[position_].kind != kind) {
      fail("Unexpected token: `" + tokens_[position_].text + "`");
    }
    return consume().text;
  }

  void body() {
    expect(Tok::Open);
    while (!probeConsume(Tok::Close)) {
      if (probe({Tok::If})) {
        ifStatement();
      } else if (probe({Tok::Return})) {
        returnStatement();
        expect(Tok::Semicolon);
      } else if (probe({Tok::Variable, Tok::Assign})) {
        assignment();
        expect(Tok::Semicolon);
      } else {
        // ImpalaJIT evaluates an expression statement and drops its value
        out_ << " local _ = ";
        expression();
        expect(Tok::Semicolon);
      }
    }
  }

  void ifStatement() {
    expect(Tok::If);
    out_ << " if ";
    boolean();
    out_ << " then ";
    body();
    if (probeConsume(Tok::Else)) {
      out_ << " else ";
      body();
    }
    out_ << " end ";
  }

  void returnStatement() {
    expect(Tok::Return);
    // `do ... end` allows a return in the middle of a block
    out_ << " do return ";
    expression();
    out_ << " end ";
  }

  void assignment() {
    const std::string variable = expect(Tok::Variable);
    expect(Tok::Assign);
    out_ << " " + luaVariable(variable) + " = ";
    expression();
    // known after its assignment, as in ImpalaJIT; an input stays one, its local copy changed
    if (inputs_.count(variable) == 0 && locals_.insert(variable).second) {
      localsInOrder_.push_back(variable);
    }
  }

  void call() {
    const std::string function = expect(Tok::Variable);
    // the functions of C's math library: min, max and abs of floats, floor and ceil to floats
    static const std::map<std::string, std::string> Renamed = {
        {"min", "fmin"}, {"max", "fmax"}, {"abs", "fabs"}};
    const auto renamed = Renamed.find(function);
    const bool rounding = function == "floor" || function == "ceil";
    if (rounding) {
      out_ << "(";
    }
    out_ << "math." + (renamed != Renamed.end() ? renamed->second : function) + "(";
    expect(Tok::BracketOpen);
    if (!probeConsume(Tok::BracketClose)) {
      expression();
      while (!probeConsume(Tok::BracketClose)) {
        expect(Tok::Comma);
        out_ << ", ";
        expression();
      }
    }
    out_ << ")";
    if (rounding) {
      out_ << " + 0.0)";
    }
  }

  void atom() {
    if (probeConsume(Tok::BracketOpen)) {
      out_ << "(";
      expression();
      expect(Tok::BracketClose);
      out_ << ")";
    } else if (probe({Tok::Variable, Tok::BracketOpen})) {
      call();
    } else if (probe({Tok::Variable})) {
      const std::string variable = tokens_[position_].text;
      if (inputs_.count(variable) == 0 && locals_.count(variable) == 0) {
        fail("Unknown variable: " + variable);
      }
      consume();
      out_ << luaVariable(variable);
    } else if (probe({Tok::Constant})) {
      out_ << luaConstant(expect(Tok::Constant));
    } else if (position_ < tokens_.size()) {
      fail("Unexpected token when parsing expression: `" + tokens_[position_].text + "`");
    } else {
      fail("Unexpected end of the function");
    }
  }

  void unary() {
    probeConsume(Tok::Add);
    if (probeConsume(Tok::Sub)) {
      out_ << "(-";
      atom();
      out_ << ")";
    } else {
      atom();
    }
  }

  /// A chain of `first` (such as -) and then `following` (such as + or -), which easi 1.5 to 1.7
  /// grouped to the right, at the line of `at`.
  void regrouped(const Token& at, const char* first, const char* following) {
    chains_.push_back(ImpalaChain{entry_.name, lineOf(at), first, following});
  }

  /// A product, grouped to the left.
  void product() {
    out_ << "(";
    unary();
    const Token* division = nullptr;
    while (probe({Tok::Mul}) || probe({Tok::Div})) {
      const Token& op = consume();
      if (division != nullptr) {
        regrouped(*division, "/", "* or /");
        division = nullptr;
      }
      if (op.kind == Tok::Div) {
        division = &op;
      }
      out_ << (op.kind == Tok::Mul ? " * " : " / ");
      unary();
    }
    out_ << ")";
  }

  /// A sum, grouped to the left.
  void sum() {
    out_ << "(";
    product();
    const Token* difference = nullptr;
    while (probe({Tok::Add}) || probe({Tok::Sub})) {
      const Token& op = consume();
      if (difference != nullptr) {
        regrouped(*difference, "-", "+ or -");
        difference = nullptr;
      }
      if (op.kind == Tok::Sub) {
        difference = &op;
      }
      out_ << (op.kind == Tok::Add ? " + " : " - ");
      product();
    }
    out_ << ")";
  }

  void expression() { sum(); }

  void comparison() {
    out_ << "(";
    expression();
    static const std::vector<std::pair<Tok, const char*>> Comparisons = {{Tok::Eq, " == "},
                                                                         {Tok::Ne, " ~= "},
                                                                         {Tok::Le, " <= "},
                                                                         {Tok::Lt, " < "},
                                                                         {Tok::Ge, " >= "},
                                                                         {Tok::Gt, " > "}};
    bool found = false;
    for (const auto& [kind, text] : Comparisons) {
      if (probeConsume(kind)) {
        out_ << text;
        found = true;
        break;
      }
    }
    if (!found) {
      fail("No comparison operator given");
    }
    expression();
    out_ << ")";
  }

  /// A comparison, or a condition in brackets: a bracket opens either, as in (a < b) && c < d and
  /// (x + 1) * 2 < y; the condition is tried first.
  void atomBoolean() {
    if (probe({Tok::BracketOpen})) {
      const std::size_t position = position_;
      const auto mark = out_.mark();
      const std::size_t chains = chains_.size();
      try {
        expect(Tok::BracketOpen);
        out_ << "(";
        boolean();
        expect(Tok::BracketClose);
        out_ << ")";
        return;
      } catch (const ImpalaError&) {
        position_ = position;
        out_.restore(mark);
        chains_.resize(chains);
      }
    }
    comparison();
  }

  void conjunction() {
    out_ << "(";
    atomBoolean();
    while (probeConsume(Tok::And)) {
      out_ << " and ";
      atomBoolean();
    }
    out_ << ")";
  }

  void boolean() {
    out_ << "(";
    conjunction();
    while (probeConsume(Tok::Or)) {
      out_ << " or ";
      conjunction();
    }
    out_ << ")";
  }

  const ImpalaEntry& entry_;
  const std::set<std::string>& inputs_;
  std::vector<ImpalaChain>& chains_;
  Layout out_;
  std::vector<Token> tokens_;
  std::size_t position_{0};
  std::set<std::string> locals_;
  std::vector<std::string> localsInOrder_;
};

} // namespace

ImpalaLua impalaToLua(const std::vector<ImpalaEntry>& entries,
                      const std::set<std::string>& inputs) {
  ImpalaLua lua;
  lua.firstLine = entries.empty() ? 1 : entries.front().line;
  Layout out(lua.firstLine);
  out << "function f(__in__main) ";
  // easi reads the map into a map by name (the last of a name wins); the functions keep their
  // variables to themselves, so the order of the parameters does not matter
  std::map<std::string, std::size_t> byName;
  for (std::size_t i = 0; i < entries.size(); ++i) {
    byName[entries[i].name] = i;
  }
  // a function by index, as a name of the map need not be a name of Lua
  const auto function = [](std::size_t index) { return "__compute" + std::to_string(index); };
  for (std::size_t i = 0; i < entries.size(); ++i) {
    if (byName.at(entries[i].name) == i) {
      Converter(entries[i], inputs, lua.chains).function(out, function(i));
    }
  }
  out << " return { ";
  for (const auto& [name, index] : byName) {
    out << "[" + luaString(name) + "] = " + function(index) + "(__in__main), ";
  }
  out << "} end";
  lua.code = out.text();
  return lua;
}

} // namespace seissol::reader::scripting
