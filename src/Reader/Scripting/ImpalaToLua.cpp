// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "ImpalaToLua.h"

#include <algorithm>
#include <cctype>
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

/// easi's tokenizer (tokenizeImpala()): the words and operators first, longest first, by prefix --
/// `iffy` is `if` and `fy` --, then numbers without an exponent, names, and space.
std::vector<Token> tokenize(const std::string& code, int firstLine, bool literal) {
  static const std::vector<std::pair<std::string, Tok>> Words = {
      {"return", Tok::Return}, {"else", Tok::Else},      {"if", Tok::If},
      {"==", Tok::Eq},         {"!=", Tok::Ne},          {"<=", Tok::Le},
      {">=", Tok::Ge},         {"&&", Tok::And},         {"||", Tok::Or},
      {"<", Tok::Lt},          {">", Tok::Gt},           {"+", Tok::Add},
      {"-", Tok::Sub},         {"*", Tok::Mul},          {"/", Tok::Div},
      {",", Tok::Comma},       {"=", Tok::Assign},       {";", Tok::Semicolon},
      {"(", Tok::BracketOpen}, {")", Tok::BracketClose}, {"{", Tok::Open},
      {"}", Tok::Close}};
  const auto isDigit = [&](std::size_t i) {
    return i < code.size() && std::isdigit(static_cast<unsigned char>(code[i])) != 0;
  };
  const auto isAsciiAlpha = [](char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); };
  std::vector<Token> tokens;
  std::size_t at = 0;
  int line = 0;
  while (at < code.size()) {
    bool found = false;
    for (const auto& [word, kind] : Words) {
      if (code.compare(at, word.size(), word) == 0) {
        tokens.push_back(Token{kind, word, line});
        at += word.size();
        found = true;
        break;
      }
    }
    if (found) {
      continue;
    }
    const char c = code[at];
    if (isDigit(at)) {
      std::size_t end = at;
      while (isDigit(end)) {
        ++end;
      }
      if (end < code.size() && code[end] == '.') {
        ++end;
        while (isDigit(end)) {
          ++end;
        }
      }
      tokens.push_back(Token{Tok::Constant, code.substr(at, end - at), line});
      at = end;
    } else if (isAsciiAlpha(c)) {
      std::size_t end = at + 1;
      while (end < code.size() && (isAsciiAlpha(code[end]) || isDigit(end) || code[end] == '_')) {
        ++end;
      }
      tokens.push_back(Token{Tok::Variable, code.substr(at, end - at), line});
      at = end;
    } else if (std::isspace(static_cast<unsigned char>(c)) != 0) {
      if (c == '\n') {
        ++line;
      }
      ++at;
    } else {
      throw ImpalaError(literal ? firstLine + line : firstLine,
                        "Tokenization error: " + code.substr(at));
    }
  }
  return tokens;
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
  [[nodiscard]] const std::string& text() const { return text_; }

  private:
  std::string text_;
  int line_;
};

/// easi's parser of a function (ParserState), which writes its Lua: every operand in parentheses,
/// the chains grouped to the right. A name is an input where it is one, else a variable, which an
/// assignment before it in the text declares; a function is that of Lua's math of its name.
class Converter {
  public:
  Converter(const ImpalaEntry& entry,
            const std::set<std::string>& inputs,
            Layout& out,
            std::vector<ImpalaChain>& chains)
      : entry_(entry), inputs_(inputs), out_(out), chains_(chains),
        tokens_(tokenize("{" + entry.code + "}", entry.line, entry.literal)) {}

  void function() {
    out_.at(entry_.line);
    out_ << " function __compute__" + entry_.name + "(__in) ";
    body();
    out_ << " end ";
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

  // easi probes several tokens at once, and fails where the code has fewer
  bool probe(std::initializer_list<Tok> kinds) {
    if (tokens_.size() < position_ + kinds.size()) {
      fail("Out of tokens");
    }
    std::size_t i = position_;
    bool match = true;
    for (const Tok kind : kinds) {
      match = match && tokens_[i++].kind == kind;
    }
    return match;
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
      fail("Out of tokens");
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
        out_ << ";";
      } else if (probe({Tok::Variable, Tok::Assign})) {
        assignment();
        expect(Tok::Semicolon);
        out_ << ";";
      } else {
        expression();
        expect(Tok::Semicolon);
        out_ << ";";
      }
    }
  }

  void ifStatement() {
    expect(Tok::If);
    out_ << "if (";
    boolean();
    out_ << ") then ";
    body();
    if (probeConsume(Tok::Else)) {
      out_ << " else ";
      body();
    }
    out_ << " end ";
  }

  void returnStatement() {
    expect(Tok::Return);
    out_ << " do return ";
    expression();
    out_ << " end ";
  }

  void assignment() {
    const std::string variable = expect(Tok::Variable);
    // declared before its value is read, as easi has it
    variables_.insert(variable);
    expect(Tok::Assign);
    out_ << variable + " = ";
    expression();
  }

  void call() {
    const std::string function = expect(Tok::Variable);
    out_ << "math." + function + "(";
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
  }

  void atom() {
    out_ << "(";
    if (probeConsume(Tok::BracketOpen)) {
      expression();
      expect(Tok::BracketClose);
    } else if (probe({Tok::Variable, Tok::BracketOpen})) {
      call();
    } else if (probe({Tok::Variable})) {
      const std::string variable = expect(Tok::Variable);
      if (inputs_.count(variable) != 0) {
        out_ << "__in[\"" + variable + "\"]";
      } else if (variables_.count(variable) != 0) {
        out_ << variable;
      } else {
        --position_;
        fail("Unknown variable: " + variable);
      }
    } else if (probe({Tok::Constant})) {
      out_ << expect(Tok::Constant);
    } else {
      fail("Unexpected token when parsing expression: `" + tokens_[position_].text + "`");
    }
    out_ << ")";
  }

  void unary() {
    out_ << "(";
    if (probeConsume(Tok::Add)) {
      out_ << "+";
    }
    if (probeConsume(Tok::Sub)) {
      out_ << "-";
    }
    atom();
    out_ << ")";
  }

  /// A product, grouped to the right; whether it has an operator.
  bool product() {
    out_ << "(";
    unary();
    bool chain = false;
    if (probeConsume(Tok::Mul)) {
      out_ << " * ";
      product();
      chain = true;
    }
    if (probe({Tok::Div})) {
      const int line = lineOf(consume());
      out_ << " / ";
      if (product()) {
        chains_.push_back(ImpalaChain{entry_.name, line, "/", "* or /"});
      }
      chain = true;
    }
    out_ << ")";
    return chain;
  }

  /// A sum, grouped to the right; whether it has an operator.
  bool sum() {
    out_ << "(";
    product();
    bool chain = false;
    if (probeConsume(Tok::Add)) {
      out_ << " + ";
      sum();
      chain = true;
    }
    if (probe({Tok::Sub})) {
      const int line = lineOf(consume());
      out_ << " - ";
      if (sum()) {
        chains_.push_back(ImpalaChain{entry_.name, line, "-", "+ or -"});
      }
      chain = true;
    }
    out_ << ")";
    return chain;
  }

  void expression() { sum(); }

  void comparison() {
    out_ << "(";
    if (probeConsume(Tok::BracketOpen)) {
      boolean();
      expect(Tok::BracketClose);
    } else {
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
    }
    out_ << ")";
  }

  void conjunction() {
    out_ << "(";
    comparison();
    if (probeConsume(Tok::And)) {
      out_ << " and ";
      conjunction();
    }
    out_ << ")";
  }

  void boolean() {
    out_ << "(";
    conjunction();
    if (probeConsume(Tok::Or)) {
      out_ << " or ";
      boolean();
    }
    out_ << ")";
  }

  const ImpalaEntry& entry_;
  const std::set<std::string>& inputs_;
  Layout& out_;
  std::vector<ImpalaChain>& chains_;
  std::vector<Token> tokens_;
  std::size_t position_{0};
  std::set<std::string> variables_;
};

} // namespace

ImpalaLua impalaToLua(const std::vector<ImpalaEntry>& entries,
                      const std::set<std::string>& inputs) {
  ImpalaLua lua;
  lua.firstLine = entries.empty() ? 1 : entries.front().line;
  Layout out(lua.firstLine);
  out << "function f(__in__main) ";
  // easi reads the map into a map by name (the last of a name wins), and returns the parameters in
  // the order of an unordered_map made of it -- which only matters where two of them share a
  // variable, a global of Lua
  std::map<std::string, std::size_t> byName;
  for (std::size_t i = 0; i < entries.size(); ++i) {
    byName[entries[i].name] = i;
  }
  for (std::size_t i = 0; i < entries.size(); ++i) {
    if (byName.at(entries[i].name) == i) {
      Converter(entries[i], inputs, out, lua.chains).function();
    }
  }
  std::map<std::string, std::string> precode;
  for (const auto& [name, index] : byName) {
    precode[name] = entries[index].code;
  }
  const std::unordered_map<std::string, std::string> code(precode.begin(), precode.end());
  out << " return { ";
  for (const auto& [name, text] : code) {
    out << name + " = __compute__" + name + "(__in__main),";
  }
  out << "} end";
  lua.code = out.text();
  return lua;
}

} // namespace seissol::reader::scripting
