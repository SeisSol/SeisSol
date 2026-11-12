// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "LuaSyntax.h"

#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace seissol::reader::scripting::lua {

namespace {

enum class Tok : std::uint8_t {
  Name,
  Number,
  String,
  // keywords
  And,
  Break,
  Do,
  Else,
  Elseif,
  End,
  False,
  For,
  Function,
  Goto,
  If,
  In,
  Local,
  Nil,
  Not,
  Or,
  Repeat,
  Return,
  Then,
  True,
  Until,
  While,
  // operators
  Plus,
  Minus,
  Star,
  Slash,
  DoubleSlash,
  Percent,
  Caret,
  Hash,
  Ampersand,
  Tilde,
  Bar,
  ShiftLeft,
  ShiftRight,
  Equal,
  NotEqual,
  LessEqual,
  GreaterEqual,
  Less,
  Greater,
  Assign,
  OpenParen,
  CloseParen,
  OpenBrace,
  CloseBrace,
  OpenBracket,
  CloseBracket,
  DoubleColon,
  Semicolon,
  Colon,
  Comma,
  Dot,
  Concat,
  Dots,
  Eof
};

struct Token {
  Tok kind{Tok::Eof};
  std::string text;
  double number{0.0};
  std::int64_t integer{0};
  bool isInteger{false};
  int line{1};
};

const std::array<std::pair<const char*, Tok>, 22> Keywords{
    {{"and", Tok::And},     {"break", Tok::Break},   {"do", Tok::Do},
     {"else", Tok::Else},   {"elseif", Tok::Elseif}, {"end", Tok::End},
     {"false", Tok::False}, {"for", Tok::For},       {"function", Tok::Function},
     {"goto", Tok::Goto},   {"if", Tok::If},         {"in", Tok::In},
     {"local", Tok::Local}, {"nil", Tok::Nil},       {"not", Tok::Not},
     {"or", Tok::Or},       {"repeat", Tok::Repeat}, {"return", Tok::Return},
     {"then", Tok::Then},   {"true", Tok::True},     {"until", Tok::Until},
     {"while", Tok::While}}};

/// Lua's lexer (llex.c).
class Lexer {
  public:
  explicit Lexer(const std::string& source) : source_(source) {
    // a first line starting with # is skipped, as lua does
    if (!source_.empty() && source_[0] == '#') {
      while (position_ < source_.size() && source_[position_] != '\n') {
        ++position_;
      }
    }
  }

  Token next() {
    skipSpaceAndComments();
    Token token;
    token.line = line_;
    if (position_ >= source_.size()) {
      return token;
    }
    const char c = source_[position_];
    if (std::isalpha(static_cast<unsigned char>(c)) != 0 || c == '_') {
      const std::size_t start = position_;
      while (position_ < source_.size() &&
             (std::isalnum(static_cast<unsigned char>(source_[position_])) != 0 ||
              source_[position_] == '_')) {
        ++position_;
      }
      token.text = source_.substr(start, position_ - start);
      token.kind = Tok::Name;
      for (const auto& [word, kind] : Keywords) {
        if (token.text == word) {
          token.kind = kind;
          break;
        }
      }
      return token;
    }
    if (std::isdigit(static_cast<unsigned char>(c)) != 0 ||
        (c == '.' && position_ + 1 < source_.size() &&
         std::isdigit(static_cast<unsigned char>(source_[position_ + 1])) != 0)) {
      return numeral();
    }
    if (c == '"' || c == '\'') {
      return shortString(c);
    }
    if (c == '[') {
      const std::size_t level = longBracketLevel();
      if (level != std::string::npos) {
        token.kind = Tok::String;
        token.text = longString(level);
        return token;
      }
    }
    const std::array<std::pair<const char*, Tok>, 33> Operators{
        {{"...", Tok::Dots},       {"..", Tok::Concat},      {"//", Tok::DoubleSlash},
         {"<<", Tok::ShiftLeft},   {">>", Tok::ShiftRight},  {"==", Tok::Equal},
         {"~=", Tok::NotEqual},    {"<=", Tok::LessEqual},   {">=", Tok::GreaterEqual},
         {"::", Tok::DoubleColon}, {"+", Tok::Plus},         {"-", Tok::Minus},
         {"*", Tok::Star},         {"/", Tok::Slash},        {"%", Tok::Percent},
         {"^", Tok::Caret},        {"#", Tok::Hash},         {"&", Tok::Ampersand},
         {"~", Tok::Tilde},        {"|", Tok::Bar},          {"<", Tok::Less},
         {">", Tok::Greater},      {"=", Tok::Assign},       {"(", Tok::OpenParen},
         {")", Tok::CloseParen},   {"{", Tok::OpenBrace},    {"}", Tok::CloseBrace},
         {"[", Tok::OpenBracket},  {"]", Tok::CloseBracket}, {";", Tok::Semicolon},
         {":", Tok::Colon},        {",", Tok::Comma},        {".", Tok::Dot}}};
    for (const auto& [spelling, kind] : Operators) {
      const std::string text(spelling);
      if (source_.compare(position_, text.size(), text) == 0) {
        position_ += text.size();
        token.kind = kind;
        token.text = text;
        return token;
      }
    }
    throw LuaError(line_, std::string("unexpected symbol near '") + c + "'");
  }

  private:
  void skipSpaceAndComments() {
    while (position_ < source_.size()) {
      const char c = source_[position_];
      if (c == '\n') {
        ++line_;
        ++position_;
      } else if (std::isspace(static_cast<unsigned char>(c)) != 0) {
        ++position_;
      } else if (source_.compare(position_, 2, "--") == 0) {
        position_ += 2;
        if (position_ < source_.size() && source_[position_] == '[') {
          const std::size_t level = longBracketLevel();
          if (level != std::string::npos) {
            static_cast<void>(longString(level));
            continue;
          }
        }
        while (position_ < source_.size() && source_[position_] != '\n') {
          ++position_;
        }
      } else {
        return;
      }
    }
  }

  /// At `[`: the level of a long bracket `[==[`, or npos where none opens here.
  std::size_t longBracketLevel() const {
    std::size_t at = position_ + 1;
    std::size_t level = 0;
    while (at < source_.size() && source_[at] == '=') {
      ++level;
      ++at;
    }
    return at < source_.size() && source_[at] == '[' ? level : std::string::npos;
  }

  std::string longString(std::size_t level) {
    const int startLine = line_;
    position_ += level + 2;
    // a newline right after the opening bracket is not part of the string
    if (position_ < source_.size() && source_[position_] == '\r') {
      ++position_;
    }
    if (position_ < source_.size() && source_[position_] == '\n') {
      ++position_;
      ++line_;
    }
    const std::string close = "]" + std::string(level, '=') + "]";
    const std::size_t end = source_.find(close, position_);
    if (end == std::string::npos) {
      throw LuaError(startLine, "unfinished long string or comment");
    }
    std::string text = source_.substr(position_, end - position_);
    for (const char c : text) {
      if (c == '\n') {
        ++line_;
      }
    }
    position_ = end + close.size();
    return text;
  }

  Token shortString(char quote) {
    Token token;
    token.kind = Tok::String;
    token.line = line_;
    ++position_;
    while (true) {
      if (position_ >= source_.size() || source_[position_] == '\n') {
        throw LuaError(line_, "unfinished string");
      }
      const char c = source_[position_++];
      if (c == quote) {
        return token;
      }
      if (c != '\\') {
        token.text.push_back(c);
        continue;
      }
      if (position_ >= source_.size()) {
        throw LuaError(line_, "unfinished string");
      }
      const char e = source_[position_++];
      switch (e) {
      case 'a':
        token.text.push_back('\a');
        break;
      case 'b':
        token.text.push_back('\b');
        break;
      case 'f':
        token.text.push_back('\f');
        break;
      case 'n':
        token.text.push_back('\n');
        break;
      case 'r':
        token.text.push_back('\r');
        break;
      case 't':
        token.text.push_back('\t');
        break;
      case 'v':
        token.text.push_back('\v');
        break;
      case '\\':
      case '"':
      case '\'':
        token.text.push_back(e);
        break;
      case '\n':
        ++line_;
        token.text.push_back('\n');
        break;
      case 'x': {
        if (position_ + 2 > source_.size()) {
          throw LuaError(line_, "hexadecimal digit expected");
        }
        token.text.push_back(
            static_cast<char>(std::stoi(source_.substr(position_, 2), nullptr, 16)));
        position_ += 2;
        break;
      }
      case 'z':
        while (position_ < source_.size() &&
               std::isspace(static_cast<unsigned char>(source_[position_])) != 0) {
          if (source_[position_] == '\n') {
            ++line_;
          }
          ++position_;
        }
        break;
      default:
        if (std::isdigit(static_cast<unsigned char>(e)) != 0) {
          int value = e - '0';
          for (int i = 0; i < 2 && position_ < source_.size() &&
                          std::isdigit(static_cast<unsigned char>(source_[position_])) != 0;
               ++i) {
            value = 10 * value + (source_[position_++] - '0');
          }
          if (value > 255) {
            throw LuaError(line_, "decimal escape too large");
          }
          token.text.push_back(static_cast<char>(value));
        } else {
          throw LuaError(line_, std::string("invalid escape sequence '\\") + e + "'");
        }
      }
    }
  }

  /// A numeral of Lua: an integer where it is one and fits (a decimal one), else a float.
  Token numeral() {
    Token token;
    token.kind = Tok::Number;
    token.line = line_;
    const std::size_t start = position_;
    bool hex = false;
    if (source_.compare(position_, 2, "0x") == 0 || source_.compare(position_, 2, "0X") == 0) {
      hex = true;
      position_ += 2;
    }
    const char exponent = hex ? 'p' : 'e';
    while (position_ < source_.size()) {
      const char c = source_[position_];
      const char lower = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      if (lower == exponent) {
        ++position_;
        if (position_ < source_.size() &&
            (source_[position_] == '+' || source_[position_] == '-')) {
          ++position_;
        }
        continue;
      }
      if (std::isxdigit(static_cast<unsigned char>(c)) != 0 || c == '.') {
        ++position_;
        continue;
      }
      break;
    }
    const std::string text = source_.substr(start, position_ - start);
    token.text = text;
    const bool fractional =
        text.find('.') != std::string::npos ||
        (hex ? text.find_first_of("pP") : text.find_first_of("eE")) != std::string::npos;
    if (!fractional) {
      if (hex) {
        // wraps around, as Lua's hexadecimal integers do
        std::uint64_t value = 0;
        for (std::size_t i = 2; i < text.size(); ++i) {
          const int digit = std::isdigit(static_cast<unsigned char>(text[i])) != 0
                                ? text[i] - '0'
                                : std::tolower(static_cast<unsigned char>(text[i])) - 'a' + 10;
          value = value * 16 + static_cast<std::uint64_t>(digit);
        }
        token.isInteger = true;
        token.integer = static_cast<std::int64_t>(value);
        token.number = static_cast<double>(token.integer);
        return token;
      }
      std::int64_t value = 0;
      const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
      if (error == std::errc() && end == text.data() + text.size()) {
        token.isInteger = true;
        token.integer = value;
        token.number = static_cast<double>(value);
        return token;
      }
      // a decimal integer that does not fit is a float
    }
    double value = 0.0;
    const char* first = text.data() + (hex ? 2 : 0);
    const auto [end, error] =
        std::from_chars(first,
                        text.data() + text.size(),
                        value,
                        hex ? std::chars_format::hex : std::chars_format::general);
    if (end != text.data() + text.size() ||
        (error != std::errc() && error != std::errc::result_out_of_range)) {
      throw LuaError(line_, "malformed number near '" + text + "'");
    }
    if (error == std::errc::result_out_of_range) {
      // beyond binary64: an infinity or zero, as strtod gives
      value = text.find_first_of("-") != std::string::npos && !hex ? 0.0 : HUGE_VAL;
    }
    token.number = value;
    return token;
  }

  const std::string& source_;
  std::size_t position_{0};
  int line_{1};
};

/// Lua's parser (lparser.c), into the tree of LuaSyntax.h.
class Parser {
  public:
  explicit Parser(const std::string& source) : lexer_(source) {
    current_ = lexer_.next();
    ahead_ = lexer_.next();
  }

  std::shared_ptr<FunctionBody> chunk() {
    auto main = std::make_shared<FunctionBody>();
    main->vararg = true;
    main->line = 0;
    main->name = "main chunk";
    block(main->block);
    if (current_.kind != Tok::Eof) {
      throw LuaError(current_.line, "'<eof>' expected near '" + current_.text + "'");
    }
    return main;
  }

  private:
  void advance() {
    current_ = std::move(ahead_);
    ahead_ = lexer_.next();
  }
  bool accept(Tok kind) {
    if (current_.kind == kind) {
      advance();
      return true;
    }
    return false;
  }
  void expect(Tok kind, const char* what) {
    if (!accept(kind)) {
      throw LuaError(current_.line,
                     std::string("'") + what + "' expected near '" + describe(current_) + "'");
    }
  }
  static std::string describe(const Token& token) {
    return token.kind == Tok::Eof ? std::string("<eof>") : token.text;
  }
  std::string name() {
    if (current_.kind != Tok::Name) {
      throw LuaError(current_.line, "<name> expected near '" + describe(current_) + "'");
    }
    std::string text = current_.text;
    advance();
    return text;
  }

  static bool blockEnds(Tok kind) {
    return kind == Tok::Else || kind == Tok::Elseif || kind == Tok::End || kind == Tok::Until ||
           kind == Tok::Eof;
  }

  void block(Block& out) {
    while (!blockEnds(current_.kind)) {
      if (current_.kind == Tok::Return) {
        out.statements.push_back(returnStatement());
        return;
      }
      out.statements.push_back(statement());
    }
  }

  Statement returnStatement() {
    Statement statement;
    statement.kind = Statement::Kind::Return;
    statement.line = current_.line;
    advance();
    if (!blockEnds(current_.kind) && current_.kind != Tok::Semicolon) {
      expressionList(statement.values);
    }
    accept(Tok::Semicolon);
    return statement;
  }

  Statement statement() {
    Statement statement;
    statement.line = current_.line;
    switch (current_.kind) {
    case Tok::Semicolon:
      advance();
      statement.kind = Statement::Kind::Empty;
      return statement;
    case Tok::If:
      return ifStatement();
    case Tok::While: {
      advance();
      statement.kind = Statement::Kind::While;
      statement.conditions.push_back(expression());
      expect(Tok::Do, "do");
      statement.body = std::make_shared<Block>();
      block(*statement.body);
      expect(Tok::End, "end");
      return statement;
    }
    case Tok::Do: {
      advance();
      statement.kind = Statement::Kind::Do;
      statement.body = std::make_shared<Block>();
      block(*statement.body);
      expect(Tok::End, "end");
      return statement;
    }
    case Tok::For:
      return forStatement();
    case Tok::Repeat: {
      advance();
      statement.kind = Statement::Kind::Repeat;
      statement.body = std::make_shared<Block>();
      block(*statement.body);
      expect(Tok::Until, "until");
      statement.conditions.push_back(expression());
      return statement;
    }
    case Tok::Function: {
      advance();
      statement.kind = Statement::Kind::Function;
      statement.path.push_back(name());
      while (accept(Tok::Dot)) {
        statement.path.push_back(name());
      }
      if (accept(Tok::Colon)) {
        statement.method = name();
      }
      std::string full = statement.path.front();
      for (std::size_t i = 1; i < statement.path.size(); ++i) {
        full += "." + statement.path[i];
      }
      if (!statement.method.empty()) {
        full += ":" + statement.method;
      }
      statement.function = functionBody(!statement.method.empty(), statement.line, full);
      return statement;
    }
    case Tok::Local: {
      advance();
      if (accept(Tok::Function)) {
        statement.kind = Statement::Kind::LocalFunction;
        statement.names.push_back(name());
        statement.function = functionBody(false, statement.line, statement.names.front());
        return statement;
      }
      statement.kind = Statement::Kind::Local;
      do {
        statement.names.push_back(name());
        std::string attribute;
        if (accept(Tok::Less)) {
          attribute = name();
          expect(Tok::Greater, ">");
        }
        statement.attributes.push_back(attribute);
      } while (accept(Tok::Comma));
      if (accept(Tok::Assign)) {
        expressionList(statement.values);
      }
      return statement;
    }
    case Tok::DoubleColon: {
      advance();
      statement.kind = Statement::Kind::Label;
      statement.names.push_back(name());
      expect(Tok::DoubleColon, "::");
      return statement;
    }
    case Tok::Return:
      return returnStatement();
    case Tok::Break:
      advance();
      statement.kind = Statement::Kind::Break;
      return statement;
    case Tok::Goto:
      advance();
      statement.kind = Statement::Kind::Goto;
      statement.names.push_back(name());
      return statement;
    default:
      return expressionStatement();
    }
  }

  Statement ifStatement() {
    Statement statement;
    statement.kind = Statement::Kind::If;
    statement.line = current_.line;
    do {
      advance(); // if or elseif
      statement.conditions.push_back(expression());
      expect(Tok::Then, "then");
      statement.blocks.push_back(std::make_shared<Block>());
      block(*statement.blocks.back());
    } while (current_.kind == Tok::Elseif);
    if (accept(Tok::Else)) {
      statement.otherwise = std::make_shared<Block>();
      block(*statement.otherwise);
    }
    expect(Tok::End, "end");
    return statement;
  }

  Statement forStatement() {
    Statement statement;
    statement.line = current_.line;
    advance();
    const std::string first = name();
    if (accept(Tok::Assign)) {
      statement.kind = Statement::Kind::NumericFor;
      statement.names.push_back(first);
      statement.start = expression();
      expect(Tok::Comma, ",");
      statement.stop = expression();
      if (accept(Tok::Comma)) {
        statement.step = expression();
      }
    } else {
      statement.kind = Statement::Kind::GenericFor;
      statement.names.push_back(first);
      while (accept(Tok::Comma)) {
        statement.names.push_back(name());
      }
      expect(Tok::In, "in");
      expressionList(statement.values);
    }
    expect(Tok::Do, "do");
    statement.body = std::make_shared<Block>();
    block(*statement.body);
    expect(Tok::End, "end");
    return statement;
  }

  Statement expressionStatement() {
    Statement statement;
    statement.line = current_.line;
    ExpressionPtr first = suffixedExpression();
    if (current_.kind == Tok::Assign || current_.kind == Tok::Comma) {
      statement.kind = Statement::Kind::Assign;
      statement.targets.push_back(std::move(first));
      while (accept(Tok::Comma)) {
        statement.targets.push_back(suffixedExpression());
      }
      expect(Tok::Assign, "=");
      expressionList(statement.values);
      for (const auto& target : statement.targets) {
        if (target->kind != Expression::Kind::Name && target->kind != Expression::Kind::Index) {
          throw LuaError(statement.line, "syntax error near '='");
        }
      }
      return statement;
    }
    if (first->kind != Expression::Kind::Call && first->kind != Expression::Kind::Method) {
      throw LuaError(statement.line, "syntax error near '" + describe(current_) + "'");
    }
    statement.kind = Statement::Kind::Call;
    statement.call = std::move(first);
    return statement;
  }

  void expressionList(std::vector<ExpressionPtr>& out) {
    out.push_back(expression());
    while (accept(Tok::Comma)) {
      out.push_back(expression());
    }
  }

  std::shared_ptr<FunctionBody> functionBody(bool method, int line, std::string name) {
    auto body = std::make_shared<FunctionBody>();
    body->line = line;
    body->name = std::move(name);
    if (method) {
      body->parameters.emplace_back("self");
    }
    expect(Tok::OpenParen, "(");
    if (current_.kind != Tok::CloseParen) {
      do {
        if (accept(Tok::Dots)) {
          body->vararg = true;
          break;
        }
        body->parameters.push_back(this->name());
      } while (accept(Tok::Comma));
    }
    expect(Tok::CloseParen, ")");
    block(body->block);
    expect(Tok::End, "end");
    return body;
  }

  ExpressionPtr make(Expression::Kind kind, int line) {
    auto expression = std::make_unique<Expression>();
    expression->kind = kind;
    expression->line = line;
    return expression;
  }

  ExpressionPtr primaryExpression() {
    const int line = current_.line;
    if (current_.kind == Tok::Name) {
      auto expression = make(Expression::Kind::Name, line);
      expression->text = name();
      return expression;
    }
    if (accept(Tok::OpenParen)) {
      auto expression = make(Expression::Kind::Paren, line);
      expression->a = this->expression();
      expect(Tok::CloseParen, ")");
      return expression;
    }
    throw LuaError(line, "unexpected symbol near '" + describe(current_) + "'");
  }

  ExpressionPtr suffixedExpression() {
    ExpressionPtr expression = primaryExpression();
    while (true) {
      const int line = current_.line;
      switch (current_.kind) {
      case Tok::Dot: {
        advance();
        auto key = make(Expression::Kind::String, line);
        key->text = name();
        auto index = make(Expression::Kind::Index, line);
        index->a = std::move(expression);
        index->b = std::move(key);
        expression = std::move(index);
        break;
      }
      case Tok::OpenBracket: {
        advance();
        auto index = make(Expression::Kind::Index, line);
        index->a = std::move(expression);
        index->b = this->expression();
        expect(Tok::CloseBracket, "]");
        expression = std::move(index);
        break;
      }
      case Tok::Colon: {
        advance();
        auto call = make(Expression::Kind::Method, line);
        call->a = std::move(expression);
        call->text = name();
        callArguments(call->args);
        expression = std::move(call);
        break;
      }
      case Tok::OpenParen:
      case Tok::String:
      case Tok::OpenBrace: {
        auto call = make(Expression::Kind::Call, line);
        call->a = std::move(expression);
        callArguments(call->args);
        expression = std::move(call);
        break;
      }
      default:
        return expression;
      }
    }
  }

  void callArguments(std::vector<ExpressionPtr>& out) {
    const int line = current_.line;
    if (current_.kind == Tok::String) {
      auto text = make(Expression::Kind::String, line);
      text->text = current_.text;
      advance();
      out.push_back(std::move(text));
      return;
    }
    if (current_.kind == Tok::OpenBrace) {
      out.push_back(tableConstructor());
      return;
    }
    expect(Tok::OpenParen, "(");
    if (current_.kind != Tok::CloseParen) {
      expressionList(out);
    }
    expect(Tok::CloseParen, ")");
  }

  ExpressionPtr tableConstructor() {
    auto table = make(Expression::Kind::Table, current_.line);
    expect(Tok::OpenBrace, "{");
    while (current_.kind != Tok::CloseBrace) {
      Expression::Field field;
      if (current_.kind == Tok::Name && ahead_.kind == Tok::Assign) {
        auto key = make(Expression::Kind::String, current_.line);
        key->text = name();
        advance(); // =
        field.key = std::move(key);
        field.value = expression();
      } else if (accept(Tok::OpenBracket)) {
        field.key = expression();
        expect(Tok::CloseBracket, "]");
        expect(Tok::Assign, "=");
        field.value = expression();
      } else {
        field.value = expression();
      }
      table->fields.push_back(std::move(field));
      if (!accept(Tok::Comma) && !accept(Tok::Semicolon)) {
        break;
      }
    }
    expect(Tok::CloseBrace, "}");
    return table;
  }

  ExpressionPtr simpleExpression() {
    const int line = current_.line;
    switch (current_.kind) {
    case Tok::Number: {
      auto number = make(Expression::Kind::Number, line);
      number->number = current_.number;
      number->integer = current_.integer;
      number->isInteger = current_.isInteger;
      number->text = current_.text;
      advance();
      return number;
    }
    case Tok::String: {
      auto text = make(Expression::Kind::String, line);
      text->text = current_.text;
      advance();
      return text;
    }
    case Tok::Nil:
      advance();
      return make(Expression::Kind::Nil, line);
    case Tok::True:
      advance();
      return make(Expression::Kind::True, line);
    case Tok::False:
      advance();
      return make(Expression::Kind::False, line);
    case Tok::Dots:
      advance();
      return make(Expression::Kind::Vararg, line);
    case Tok::OpenBrace:
      return tableConstructor();
    case Tok::Function: {
      advance();
      auto function = make(Expression::Kind::Function, line);
      function->function = functionBody(false, line, "anonymous function");
      return function;
    }
    default:
      return suffixedExpression();
    }
  }

  struct Priority {
    int left;
    int right;
  };

  /// The binary operator at the current token, and its priorities (lparser.c).
  static bool binary(Tok kind, Expression::Op& op, Priority& priority) {
    using Op = Expression::Op;
    switch (kind) {
    case Tok::Plus:
      op = Op::Add;
      priority = {10, 10};
      return true;
    case Tok::Minus:
      op = Op::Sub;
      priority = {10, 10};
      return true;
    case Tok::Star:
      op = Op::Mul;
      priority = {11, 11};
      return true;
    case Tok::Percent:
      op = Op::Mod;
      priority = {11, 11};
      return true;
    case Tok::Caret:
      op = Op::Pow;
      priority = {14, 13};
      return true;
    case Tok::Slash:
      op = Op::Div;
      priority = {11, 11};
      return true;
    case Tok::DoubleSlash:
      op = Op::IntDiv;
      priority = {11, 11};
      return true;
    case Tok::Ampersand:
      op = Op::BitAnd;
      priority = {6, 6};
      return true;
    case Tok::Bar:
      op = Op::BitOr;
      priority = {4, 4};
      return true;
    case Tok::Tilde:
      op = Op::BitXor;
      priority = {5, 5};
      return true;
    case Tok::ShiftLeft:
      op = Op::Shl;
      priority = {7, 7};
      return true;
    case Tok::ShiftRight:
      op = Op::Shr;
      priority = {7, 7};
      return true;
    case Tok::Concat:
      op = Op::Concat;
      priority = {9, 8};
      return true;
    case Tok::NotEqual:
      op = Op::Ne;
      priority = {3, 3};
      return true;
    case Tok::Equal:
      op = Op::Eq;
      priority = {3, 3};
      return true;
    case Tok::Less:
      op = Op::Lt;
      priority = {3, 3};
      return true;
    case Tok::LessEqual:
      op = Op::Le;
      priority = {3, 3};
      return true;
    case Tok::Greater:
      op = Op::Gt;
      priority = {3, 3};
      return true;
    case Tok::GreaterEqual:
      op = Op::Ge;
      priority = {3, 3};
      return true;
    case Tok::And:
      op = Op::And;
      priority = {2, 2};
      return true;
    case Tok::Or:
      op = Op::Or;
      priority = {1, 1};
      return true;
    default:
      return false;
    }
  }

  static constexpr int UnaryPriority = 12;

  /// subexpr of lparser.c: an expression whose binary operators bind tighter than `limit`.
  ExpressionPtr subExpression(int limit) {
    ExpressionPtr left;
    const int line = current_.line;
    Expression::Op unary{};
    bool isUnary = true;
    switch (current_.kind) {
    case Tok::Not:
      unary = Expression::Op::Not;
      break;
    case Tok::Minus:
      unary = Expression::Op::Negate;
      break;
    case Tok::Hash:
      unary = Expression::Op::Length;
      break;
    case Tok::Tilde:
      unary = Expression::Op::BitNot;
      break;
    default:
      isUnary = false;
    }
    if (isUnary) {
      advance();
      left = make(Expression::Kind::Unary, line);
      left->op = unary;
      left->a = subExpression(UnaryPriority);
    } else {
      left = simpleExpression();
    }
    Expression::Op op{};
    Priority priority{};
    while (binary(current_.kind, op, priority) && priority.left > limit) {
      const int opLine = current_.line;
      advance();
      auto node = make(Expression::Kind::Binary, opLine);
      node->op = op;
      node->a = std::move(left);
      node->b = subExpression(priority.right);
      left = std::move(node);
    }
    return left;
  }

  ExpressionPtr expression() { return subExpression(0); }

  Lexer lexer_;
  Token current_;
  Token ahead_;
};

} // namespace

std::shared_ptr<FunctionBody> parseChunk(const std::string& source) {
  return Parser(source).chunk();
}

} // namespace seissol::reader::scripting::lua
