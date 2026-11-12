// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_LUASYNTAX_H_
#define SEISSOL_SRC_READER_SCRIPTING_LUASYNTAX_H_

// The syntax of Lua 5.4, as its reference manual has it (§9), read into a tree: the frontend of
// Lua (LuaFrontend.h) evaluates the tree, without a Lua VM. Every construct of the language is
// read -- what the frontend does not evaluate it refuses with its line, where it would be
// evaluated -- so that a script is either read whole or not at all.

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::reader::scripting::lua {

/// A script that is no Lua, or that fails at every point it is evaluated at, at a line of it.
class LuaError : public std::runtime_error {
  public:
  /// At `line`, or at none (0).
  LuaError(int line, const std::string& what)
      : std::runtime_error(line > 0 ? "line " + std::to_string(line) + ": " + what : what),
        line_(line), message_(what) {}
  [[nodiscard]] int line() const { return line_; }
  /// The message, without the line.
  [[nodiscard]] const std::string& message() const { return message_; }

  private:
  int line_;
  std::string message_;
};

struct Block;
struct FunctionBody;

/// An expression.
struct Expression {
  enum class Kind : std::uint8_t {
    Nil,
    True,
    False,
    Number,
    String,
    Vararg,
    Function,
    Table,
    Binary,
    Unary,
    Name,
    Index,
    Call,
    Method,
    Paren
  };
  /// The binary operators, and the unary ones (Negate, Not, Length, BitNot).
  enum class Op : std::uint8_t {
    Add,
    Sub,
    Mul,
    Div,
    IntDiv,
    Mod,
    Pow,
    Concat,
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    And,
    Or,
    BitAnd,
    BitOr,
    BitXor,
    Shl,
    Shr,
    Negate,
    Not,
    Length,
    BitNot
  };
  /// A field of a table constructor: `[key] = value`, `name = value` (a String key), or
  /// `value` (positional, key null).
  struct Field {
    std::unique_ptr<Expression> key;
    std::unique_ptr<Expression> value;
  };

  Kind kind{Kind::Nil};
  int line{0};
  /// Number: its value, and whether the literal is an integer of Lua
  double number{0.0};
  std::int64_t integer{0};
  bool isInteger{false};
  /// String, Name, and the name of a Method
  std::string text;
  Op op{Op::Add};
  /// Binary: a and b; Unary, Paren: a; Index: a[b]; Call: a(args); Method: a:text(args)
  std::unique_ptr<Expression> a;
  std::unique_ptr<Expression> b;
  std::vector<std::unique_ptr<Expression>> args;
  std::vector<Field> fields;
  std::shared_ptr<FunctionBody> function;
};

using ExpressionPtr = std::unique_ptr<Expression>;

/// A statement.
struct Statement {
  enum class Kind : std::uint8_t {
    Local,
    Assign,
    Call,
    Do,
    While,
    Repeat,
    If,
    NumericFor,
    GenericFor,
    Function,
    LocalFunction,
    Return,
    Break,
    Goto,
    Label,
    Empty
  };
  Kind kind{Kind::Empty};
  int line{0};
  /// Local: names (with their attributes `<const>` or `<close>` in `attributes`); NumericFor:
  /// the variable; GenericFor: the variables; LocalFunction: the name; Goto, Label: the label
  std::vector<std::string> names;
  std::vector<std::string> attributes;
  /// Assign: the targets (names, indexings); Local, Assign, Return, GenericFor: the expressions;
  /// Call: the call; NumericFor: start, stop, step (null without)
  std::vector<ExpressionPtr> targets;
  std::vector<ExpressionPtr> values;
  ExpressionPtr call;
  /// If: the conditions and their blocks, and the block of `else` (null without); While,
  /// Repeat: the condition
  std::vector<ExpressionPtr> conditions;
  std::vector<std::shared_ptr<Block>> blocks;
  std::shared_ptr<Block> otherwise;
  ExpressionPtr start;
  ExpressionPtr stop;
  ExpressionPtr step;
  /// Do, While, Repeat, NumericFor, GenericFor: the body
  std::shared_ptr<Block> body;
  /// Function: the name `a.b.c` (and `:m` as `method`), the function
  std::vector<std::string> path;
  std::string method;
  std::shared_ptr<FunctionBody> function;
};

struct Block {
  std::vector<Statement> statements;
};

struct FunctionBody {
  std::vector<std::string> parameters;
  bool vararg{false};
  Block block;
  int line{0};
  std::string name;
};

/// The chunk `source`: its block, as the body of a vararg function. Throws LuaError.
[[nodiscard]] std::shared_ptr<FunctionBody> parseChunk(const std::string& source);

} // namespace seissol::reader::scripting::lua

#endif // SEISSOL_SRC_READER_SCRIPTING_LUASYNTAX_H_
