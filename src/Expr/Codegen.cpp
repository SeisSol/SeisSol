// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Expr/Codegen.h"

#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr::codegen {

namespace {

bool identifierChar(char c) {
  return (std::isalnum(static_cast<unsigned char>(c)) != 0) || c == '_';
}

/// The type of transient `slot` of `stage`.
ComputeType slotType(const StageCode& stage, std::int32_t slot) {
  return static_cast<std::size_t>(slot) < stage.slotTypes.size()
             ? stage.slotTypes[static_cast<std::size_t>(slot)]
             : ComputeType::F64;
}

/// `value`, of type `from`, as an expression of type `to`: itself where the two agree, so that a
/// program in one type gets the same source as before there were types.
std::string converted(const std::string& value, ComputeType from, ComputeType to) {
  return from == to ? value : "(" + std::string(typeName(to)) + ")(" + value + ")";
}

/// The text of a pointwise operation: that of the interpreter's table, or, computed as written
/// (Instruction::strict), the macro of its operation (emitStrictMacros()); with the operand `a0`
/// of type `first`. The operations of the lowering's types are spelled out: a comparison of type
/// Bool as the predicate of its operands, `and` and `or` of Bools, a select with a Bool condition,
/// and the minimum and maximum of whole numbers, which the table spells as fmin and fmax -- in
/// binary64.
std::string pwText(const Instruction& inst,
                   const std::string& a0,
                   const std::string& a1,
                   const std::string& a2,
                   const std::string& type,
                   ComputeType first,
                   MathStyle style) {
  if (inst.type == ComputeType::Bool) {
    switch (inst.fn) {
    case Fn::Lt:
      return "(" + a0 + " < " + a1 + ")";
    case Fn::Le:
      return "(" + a0 + " <= " + a1 + ")";
    case Fn::Eq:
      return "(" + a0 + " == " + a1 + ")";
    case Fn::And:
      return "(" + a0 + " && " + a1 + ")";
    case Fn::Or:
      return "(" + a0 + " || " + a1 + ")";
    default:
      break;
    }
  }
  if (inst.fn == Fn::Select && first == ComputeType::Bool) {
    return "(" + a0 + " ? " + a1 + " : " + a2 + ")";
  }
  if (inst.fn == Fn::Select && inst.type == ComputeType::Bool) {
    // truth values chosen by a condition of another type: where the check of a lookup counts
    // (Branches.h)
    return "((" + a0 + " != 0) ? " + a1 + " : " + a2 + ")";
  }
  if (inst.type == ComputeType::I64 && (inst.fn == Fn::Min || inst.fn == Fn::Max)) {
    const char* compare = inst.fn == Fn::Min ? " < " : " > ";
    return "(" + a0 + compare + a1 + " ? " + a0 + " : " + a1 + ")";
  }
  if (inst.strict && inst.type == ComputeType::F64) {
    const char* macro = inst.fn == Fn::Add   ? "SEISSOL_EXPR_DADD"
                        : inst.fn == Fn::Sub ? "SEISSOL_EXPR_DSUB"
                        : inst.fn == Fn::Mul ? "SEISSOL_EXPR_DMUL"
                                             : nullptr;
    if (macro != nullptr) {
      return std::string(macro) + "(" + a0 + ", " + a1 + ")";
    }
  }
  return substitute(expressionText(inst.fn), a0, a1, a2, type, style);
}

/// The value of a GridParam or GridLoad instruction, in the type of the instruction, with `operand`
/// the text of the operand of a GridLoad.
std::string gridText(const Instruction& inst,
                     const StageAddressing& addressing,
                     const std::string& operand) {
  if (addressing.gridParam == nullptr || addressing.gridLoad == nullptr) {
    throw std::invalid_argument(
        "expr: a grid lookup reached a code generator that does not sample grids");
  }
  if (inst.op == Opcode::GridParam) {
    // the field as the descriptor holds it, in the type of the instruction
    return "(" + std::string(typeName(inst.type)) + ")(" +
           addressing.gridParam(inst.grid, static_cast<GridField>(inst.comp), inst.imm) + ")";
  }
  return converted(addressing.gridLoad(inst.grid, inst.comp, operand), ComputeType::F64, inst.type);
}

/// Whether the condition of `branch` takes it, of the slot `slot` (a lane of it, in the tiled
/// code): a condition holds where it is not zero, a NaN included.
std::string takes(const Branch& branch, const std::string& slot) {
  return "(" + slot + (branch.side ? " != 0" : " == 0") + ")";
}

void emitContraction(std::ostringstream& out,
                     const Instruction& inst,
                     const std::string& dst,
                     const std::string& computeType,
                     const ContractAddressing& contract,
                     const char* indent) {
  const auto& shape = contract.shapes.at(static_cast<std::size_t>(inst.matrix));
  const auto& matrixType = contract.matrixTypes.at(static_cast<std::size_t>(inst.matrix));
  const auto& blockType = contract.blockTypes.at(static_cast<std::size_t>(inst.block));
  const auto& index = contract.indexType;
  const auto& global = contract.global;
  const std::string inner = std::string(indent) + "  ";

  // C-style casts and no initializer syntax, so that the text is valid C++ and OpenCL C alike.
  out << indent << "{\n"
      << inner << "const " << index << " c_cell = (" << contract.point << ") / " << shape.rows
      << ";\n"
      << inner << "const " << index << " c_row = (" << contract.point << ") - c_cell * "
      << shape.rows << ";\n"
      << inner << global << "const unsigned int* c_index = " << contract.cellIndex(inst.block)
      << ";\n"
      << inner << global << "const char* c_block = (" << global << "const char*)("
      << contract.blockBase(inst.block) << ") + (c_index != 0 ? (" << index
      << ")c_index[c_cell] : c_cell) * (" << contract.cellStride(inst.block) << ");\n"
      << inner << "const " << index << " c_step = " << contract.modeStride(inst.block) << ";\n"
      << inner << global << "const " << matrixType << "* c_matrix = (" << global << "const "
      << matrixType << "*)(" << contract.matrix(inst.matrix) << ") + c_row;\n"
      << inner << computeType << " c_sum = (" << computeType << ")0;\n"
      << inner << "for (" << index << " c_m = 0; c_m < " << shape.cols << "; ++c_m) {\n"
      << inner << "  c_sum += (" << computeType << ")c_matrix[c_m * " << shape.leadingDimension
      << "] * (" << computeType << ")(*(" << global << "const " << blockType
      << "*)(c_block + c_m * c_step));\n"
      << inner << "}\n"
      << inner << dst << " = c_sum;\n"
      << indent << "}\n";
}

/// A percentile: its values gathered into an array, which the function of its type
/// (emitPercentileFunctions()) reorders.
void emitPercentile(std::ostringstream& out,
                    const Instruction& inst,
                    const std::vector<std::int32_t>& operands,
                    const std::string& dst,
                    const std::function<std::string(std::int32_t)>& slot,
                    const std::string& indent) {
  const std::string type = typeName(inst.type);
  const std::int32_t count = inst.operandCount - 1;
  out << indent << "{\n" << indent << "  " << type << " values[" << count << "] = {";
  for (std::int32_t k = 0; k < count; ++k) {
    if (k > 0) {
      out << (k % 8 == 0 ? ",\n" + indent + "      " : std::string(", "));
    }
    out << slot(operands[inst.operandBegin + 1 + k]);
  }
  out << "};\n"
      << indent << "  " << dst << " = " << percentileFunction(inst.type) << "("
      << slot(operands[inst.operandBegin]) << ", values, " << count << ");\n"
      << indent << "}\n";
}

} // namespace

std::string percentileFunction(ComputeType type) {
  return std::string("seissol_expr_percentile_") + (type == ComputeType::F32 ? "f32" : "f64");
}

void emitPercentileFunctions(std::ostringstream& out,
                             const LoweredProgram& lowered,
                             const PercentileArithmetic& arithmetic) {
  std::vector<ComputeType> types;
  for (const StageCode* stage : {&lowered.precompute(), &lowered.run()}) {
    for (const Instruction& inst : stage->code) {
      if (inst.op == Opcode::Percentile &&
          std::find(types.begin(), types.end(), inst.type) == types.end()) {
        types.push_back(inst.type);
      }
    }
  }
  const auto op = [](const std::string& text, const std::string& x, const std::string& y) {
    return substitute(text, x, y, "", "double", MathStyle::Namespaced);
  };
  const std::string add = arithmetic.add;
  const std::string sub = arithmetic.sub;
  const std::string mul = arithmetic.mul;
  for (const ComputeType type : types) {
    const std::string t = typeName(type);
    // a negative zero before a positive one, which compare equal: by their reciprocals
    const auto before = [&](const std::string& x, const std::string& y) {
      return "(" + x + " < " + y + " || (" + x + " == " + y + " && (" + t + ")1 / " + x + " < (" +
             t + ")1 / " + y + "))";
    };
    const std::string infinite = " > 1.7976931348623157e308";
    out << "/* the percentile p of the n values at v, which it reorders (Kind::Percentile) */\n"
        << arithmetic.qualifier << " " << t << " " << percentileFunction(type) << "(" << t << " p, "
        << t << "* v, int n) {\n";
    if (!arithmetic.prologue.empty()) {
      out << arithmetic.prologue << "\n";
    }
    out << "  for (int i = 0; i < n; ++i) {\n"
        << "    if (v[i] != v[i]) {\n"
        << "      return (" << t << ")" << arithmetic.nan << ";\n"
        << "    }\n"
        << "  }\n"
        << "  const double at = (double)p;\n"
        << "  if (!(at >= 0.0 && at <= 100.0)) {\n"
        << "    return (" << t << ")" << arithmetic.nan << ";\n"
        << "  }\n"
        << "  const double q = " << op(arithmetic.div, "at", "100.0") << ";\n"
        << "  const double position = " << op(mul, "(double)(n - 1)", "q") << ";\n"
        << "  const double below = " << arithmetic.floor << "(position);\n"
        << "  const int lower = (int)below;\n"
        << "  const double t = " << op(sub, "position", "below") << ";\n"
        << "  int left = 0;\n"
        << "  int right = n - 1;\n"
        << "  while (left < right) {\n"
        << "    const " << t << " pivot = v[left + (right - left) / 2];\n"
        << "    int i = left;\n"
        << "    int j = right;\n"
        << "    while (i <= j) {\n"
        << "      while (" << before("v[i]", "pivot") << ") {\n"
        << "        ++i;\n"
        << "      }\n"
        << "      while (" << before("pivot", "v[j]") << ") {\n"
        << "        --j;\n"
        << "      }\n"
        << "      if (i <= j) {\n"
        << "        const " << t << " swap = v[i];\n"
        << "        v[i] = v[j];\n"
        << "        v[j] = swap;\n"
        << "        ++i;\n"
        << "        --j;\n"
        << "      }\n"
        << "    }\n"
        << "    if (lower <= j) {\n"
        << "      right = j;\n"
        << "    } else if (lower >= i) {\n"
        << "      left = i;\n"
        << "    } else {\n"
        << "      break;\n"
        << "    }\n"
        << "  }\n"
        << "  const " << t << " a = v[lower];\n"
        << "  " << t << " b = a;\n"
        << "  if (lower + 1 < n) {\n"
        << "    b = v[lower + 1];\n"
        << "    for (int k = lower + 2; k < n; ++k) {\n"
        << "      if (" << before("v[k]", "b") << ") {\n"
        << "        b = v[k];\n"
        << "      }\n"
        << "    }\n"
        << "  }\n"
        << "  const double x = (double)a;\n"
        << "  const double y = (double)b;\n"
        << "  double r;\n"
        << "  if (" << arithmetic.fabs << "(x)" << infinite << " || " << arithmetic.fabs << "(y)"
        << infinite << ") {\n"
        << "    r = (t == 0.0 || x == y) ? x : " << op(add, "x", "y") << ";\n"
        << "  } else {\n"
        << "    const double d = " << op(sub, "y", "x") << ";\n"
        << "    r = t < 0.5 ? " << op(add, "x", op(mul, "d", "t")) << " : "
        << op(sub, "y", op(mul, "d", op(sub, "1.0", "t"))) << ";\n"
        << "  }\n"
        << "  return (" << t << ")r;\n"
        << "}\n\n";
  }
}

const char* expressionText(Fn fn) {
  switch (fn) {
#define SEISSOL_EXPR_TEXT_UNARY(NAME, EXPR)                                                        \
  case Fn::NAME:                                                                                   \
    return #EXPR;
#define SEISSOL_EXPR_TEXT_BINARY(NAME, EXPR)                                                       \
  case Fn::NAME:                                                                                   \
    return #EXPR;
#define SEISSOL_EXPR_TEXT_TERNARY(NAME, EXPR)                                                      \
  case Fn::NAME:                                                                                   \
    return #EXPR;
    SEISSOL_EXPR_PW_LIST(
        SEISSOL_EXPR_TEXT_UNARY, SEISSOL_EXPR_TEXT_BINARY, SEISSOL_EXPR_TEXT_TERNARY)
#undef SEISSOL_EXPR_TEXT_UNARY
#undef SEISSOL_EXPR_TEXT_BINARY
#undef SEISSOL_EXPR_TEXT_TERNARY
  }
  return nullptr;
}

std::string substitute(const std::string& text,
                       const std::string& x,
                       const std::string& y,
                       const std::string& z,
                       const std::string& computeType,
                       MathStyle style) {
  std::string out;
  out.reserve(text.size() * 2);
  std::size_t i = 0;
  while (i < text.size()) {
    if (!identifierChar(text[i])) {
      out.push_back(text[i]);
      ++i;
      continue;
    }
    const std::size_t begin = i;
    while (i < text.size() && identifierChar(text[i])) {
      ++i;
    }
    const std::string token = text.substr(begin, i - begin);

    if (token == "std" && style != MathStyle::Namespaced && text.compare(i, 2, "::") == 0) {
      // Device code has the maths built-ins unqualified, and NVRTC has no
      // <cmath> to bring `std` into scope. Dropping the qualifier here rather
      // than keeping a second expression table is the entire point of routing
      // both backends through one substitution.
      i += 2;
      // the function after the qualifier: the absolute value of a floating-point number is fabs,
      // in C, and in OpenCL C, whose abs takes integers only; by its C name for the compute type
      std::size_t end = i;
      while (end < text.size() && identifierChar(text[end])) {
        ++end;
      }
      std::string function = text.substr(i, end - i);
      if (function == "abs") {
        function = "fabs";
      }
      if (style == MathStyle::CNames && computeType == "float") {
        function += 'f';
      }
      out += function;
      i = end;
      continue;
    }
    if (token == "x") {
      out += x;
    } else if (token == "y") {
      out += y;
    } else if (token == "z") {
      out += z;
    } else if (token == "T") {
      // T(v) is a functional cast, which C++ has and OpenCL C does not; (T)(v) is the same cast in
      // both
      if (i < text.size() && text[i] == '(') {
        out += "(" + computeType + ")";
      } else {
        out += computeType;
      }
    } else {
      out += token;
    }
  }
  return out;
}

std::string slotName(std::int32_t slot) { return "s" + std::to_string(slot); }

void dump(const std::string& name, const std::string& text) {
  static const char* directory = std::getenv("SEISSOL_EXPR_DUMP");
  if (directory == nullptr || *directory == '\0') {
    return;
  }
  const std::string path = std::string(directory) + "/" + name;
  std::FILE* file = std::fopen(path.c_str(), "w");
  if (file == nullptr) {
    logWarning() << "expr: could not write" << path.c_str() << "(SEISSOL_EXPR_DUMP).";
    return;
  }
  std::fwrite(text.data(), 1, text.size(), file);
  std::fclose(file);
}

std::string literal(double value, const std::string& computeType) {
  // what %.17g writes of these is no C: `nan` and `inf` no numbers, `-0` the whole number 0
  if (std::isnan(value)) {
    return "((" + computeType + ")(0.0 / 0.0))";
  }
  if (std::isinf(value)) {
    return "((" + computeType + ")(" + (value < 0 ? "-" : "") + "(1.0 / 0.0)))";
  }
  if (value == 0.0 && std::signbit(value)) {
    return "((" + computeType + ")-0.0)";
  }
  std::array<char, 48> buffer{};
  std::snprintf(buffer.data(), buffer.size(), "%.17g", value);
  return "((" + computeType + ")" + buffer.data() + ")";
}

void emitStageBody(std::ostringstream& out,
                   const StageCode& stage,
                   const std::vector<std::int32_t>& operands,
                   ComputeType dataType,
                   MathStyle style,
                   const StageAddressing& addressing,
                   const char* indent) {
  for (std::int32_t slot = 0; slot < stage.slotCount; ++slot) {
    const char* type = typeName(slotType(stage, slot));
    out << indent << type << " " << slotName(slot) << " = (" << type << ")0;\n";
  }
  if (stage.slotCount > 0) {
    // Silences the unused-variable warning for a slot that is written but never
    // read, which liveness can legitimately produce.
    out << indent << "(void)" << slotName(0) << ";\n";
  }

  // The writes run before the instruction at their positions, never in a branch (Store::position).
  const auto value = [&](const Store& store) {
    return converted(slotName(store.source), slotType(stage, store.source), dataType);
  };
  const auto writesAt = [&](std::size_t index) {
    return std::any_of(stage.outputs.begin(),
                       stage.outputs.end(),
                       [&](const Store& store) {
                         return store.position == static_cast<std::int32_t>(index);
                       }) ||
           std::any_of(stage.persistent.begin(), stage.persistent.end(), [&](const Store& store) {
             return store.position == static_cast<std::int32_t>(index);
           });
  };
  const auto writeAt = [&](std::size_t index) {
    for (const Store& store : stage.outputs) {
      if (store.position == static_cast<std::int32_t>(index)) {
        out << indent << addressing.storeOutput(store.target, value(store)) << ";\n";
      }
    }
    for (const Store& store : stage.persistent) {
      if (store.position != static_cast<std::int32_t>(index)) {
        continue;
      }
      if (store.target < addressing.stateSlots) {
        if (addressing.storeState == nullptr) {
          throw std::invalid_argument(
              "expr: a state reached a code generator that cannot address states");
        }
        out << indent << addressing.storeState(store.target, value(store)) << ";\n";
      } else {
        out << indent
            << addressing.storePersistent(store.target - addressing.stateSlots, value(store))
            << ";\n";
      }
    }
  };

  // A branch is an `if` around its part of the code: the points of a GPU warp that do not take it
  // wait, those that do run it, and a CPU skips it where the point does not take it. Two parts of
  // the same condition, one for either side, become an `if` and its `else`, where no write lies
  // between them.
  std::string current = indent;
  std::vector<Branch> open;
  std::size_t nextBranch = 0;
  const auto& branches = stage.branches;
  for (std::size_t index = 0; index < stage.code.size(); ++index) {
    writeAt(index);
    while (nextBranch < branches.size() &&
           branches[nextBranch].begin == static_cast<std::int32_t>(index)) {
      const Branch& branch = branches[nextBranch++];
      if (branch.begin == branch.end) {
        continue;
      }
      out << current << "if " << takes(branch, slotName(branch.condition)) << " {\n";
      current += "  ";
      open.push_back(branch);
    }
    const Instruction& inst = stage.code[index];
    const char* const pad = current.c_str();
    const std::string dst = slotName(inst.dst);
    const std::string type = typeName(inst.type);
    switch (inst.op) {
    case Opcode::Const:
      out << pad << dst << " = " << literal(inst.value, type) << ";\n";
      break;
    case Opcode::LoadInput:
      out << pad << dst << " = " << converted(addressing.loadInput(inst.imm), dataType, inst.type)
          << ";\n";
      break;
    case Opcode::LoadPersistent:
      if (inst.imm < addressing.stateSlots) {
        if (addressing.loadState == nullptr) {
          throw std::invalid_argument(
              "expr: a state reached a code generator that cannot address states");
        }
        out << pad << dst << " = " << converted(addressing.loadState(inst.imm), dataType, inst.type)
            << ";\n";
      } else {
        out << pad << dst << " = "
            << converted(
                   addressing.loadPersistent(inst.imm - addressing.stateSlots), dataType, inst.type)
            << ";\n";
      }
      break;
    case Opcode::Pw: {
      const std::string a0 = slotName(operands[inst.operandBegin]);
      // Unary and binary ops never read y or z, but the interpreter aliases
      // them onto a0, so mirroring that keeps the emitted text well formed for
      // every arity without a special case per arity.
      const std::string a1 = inst.operandCount > 1 ? slotName(operands[inst.operandBegin + 1]) : a0;
      const std::string a2 = inst.operandCount > 2 ? slotName(operands[inst.operandBegin + 2]) : a0;
      out << pad << dst << " = "
          << pwText(inst, a0, a1, a2, type, slotType(stage, operands[inst.operandBegin]), style)
          << ";\n";
      break;
    }
    case Opcode::Convert: {
      const std::int32_t source = operands[inst.operandBegin];
      out << pad << dst << " = " << converted(slotName(source), slotType(stage, source), inst.type)
          << ";\n";
      break;
    }
    case Opcode::GridParam:
    case Opcode::GridLoad: {
      std::string operand;
      if (inst.op == Opcode::GridLoad) {
        const std::int32_t source = operands[inst.operandBegin];
        operand = converted(slotName(source), slotType(stage, source), ComputeType::I64);
      }
      out << pad << dst << " = " << gridText(inst, addressing, operand) << ";\n";
      break;
    }
    case Opcode::Contract:
      if (addressing.contract == nullptr) {
        throw std::invalid_argument(
            "expr: a contraction reached a code generator that cannot address blocks");
      }
      emitContraction(out, inst, dst, type, *addressing.contract, pad);
      break;
    case Opcode::Percentile:
      emitPercentile(out, inst, operands, dst, slotName, current);
      break;
    }
    const auto after = static_cast<std::int32_t>(index) + 1;
    while (!open.empty() && open.back().end == after) {
      const Branch closed = open.back();
      open.pop_back();
      current.resize(current.size() - 2);
      // the other side of the same condition next, at the same depth: its `else`
      const bool outermost = open.empty() || open.back().end != after;
      if (outermost && nextBranch < branches.size() && branches[nextBranch].begin == after &&
          branches[nextBranch].end > after && branches[nextBranch].condition == closed.condition &&
          branches[nextBranch].side != closed.side && !writesAt(static_cast<std::size_t>(after))) {
        out << current << "} else {\n";
        current += "  ";
        open.push_back(branches[nextBranch++]);
        break;
      }
      out << current << "}\n";
    }
  }

  writeAt(stage.code.size());
  for (const CheckStore& store : stage.checks) {
    if (addressing.storeCheck == nullptr) {
      throw std::invalid_argument(
          "expr: a grid lookup reached a code generator that does not sample grids");
    }
    out << indent << addressing.storeCheck(store.grid, store.kind, slotName(store.source)) << ";\n";
  }
}

void emitStageBodyTiled(std::ostringstream& out,
                        const StageCode& stage,
                        const std::vector<std::int32_t>& operands,
                        ComputeType dataType,
                        MathStyle style,
                        const StageAddressing& addressing,
                        std::size_t block,
                        const char* indent) {
  const std::string blockSize = std::to_string(block) + "ul";
  const std::string inner = std::string(indent) + "  ";
  std::string current = inner;
  const auto slot = [](std::int32_t index) { return slotName(index) + "[j]"; };
  // one instruction over the lanes of the block; `l` where the addressing reads the point
  const auto laneLoop = [&](const std::string& statement, bool readsPoint) {
    const std::string lane = current + "  ";
    out << "#pragma omp simd\n" << current << "for (unsigned long j = 0; j < n; ++j) {\n";
    if (readsPoint) {
      out << lane << "const unsigned long l = base + j;\n";
    }
    out << lane << statement << ";\n" << current << "}\n";
  };

  // With branches, the transients outlive a block: a select reads the lanes of the operand it does
  // not choose also where the block skipped the branch that computes it, which then hold the
  // values of an earlier block, or zero.
  const auto& branches = stage.branches;
  if (!branches.empty()) {
    for (std::int32_t index = 0; index < stage.slotCount; ++index) {
      out << indent << typeName(slotType(stage, index)) << " " << slotName(index) << "[" << block
          << "] = {0};\n";
    }
  }
  out << indent << "for (unsigned long base = 0; base < count; base += " << blockSize << ") {\n"
      << inner << "const unsigned long n = count - base < " << blockSize
      << " ? count - base : " << blockSize << ";\n";
  if (branches.empty()) {
    for (std::int32_t index = 0; index < stage.slotCount; ++index) {
      out << inner << typeName(slotType(stage, index)) << " " << slotName(index) << "[" << block
          << "];\n";
    }
  }

  // The writes run before the instruction at their positions, never in a branch (Store::position).
  const auto value = [&](const Store& store) {
    return converted(slot(store.source), slotType(stage, store.source), dataType);
  };
  const auto writeAt = [&](std::size_t index) {
    for (const Store& store : stage.outputs) {
      if (store.position == static_cast<std::int32_t>(index)) {
        laneLoop(addressing.storeOutput(store.target, value(store)), true);
      }
    }
    for (const Store& store : stage.persistent) {
      if (store.position != static_cast<std::int32_t>(index)) {
        continue;
      }
      if (store.target < addressing.stateSlots) {
        if (addressing.storeState == nullptr) {
          throw std::invalid_argument(
              "expr: a state reached a code generator that cannot address states");
        }
        laneLoop(addressing.storeState(store.target, value(store)), true);
      } else {
        laneLoop(addressing.storePersistent(store.target - addressing.stateSlots, value(store)),
                 true);
      }
    }
  };

  // A branch runs for the block where a lane of it takes the branch.
  std::vector<std::int32_t> ends;
  std::size_t nextBranch = 0;
  for (std::size_t index = 0; index < stage.code.size(); ++index) {
    writeAt(index);
    while (nextBranch < branches.size() &&
           branches[nextBranch].begin == static_cast<std::int32_t>(index)) {
      const std::string flag = "b" + std::to_string(nextBranch);
      const Branch& branch = branches[nextBranch++];
      if (branch.begin == branch.end) {
        continue;
      }
      out << current << "int " << flag << " = 0;\n"
          << "#pragma omp simd reduction(|:" << flag << ")\n"
          << current << "for (unsigned long j = 0; j < n; ++j) {\n"
          << current << "  " << flag << " |= " << takes(branch, slot(branch.condition)) << ";\n"
          << current << "}\n"
          << current << "if (" << flag << " != 0) {\n";
      current += "  ";
      ends.push_back(branch.end);
    }
    const Instruction& inst = stage.code[index];
    const std::string dst = slot(inst.dst);
    const std::string type = typeName(inst.type);
    switch (inst.op) {
    case Opcode::Const:
      laneLoop(dst + " = " + literal(inst.value, type), false);
      break;
    case Opcode::LoadInput:
      laneLoop(dst + " = " + converted(addressing.loadInput(inst.imm), dataType, inst.type), true);
      break;
    case Opcode::LoadPersistent:
      if (inst.imm < addressing.stateSlots) {
        if (addressing.loadState == nullptr) {
          throw std::invalid_argument(
              "expr: a state reached a code generator that cannot address states");
        }
        laneLoop(dst + " = " + converted(addressing.loadState(inst.imm), dataType, inst.type),
                 true);
      } else {
        laneLoop(dst + " = " +
                     converted(addressing.loadPersistent(inst.imm - addressing.stateSlots),
                               dataType,
                               inst.type),
                 true);
      }
      break;
    case Opcode::Pw: {
      const std::string a0 = slot(operands[inst.operandBegin]);
      const std::string a1 = inst.operandCount > 1 ? slot(operands[inst.operandBegin + 1]) : a0;
      const std::string a2 = inst.operandCount > 2 ? slot(operands[inst.operandBegin + 2]) : a0;
      laneLoop(
          dst + " = " +
              pwText(inst, a0, a1, a2, type, slotType(stage, operands[inst.operandBegin]), style),
          false);
      break;
    }
    case Opcode::Convert: {
      const std::int32_t source = operands[inst.operandBegin];
      laneLoop(dst + " = " + converted(slot(source), slotType(stage, source), inst.type), false);
      break;
    }
    case Opcode::GridParam:
    case Opcode::GridLoad: {
      std::string operand;
      if (inst.op == Opcode::GridLoad) {
        const std::int32_t source = operands[inst.operandBegin];
        operand = converted(slot(source), slotType(stage, source), ComputeType::I64);
      }
      laneLoop(dst + " = " + gridText(inst, addressing, operand), false);
      break;
    }
    case Opcode::Contract:
      throw std::invalid_argument(
          "expr: a contraction reached the tiled code generator; emit one point after the other");
    case Opcode::Percentile:
      // one lane after the other: a selection does not vectorize
      out << current << "for (unsigned long j = 0; j < n; ++j) {\n";
      emitPercentile(out, inst, operands, dst, slot, current + "  ");
      out << current << "}\n";
      break;
    }
    while (!ends.empty() && ends.back() == static_cast<std::int32_t>(index) + 1) {
      ends.pop_back();
      current.resize(current.size() - 2);
      out << current << "}\n";
    }
  }

  writeAt(stage.code.size());
  for (const CheckStore& store : stage.checks) {
    if (addressing.storeCheck == nullptr) {
      throw std::invalid_argument(
          "expr: a grid lookup reached a code generator that does not sample grids");
    }
    laneLoop(addressing.storeCheck(store.grid, store.kind, slot(store.source)), false);
  }
  out << indent << "}\n";
}

bool callsLibrary(Fn fn) {
  switch (fn) {
  case Fn::Exp:
  case Fn::Log:
  case Fn::Log2:
  case Fn::Log10:
  case Fn::Sin:
  case Fn::Cos:
  case Fn::Tan:
  case Fn::Asin:
  case Fn::Acos:
  case Fn::Atan:
  case Fn::Sinh:
  case Fn::Cosh:
  case Fn::Tanh:
  case Fn::Asinh:
  case Fn::Acosh:
  case Fn::Atanh:
  case Fn::Erf:
  case Fn::Pow:
  case Fn::Atan2:
  case Fn::Mod:
    return true;
  default:
    return false;
  }
}

std::int32_t criticalPath(const StageCode& stage, const std::vector<std::int32_t>& operands) {
  // the depth of the value each transient holds as the instructions go
  std::vector<std::int32_t> depth(static_cast<std::size_t>(std::max(stage.slotCount, 1)), 0);
  std::int32_t longest = 0;
  for (const Instruction& inst : stage.code) {
    std::int32_t before = 0;
    for (std::int32_t k = 0; k < inst.operandCount; ++k) {
      before = std::max(before, depth[static_cast<std::size_t>(operands[inst.operandBegin + k])]);
    }
    depth[static_cast<std::size_t>(inst.dst)] = before + 1;
    longest = std::max(longest, before + 1);
  }
  return longest;
}

GridDescriptor gridDescriptor(const reader::datafield::ArrayView& view, const void* data) {
  GridDescriptor descriptor;
  descriptor.data = static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(data));
  for (unsigned d = 0; d < view.dimensions && d < reader::datafield::MaxGridDimensions; ++d) {
    descriptor.min[d] = view.min[d];
    descriptor.deltaInv[d] = view.deltaInv[d];
    descriptor.num[d] = view.num[d];
    descriptor.stride[d] = view.stride[d];
    descriptor.before[d] = view.before[d];
    descriptor.after[d] = view.after[d];
    descriptor.knots[d] =
        static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(view.knots[d]));
    descriptor.inverseSpacings[d] =
        static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(view.inverseSpacings[d]));
  }
  descriptor.componentStride = view.componentStride;
  return descriptor;
}

std::size_t gridBytes(const reader::datafield::ArrayView& view) {
  // up to the last element: the far corner of the last component
  std::size_t last = (view.components > 0 ? view.components - 1 : 0) * view.componentStride;
  for (unsigned d = 0; d < view.dimensions && d < reader::datafield::MaxGridDimensions; ++d) {
    if (view.num[d] == 0) {
      return 0;
    }
    last += static_cast<std::size_t>(view.num[d] - 1) * view.stride[d];
  }
  const std::size_t element =
      view.type == reader::datafield::ElementType::Float ? sizeof(float) : sizeof(double);
  return (last + 1) * element;
}

void emitGridDescriptor(std::ostringstream& out,
                        const std::vector<reader::datafield::ElementType>& elements,
                        const std::vector<GridUse>& uses,
                        const std::string& u64,
                        const std::string& qualifier,
                        const std::string& global) {
  constexpr auto Axes = reader::datafield::MaxGridDimensions;
  out << "// the byte image of codegen::GridDescriptor\n"
      << "typedef struct {\n"
      << "  " << u64 << " data;\n"
      << "  double min[" << Axes << "];\n"
      << "  double deltaInv[" << Axes << "];\n"
      << "  " << u64 << " num[" << Axes << "];\n"
      << "  " << u64 << " stride[" << Axes << "];\n"
      << "  " << u64 << " before[" << Axes << "];\n"
      << "  " << u64 << " after[" << Axes << "];\n"
      << "  " << u64 << " componentStride;\n"
      << "  " << u64 << " knots[" << Axes << "];\n"
      << "  " << u64 << " inverseSpacings[" << Axes << "];\n"
      << "  " << u64 << " outside;\n"
      << "  " << u64 << " stale;\n"
      << "} SeissolExprGrid;\n\n";
  for (std::size_t grid = 0; grid < uses.size() && grid < elements.size(); ++grid) {
    if (uses[grid].dimensions == 0) {
      continue;
    }
    const char* type = elements[grid] == reader::datafield::ElementType::Float ? "float" : "double";
    out << qualifier << " double seissol_expr_load" << grid << "(" << global
        << "const SeissolExprGrid* g, SeissolI64 offset, SeissolI64 component) {\n"
        << "  return (double)((" << global << "const " << type
        << "*)(g->data))[offset + component * ("
        << "SeissolI64)g->componentStride];\n"
        << "}\n\n";
    // a knot, or the reciprocal of a spacing, of an axis of knots (expr::knotsOf())
    out << qualifier << " double seissol_expr_knot" << grid << "(" << global
        << "const SeissolExprGrid* g, SeissolI64 offset, SeissolI64 axis, SeissolI64 spacing) {\n"
        << "  return ((" << global << "const double*)(spacing ? g->inverseSpacings[axis] : "
        << "g->knots[axis]))[offset];\n"
        << "}\n\n";
  }
}

void emitStrictMacros(std::ostringstream& out,
                      const std::string& add,
                      const std::string& sub,
                      const std::string& mul) {
  out << "// each operation rounded on its own, never contracted (Instruction::strict)\n"
      << "#ifndef SEISSOL_EXPR_DADD\n"
      << "#define SEISSOL_EXPR_DADD(x, y) " << add << "\n"
      << "#define SEISSOL_EXPR_DSUB(x, y) " << sub << "\n"
      << "#define SEISSOL_EXPR_DMUL(x, y) " << mul << "\n"
      << "#endif\n\n";
}

bool containsContraction(const LoweredProgram& lowered) {
  for (const auto* stage : {&lowered.precompute(), &lowered.run()}) {
    for (const Instruction& inst : stage->code) {
      if (inst.op == Opcode::Contract) {
        return true;
      }
    }
  }
  return false;
}

const char* typeName(ComputeType type) {
  switch (type) {
  case ComputeType::F32:
    return "float";
  case ComputeType::F64:
    return "double";
  case ComputeType::I64:
    // defined by every source: `long long` in C++, `long` in OpenCL C (emitIntegerType())
    return "SeissolI64";
  case ComputeType::Bool:
    return "bool";
  }
  return "double";
}

void emitIntegerType(std::ostringstream& out, bool openCl) {
  out << "#ifndef SEISSOL_EXPR_I64\n"
      << "#define SEISSOL_EXPR_I64\n"
      << "typedef " << (openCl ? "long" : "long long") << " SeissolI64;\n"
      << "#endif\n\n";
}

const char* elementTypeName(reader::scripting::DataType type) {
  switch (type) {
  case reader::scripting::DataType::F32:
    return "float";
  case reader::scripting::DataType::F64:
    return "double";
  case reader::scripting::DataType::I32:
    return "int";
  case reader::scripting::DataType::I64:
    return "long long";
  }
  return "double";
}

} // namespace seissol::expr::codegen
