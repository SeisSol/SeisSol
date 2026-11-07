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
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr::codegen {

namespace {

bool identifierChar(char c) {
  return (std::isalnum(static_cast<unsigned char>(c)) != 0) || c == '_';
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

} // namespace

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
      if (style == MathStyle::CNames) {
        // the function after the qualifier, by its C name for the compute type
        std::size_t end = i;
        while (end < text.size() && identifierChar(text[end])) {
          ++end;
        }
        std::string function = text.substr(i, end - i);
        if (function == "abs") {
          function = "fabs";
        }
        if (computeType == "float") {
          function += "f";
        }
        out += function;
        i = end;
      }
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
  std::array<char, 48> buffer{};
  std::snprintf(buffer.data(), buffer.size(), "%.17g", value);
  return "((" + computeType + ")" + buffer.data() + ")";
}

void emitStageBody(std::ostringstream& out,
                   const StageCode& stage,
                   const std::vector<std::int32_t>& operands,
                   const std::string& computeType,
                   MathStyle style,
                   const StageAddressing& addressing,
                   const char* indent) {
  for (std::int32_t slot = 0; slot < stage.slotCount; ++slot) {
    out << indent << computeType << " " << slotName(slot) << " = (" << computeType << ")0;\n";
  }
  if (stage.slotCount > 0) {
    // Silences the unused-variable warning for a slot that is written but never
    // read, which liveness can legitimately produce.
    out << indent << "(void)" << slotName(0) << ";\n";
  }

  for (const Instruction& inst : stage.code) {
    const std::string dst = slotName(inst.dst);
    switch (inst.op) {
    case Opcode::Const:
      out << indent << dst << " = " << literal(inst.value, computeType) << ";\n";
      break;
    case Opcode::LoadInput:
      out << indent << dst << " = " << addressing.loadInput(inst.imm) << ";\n";
      break;
    case Opcode::LoadPersistent:
      if (inst.imm < addressing.stateSlots) {
        if (addressing.loadState == nullptr) {
          throw std::invalid_argument(
              "expr: a state reached a code generator that cannot address states");
        }
        out << indent << dst << " = " << addressing.loadState(inst.imm) << ";\n";
      } else {
        out << indent << dst << " = " << addressing.loadPersistent(inst.imm - addressing.stateSlots)
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
      out << indent << dst << " = "
          << substitute(expressionText(inst.fn), a0, a1, a2, computeType, style) << ";\n";
      break;
    }
    case Opcode::Lookup:
      throw std::invalid_argument(
          "expr: a grid lookup reached the code generator; gate on containsLookup() first");
    case Opcode::Contract:
      if (addressing.contract == nullptr) {
        throw std::invalid_argument(
            "expr: a contraction reached a code generator that cannot address blocks");
      }
      emitContraction(out, inst, dst, computeType, *addressing.contract, indent);
      break;
    }
  }

  for (const Store& store : stage.outputs) {
    out << indent << addressing.storeOutput(store.target, slotName(store.source)) << ";\n";
  }
  for (const Store& store : stage.persistent) {
    if (store.target < addressing.stateSlots) {
      if (addressing.storeState == nullptr) {
        throw std::invalid_argument(
            "expr: a state reached a code generator that cannot address states");
      }
      out << indent << addressing.storeState(store.target, slotName(store.source)) << ";\n";
    } else {
      out << indent
          << addressing.storePersistent(store.target - addressing.stateSlots,
                                        slotName(store.source))
          << ";\n";
    }
  }
}

void emitStageBodyTiled(std::ostringstream& out,
                        const StageCode& stage,
                        const std::vector<std::int32_t>& operands,
                        const std::string& computeType,
                        MathStyle style,
                        const StageAddressing& addressing,
                        std::size_t block,
                        const char* indent) {
  const std::string blockSize = std::to_string(block) + "ul";
  const std::string inner = std::string(indent) + "  ";
  const std::string lane = inner + "  ";
  const auto slot = [](std::int32_t index) { return slotName(index) + "[j]"; };
  // one instruction over the lanes of the block; `l` where the addressing reads the point
  const auto laneLoop = [&](const std::string& statement, bool readsPoint) {
    out << "#pragma omp simd\n" << inner << "for (unsigned long j = 0; j < n; ++j) {\n";
    if (readsPoint) {
      out << lane << "const unsigned long l = base + j;\n";
    }
    out << lane << statement << ";\n" << inner << "}\n";
  };

  out << indent << "for (unsigned long base = 0; base < count; base += " << blockSize << ") {\n"
      << inner << "const unsigned long n = count - base < " << blockSize
      << " ? count - base : " << blockSize << ";\n";
  for (std::int32_t index = 0; index < stage.slotCount; ++index) {
    out << inner << computeType << " " << slotName(index) << "[" << block << "];\n";
  }

  for (const Instruction& inst : stage.code) {
    const std::string dst = slot(inst.dst);
    switch (inst.op) {
    case Opcode::Const:
      laneLoop(dst + " = " + literal(inst.value, computeType), false);
      break;
    case Opcode::LoadInput:
      laneLoop(dst + " = " + addressing.loadInput(inst.imm), true);
      break;
    case Opcode::LoadPersistent:
      if (inst.imm < addressing.stateSlots) {
        if (addressing.loadState == nullptr) {
          throw std::invalid_argument(
              "expr: a state reached a code generator that cannot address states");
        }
        laneLoop(dst + " = " + addressing.loadState(inst.imm), true);
      } else {
        laneLoop(dst + " = " + addressing.loadPersistent(inst.imm - addressing.stateSlots), true);
      }
      break;
    case Opcode::Pw: {
      const std::string a0 = slot(operands[inst.operandBegin]);
      const std::string a1 = inst.operandCount > 1 ? slot(operands[inst.operandBegin + 1]) : a0;
      const std::string a2 = inst.operandCount > 2 ? slot(operands[inst.operandBegin + 2]) : a0;
      laneLoop(dst + " = " + substitute(expressionText(inst.fn), a0, a1, a2, computeType, style),
               false);
      break;
    }
    case Opcode::Lookup:
      throw std::invalid_argument(
          "expr: a grid lookup reached the code generator; gate on containsLookup() first");
    case Opcode::Contract:
      throw std::invalid_argument(
          "expr: a contraction reached the tiled code generator; emit one point after the other");
    }
  }

  for (const Store& store : stage.outputs) {
    laneLoop(addressing.storeOutput(store.target, slot(store.source)), true);
  }
  for (const Store& store : stage.persistent) {
    if (store.target < addressing.stateSlots) {
      if (addressing.storeState == nullptr) {
        throw std::invalid_argument(
            "expr: a state reached a code generator that cannot address states");
      }
      laneLoop(addressing.storeState(store.target, slot(store.source)), true);
    } else {
      laneLoop(addressing.storePersistent(store.target - addressing.stateSlots, slot(store.source)),
               true);
    }
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

bool containsLookup(const LoweredProgram& lowered) {
  for (const auto* stage : {&lowered.precompute(), &lowered.run()}) {
    for (const Instruction& inst : stage->code) {
      if (inst.op == Opcode::Lookup) {
        return true;
      }
    }
  }
  return false;
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
