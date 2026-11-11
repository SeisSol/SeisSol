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
#include <cstddef>
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

  for (const Instruction& inst : stage.code) {
    const std::string dst = slotName(inst.dst);
    const std::string type = typeName(inst.type);
    switch (inst.op) {
    case Opcode::Const:
      out << indent << dst << " = " << literal(inst.value, type) << ";\n";
      break;
    case Opcode::LoadInput:
      out << indent << dst << " = "
          << converted(addressing.loadInput(inst.imm), dataType, inst.type) << ";\n";
      break;
    case Opcode::LoadPersistent:
      if (inst.imm < addressing.stateSlots) {
        if (addressing.loadState == nullptr) {
          throw std::invalid_argument(
              "expr: a state reached a code generator that cannot address states");
        }
        out << indent << dst << " = "
            << converted(addressing.loadState(inst.imm), dataType, inst.type) << ";\n";
      } else {
        out << indent << dst << " = "
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
      out << indent << dst << " = "
          << pwText(inst, a0, a1, a2, type, slotType(stage, operands[inst.operandBegin]), style)
          << ";\n";
      break;
    }
    case Opcode::Convert: {
      const std::int32_t source = operands[inst.operandBegin];
      out << indent << dst << " = "
          << converted(slotName(source), slotType(stage, source), inst.type) << ";\n";
      break;
    }
    case Opcode::GridParam:
    case Opcode::GridLoad: {
      std::string operand;
      if (inst.op == Opcode::GridLoad) {
        const std::int32_t source = operands[inst.operandBegin];
        operand = converted(slotName(source), slotType(stage, source), ComputeType::I64);
      }
      out << indent << dst << " = " << gridText(inst, addressing, operand) << ";\n";
      break;
    }
    case Opcode::Contract:
      if (addressing.contract == nullptr) {
        throw std::invalid_argument(
            "expr: a contraction reached a code generator that cannot address blocks");
      }
      emitContraction(out, inst, dst, type, *addressing.contract, indent);
      break;
    }
  }

  const auto value = [&](const Store& store) {
    return converted(slotName(store.source), slotType(stage, store.source), dataType);
  };
  for (const Store& store : stage.outputs) {
    out << indent << addressing.storeOutput(store.target, value(store)) << ";\n";
  }
  for (const Store& store : stage.persistent) {
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
    out << inner << typeName(slotType(stage, index)) << " " << slotName(index) << "[" << block
        << "];\n";
  }

  for (const Instruction& inst : stage.code) {
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
    }
  }

  const auto value = [&](const Store& store) {
    return converted(slot(store.source), slotType(stage, store.source), dataType);
  };
  for (const Store& store : stage.outputs) {
    laneLoop(addressing.storeOutput(store.target, value(store)), true);
  }
  for (const Store& store : stage.persistent) {
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
