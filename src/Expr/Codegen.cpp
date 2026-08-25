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

#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr::codegen {

namespace {

bool identifierChar(char c) {
  return (std::isalnum(static_cast<unsigned char>(c)) != 0) || c == '_';
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

    if (token == "std" && style == MathStyle::Unqualified && text.compare(i, 2, "::") == 0) {
      // Device code has the maths built-ins unqualified, and NVRTC has no
      // <cmath> to bring `std` into scope. Dropping the qualifier here rather
      // than keeping a second expression table is the entire point of routing
      // both backends through one substitution.
      i += 2;
      continue;
    }
    if (token == "x") {
      out += x;
    } else if (token == "y") {
      out += y;
    } else if (token == "z") {
      out += z;
    } else if (token == "T") {
      out += computeType;
    } else {
      out += token;
    }
  }
  return out;
}

std::string slotName(std::int32_t slot) { return "s" + std::to_string(slot); }

std::string literal(double value, const std::string& computeType) {
  std::array<char, 48> buffer{};
  std::snprintf(buffer.data(), buffer.size(), "%.17g", value);
  return computeType + "(" + buffer.data() + ")";
}

void emitStageBody(std::ostringstream& out,
                   const StageCode& stage,
                   const std::vector<std::int32_t>& operands,
                   const std::string& computeType,
                   MathStyle style,
                   const StageAddressing& addressing,
                   const char* indent) {
  for (std::int32_t slot = 0; slot < stage.slotCount; ++slot) {
    out << indent << computeType << " " << slotName(slot) << " = " << computeType << "(0);\n";
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
      out << indent << dst << " = " << addressing.loadPersistent(inst.imm) << ";\n";
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
    }
  }

  for (const Store& store : stage.outputs) {
    out << indent << addressing.storeOutput(store.target) << " = " << slotName(store.source)
        << ";\n";
  }
  for (const Store& store : stage.persistent) {
    out << indent << addressing.storePersistent(store.target) << " = " << slotName(store.source)
        << ";\n";
  }
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

} // namespace seissol::expr::codegen
