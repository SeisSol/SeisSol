// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Expr/Cost.h"

#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"

#include <cstdint>
#include <sstream>
#include <string>

namespace seissol::expr {

namespace {

enum class Class : std::uint8_t { Additive, Multiplicative, Division, Transcendental };

/// Which cost class each Fn falls into.
///
/// An explicit switch with NO default, deliberately. It is a second table
/// alongside SEISSOL_EXPR_PW_LIST, which is exactly the kind of duplication the
/// rest of this codebase avoids -- but here it cannot be derived: the table says
/// what an operation COMPUTES, and nothing in `(std::fmod(x, y))` says whether
/// fmod is one instruction or thirty. Without a default, adding an Fn fails to
/// compile until someone classifies it, which is the only guarantee available.
Class classify(Fn fn) {
  switch (fn) {
  case Fn::Add:
  case Fn::Sub:
  case Fn::Neg:
  case Fn::Abs:
  case Fn::Sign:
  case Fn::Min:
  case Fn::Max:
  case Fn::Floor:
  case Fn::Ceil:
  case Fn::Round:
  case Fn::Lt:
  case Fn::Le:
  case Fn::Eq:
  case Fn::And:
  case Fn::Or:
  case Fn::Select:
    return Class::Additive;
  case Fn::Mul:
    return Class::Multiplicative;
  case Fn::Div:
  case Fn::Rcp:
  case Fn::Mod:
    return Class::Division;
  case Fn::Sqrt:
  case Fn::Exp:
  case Fn::Log:
  case Fn::Log2:
  case Fn::Log10:
  case Fn::Pow:
  case Fn::Sin:
  case Fn::Cos:
  case Fn::Tan:
  case Fn::Asin:
  case Fn::Acos:
  case Fn::Atan:
  case Fn::Atan2:
  case Fn::Sinh:
  case Fn::Cosh:
  case Fn::Tanh:
  case Fn::Asinh:
  case Fn::Acosh:
  case Fn::Atanh:
  case Fn::Erf:
    return Class::Transcendental;
  }
  return Class::Additive;
}

StageCost countStage(const StageCode& stage) {
  StageCost costs;
  for (const Instruction& inst : stage.code) {
    switch (inst.op) {
    case Opcode::Const:
      // Materialised into an immediate by every backend, and by the interpreter
      // into a slot write. Neither is arithmetic.
      break;
    case Opcode::LoadInput:
    case Opcode::LoadPersistent:
      ++costs.loads;
      break;
    case Opcode::Lookup:
      ++costs.lookups;
      break;
    case Opcode::Pw:
      switch (classify(inst.fn)) {
      case Class::Additive:
        ++costs.additive;
        break;
      case Class::Multiplicative:
        ++costs.multiplicative;
        break;
      case Class::Division:
        ++costs.divisions;
        break;
      case Class::Transcendental:
        ++costs.transcendentals;
        break;
      }
      break;
    }
  }
  costs.stores = stage.outputs.size() + stage.persistent.size();
  return costs;
}

const char* typeName(ComputeType type) { return type == ComputeType::F32 ? "f32" : "f64"; }

} // namespace

std::uint64_t StageCost::bytes(ComputeType type) const {
  const std::uint64_t width = type == ComputeType::F32 ? 4 : 8;
  return (loads + stores) * width;
}

double StageCost::intensity(ComputeType type) const {
  const std::uint64_t moved = bytes(type);
  return moved == 0 ? 0.0 : static_cast<double>(weighted()) / static_cast<double>(moved);
}

std::string ProgramCost::summary(ComputeType type) const {
  std::ostringstream out;
  const auto line = [&out, type](const char* label, const StageCost& costs) {
    out << label << ": " << costs.operations() << " ops (" << costs.additive << "+/-, "
        << costs.multiplicative << "*, " << costs.divisions << "/, " << costs.transcendentals
        << " transc.";
    if (costs.lookups > 0) {
      out << ", " << costs.lookups << " lookups";
    }
    out << "), weighted " << costs.weighted() << ", " << costs.bytes(type) << " B/point, "
        << costs.intensity(type) << " ops/B";
  };
  if (precompute.operations() > 0 || precompute.lookups > 0) {
    line("precompute", precompute);
    out << "; ";
  }
  line("run", run);
  out << " [" << typeName(type) << "]";
  return out.str();
}

ProgramCost cost(const LoweredProgram& lowered, ComputeType /*type*/) {
  ProgramCost costs;
  costs.precompute = countStage(lowered.precompute());
  costs.run = countStage(lowered.run());
  return costs;
}

} // namespace seissol::expr
