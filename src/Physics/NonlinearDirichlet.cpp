// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "NonlinearDirichlet.h"

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Physics/GhostState.h"
#include "Physics/ScriptScope.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/ReaderBuilder.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::physics {

namespace {

using reader::scripting::DataTable;
using reader::scripting::DataType;
using reader::scripting::Direction;

constexpr const char* Frame = "frame";

std::string listOf(const std::vector<std::string>& names) {
  std::string list;
  for (std::size_t i = 0; i < names.size(); ++i) {
    list.append(i == 0 ? "" : ", ").append(names[i]);
  }
  return list;
}

/// The value of `root`, which reads no channel.
double constantOf(const expr::Program& program, expr::NodeId root) {
  expr::Program constant;
  constant.arena() = program.arena();
  constant.setComputeType(program.computeType());
  for (const auto& grid : program.grids()) {
    constant.internGrid(grid);
  }
  constant.addOutput(Frame, DataType::F64, root);
  double value = std::numeric_limits<double>::quiet_NaN();
  DataTable table(1);
  table.bindView<double>(Frame, Direction::Out, &value);
  auto binding = expr::Binding::bind(constant, table);
  expr::BackendOptions options;
  options.quiet = true;
  const auto kernel =
      expr::makeKernel(constant, binding, reader::datafield::sharedGridStore(), options);
  kernel->precompute(table);
  kernel->run(table);
  return value;
}

} // namespace

NonlinearDirichlet::NonlinearDirichlet(const std::string& path,
                                       std::vector<std::string> quantities,
                                       std::vector<std::string> parameters,
                                       std::size_t timePoints,
                                       const std::vector<std::string>& inputs,
                                       const ProgramFilter& filter)
    : name_(path), quantities_(std::move(quantities)), parameters_(std::move(parameters)) {
  // in sderiv, a quantity reads the inner state also next to the definition of its ghost value
  expr::SderivOptions options;
  const auto& innerState = inputs.empty() ? quantities_ : inputs;
  options.inputs.insert(innerState.begin(), innerState.end());
  std::string reason;
  const auto program = reader::scripting::buildProgram(path, &reason, options);
  if (!program.has_value()) {
    logError() << "The nonlinear Dirichlet boundary" << path
               << "has to be a script that compiles, an sderiv module or a Lua model that "
                  "traces --"
               << reason;
    return;
  }
  try {
    setUp(filter ? filter(*program) : *program, timePoints);
  } catch (const std::invalid_argument& error) {
    logError() << "Nonlinear Dirichlet boundary:" << error.what();
  }
  logInfo() << "The nonlinear Dirichlet boundary" << path << "gives the ghost state in"
            << (faceAligned_ ? "the face-aligned basis." : "global coordinates.");
}

NonlinearDirichlet::NonlinearDirichlet(const expr::Program& program,
                                       const std::string& name,
                                       std::vector<std::string> quantities,
                                       std::vector<std::string> parameters,
                                       std::size_t timePoints)
    : name_(name), quantities_(std::move(quantities)), parameters_(std::move(parameters)) {
  setUp(program, timePoints);
}

NonlinearDirichlet::~NonlinearDirichlet() = default;

void NonlinearDirichlet::setUp(const expr::Program& program, std::size_t timePoints) {
  if (!program.state().empty()) {
    throw std::invalid_argument(name_ + " keeps state, which a boundary condition cannot");
  }

  // the ghost values and the frame
  expr::Program ghost;
  ghost.arena() = program.arena();
  ghost.setComputeType(program.computeType());
  for (const auto& grid : program.grids()) {
    ghost.internGrid(grid);
  }
  for (const auto& input : program.inputs()) {
    ghost.addInput(input.name, input.type);
  }
  for (std::size_t i = 0; i < program.outputs().size(); ++i) {
    const auto& output = program.outputs()[i];
    if (output.name == Frame) {
      const auto frame = program.roots()[i];
      if (!expr::channelsRead(program, {frame}).empty()) {
        throw std::invalid_argument("the frame of " + name_ +
                                    " is the same on all faces, so it may read nothing");
      }
      const double value = constantOf(program, frame);
      if (value != 0.0 && value != 1.0) {
        throw std::invalid_argument("the frame of " + name_ +
                                    " has to be 0, for global coordinates, or 1, for the "
                                    "face-aligned basis");
      }
      faceAligned_ = value == 1.0;
      continue;
    }
    if (std::find(quantities_.begin(), quantities_.end(), output.name) == quantities_.end()) {
      std::string message = name_ + " gives `" + output.name;
      message += "`, which is no quantity of the material, nor the frame; the quantities are ";
      message += listOf(quantities_);
      throw std::invalid_argument(message);
    }
    ghost.addOutput(output.name, output.type, program.roots()[i]);
  }

  integral_ = std::make_unique<GhostStateIntegral>(ghost,
                                                   name_,
                                                   quantities_,
                                                   quantities_,
                                                   parameters_,
                                                   GhostStateIntegral::Missing::Inner,
                                                   timePoints);
#ifdef ACL_DEVICE
  deviceIntegral_ = std::make_unique<GhostStateIntegral>(ghost,
                                                         name_,
                                                         quantities_,
                                                         quantities_,
                                                         parameters_,
                                                         GhostStateIntegral::Missing::Zero,
                                                         timePoints,
                                                         false);
#endif
}

bool NonlinearDirichlet::faceAligned() const { return faceAligned_; }

std::size_t NonlinearDirichlet::quantityCount() const { return quantities_.size(); }

const GhostStateIntegral& NonlinearDirichlet::integral() const { return *integral_; }

const GhostStateIntegral* NonlinearDirichlet::deviceIntegral() const {
  return deviceIntegral_.get();
}

} // namespace seissol::physics
