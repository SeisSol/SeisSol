// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GhostState.h"

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Expr/Rewrite.h"
#include "Model/ScriptParameters.h"
#include "Numerical/Quadrature.h"
#include "Parallel/OpenMP.h"
#include "Physics/InitialField.h"
#include "Physics/ScriptField.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::physics {

namespace {

using reader::scripting::DataTable;
using reader::scripting::DataType;
using reader::scripting::Direction;

/// Points a kernel is bound for; a call with more runs in pieces. A face has fewer nodes.
constexpr std::size_t Capacity = 64;

/// The inputs of the integral besides those of the script: the step, and the state inside per
/// time point.
constexpr const char* Start = "$t0";
constexpr const char* Width = "$dt";

/// What an input of the integral reads.
enum class Input : std::uint8_t { X, Y, Z, Start, Width, Simulation, Parameter, Inner };

/// What an input of the script reads, besides the state inside, the time and the material.
std::optional<Input> scriptInputOf(const std::string& name) {
  const std::array<std::pair<const char*, Input>, 4> inputs{
      {{"x", Input::X}, {"y", Input::Y}, {"z", Input::Z}, {"sim", Input::Simulation}}};
  for (const auto& [channel, input] : inputs) {
    if (name == channel) {
      return input;
    }
  }
  return std::nullopt;
}

/// The input of inner component `component` at time point `point`.
std::string innerName(std::size_t point, const std::string& component) {
  std::string name = "$q";
  name.append(std::to_string(point)).append("_").append(component);
  return name;
}

/// The output of the integral of component `component`: a table has one column per name, and an
/// inner component takes the name of the component.
std::string integralName(const std::string& component) { return "$ghost_" + component; }

} // namespace

struct GhostStateIntegral::Impl {
  std::string name;
  std::vector<std::string> components;
  // the inner components the integral reads per time point; empty if none
  std::vector<std::string> inner;
  std::vector<double> quadPoints;
  std::vector<double> quadWeights;

  expr::Program program;
  std::vector<Input> inputs;
  // per input: for an inner one, its index in the inner state, i * inner.size() + k; for a
  // parameter of the material, its index among them
  std::vector<std::size_t> inputIndex;
  // per output, its component
  std::vector<std::size_t> outputComponent;

  // the storage the tables are bound to; every call moves the bases
  std::vector<std::array<double, 3>> boundPoints;
  std::vector<double> boundValues;
  double boundScalar{0};

  struct Thread {
    std::unique_ptr<DataTable> table;
    std::unique_ptr<expr::Binding> binding;
    std::unique_ptr<expr::Kernel> kernel;
    std::vector<const void*> inputBases;
    std::vector<void*> outputBases;
  };
  std::vector<Thread> threads;

  Thread& thread() {
    const auto id = seissol::OpenMP::threadId();
    if (id >= threads.size()) {
      logError() << "The ghost state of" << name << "was evaluated by thread" << id << "of"
                 << threads.size() << ".";
    }
    return threads[id];
  }
};

GhostStateIntegral::GhostStateIntegral(const expr::Program& ghost,
                                       const std::string& name,
                                       std::vector<std::string> components,
                                       std::vector<std::string> inner,
                                       const std::vector<std::string>& parameters,
                                       Missing missing,
                                       std::size_t timePoints)
    : impl_(std::make_unique<Impl>()) {
  auto& impl = *impl_;
  impl.name = name;
  impl.components = std::move(components);
  const auto& names = impl.components;

  if (!ghost.state().empty()) {
    throw std::invalid_argument(name + " keeps state, which a boundary condition cannot");
  }
  for (const auto& input : ghost.inputs()) {
    if (input.name != "t" && !scriptInputOf(input.name).has_value() &&
        std::find(inner.begin(), inner.end(), input.name) == inner.end() &&
        !model::parameterIndexOf(parameters, input.name).has_value()) {
      throw std::invalid_argument(
          std::string(name)
              .append(" reads `")
              .append(input.name)
              .append("`, which a boundary condition does not have; it reads x, y, z, t, sim")
              .append(inner.empty() ? "" : ", the state inside")
              .append(" and the parameters of the material: ")
              .append(model::parameterList(parameters)));
    }
  }

  // the script per component; what else it gives is not asked for
  std::vector<expr::NodeId> roots(names.size(), expr::NoNode);
  for (std::size_t i = 0; i < ghost.outputs().size(); ++i) {
    const auto found = std::find(names.begin(), names.end(), ghost.outputs()[i].name);
    if (found != names.end()) {
      roots[static_cast<std::size_t>(found - names.begin())] = ghost.roots()[i];
    }
  }

  impl.quadPoints.resize(timePoints);
  impl.quadWeights.resize(timePoints);
  quadrature::GaussLegendre(
      impl.quadPoints.data(), impl.quadWeights.data(), static_cast<unsigned>(timePoints));

  auto& program = impl.program;
  program.arena() = ghost.arena();
  program.setComputeType(ghost.computeType());
  for (const auto& grid : ghost.grids()) {
    program.internGrid(grid);
  }
  auto& arena = program.arena();
  const auto start = arena.field(Start);
  const auto width = arena.field(Width);
  const auto timeChannel = arena.findChannel("t");
  const auto isInner = [&](const std::string& component) {
    return std::find(inner.begin(), inner.end(), component) != inner.end();
  };

  std::vector<expr::NodeId> integrals(names.size(), expr::NoNode);
  for (std::size_t point = 0; point < timePoints; ++point) {
    // the time of the point and its weight, as kernels::AnalyticalBoundary computes them
    const auto time = arena.pw(
        expr::Fn::Div,
        arena.pw(expr::Fn::Add,
                 arena.pw(expr::Fn::Add,
                          arena.pw(expr::Fn::Mul, width, arena.konst(impl.quadPoints[point])),
                          arena.pw(expr::Fn::Mul, arena.konst(2.0), start)),
                 width),
        arena.konst(2.0));
    const auto weight = arena.pw(expr::Fn::Mul,
                                 arena.pw(expr::Fn::Mul, arena.konst(0.5), width),
                                 arena.konst(impl.quadWeights[point]));

    // the script at that time, reading the state inside at that time
    std::unordered_map<int, expr::NodeId> replacements;
    if (timeChannel >= 0) {
      replacements.emplace(timeChannel, time);
    }
    for (const auto& component : inner) {
      const auto channel = arena.findChannel(component);
      if (channel >= 0) {
        replacements.emplace(channel, arena.field(innerName(point, component)));
      }
    }
    std::unordered_map<expr::NodeId, expr::NodeId> memo;
    for (std::size_t j = 0; j < names.size(); ++j) {
      expr::NodeId value = expr::NoNode;
      if (roots[j] != expr::NoNode) {
        value = expr::replaceChannels(arena, roots[j], replacements, memo);
      } else if (missing == Missing::Inner && isInner(names[j])) {
        value = arena.field(innerName(point, names[j]));
      }
      if (value == expr::NoNode) {
        continue;
      }
      // added up in the order of the time points, as the update of the nodal state does
      const auto term = arena.pw(expr::Fn::Mul, weight, value);
      integrals[j] =
          integrals[j] == expr::NoNode ? term : arena.pw(expr::Fn::Add, integrals[j], term);
    }
  }

  std::vector<expr::NodeId> outputRoots;
  for (std::size_t j = 0; j < names.size(); ++j) {
    if (integrals[j] != expr::NoNode) {
      program.addOutput(integralName(names[j]), DataType::F64, integrals[j]);
      impl.outputComponent.push_back(j);
      outputRoots.push_back(integrals[j]);
    }
  }

  // the inputs the integral reads
  const auto read = expr::channelsRead(program, outputRoots);
  bool readsInner = false;
  for (std::size_t point = 0; point < timePoints && !readsInner; ++point) {
    for (const auto& component : inner) {
      readsInner = readsInner || read.count(innerName(point, component)) > 0;
    }
  }
  if (readsInner) {
    impl.inner = inner;
  }
  for (const auto& channel : read) {
    std::optional<std::pair<Input, std::size_t>> meaning;
    if (channel == Start) {
      meaning = std::pair{Input::Start, std::size_t{0}};
    } else if (channel == Width) {
      meaning = std::pair{Input::Width, std::size_t{0}};
    } else if (const auto scriptInput = scriptInputOf(channel)) {
      meaning = std::pair{*scriptInput, std::size_t{0}};
    } else if (const auto parameter = model::parameterIndexOf(parameters, channel)) {
      meaning = std::pair{Input::Parameter, *parameter};
    } else {
      for (std::size_t point = 0; point < timePoints && !meaning.has_value(); ++point) {
        for (std::size_t k = 0; k < impl.inner.size(); ++k) {
          if (channel == innerName(point, impl.inner[k])) {
            meaning = std::pair{Input::Inner, point * impl.inner.size() + k};
            break;
          }
        }
      }
    }
    if (!meaning.has_value()) {
      throw std::invalid_argument(std::string(name)
                                      .append(" reads `")
                                      .append(channel)
                                      .append("` at the time points of a step, which it cannot"));
    }
    program.addInput(channel, DataType::F64);
    impl.inputs.push_back(meaning->first);
    impl.inputIndex.push_back(meaning->second);
  }

  // one kernel per thread, bound to storage of its own; a call moves the bases
  impl.threads.resize(std::max<std::size_t>(1, seissol::OpenMP::threadCount()));
  if (program.outputs().empty()) {
    return;
  }
  const auto innerInputs = impl.inner.size() * timePoints;
  impl.boundPoints.assign(Capacity, {0, 0, 0});
  impl.boundValues.assign(Capacity * (innerInputs + program.outputs().size()), 0);
  expr::BackendOptions options;
  options.preferred = expr::BackendKind::RtcCpu;
  for (auto& thread : impl.threads) {
    thread.table = std::make_unique<DataTable>(Capacity);
    auto& table = *thread.table;
    for (std::size_t i = 0; i < program.inputs().size(); ++i) {
      const auto& input = program.inputs()[i].name;
      switch (impl.inputs[i]) {
      case Input::X:
        table.bindViewConst<double>(input, Direction::In, impl.boundPoints.data()->data(), 3, 0);
        break;
      case Input::Y:
        table.bindViewConst<double>(input, Direction::In, impl.boundPoints.data()->data(), 3, 1);
        break;
      case Input::Z:
        table.bindViewConst<double>(input, Direction::In, impl.boundPoints.data()->data(), 3, 2);
        break;
      case Input::Inner:
        table.bindViewConst<double>(
            input, Direction::In, impl.boundValues.data() + impl.inputIndex[i] * Capacity);
        break;
      default:
        table.bindViewConst<double>(input, Direction::In, &impl.boundScalar, 0);
        break;
      }
    }
    for (std::size_t j = 0; j < program.outputs().size(); ++j) {
      table.bindView<double>(program.outputs()[j].name,
                             Direction::Out,
                             impl.boundValues.data() + (innerInputs + j) * Capacity);
    }
    thread.binding = std::make_unique<expr::Binding>(expr::Binding::bind(program, table));
    thread.kernel =
        expr::makeKernel(program, *thread.binding, reader::datafield::sharedGridStore(), options);
    thread.kernel->precompute(table);
    // the first kernel says how it went; the others follow it quietly
    options.preferred = thread.kernel->kind();
    options.quiet = true;
    thread.inputBases.resize(program.inputs().size());
    thread.outputBases.resize(program.outputs().size());
  }
}

GhostStateIntegral::~GhostStateIntegral() = default;

std::size_t GhostStateIntegral::timePointCount() const { return impl_->quadPoints.size(); }

std::size_t GhostStateIntegral::componentCount() const { return impl_->components.size(); }

std::size_t GhostStateIntegral::innerCount() const { return impl_->inner.size(); }

void GhostStateIntegral::timePoints(double start, double width, double* times) const {
  for (std::size_t point = 0; point < impl_->quadPoints.size(); ++point) {
    times[point] = (width * impl_->quadPoints[point] + 2 * start + width) / 2;
  }
}

void GhostStateIntegral::evaluate(double start,
                                  double width,
                                  std::size_t simulation,
                                  const std::array<double, 3>* points,
                                  std::size_t count,
                                  const double* parameters,
                                  const double* inner,
                                  double* ghost) const {
  const auto& impl = *impl_;
  // what the integral does not give is zero
  std::fill_n(ghost, impl.components.size() * count, 0.0);
  if (impl.outputComponent.empty()) {
    return;
  }
  auto& thread = impl_->thread();
  const auto simulationValue = static_cast<double>(simulation);

  for (std::size_t first = 0; first < count; first += Capacity) {
    const std::size_t pieceCount = std::min(Capacity, count - first);
    for (std::size_t i = 0; i < impl.inputs.size(); ++i) {
      switch (impl.inputs[i]) {
      case Input::X:
      case Input::Y:
      case Input::Z:
        // the base of the array; the binding adds the offset of the coordinate
        thread.inputBases[i] = points[first].data();
        break;
      case Input::Start:
        thread.inputBases[i] = &start;
        break;
      case Input::Width:
        thread.inputBases[i] = &width;
        break;
      case Input::Simulation:
        thread.inputBases[i] = &simulationValue;
        break;
      case Input::Parameter:
        thread.inputBases[i] = parameters + impl.inputIndex[i];
        break;
      case Input::Inner:
        thread.inputBases[i] = inner + impl.inputIndex[i] * count + first;
        break;
      }
    }
    for (std::size_t j = 0; j < impl.outputComponent.size(); ++j) {
      thread.outputBases[j] = ghost + impl.outputComponent[j] * count + first;
    }
    expr::KernelArgs args;
    args.inputs = thread.inputBases.data();
    args.inputCount = thread.inputBases.size();
    args.outputs = thread.outputBases.data();
    args.outputCount = thread.outputBases.size();
    args.first = 0;
    args.count = pieceCount;
    thread.kernel->run(args);
  }
}

namespace {

/// Whether two programs compute the same outputs from the same inputs: the canonical form leaves
/// out the names.
bool sameProgram(const expr::Program& first, const expr::Program& second) {
  const auto namesOf = [](const std::vector<expr::VarSpec>& variables) {
    std::vector<std::string> names;
    names.reserve(variables.size());
    for (const auto& variable : variables) {
      names.push_back(variable.name);
    }
    return names;
  };
  return namesOf(first.inputs()) == namesOf(second.inputs()) &&
         namesOf(first.outputs()) == namesOf(second.outputs()) &&
         first.canonicalForm() == second.canonicalForm();
}

} // namespace

std::unique_ptr<AnalyticalGhostState>
    AnalyticalGhostState::fromFields(const std::vector<std::unique_ptr<InitialField>>& fields,
                                     const std::vector<std::string>& components,
                                     const std::vector<std::string>& parameters,
                                     std::size_t simulations,
                                     std::size_t timePoints) {
  std::vector<const ScriptField*> scripts;
  for (const auto& field : fields) {
    const auto* script = dynamic_cast<const ScriptField*>(field.get());
    if (script == nullptr || !script->compiled()) {
      return nullptr;
    }
    if (!script->error().empty()) {
      logError() << "The analytical boundary takes the initial condition at its faces, but"
                 << script->error();
    }
    scripts.push_back(script);
  }
  if (scripts.empty()) {
    return nullptr;
  }

  auto state = std::make_unique<AnalyticalGhostState>();
  // the program of each integral, and the integral of each field
  std::vector<const expr::Program*> programs;
  std::vector<std::size_t> integralOfField;
  for (const auto* script : scripts) {
    const auto& program = *script->program();
    const auto found = std::find_if(programs.begin(), programs.end(), [&](const auto* other) {
      return sameProgram(*other, program);
    });
    integralOfField.push_back(static_cast<std::size_t>(found - programs.begin()));
    if (found != programs.end()) {
      continue;
    }
    programs.push_back(&program);
    try {
      state->integrals_.push_back(
          std::make_unique<GhostStateIntegral>(program,
                                               "the initial condition",
                                               components,
                                               std::vector<std::string>{},
                                               parameters,
                                               GhostStateIntegral::Missing::Zero,
                                               timePoints));
    } catch (const std::invalid_argument& error) {
      logError() << "Analytical boundary:" << error.what();
    }
  }
  for (std::size_t simulation = 0; simulation < simulations; ++simulation) {
    const auto field = simulation % scripts.size();
    state->integralOf_.push_back(integralOfField[field]);
    state->scriptSimulation_.push_back(scripts[field]->simulation());
  }
  logInfo() << "The analytical boundary integrates the initial condition over a time step in"
            << programs.size() << "program(s) for" << simulations << "simulation(s).";
  return state;
}

const GhostStateIntegral& AnalyticalGhostState::integral(std::size_t simulation) const {
  return *integrals_[integralOf_[simulation]];
}

std::size_t AnalyticalGhostState::scriptSimulation(std::size_t simulation) const {
  return scriptSimulation_[simulation];
}

} // namespace seissol::physics
