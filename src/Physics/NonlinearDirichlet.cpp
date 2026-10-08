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
#include "Initializer/Typedefs.h"
#include "Parallel/OpenMP.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/ReaderBuilder.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
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

/// Points a kernel is bound for; a call with more runs in pieces. A face has fewer nodes.
constexpr std::size_t Capacity = 64;

constexpr const char* Frame = "frame";

/// What an input of the script reads.
enum class Input : std::uint8_t { X, Y, Z, Time, Simulation, Density, Mu, Lambda, Quantity };

std::optional<Input> inputOf(const std::string& name) {
  const std::array<std::pair<const char*, Input>, 8> inputs{{{"x", Input::X},
                                                             {"y", Input::Y},
                                                             {"z", Input::Z},
                                                             {"t", Input::Time},
                                                             {"sim", Input::Simulation},
                                                             {"rho", Input::Density},
                                                             {"mu", Input::Mu},
                                                             {"lambda", Input::Lambda}}};
  for (const auto& [channel, input] : inputs) {
    if (name == channel) {
      return input;
    }
  }
  return std::nullopt;
}

/// The name of the ghost value of `quantity` in the bound program: a table has one column per
/// name, and the inner state takes the name of the quantity.
std::string ghostName(const std::string& quantity) { return "$ghost_" + quantity; }

/// The material of a cell as the script sees it.
std::array<double, 3> materialOf(const CellMaterialData& materialData) {
  if (materialData.local == nullptr) {
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    return {NaN, NaN, NaN};
  }
  return {materialData.local->getDensity(),
          materialData.local->getMuBar(),
          materialData.local->getLambdaBar()};
}

std::string listOf(const std::vector<std::string>& names) {
  std::string list;
  for (std::size_t i = 0; i < names.size(); ++i) {
    list += (i == 0 ? "" : ", ") + names[i];
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

struct NonlinearDirichlet::Impl {
  std::string name;
  std::vector<std::string> quantities;
  bool faceAligned{false};

  // the program as bound: the inputs it reads, the ghost values under ghostName()
  expr::Program program;
  std::vector<Input> inputs;
  // per input of program() that is a quantity, which one
  std::vector<std::size_t> inputQuantity;
  // per output, its quantity
  std::vector<std::size_t> outputQuantity;

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
      logError() << "The nonlinear Dirichlet boundary was evaluated by thread" << id << "of"
                 << threads.size() << ".";
    }
    return threads[id];
  }
};

NonlinearDirichlet::NonlinearDirichlet(const std::string& path, std::vector<std::string> quantities)
    : impl_(std::make_unique<Impl>()) {
  impl_->name = path;
  impl_->quantities = std::move(quantities);

  // in sderiv, a quantity reads the inner state also next to the definition of its ghost value
  expr::SderivOptions options;
  options.inputs.insert(impl_->quantities.begin(), impl_->quantities.end());
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
    setUp(*program);
  } catch (const std::invalid_argument& error) {
    logError() << "Nonlinear Dirichlet boundary:" << error.what();
  }
  logInfo() << "The nonlinear Dirichlet boundary" << path << "gives the ghost state in"
            << (impl_->faceAligned ? "the face-aligned basis." : "global coordinates.");
}

NonlinearDirichlet::NonlinearDirichlet(const expr::Program& program,
                                       const std::string& name,
                                       std::vector<std::string> quantities)
    : impl_(std::make_unique<Impl>()) {
  impl_->name = name;
  impl_->quantities = std::move(quantities);
  setUp(program);
}

NonlinearDirichlet::~NonlinearDirichlet() = default;

void NonlinearDirichlet::setUp(const expr::Program& program) {
  auto& impl = *impl_;
  const auto& quantities = impl.quantities;
  if (!program.state().empty()) {
    throw std::invalid_argument(impl.name + " keeps state, which a boundary condition cannot");
  }

  // the ghost values and the frame
  std::vector<std::pair<std::size_t, expr::NodeId>> ghost;
  std::optional<expr::NodeId> frame;
  for (std::size_t i = 0; i < program.outputs().size(); ++i) {
    const auto& output = program.outputs()[i].name;
    if (output == Frame) {
      frame = program.roots()[i];
      continue;
    }
    const auto found = std::find(quantities.begin(), quantities.end(), output);
    if (found == quantities.end()) {
      std::string message = impl.name + " gives `" + output;
      message += "`, which is no quantity of the material, nor the frame; the quantities are ";
      message += listOf(quantities);
      throw std::invalid_argument(message);
    }
    ghost.emplace_back(static_cast<std::size_t>(found - quantities.begin()), program.roots()[i]);
  }
  if (frame.has_value()) {
    if (!expr::channelsRead(program, {*frame}).empty()) {
      throw std::invalid_argument("the frame of " + impl.name +
                                  " is the same on all faces, so it may read nothing");
    }
    const double value = constantOf(program, *frame);
    if (value != 0.0 && value != 1.0) {
      throw std::invalid_argument("the frame of " + impl.name +
                                  " has to be 0, for global coordinates, or 1, for the "
                                  "face-aligned basis");
    }
    impl.faceAligned = value == 1.0;
  }

  // the program of the ghost values, with the inputs they read
  std::vector<expr::NodeId> roots;
  roots.reserve(ghost.size());
  for (const auto& entry : ghost) {
    roots.push_back(entry.second);
  }
  const auto read = expr::channelsRead(program, roots);
  auto& bound = impl.program;
  bound.arena() = program.arena();
  bound.setComputeType(program.computeType());
  for (const auto& grid : program.grids()) {
    bound.internGrid(grid);
  }
  for (const auto& input : program.inputs()) {
    if (read.count(input.name) == 0) {
      continue;
    }
    const auto found = std::find(quantities.begin(), quantities.end(), input.name);
    if (found != quantities.end()) {
      impl.inputs.push_back(Input::Quantity);
      impl.inputQuantity.push_back(static_cast<std::size_t>(found - quantities.begin()));
    } else {
      const auto meaning = inputOf(input.name);
      if (!meaning.has_value()) {
        std::string message = impl.name + " reads `" + input.name;
        message += "`; the boundary reads the quantities of the material (";
        message += listOf(quantities);
        message += "), x, y, z, t, sim, rho, mu and lambda";
        throw std::invalid_argument(message);
      }
      impl.inputs.push_back(*meaning);
      impl.inputQuantity.push_back(0);
    }
    bound.addInput(input.name, input.type);
  }
  for (const auto& [quantity, root] : ghost) {
    bound.addOutput(ghostName(quantities[quantity]), DataType::F64, root);
    impl.outputQuantity.push_back(quantity);
  }

  // one kernel per thread, bound to storage of its own; a call moves the bases
  impl.threads.resize(std::max<std::size_t>(1, seissol::OpenMP::threadCount()));
  if (bound.outputs().empty()) {
    // the ghost state is the inner one
    return;
  }
  impl.boundPoints.assign(Capacity, {0, 0, 0});
  impl.boundValues.assign(Capacity * (quantities.size() + bound.outputs().size()), 0);
  expr::BackendOptions options;
  options.preferred = expr::BackendKind::RtcCpu;
  for (auto& thread : impl.threads) {
    thread.table = std::make_unique<DataTable>(Capacity);
    auto& table = *thread.table;
    for (std::size_t i = 0; i < bound.inputs().size(); ++i) {
      const auto& name = bound.inputs()[i].name;
      switch (impl.inputs[i]) {
      case Input::X:
        table.bindViewConst<double>(name, Direction::In, impl.boundPoints.data()->data(), 3, 0);
        break;
      case Input::Y:
        table.bindViewConst<double>(name, Direction::In, impl.boundPoints.data()->data(), 3, 1);
        break;
      case Input::Z:
        table.bindViewConst<double>(name, Direction::In, impl.boundPoints.data()->data(), 3, 2);
        break;
      case Input::Quantity:
        table.bindViewConst<double>(
            name, Direction::In, impl.boundValues.data() + impl.inputQuantity[i] * Capacity);
        break;
      default:
        table.bindViewConst<double>(name, Direction::In, &impl.boundScalar, 0);
        break;
      }
    }
    for (std::size_t j = 0; j < bound.outputs().size(); ++j) {
      table.bindView<double>(bound.outputs()[j].name,
                             Direction::Out,
                             impl.boundValues.data() + (quantities.size() + j) * Capacity);
    }
    thread.binding = std::make_unique<expr::Binding>(expr::Binding::bind(bound, table));
    thread.kernel =
        expr::makeKernel(bound, *thread.binding, reader::datafield::sharedGridStore(), options);
    thread.kernel->precompute(table);
    // the first kernel says how it went; the others follow it quietly
    options.preferred = thread.kernel->kind();
    options.quiet = true;
    thread.inputBases.resize(bound.inputs().size());
    thread.outputBases.resize(bound.outputs().size());
  }
}

bool NonlinearDirichlet::faceAligned() const { return impl_->faceAligned; }

std::size_t NonlinearDirichlet::quantityCount() const { return impl_->quantities.size(); }

void NonlinearDirichlet::evaluate(double time,
                                  std::size_t simulation,
                                  const std::array<double, 3>* points,
                                  std::size_t count,
                                  const CellMaterialData& materialData,
                                  const double* inner,
                                  double* ghost) const {
  const auto& impl = *impl_;
  auto& thread = impl_->thread();
  const auto material = materialOf(materialData);
  const auto simulationValue = static_cast<double>(simulation);

  // what the script does not give is the inner state
  std::copy_n(inner, impl.quantities.size() * count, ghost);
  if (impl.outputQuantity.empty()) {
    return;
  }

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
      case Input::Time:
        thread.inputBases[i] = &time;
        break;
      case Input::Simulation:
        thread.inputBases[i] = &simulationValue;
        break;
      case Input::Density:
        thread.inputBases[i] = material.data();
        break;
      case Input::Mu:
        thread.inputBases[i] = material.data() + 1;
        break;
      case Input::Lambda:
        thread.inputBases[i] = material.data() + 2;
        break;
      case Input::Quantity:
        thread.inputBases[i] = inner + impl.inputQuantity[i] * count + first;
        break;
      }
    }
    for (std::size_t j = 0; j < impl.outputQuantity.size(); ++j) {
      thread.outputBases[j] = ghost + impl.outputQuantity[j] * count + first;
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

} // namespace seissol::physics
