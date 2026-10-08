// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ScriptField.h"

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Program.h"
#include "Initializer/Typedefs.h"
#include "Parallel/OpenMP.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataReader.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/ReaderBuilder.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::physics {

namespace {

using reader::scripting::DataTable;
using reader::scripting::Direction;

/// Points a kernel is bound for; a call with more runs in pieces.
constexpr std::size_t Capacity = 1024;

/// What a field reads besides the point.
enum class Input : std::uint8_t { X, Y, Z, Time, Simulation, Density, Mu, Lambda };

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

/// The material of a cell as a field sees it; NaN where a caller has none.
std::array<double, 3> materialOf(const CellMaterialData& materialData) {
  if (materialData.local == nullptr) {
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    return {NaN, NaN, NaN};
  }
  return {materialData.local->getDensity(),
          materialData.local->getMuBar(),
          materialData.local->getLambdaBar()};
}

} // namespace

struct ScriptField::Impl {
  std::string path;
  // why the script cannot be evaluated as a field, if it cannot
  std::string error;
  std::vector<std::string> quantities;
  double simulation{0};

  // a compiled script: per program input what it reads, per program output its quantity
  std::optional<expr::Program> program;
  std::vector<Input> inputs;
  std::vector<std::size_t> outputQuantity;
  // the storage the tables are bound to; every call moves the bases
  std::vector<std::array<double, 3>> boundPoints;
  std::vector<double> boundValues;
  double boundScalar{0};

  struct Thread {
    std::unique_ptr<DataTable> table;
    std::unique_ptr<expr::Binding> binding;
    std::unique_ptr<expr::Kernel> kernel;
    std::unique_ptr<reader::scripting::DataReader> reader;
    std::vector<const void*> inputBases;
    std::vector<void*> outputBases;
    std::vector<double> scratch;
  };
  std::vector<Thread> threads;

  Thread& thread() {
    const auto id = seissol::OpenMP::threadId();
    if (id >= threads.size()) {
      logError() << "A scripted field was evaluated by thread" << id << "of" << threads.size()
                 << ".";
    }
    return threads[id];
  }
};

ScriptField::ScriptField(const std::string& path,
                         std::vector<std::string> quantities,
                         std::size_t simulation,
                         bool hasTime)
    : impl_(std::make_unique<Impl>()) {
  impl_->path = path;
  impl_->quantities = std::move(quantities);
  impl_->simulation = static_cast<double>(simulation);
  impl_->threads.resize(std::max<std::size_t>(1, seissol::OpenMP::threadCount()));

  std::string reason;
  impl_->program = reader::scripting::buildProgram(path, &reason);
  if (!impl_->program.has_value()) {
    logInfo() << "The scripted field" << path
              << "is evaluated through its interpreted reader, which is slow --" << reason;
    for (auto& thread : impl_->threads) {
      thread.reader = reader::scripting::buildInterpretedReader(
          path,
          hasTime ? std::vector<std::string>{"t", "x", "y", "z"}
                  : std::vector<std::string>{"x", "y", "z"});
    }
    return;
  }

  const auto& program = *impl_->program;
  if (!program.state().empty()) {
    impl_->error =
        "the scripted field " + path + " keeps state, which a field of position and time cannot";
    return;
  }
  for (const auto& input : program.inputs()) {
    const auto meaning = inputOf(input.name);
    if (!meaning.has_value()) {
      impl_->error = "the scripted field " + path + " reads `" + input.name +
                     "`; a field evaluated at arbitrary points reads x, y, z, t, sim, rho, mu "
                     "and lambda";
      return;
    }
    impl_->inputs.push_back(*meaning);
  }
  for (const auto& output : program.outputs()) {
    const auto found = std::find(impl_->quantities.begin(), impl_->quantities.end(), output.name);
    if (found == impl_->quantities.end()) {
      impl_->error = "the scripted field " + path + " gives `" + output.name +
                     "`, which is not a quantity of the material";
      return;
    }
    impl_->outputQuantity.push_back(static_cast<std::size_t>(found - impl_->quantities.begin()));
  }

  impl_->boundPoints.assign(Capacity, {0, 0, 0});
  impl_->boundValues.assign(Capacity * program.outputs().size(), 0);

  expr::BackendOptions options;
  options.preferred = expr::BackendKind::RtcCpu;
  for (auto& thread : impl_->threads) {
    thread.table = std::make_unique<DataTable>(Capacity);
    auto& table = *thread.table;
    for (std::size_t i = 0; i < program.inputs().size(); ++i) {
      const auto& name = program.inputs()[i].name;
      switch (impl_->inputs[i]) {
      case Input::X:
        table.bindViewConst<double>(name, Direction::In, impl_->boundPoints.data()->data(), 3, 0);
        break;
      case Input::Y:
        table.bindViewConst<double>(name, Direction::In, impl_->boundPoints.data()->data(), 3, 1);
        break;
      case Input::Z:
        table.bindViewConst<double>(name, Direction::In, impl_->boundPoints.data()->data(), 3, 2);
        break;
      default:
        table.bindViewConst<double>(name, Direction::In, &impl_->boundScalar, 0);
        break;
      }
    }
    for (std::size_t j = 0; j < program.outputs().size(); ++j) {
      table.bindView<double>(
          program.outputs()[j].name, Direction::Out, impl_->boundValues.data() + j * Capacity);
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

ScriptField::~ScriptField() = default;

std::size_t ScriptField::quantityCount() const { return impl_->quantities.size(); }

bool ScriptField::compiled() const { return impl_->program.has_value(); }

std::vector<double>& ScriptField::scratch(std::size_t count) const {
  auto& values = impl_->thread().scratch;
  values.resize(count * impl_->quantities.size());
  return values;
}

void ScriptField::evaluateValues(double time,
                                 const std::array<double, 3>* points,
                                 std::size_t count,
                                 const CellMaterialData& materialData,
                                 double* values) const {
  if (!impl_->error.empty()) {
    logError() << impl_->error;
  }
  auto& thread = impl_->thread();
  const auto material = materialOf(materialData);
  std::fill_n(values, count * impl_->quantities.size(), 0.0);

  if (!impl_->program.has_value()) {
    DataTable table(count);
    table.bindViewConst<double>("x", Direction::In, points->data(), 3, 0);
    table.bindViewConst<double>("y", Direction::In, points->data(), 3, 1);
    table.bindViewConst<double>("z", Direction::In, points->data(), 3, 2);
    table.bindConstant<double>("t", time);
    table.bindConstant<double>("sim", impl_->simulation);
    table.bindConstant<double>("rho", material[0]);
    table.bindConstant<double>("mu", material[1]);
    table.bindConstant<double>("lambda", material[2]);
    // every output of the reader needs a column; those that are no quantity go to the scratch
    std::vector<double> discard(count);
    for (const auto& name : thread.reader->outputVars()) {
      const auto found = std::find(impl_->quantities.begin(), impl_->quantities.end(), name);
      table.bindView<double>(name,
                             Direction::Out,
                             found == impl_->quantities.end()
                                 ? discard.data()
                                 : values + (found - impl_->quantities.begin()) * count);
    }
    thread.reader->call(table);
    return;
  }

  for (std::size_t first = 0; first < count; first += Capacity) {
    const std::size_t pieceCount = std::min(Capacity, count - first);
    for (std::size_t i = 0; i < impl_->inputs.size(); ++i) {
      switch (impl_->inputs[i]) {
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
        thread.inputBases[i] = &impl_->simulation;
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
      }
    }
    for (std::size_t j = 0; j < impl_->outputQuantity.size(); ++j) {
      thread.outputBases[j] = values + impl_->outputQuantity[j] * count + first;
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
