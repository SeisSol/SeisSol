// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SlipRateScript.h"

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Ir.h"
#include "Expr/Program.h"
#include "Parallel/OpenMP.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/ReaderBuilder.h"

#include <algorithm>
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

namespace seissol::dr::friction_law {

namespace {

using expr::Arena;
using expr::Fn;
using expr::Kind;
using expr::NodeId;
using reader::scripting::DataTable;
using reader::scripting::DataType;
using reader::scripting::Direction;

constexpr const char* SlipStrike = "slip_strike";
constexpr const char* SlipDip = "slip_dip";
constexpr const char* SlipRateStrike = "slip_rate_strike";
constexpr const char* SlipRateDip = "slip_rate_dip";

/// `root` with every read of the channel `from` replaced by the node `to`, in the same arena: a
/// subtree that does not read `from` stays the node it is, so the program computes it once for
/// both uses.
NodeId replaceChannel(
    Arena& arena, NodeId root, int from, NodeId to, std::unordered_map<NodeId, NodeId>& memo) {
  const auto found = memo.find(root);
  if (found != memo.end()) {
    return found->second;
  }
  // a copy: interning a node below may move the nodes of the arena
  const expr::Node node = arena[root];
  NodeId result = root;
  switch (node.kind) {
  case Kind::Const:
  case Kind::Contract:
    break;
  case Kind::Field:
    if (node.ch == from) {
      result = to;
    }
    break;
  case Kind::PW: {
    const int count = expr::arity(node.fn);
    const NodeId a = replaceChannel(arena, node.a, from, to, memo);
    const NodeId b = count > 1 ? replaceChannel(arena, node.b, from, to, memo) : expr::NoNode;
    const NodeId c = count > 2 ? replaceChannel(arena, node.c, from, to, memo) : expr::NoNode;
    if (a != node.a || b != node.b || c != node.c) {
      if (count == 1) {
        result = arena.pw(node.fn, a);
      } else if (count == 2) {
        result = arena.pw(node.fn, a, b);
      } else {
        result = arena.pw(node.fn, a, b, c);
      }
    }
    break;
  }
  case Kind::Lookup: {
    std::vector<NodeId> coordinates(arena.args(node), arena.args(node) + node.argCount);
    bool changed = false;
    for (auto& coordinate : coordinates) {
      const NodeId replaced = replaceChannel(arena, coordinate, from, to, memo);
      changed |= replaced != coordinate;
      coordinate = replaced;
    }
    if (changed) {
      result = arena.lookup(node.grid, node.comp, coordinates);
    }
    break;
  }
  default:
    throw std::invalid_argument(std::string("a slip rate script cannot contain a node of kind ") +
                                expr::name(node.kind));
  }
  memo.emplace(root, result);
  return result;
}

bool readsInput(const expr::Program& program, const std::string& name) {
  return std::any_of(program.inputs().begin(),
                     program.inputs().end(),
                     [&](const expr::VarSpec& input) { return input.name == name; });
}

/// The number of points a host kernel evaluates per call.
constexpr std::size_t Capacity = 1024;

std::size_t elementSize(DataType type) {
  return type == DataType::F32 ? sizeof(float) : sizeof(double);
}

/// The backend of the device of this build.
expr::BackendKind deviceBackend() {
#if defined(SEISSOL_KERNELS_CUDA)
  return expr::BackendKind::RtcCuda;
#elif defined(SEISSOL_KERNELS_HIP)
  return expr::BackendKind::RtcHip;
#elif defined(SEISSOL_KERNELS_SYCL)
  return expr::BackendKind::RtcOpenCl;
#else
  return expr::BackendKind::RtcCpu;
#endif
}

} // namespace

SlipRateScript::SlipRateScript(const std::string& path) : path_(path) {
  std::string reason;
  auto script = reader::scripting::buildProgram(path, &reason);
  if (!script.has_value()) {
    logError() << "imposed slip rates: the script" << path << "cannot be used --" << reason
               << "; it is an sderiv module or a Lua model that traces.";
  }
  try {
    *this = SlipRateScript(*script, path);
  } catch (const std::invalid_argument& error) {
    logError() << "imposed slip rates:" << error.what();
  }
  if (cumulative_ && !readsInput(*script, Time)) {
    logWarning() << "imposed slip rates: the slip of" << path
                 << "does not depend on the time t, so nothing slips.";
  }
}

SlipRateScript::SlipRateScript(const expr::Program& script, const std::string& path) : path_(path) {
  const auto rootOf = [&](const char* name) -> std::optional<NodeId> {
    for (std::size_t i = 0; i < script.outputs().size(); ++i) {
      if (script.outputs()[i].name == name) {
        return script.roots()[i];
      }
    }
    return std::nullopt;
  };
  for (const auto& output : script.outputs()) {
    const auto& name = output.name;
    if (name != SlipStrike && name != SlipDip && name != SlipRateStrike && name != SlipRateDip) {
      std::string message = "the script ";
      message += path;
      message += " gives `";
      message += name;
      message += "`; a slip rate script gives slip_strike and slip_dip, or slip_rate_strike and "
                 "slip_rate_dip";
      throw std::invalid_argument(message);
    }
  }
  const auto slipStrike = rootOf(SlipStrike);
  const auto slipDip = rootOf(SlipDip);
  const auto rateStrike = rootOf(SlipRateStrike);
  const auto rateDip = rootOf(SlipRateDip);
  if (slipStrike.has_value() && slipDip.has_value() && !rateStrike.has_value() &&
      !rateDip.has_value()) {
    cumulative_ = true;
  } else if (rateStrike.has_value() && rateDip.has_value() && !slipStrike.has_value() &&
             !slipDip.has_value()) {
    cumulative_ = false;
  } else {
    throw std::invalid_argument("the script " + path +
                                " gives neither both of slip_strike and slip_dip nor both of "
                                "slip_rate_strike and slip_rate_dip");
  }
  if (!script.state().empty()) {
    throw std::invalid_argument("the script " + path +
                                " keeps state, which the slip at a point and time cannot");
  }
  if (!script.matrices().empty() || !script.blocks().empty()) {
    throw std::invalid_argument("the script " + path + " reads contractions");
  }

  expr::Program program;
  program.arena() = script.arena();
  program.setComputeType(script.computeType());
  for (const auto& grid : script.grids()) {
    program.internGrid(grid);
  }
  for (const auto& input : script.inputs()) {
    program.addInput(input.name, input.type);
  }
  auto& arena = program.arena();

  NodeId strike = cumulative_ ? *slipStrike : *rateStrike;
  NodeId dip = cumulative_ ? *slipDip : *rateDip;
  if (cumulative_) {
    // the increment over the sub-step, from its start to its end t, divided by its length
    NodeId strikeStart = strike;
    NodeId dipStart = dip;
    if (readsInput(script, Time)) {
      program.addInput(StartTime, DataType::F64);
      const NodeId start = arena.field(StartTime);
      const int time = arena.findChannel(Time);
      std::unordered_map<NodeId, NodeId> memo;
      strikeStart = replaceChannel(arena, strike, time, start, memo);
      dipStart = replaceChannel(arena, dip, time, start, memo);
    }
    if (!readsInput(script, TimeStep)) {
      program.addInput(TimeStep, DataType::F64);
    }
    const NodeId length = arena.field(TimeStep);
    strike = arena.pw(Fn::Div, arena.pw(Fn::Sub, strike, strikeStart), length);
    dip = arena.pw(Fn::Div, arena.pw(Fn::Sub, dip, dipStart), length);
  }

  // from strike and dip into the directions of the face, as FL 33 rotates its slip
  program.addInput(RotationCos, DataType::F64);
  program.addInput(RotationSin, DataType::F64);
  const NodeId cos = arena.field(RotationCos);
  const NodeId sin = arena.field(RotationSin);
  const NodeId rate1 =
      arena.pw(Fn::Add, arena.pw(Fn::Mul, cos, strike), arena.pw(Fn::Mul, sin, dip));
  const NodeId rate2 =
      arena.pw(Fn::Sub, arena.pw(Fn::Mul, cos, dip), arena.pw(Fn::Mul, sin, strike));
  program.addOutput(Output1, DataType::F64, rate1);
  program.addOutput(Output2, DataType::F64, rate2);
  expr::validate(program);

  for (const auto& input : program.inputs()) {
    const auto& name = input.name;
    if (name == Time || name == TimeStep || name == StartTime) {
      continue;
    }
    rowNames_.push_back(name);
    if (name == "x") {
      rowKinds_.push_back(Row::X);
    } else if (name == "y") {
      rowKinds_.push_back(Row::Y);
    } else if (name == "z") {
      rowKinds_.push_back(Row::Z);
    } else if (name == "sim") {
      rowKinds_.push_back(Row::Simulation);
    } else if (name == RotationCos) {
      rowKinds_.push_back(Row::Cos);
    } else if (name == RotationSin) {
      rowKinds_.push_back(Row::Sin);
    } else {
      rowKinds_.push_back(Row::Parameter);
    }
  }
  program_ = std::move(program);
}

struct SlipRateEvaluator::Impl {
  /// What an input of the program reads: a row, or a time of the sub-step.
  enum class Input : std::uint8_t { Row, End, Length, Start };

  struct Kernel {
    std::unique_ptr<DataTable> table;
    std::unique_ptr<expr::Binding> binding;
    std::unique_ptr<expr::Kernel> kernel;
    std::vector<const void*> inputs;
    std::vector<void*> outputs;
  };

  std::shared_ptr<const SlipRateScript> script;
  DataType slipRateType{DataType::F64};
  std::size_t timeSteps{0};
  std::size_t points{0};
  LayerData host;
  LayerData device;

  std::vector<Input> inputKinds;
  std::vector<std::size_t> inputRows;

  /// End, length and start of the sub-steps of the current step.
  std::vector<double> ends;
  std::vector<double> lengths;
  std::vector<double> starts;

  std::vector<Kernel> hostKernels;
  Kernel deviceKernel;
  bool deviceTried{false};

  /// A table of `numPoints` points over the arrays of `data`, through which a kernel is bound;
  /// a call moves the bases.
  [[nodiscard]] std::unique_ptr<DataTable> makeTable(std::size_t numPoints,
                                                     const LayerData& data) const {
    auto table = std::make_unique<DataTable>(numPoints);
    const auto& program = script->program();
    for (std::size_t i = 0; i < program.inputs().size(); ++i) {
      const auto& name = program.inputs()[i].name;
      if (inputKinds[i] == Input::Row) {
        table->bindViewConst<double>(
            name, Direction::In, data.parameters + inputRows[i] * points, 1);
      } else {
        table->bindConstant<double>(name, 0.0);
      }
    }
    for (std::size_t d = 0; d < 2; ++d) {
      const char* name = d == 0 ? SlipRateScript::Output1 : SlipRateScript::Output2;
      if (slipRateType == DataType::F32) {
        table->bindView<float>(
            name, Direction::Out, static_cast<float*>(data.slipRates) + d * points);
      } else {
        table->bindView<double>(
            name, Direction::Out, static_cast<double*>(data.slipRates) + d * points);
      }
    }
    return table;
  }

  /// The bases of the call for sub-step `step` over the points [first, first + count) of `data`.
  void setBases(Kernel& kernel, const LayerData& data, std::size_t step, std::size_t first) const {
    for (std::size_t i = 0; i < inputKinds.size(); ++i) {
      switch (inputKinds[i]) {
      case Input::Row:
        kernel.inputs[i] = data.parameters + inputRows[i] * points + first;
        break;
      case Input::End:
        kernel.inputs[i] = &ends[step];
        break;
      case Input::Length:
        kernel.inputs[i] = &lengths[step];
        break;
      case Input::Start:
        kernel.inputs[i] = &starts[step];
        break;
      }
    }
    const std::size_t width = elementSize(slipRateType);
    for (std::size_t d = 0; d < 2; ++d) {
      kernel.outputs[d] =
          static_cast<char*>(data.slipRates) + ((2 * step + d) * points + first) * width;
    }
  }

  void run(Kernel& kernel, std::size_t count, void* stream) const {
    expr::KernelArgs args;
    args.inputs = kernel.inputs.data();
    args.inputCount = kernel.inputs.size();
    args.outputs = kernel.outputs.data();
    args.outputCount = kernel.outputs.size();
    args.first = 0;
    args.count = count;
    args.stream = stream;
    kernel.kernel->run(args);
  }

  /// One kernel per thread.
  void prepareHost() {
    const std::size_t threads = std::max<std::size_t>(1, seissol::OpenMP::threadCount());
    if (hostKernels.size() >= threads) {
      return;
    }
    expr::BackendOptions options;
    options.preferred =
        hostKernels.empty() ? expr::BackendKind::RtcCpu : hostKernels.front().kernel->kind();
    options.quiet = !hostKernels.empty();
    const std::size_t made = hostKernels.size();
    hostKernels.resize(threads);
    for (std::size_t i = made; i < threads; ++i) {
      auto& kernel = hostKernels[i];
      kernel.table = makeTable(Capacity, host);
      kernel.binding =
          std::make_unique<expr::Binding>(expr::Binding::bind(script->program(), *kernel.table));
      kernel.kernel = expr::makeKernel(
          script->program(), *kernel.binding, reader::datafield::sharedGridStore(), options);
      kernel.kernel->precompute(*kernel.table);
      kernel.inputs.resize(inputKinds.size());
      kernel.outputs.resize(2);
      // the first kernel says how it went; the others follow it quietly
      options.preferred = kernel.kernel->kind();
      options.quiet = true;
    }
  }

  void evaluateHost() {
    prepareHost();
    const std::size_t chunks = (points + Capacity - 1) / Capacity;
#pragma omp parallel for schedule(static)
    for (std::size_t chunk = 0; chunk < chunks; ++chunk) {
      auto& kernel = hostKernels[seissol::OpenMP::threadId()];
      const std::size_t first = chunk * Capacity;
      const std::size_t count = std::min(Capacity, points - first);
      for (std::size_t step = 0; step < timeSteps; ++step) {
        setBases(kernel, host, step, first);
        run(kernel, count, nullptr);
      }
    }
  }

  /// Whether a kernel runs on the device; tried once.
  bool prepareDevice(void* stream) {
    if (deviceTried) {
      return deviceKernel.kernel != nullptr;
    }
    deviceTried = true;
    const auto backend = deviceBackend();
    if (backend == expr::BackendKind::RtcCpu) {
      return false;
    }
    auto table = makeTable(points, device);
    auto binding = std::make_unique<expr::Binding>(expr::Binding::bind(script->program(), *table));
    expr::BackendOptions options;
    options.preferred = backend;
    options.deviceQueue = stream;
    auto kernel = expr::makeKernel(
        script->program(), *binding, reader::datafield::sharedGridStore(), options);
    if (kernel->kind() != backend) {
      logWarning() << "imposed slip rates: the script" << script->path()
                   << "is evaluated on the host, and its slip rates are copied to the device.";
      return false;
    }
    kernel->precompute(*table);
    deviceKernel.table = std::move(table);
    deviceKernel.binding = std::move(binding);
    deviceKernel.kernel = std::move(kernel);
    deviceKernel.inputs.resize(inputKinds.size());
    deviceKernel.outputs.resize(2);
    return true;
  }
};

SlipRateEvaluator::SlipRateEvaluator(std::shared_ptr<const SlipRateScript> script,
                                     DataType slipRateType,
                                     std::size_t timeSteps)
    : impl_(std::make_unique<Impl>()) {
  impl_->script = std::move(script);
  impl_->slipRateType = slipRateType;
  impl_->timeSteps = timeSteps;
  const auto& program = impl_->script->program();
  const auto& rowNames = impl_->script->rowNames();
  for (const auto& input : program.inputs()) {
    const auto& name = input.name;
    if (name == SlipRateScript::Time) {
      impl_->inputKinds.push_back(Impl::Input::End);
    } else if (name == SlipRateScript::TimeStep) {
      impl_->inputKinds.push_back(Impl::Input::Length);
    } else if (name == SlipRateScript::StartTime) {
      impl_->inputKinds.push_back(Impl::Input::Start);
    } else {
      impl_->inputKinds.push_back(Impl::Input::Row);
    }
    const auto row = std::find(rowNames.begin(), rowNames.end(), name);
    impl_->inputRows.push_back(static_cast<std::size_t>(row - rowNames.begin()));
  }
  impl_->ends.resize(timeSteps);
  impl_->lengths.resize(timeSteps);
  impl_->starts.resize(timeSteps);
}

SlipRateEvaluator::~SlipRateEvaluator() = default;

void SlipRateEvaluator::setLayer(std::size_t points,
                                 const LayerData& host,
                                 const LayerData& device) {
  impl_->points = points;
  impl_->host = host;
  impl_->device = device;
  impl_->hostKernels.clear();
  impl_->deviceKernel = Impl::Kernel{};
  impl_->deviceTried = false;
}

bool SlipRateEvaluator::evaluate(double startTime,
                                 const std::vector<double>& deltaT,
                                 void* stream,
                                 bool onDevice) {
  auto& impl = *impl_;
  if (impl.points == 0) {
    return true;
  }
  if (deltaT.size() < impl.timeSteps) {
    logError() << "imposed slip rates: a step of" << deltaT.size() << "sub-steps, expected"
               << impl.timeSteps;
  }
  // the start of a sub-step is the very end of the previous one, so that the increments of the
  // slip add up to it exactly
  double time = startTime;
  for (std::size_t step = 0; step < impl.timeSteps; ++step) {
    impl.starts[step] = time;
    impl.lengths[step] = deltaT[step];
    time += deltaT[step];
    impl.ends[step] = time;
  }

  if (onDevice && impl.prepareDevice(stream)) {
    for (std::size_t step = 0; step < impl.timeSteps; ++step) {
      impl.setBases(impl.deviceKernel, impl.device, step, 0);
      impl.run(impl.deviceKernel, impl.points, stream);
    }
    return true;
  }
  impl.evaluateHost();
  return false;
}

} // namespace seissol::dr::friction_law
