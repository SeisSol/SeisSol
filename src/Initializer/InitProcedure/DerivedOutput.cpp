// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DerivedOutput.h"

#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Ir.h"
#include "Expr/Rewrite.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include "Parallel/OpenMP.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/ReaderBuilder.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

namespace {

namespace projection = seissol::numerical::projection;
using expr::Arena;
using expr::Fn;
using expr::NodeId;
using reader::scripting::DataTable;
using reader::scripting::DataType;
using reader::scripting::Direction;

constexpr std::array<char, 3> AxisNames{'x', 'y', 'z'};

/// The prefix of a quantity that names its time integral.
const std::string IntegralPrefix = "int_";

std::string jacobianName(std::size_t row, std::size_t column) {
  return "jinv" + std::to_string(row) + std::to_string(column);
}

std::string physicalName(std::size_t axis, const std::string& quantity) {
  return std::string("d") + AxisNames[axis] + "_" + quantity;
}

/// What an input channel of a derived program reads.
struct Channel {
  enum class Kind : std::uint8_t {
    Value,
    Reference,
    Physical,
    Jacobian,
    Coordinate,
    Time,
    TimeStep
  };
  Kind kind{Kind::Value};
  std::size_t source{0};
  /// The reference direction (Reference), the axis (Physical, Coordinate), or the entry
  /// row * 3 + column (Jacobian).
  std::size_t direction{0};
};

std::optional<std::size_t> findSource(const std::vector<DerivedSource>& sources,
                                      const std::string& name) {
  for (std::size_t i = 0; i < sources.size(); ++i) {
    if (sources[i].name == name) {
      return i;
    }
  }
  return std::nullopt;
}

std::optional<std::size_t> axisOf(char letter) {
  for (std::size_t axis = 0; axis < AxisNames.size(); ++axis) {
    if (AxisNames[axis] == letter) {
      return axis;
    }
  }
  return std::nullopt;
}

std::optional<std::size_t> directionOf(char digit) {
  if (digit >= '0' && digit <= '2') {
    return static_cast<std::size_t>(digit - '0');
  }
  return std::nullopt;
}

/// The meaning of the channel `name`, or nothing if it is outside the vocabulary. A source name
/// is looked up first, so a quantity of the material is never shadowed by a derived spelling.
std::optional<Channel> classify(const std::string& name,
                                const std::vector<DerivedSource>& sources) {
  using Kind = Channel::Kind;
  if (const auto source = findSource(sources, name)) {
    return Channel{Kind::Value, *source, 0};
  }
  if (name.size() == 1) {
    if (const auto axis = axisOf(name[0])) {
      return Channel{Kind::Coordinate, 0, *axis};
    }
  }
  if (name == "t") {
    return Channel{Kind::Time, 0, 0};
  }
  if (name == "dt") {
    return Channel{Kind::TimeStep, 0, 0};
  }
  if (name.size() == 6 && name.compare(0, 4, "jinv") == 0) {
    const auto row = directionOf(name[4]);
    const auto column = directionOf(name[5]);
    if (row.has_value() && column.has_value()) {
      return Channel{Kind::Jacobian, 0, *row * 3 + *column};
    }
  }
  if (name.size() > 3 && name.compare(name.size() - 3, 2, "_r") == 0) {
    const auto direction = directionOf(name.back());
    const auto source = findSource(sources, name.substr(0, name.size() - 3));
    if (direction.has_value() && source.has_value()) {
      return Channel{Kind::Reference, *source, *direction};
    }
  }
  if (name.size() > 3 && name[0] == 'd' && name[2] == '_') {
    const auto axis = axisOf(name[1]);
    const auto source = findSource(sources, name.substr(3));
    if (axis.has_value() && source.has_value()) {
      return Channel{Kind::Physical, *source, *axis};
    }
  }
  return std::nullopt;
}

std::string unknownChannel(const std::string& name, const std::vector<DerivedSource>& sources) {
  std::ostringstream message;
  message << "derived output: the program reads '" << name
          << "', which is neither a quantity nor derived from one. Available are the quantities";
  for (const auto& source : sources) {
    message << " " << source.name;
  }
  message << "; for a quantity q also q_r0, q_r1, q_r2 (reference derivatives) and dx_q, dy_q, "
             "dz_q; and jinv00 ... jinv22, x, y, z, t, dt.";
  return message.str();
}

/// Matrices are keyed by their source representation and derivative: 0 is the value, 1 + k the
/// derivative along reference direction k; nodal sources add 4.
constexpr std::size_t MatrixKinds = 8;

std::size_t matrixKey(bool nodal, std::optional<std::size_t> derivative) {
  return (nodal ? 4 : 0) + (derivative.has_value() ? 1 + *derivative : 0);
}

std::string matrixName(std::size_t key) {
  const std::string base = key >= 4 ? "nodal" : "modal";
  const std::size_t derivative = key % 4;
  return derivative == 0 ? base : base + "_r" + std::to_string(derivative - 1);
}

/// The projection matrices of all subcells stacked into one: row s * pointsPerSubcell + p is
/// output point p of subcell s, entry (row, m) at row + m * rows.
std::vector<double> stackedMatrix(const DerivedGeometry& geometry,
                                  std::size_t key,
                                  std::size_t cols,
                                  std::size_t rows) {
  projection::Spec spec;
  spec.order = geometry.order;
  spec.source = key >= 4 ? projection::Source::Nodal : projection::Source::Modal;
  spec.target = geometry.target;
  spec.nodalSet = geometry.nodalSet;
  if (key % 4 != 0) {
    spec.derivative = key % 4 - 1;
  }

  const std::size_t pointsPerSubcell = geometry.dataBase.size();
  std::vector<double> values(rows * cols);
  for (std::size_t subcell = 0; subcell < geometry.subcells.size(); ++subcell) {
    const auto local = projection::build<3, 3>(
        geometry.dataBase, geometry.dataOrder, geometry.subcells[subcell], spec);
    if (local.rows() != pointsPerSubcell || local.cols() != cols) {
      throw std::logic_error("derived output: a projection matrix has an unexpected shape");
    }
    for (std::size_t point = 0; point < pointsPerSubcell; ++point) {
      const std::size_t row = subcell * pointsPerSubcell + point;
      for (std::size_t m = 0; m < cols; ++m) {
        values[row + m * rows] = local(point, m);
      }
    }
  }
  return values;
}

bool startsWith(const std::string& text, const std::string& prefix) {
  return text.compare(0, prefix.size(), prefix) == 0;
}

} // namespace

DerivedProgram::DerivedProgram(expr::Program program,
                               const std::vector<DerivedSource>& sources,
                               const DerivedGeometry& geometry)
    : program_(std::move(program)), pointsPerSubcell_(geometry.dataBase.size()),
      pointsPerCell_(geometry.subcells.size() * geometry.dataBase.size()) {
  if (pointsPerCell_ == 0) {
    throw std::invalid_argument("derived output: there are no output points");
  }
  for (const auto& output : program_.outputs()) {
    if (output.type != DataType::F64) {
      throw std::invalid_argument("derived output: the output " + output.name +
                                  " is not a 64-bit floating-point value");
    }
  }

  // The channels the program reads, and what they need.
  using Kind = Channel::Kind;
  std::vector<std::pair<std::string, Channel>> channels;
  std::vector<bool> sourceUsed(sources.size(), false);
  std::array<bool, MatrixKinds> matrixUsed{};
  for (const auto& input : program_.inputs()) {
    const auto channel = classify(input.name, sources);
    if (!channel.has_value()) {
      throw std::invalid_argument(unknownChannel(input.name, sources));
    }
    const bool readsSource = channel->kind == Kind::Value || channel->kind == Kind::Reference ||
                             channel->kind == Kind::Physical;
    const bool nodal = readsSource && sources[channel->source].nodal;
    switch (channel->kind) {
    case Kind::Value:
      sourceUsed[channel->source] = true;
      matrixUsed[matrixKey(nodal, std::nullopt)] = true;
      break;
    case Kind::Reference:
      sourceUsed[channel->source] = true;
      matrixUsed[matrixKey(nodal, channel->direction)] = true;
      break;
    case Kind::Physical:
      sourceUsed[channel->source] = true;
      for (std::size_t k = 0; k < 3; ++k) {
        matrixUsed[matrixKey(nodal, k)] = true;
      }
      readsJacobian_ = true;
      break;
    case Kind::Jacobian:
      readsJacobian_ = true;
      break;
    case Kind::Coordinate:
      readsCoordinates_ = true;
      break;
    case Kind::Time:
      readsTime_ = true;
      break;
    case Kind::TimeStep:
      readsTimeStep_ = true;
      break;
    }
    channels.emplace_back(input.name, *channel);
  }

  // One matrix per representation and derivative, one block per quantity.
  const std::size_t modalCols = projection::modalSize(3, geometry.order);
  const std::size_t nodalCols = projection::nodalSize(3, geometry.order, geometry.nodalSet);
  std::array<expr::MatrixId, MatrixKinds> matrixIds{};
  matrixIds.fill(expr::NoMatrix);
  for (std::size_t key = 0; key < MatrixKinds; ++key) {
    if (!matrixUsed[key]) {
      continue;
    }
    const std::size_t cols = key >= 4 ? nodalCols : modalCols;
    matrixIds[key] = program_.internMatrix(matrixName(key),
                                           expr::MatrixShape{pointsPerCell_, cols, pointsPerCell_});
    matrices_.push_back(
        Matrix{matrixName(key), stackedMatrix(geometry, key, cols, pointsPerCell_), cols});
  }
  std::vector<expr::BlockId> blockIds(sources.size(), expr::NoBlock);
  for (std::size_t i = 0; i < sources.size(); ++i) {
    if (sourceUsed[i]) {
      blockIds[i] = program_.internBlock(sources[i].name, sources[i].nodal ? nodalCols : modalCols);
      usedSources_.push_back(i);
    }
  }

  // Values and reference derivatives become contractions; a physical derivative becomes the
  // chain rule over the reference derivatives, summed in the order of the reference directions
  // and with the derivative on the left, as the hand-written outputs did.
  std::map<std::string, expr::ChannelBuilder> builders;
  for (const auto& [name, channel] : channels) {
    const auto block = channel.kind == Kind::Value || channel.kind == Kind::Reference ||
                               channel.kind == Kind::Physical
                           ? blockIds[channel.source]
                           : expr::NoBlock;
    const bool nodal = block != expr::NoBlock && sources[channel.source].nodal;
    switch (channel.kind) {
    case Kind::Value:
    case Kind::Reference: {
      const auto matrix = matrixIds[matrixKey(nodal,
                                              channel.kind == Kind::Value
                                                  ? std::nullopt
                                                  : std::optional<std::size_t>(channel.direction))];
      builders.emplace(name,
                       [matrix, block](Arena& arena) { return arena.contract(matrix, block); });
      break;
    }
    case Kind::Physical: {
      const std::array<expr::MatrixId, 3> matrices{matrixIds[matrixKey(nodal, 0)],
                                                   matrixIds[matrixKey(nodal, 1)],
                                                   matrixIds[matrixKey(nodal, 2)]};
      const std::size_t axis = channel.direction;
      builders.emplace(name, [matrices, block, axis](Arena& arena) {
        NodeId sum = expr::NoNode;
        for (std::size_t k = 0; k < 3; ++k) {
          const NodeId reference = arena.contract(matrices[k], block);
          const NodeId jacobian = arena.field(jacobianName(k, axis));
          const NodeId term = arena.pw(Fn::Mul, reference, jacobian);
          sum = sum == expr::NoNode ? term : arena.pw(Fn::Add, sum, term);
        }
        return sum;
      });
      break;
    }
    default:
      break;
    }
  }
  if (!builders.empty()) {
    expr::substituteChannels(program_, builders);
  }

  referencePoints_.reserve(pointsPerCell_);
  for (const auto& subcell : geometry.subcells) {
    for (const auto& point : geometry.dataBase) {
      referencePoints_.push_back(subcell(point));
    }
  }
}

void DerivedProgram::bindMatrices(DataTable& table) const {
  for (const auto& matrix : matrices_) {
    table.bindMatrix<double>(
        matrix.name, matrix.values.data(), pointsPerCell_, matrix.cols, pointsPerCell_);
  }
}

bool readsTimeIntegral(const std::string& name) {
  if (startsWith(name, IntegralPrefix)) {
    return true;
  }
  return name.size() > 3 && name[0] == 'd' && axisOf(name[1]).has_value() && name[2] == '_' &&
         name.compare(3, IntegralPrefix.size(), IntegralPrefix) == 0;
}

expr::Program waveFieldProgram(const WaveFieldSelection& selection) {
  expr::Program program;
  Arena& arena = program.arena();
  const auto output = [&](const std::string& name, NodeId root) {
    program.addOutput(name, DataType::F64, root);
  };
  const auto velocity = [&](std::size_t i) {
    return selection.quantities.at(selection.velocityOffset + i);
  };

  // The inputs in the order of first use, which is what the binding sees.
  std::vector<std::string> inputs;
  const auto read = [&](const std::string& name) {
    if (std::find(inputs.begin(), inputs.end(), name) == inputs.end()) {
      inputs.push_back(name);
    }
    return arena.field(name);
  };

  for (std::size_t q = 0; q < selection.quantities.size(); ++q) {
    if (q < selection.outputMask.size() && selection.outputMask[q]) {
      output(selection.quantities[q], read(selection.quantities[q]));
    }
    if (q < selection.integrationMask.size() && selection.integrationMask[q]) {
      output("int-" + selection.quantities[q], read(IntegralPrefix + selection.quantities[q]));
    }
  }

  using Index = std::tuple<std::string, std::size_t, std::size_t>;
  if (selection.strain) {
    // (d_j v_i + d_i v_j) / 2 of the time integral of the velocity; d_i v_i on the diagonal
    for (const auto& [name, i, j] : {Index{"xx", 0, 0},
                                     Index{"yy", 1, 1},
                                     Index{"zz", 2, 2},
                                     Index{"xy", 0, 1},
                                     Index{"yz", 1, 2},
                                     Index{"xz", 0, 2}}) {
      const NodeId first = read(physicalName(j, IntegralPrefix + velocity(i)));
      if (i == j) {
        output("eps" + name, first);
      } else {
        const NodeId second = read(physicalName(i, IntegralPrefix + velocity(j)));
        output("eps" + name, arena.pw(Fn::Div, arena.pw(Fn::Add, first, second), arena.konst(2.0)));
      }
    }
  }

  if (selection.rotation) {
    // d_j v_i - d_i v_j
    for (const auto& [name, i, j] : {Index{"1", 2, 1}, Index{"2", 0, 2}, Index{"3", 1, 0}}) {
      const NodeId first = read(physicalName(j, velocity(i)));
      const NodeId second = read(physicalName(i, velocity(j)));
      output("rot" + name, arena.pw(Fn::Sub, first, second));
    }
  }

  for (std::size_t p = 0; p < selection.plasticQuantities.size(); ++p) {
    if (p < selection.plasticityMask.size() && selection.plasticityMask[p]) {
      output(selection.plasticQuantities[p], read(selection.plasticQuantities[p]));
    }
  }

  for (const auto& name : inputs) {
    program.addInput(name, DataType::F64);
  }
  expr::validate(program);
  return program;
}

expr::Program loadDerivedProgram(const std::string& path) {
  std::string reason;
  auto program = reader::scripting::buildProgram(path, &reason);
  if (!program.has_value()) {
    // there is no interpreted fallback for a program that reads contractions
    logError() << "derived output: the program" << path << "cannot be used --" << reason
               << "; a derived output is an sderiv module or a Lua model that traces.";
  }
  return std::move(*program);
}

bool waveFieldOutputReadsIntegrals(seissol::SeisSol& seissolInstance) {
  const auto& parameters = seissolInstance.parameters().output.waveFieldParameters;
  const bool masked = std::any_of(parameters.integrationMask.begin(),
                                  parameters.integrationMask.end(),
                                  [](bool integrate) { return integrate; });
  if (masked || parameters.computeStrain) {
    return true;
  }
  if (parameters.enabled && !parameters.script.empty()) {
    const auto program = loadDerivedProgram(parameters.script);
    return std::any_of(program.inputs().begin(),
                       program.inputs().end(),
                       [](const expr::VarSpec& input) { return readsTimeIntegral(input.name); });
  }
  return false;
}

namespace {

/// Which array of a cell a source lives in, and where in it.
enum class Storage : std::uint8_t { Dofs, Integrals, PlasticStrain };

struct SourceLocation {
  Storage storage{Storage::Dofs};
  /// In reals, from the start of the array of the cell.
  std::size_t offset{0};
};

/// A written cell of the configuration, where it is stored, and which one it is.
struct WrittenCell {
  std::size_t color{0};
  std::size_t cell{0};
  std::size_t index{0};
};

/// Cells per work item of an evaluation: enough to amortise a call, few enough to balance.
constexpr std::size_t ChunkCells = 64;

template <typename Cfg>
class DerivedVolumeOutput final : public DerivedOutput {
  public:
  using RealT = Real<Cfg>;

  DerivedVolumeOutput(seissol::SeisSol& seissolInstance,
                      const std::vector<std::size_t>& cells,
                      std::vector<WrittenCell> written,
                      const DerivedGeometry& geometry,
                      const expr::Program& program);

  [[nodiscard]] const std::vector<std::string>& names() const override { return names_; }

  [[nodiscard]] bool accumulates() const override { return !derived_.program().state().empty(); }

  void write(double time) override {
    const std::lock_guard lock(mutex_);
    std::vector<std::size_t> pending;
    for (std::size_t i = 0; i < runs_.size(); ++i) {
      if (!accumulates() || !runs_[i].stepped) {
        pending.push_back(i);
      }
      runs_[i].stepped = false;
    }
    evaluate(pending, time);
  }

  void step(std::size_t layerId, double time) override {
    const std::lock_guard lock(mutex_);
    const auto found = runOfLayer_.find(layerId);
    if (found == runOfLayer_.end()) {
      return;
    }
    runs_[found->second].stepped = true;
    evaluate({found->second}, time);
  }

  void copy(std::size_t output,
            double* target,
            std::size_t index,
            std::size_t subcell) const override {
    const auto& values = simulations_[output / outputCount_].values;
    const std::size_t pointsPerSubcell = derived_.pointsPerSubcell();
    const std::size_t first = (output % outputCount_) * numPoints_ +
                              location_.at(index) * derived_.pointsPerCell() +
                              subcell * pointsPerSubcell;
    std::copy_n(values.data() + first, pointsPerSubcell, target);
  }

  private:
  /// The written cells of one layer: a contiguous range of the point set.
  struct LayerRun {
    std::size_t layer{0};
    std::size_t firstCell{0};
    std::size_t cellCount{0};
    /// Per simulation, the bases of the blocks of the program in this layer.
    std::vector<std::vector<const void*>> blocks;
    /// The time of the previous evaluation, and the time and time step of the current one.
    std::optional<double> lastTime;
    double time{0};
    double timeStep{0};
    bool stepped{false};
  };

  struct Simulation {
    std::unique_ptr<DataTable> table;
    std::unique_ptr<expr::Binding> binding;
    /// One per thread: a kernel is not to be run from two threads at once.
    std::vector<std::unique_ptr<expr::Kernel>> kernels;
    /// Output j at point p: values[j * numPoints + p].
    std::vector<double> values;
  };

  /// Every quantity a program of this configuration can read; the time integral and the
  /// plastic strain are checked to be stored only once a program reads them.
  static std::vector<DerivedSource> sourcesOf(bool plasticity,
                                              std::vector<SourceLocation>& locations);

  const void* blockBase(const SourceLocation& location, std::size_t color, std::size_t simulation);

  void evaluate(const std::vector<std::size_t>& runs, double time);

  seissol::SeisSol& seissolInstance_;
  std::vector<SourceLocation> locations_;
  DerivedProgram derived_;
  std::vector<std::string> names_;
  std::size_t outputCount_{0};
  std::size_t numPoints_{0};
  /// Written cell -> cell of the point set; only read for the cells of this configuration.
  std::vector<std::size_t> location_;
  /// Cell of the point set -> cell in its layer.
  std::vector<std::uint32_t> cellIndex_;
  /// Per cell of the point set: the inverse Jacobian, row-major, and the cell transform as its
  /// offset followed by its matrix, column by column.
  std::vector<double> jacobians_;
  std::vector<double> transforms_;
  std::vector<LayerRun> runs_;
  std::map<std::size_t, std::size_t> runOfLayer_;
  std::vector<Simulation> simulations_;
  std::optional<std::size_t> timeInput_;
  std::optional<std::size_t> timeStepInput_;
  // the bases the time columns are bound with; every call moves them to its run
  double boundTime_{0};
  double boundTimeStep_{0};
  std::mutex mutex_;
};

template <typename Cfg>
std::vector<DerivedSource>
    DerivedVolumeOutput<Cfg>::sourcesOf(bool plasticity, std::vector<SourceLocation>& locations) {
  using MaterialT = model::MaterialOf<Cfg>;
  // the stride between two quantities, i.e. the padded coefficients of one
  constexpr auto ModalStride =
      tensor::Q<Cfg>::Size / tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];
  constexpr auto NodalStride = tensor::QStressNodal<Cfg>::Size /
                               tensor::QStressNodal<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];

  std::vector<DerivedSource> sources;
  for (std::size_t q = 0; q < MaterialT::Quantities.size(); ++q) {
    sources.push_back(DerivedSource{MaterialT::Quantities[q], false});
    locations.push_back(SourceLocation{Storage::Dofs, q * ModalStride});
  }
  for (std::size_t q = 0; q < MaterialT::Quantities.size(); ++q) {
    sources.push_back(DerivedSource{IntegralPrefix + MaterialT::Quantities[q], false});
    locations.push_back(SourceLocation{Storage::Integrals, q * ModalStride});
  }
  if (plasticity) {
    for (std::size_t p = 0; p < model::PlasticityData<Cfg>::Quantities.size(); ++p) {
      sources.push_back(DerivedSource{model::PlasticityData<Cfg>::Quantities[p], true});
      locations.push_back(SourceLocation{Storage::PlasticStrain, p * NodalStride});
    }
  }
  return sources;
}

template <typename Cfg>
const void* DerivedVolumeOutput<Cfg>::blockBase(const SourceLocation& location,
                                                std::size_t color,
                                                std::size_t simulation) {
  auto& layer = seissolInstance_.memoryManager().ltsStorage().layer(color);
  const RealT* base = nullptr;
  switch (location.storage) {
  case Storage::Dofs:
    base = reinterpret_cast<const RealT*>(layer.template var<LTS::Dofs>(Cfg()));
    break;
  case Storage::Integrals:
    base = reinterpret_cast<const RealT*>(layer.template var<LTS::Integrals>(Cfg()));
    break;
  case Storage::PlasticStrain:
    base = reinterpret_cast<const RealT*>(layer.template var<LTS::PStrain>(Cfg()));
    break;
  }
  // the simulation is the fastest dimension of the fused layout
  return base == nullptr ? nullptr : base + location.offset + simulation;
}

template <typename Cfg>
DerivedVolumeOutput<Cfg>::DerivedVolumeOutput(seissol::SeisSol& seissolInstance,
                                              const std::vector<std::size_t>& cells,
                                              std::vector<WrittenCell> written,
                                              const DerivedGeometry& geometry,
                                              const expr::Program& program)
    : seissolInstance_(seissolInstance),
      derived_(
          program, sourcesOf(seissolInstance.parameters().model.plasticity, locations_), geometry) {
  constexpr std::size_t Simulations = Cfg::NumSimulations;
  const auto& meshReader = seissolInstance.meshReader();
  const std::size_t pointsPerCell = derived_.pointsPerCell();

  for (std::size_t sim = 0; sim < Simulations; ++sim) {
    for (const auto& output : derived_.program().outputs()) {
      names_.push_back(multisim::MultisimHelperWrapper<Cfg>::MultisimEnabled
                           ? output.name + "-" + std::to_string(sim + 1)
                           : output.name);
    }
  }
  outputCount_ = derived_.program().outputs().size();

  // The point set: the written cells of this configuration, layer by layer, in memory order.
  std::sort(written.begin(), written.end(), [](const WrittenCell& a, const WrittenCell& b) {
    return std::tie(a.color, a.cell) < std::tie(b.color, b.cell);
  });
  numPoints_ = written.size() * pointsPerCell;
  location_.assign(cells.size(), std::numeric_limits<std::size_t>::max());
  cellIndex_.reserve(written.size());
  for (std::size_t i = 0; i < written.size(); ++i) {
    location_[written[i].index] = i;
    cellIndex_.push_back(static_cast<std::uint32_t>(written[i].cell));
    if (runs_.empty() || runs_.back().layer != written[i].color) {
      runOfLayer_[written[i].color] = runs_.size();
      runs_.emplace_back();
      runs_.back().layer = written[i].color;
      runs_.back().firstCell = i;
    }
    ++runs_.back().cellCount;
  }

  // The bases of the blocks, per layer and simulation.
  const auto& sources = derived_.usedSources();
  for (auto& run : runs_) {
    run.blocks.resize(Simulations);
    for (std::size_t sim = 0; sim < Simulations; ++sim) {
      for (std::size_t b = 0; b < sources.size(); ++b) {
        const void* base = blockBase(locations_[sources[b]], run.layer, sim);
        if (base == nullptr) {
          logError() << "derived output: the program reads"
                     << derived_.program().blocks()[b].name.c_str()
                     << ", which this run does not keep (the time integral of the solution is "
                        "only kept when an output reads it, the plastic strain only with "
                        "plasticity).";
        }
        run.blocks[sim].push_back(base);
      }
    }
  }

  if (derived_.readsJacobian() || derived_.readsCoordinates()) {
    const auto barycenter =
        seissol::geometry::CellTransform::VectorEigenT(Cell::ReferenceBarycenter.data());
    jacobians_.resize(written.size() * 9);
    transforms_.resize(written.size() * 12);
    for (std::size_t i = 0; i < written.size(); ++i) {
      const auto transform =
          seissol::geometry::AffineTransform::fromMeshCell(cells[written[i].index], meshReader);
      // the transform is affine, so its inverse Jacobian is the one at the barycenter
      const auto inverse = transform.refToSpaceJacobianInverse(barycenter);
      for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
          jacobians_[i * 9 + row * 3 + column] = inverse(row, column);
        }
      }
      const auto origin = transform.refToSpace(std::array<double, 3>{0, 0, 0});
      for (std::size_t axis = 0; axis < 3; ++axis) {
        transforms_[i * 12 + axis] = origin[axis];
      }
      for (std::size_t direction = 0; direction < 3; ++direction) {
        std::array<double, 3> unit{0, 0, 0};
        unit[direction] = 1;
        const auto image = transform.refToSpace(unit);
        for (std::size_t axis = 0; axis < 3; ++axis) {
          transforms_[i * 12 + 3 + direction * 3 + axis] = image[axis] - origin[axis];
        }
      }
    }
  }

  const auto& inputs = derived_.program().inputs();
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    if (inputs[i].name == "t") {
      timeInput_ = i;
    }
    if (inputs[i].name == "dt") {
      timeStepInput_ = i;
    }
  }

  // One table, binding and set of kernels per simulation: the state of a point is its own in
  // every simulation.
  const std::size_t threads = std::max<std::size_t>(1, seissol::OpenMP::threadCount());
  auto& grids = reader::datafield::sharedGridStore();
  simulations_.resize(Simulations);
  for (std::size_t sim = 0; sim < Simulations; ++sim) {
    auto& simulation = simulations_[sim];
    simulation.values.assign(outputCount_ * numPoints_, 0.0);
    simulation.table = std::make_unique<DataTable>(numPoints_);
    DataTable& table = *simulation.table;

    derived_.bindMatrices(table);
    for (std::size_t b = 0; b < sources.size(); ++b) {
      const auto& location = locations_[sources[b]];
      const std::size_t cellStride =
          location.storage == Storage::PlasticStrain
              ? tensor::QStressNodal<Cfg>::size() + tensor::QEtaNodal<Cfg>::size()
              : tensor::Q<Cfg>::size();
      reader::scripting::BlockInput block;
      block.base = runs_.front().blocks[sim][b];
      block.cellStride = cellStride * sizeof(RealT);
      block.modeStride = Simulations * sizeof(RealT);
      block.cellIndex = cellIndex_.data();
      block.type = reader::scripting::DataTypeTraits<RealT>::Type;
      block.length = derived_.program().blocks()[b].length;
      table.bindBlock(derived_.program().blocks()[b].name, block);
    }
    if (derived_.readsJacobian()) {
      for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
          table.bindCellView<double>(
              jacobianName(row, column), jacobians_.data(), pointsPerCell, 9, row * 3 + column);
        }
      }
    }
    if (derived_.readsCoordinates()) {
      for (std::size_t axis = 0; axis < 3; ++axis) {
        table.bindComputedBatch<double>(
            std::string(1, AxisNames[axis]),
            [this, axis, pointsPerCell](std::size_t first, std::size_t count, double* out) {
              const auto& reference = derived_.referencePoints();
              for (std::size_t i = 0; i < count; ++i) {
                const std::size_t point = first + i;
                const double* transform = transforms_.data() + (point / pointsPerCell) * 12;
                const auto& xi = reference[point % pointsPerCell];
                out[i] = transform[axis] + transform[3 + axis] * xi[0] +
                         transform[6 + axis] * xi[1] + transform[9 + axis] * xi[2];
              }
            });
      }
    }
    if (timeInput_.has_value()) {
      table.bindViewConst<double>("t", Direction::In, &boundTime_, 0);
    }
    if (timeStepInput_.has_value()) {
      table.bindViewConst<double>("dt", Direction::In, &boundTimeStep_, 0);
    }
    for (std::size_t j = 0; j < outputCount_; ++j) {
      table.bindView<double>(derived_.program().outputs()[j].name,
                             Direction::Out,
                             simulation.values.data() + j * numPoints_);
    }

    simulation.binding =
        std::make_unique<expr::Binding>(expr::Binding::bind(derived_.program(), table));

    // The first kernel says which backend it got and how it went; the others, made for the same
    // program, follow it quietly.
    expr::BackendOptions options;
    options.preferred = expr::BackendKind::RtcCpu;
    for (std::size_t thread = 0; thread < threads; ++thread) {
      auto kernel = expr::makeKernel(derived_.program(), *simulation.binding, grids, options);
      kernel->precompute(table);
      options.preferred = kernel->kind();
      options.quiet = true;
      simulation.kernels.push_back(std::move(kernel));
    }
  }
}

template <typename Cfg>
void DerivedVolumeOutput<Cfg>::evaluate(const std::vector<std::size_t>& runs, double time) {
  struct WorkItem {
    std::size_t run;
    std::size_t simulation;
    std::size_t firstCell;
    std::size_t cellCount;
  };
  std::vector<WorkItem> items;
  for (const auto r : runs) {
    auto& run = runs_[r];
    run.time = time;
    run.timeStep = run.lastTime.has_value() ? time - *run.lastTime : 0.0;
    run.lastTime = time;
    for (std::size_t sim = 0; sim < simulations_.size(); ++sim) {
      for (std::size_t cell = 0; cell < run.cellCount; cell += ChunkCells) {
        items.push_back(
            WorkItem{r, sim, run.firstCell + cell, std::min(ChunkCells, run.cellCount - cell)});
      }
    }
  }

  const std::size_t pointsPerCell = derived_.pointsPerCell();
  const std::size_t inputCount = derived_.program().inputs().size();
  // at most one thread per kernel
  const int threads = static_cast<int>(simulations_.front().kernels.size());
#pragma omp parallel num_threads(threads)
  {
    std::vector<const void*> inputs(inputCount, nullptr);
    const std::size_t thread = seissol::OpenMP::threadId();
#pragma omp for schedule(dynamic)
    for (std::size_t i = 0; i < items.size(); ++i) {
      const auto& item = items[i];
      const auto& run = runs_[item.run];
      auto& simulation = simulations_[item.simulation];
      if (timeInput_.has_value()) {
        inputs[*timeInput_] = &run.time;
      }
      if (timeStepInput_.has_value()) {
        inputs[*timeStepInput_] = &run.timeStep;
      }
      expr::KernelArgs args;
      args.inputs = inputs.data();
      args.inputCount = inputs.size();
      args.blocks = run.blocks[item.simulation].data();
      args.blockCount = run.blocks[item.simulation].size();
      args.first = item.firstCell * pointsPerCell;
      args.count = item.cellCount * pointsPerCell;
      simulation.kernels[thread]->run(args);
    }
  }
}

} // namespace

template <typename Cfg>
std::shared_ptr<DerivedOutput>
    makeDerivedVolumeOutput(seissol::SeisSol& seissolInstance,
                            const std::shared_ptr<const std::vector<std::size_t>>& cells,
                            const DerivedGeometry& geometry,
                            const expr::Program& program) {
  auto& storage = seissolInstance.memoryManager().ltsStorage();
  auto& backmap = seissolInstance.memoryManager().backmap();
  std::vector<WrittenCell> written;
  for (std::size_t index = 0; index < cells->size(); ++index) {
    const auto position = backmap.get(cells->at(index));
    if (storage.layer(position.color).getIdentifier().config == configIdOf<Cfg>()) {
      written.push_back(WrittenCell{position.color, position.cell, index});
    }
  }
  if (written.empty()) {
    return nullptr;
  }
  try {
    return std::make_shared<DerivedVolumeOutput<Cfg>>(
        seissolInstance, *cells, std::move(written), geometry, program);
  } catch (const std::invalid_argument& error) {
    logError() << error.what();
  }
  return nullptr;
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template std::shared_ptr<DerivedOutput> makeDerivedVolumeOutput<Cfg>(                            \
      seissol::SeisSol&,                                                                           \
      const std::shared_ptr<const std::vector<std::size_t>>&,                                      \
      const DerivedGeometry&,                                                                      \
      const expr::Program&);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::initializer
