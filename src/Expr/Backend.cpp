// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Expr/Backend.h"

#include "Expr/Binding.h"
#include "Expr/Cost.h"
#include "Expr/Expand.h"
#include "Expr/Interp.h"
#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/RtcCpu.h"
#include "Expr/RtcGpuDriver.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

using reader::datafield::Grid;
using reader::datafield::GridStore;
using reader::scripting::DataTable;

// --- grid sampling seam -----------------------------------------------------

/// Interp.h's GridArrays over a GridStore: the array of each grid of a program (by GridId), read
/// where the instructions of a written-out lookup ask for it, which is where its window currently
/// lies.
class StoreArrays final : public GridArrays {
  public:
  StoreArrays(const Program& program, std::vector<std::size_t> storeIds, GridStore& store)
      : program_(&program), store_(&store), storeIds_(std::move(storeIds)) {}

  [[nodiscard]] const reader::datafield::ArrayView& array(GridId grid) const override {
    const reader::datafield::ArrayView* view =
        store_->get(storeIds_.at(static_cast<std::size_t>(grid))).arrayView();
    if (view == nullptr) {
      logError() << "expr: the grid" << pathOf(grid) << "keeps no array to sample.";
    }
    return *view;
  }

  void report(GridId grid, CheckKind kind) const override { reportCheck(kind, pathOf(grid)); }

  void failed(std::int32_t check,
              std::size_t point,
              const std::vector<double>& inputs) const override {
    std::string where;
    for (std::size_t i = 0; i < inputs.size() && i < program_->inputs().size(); ++i) {
      where +=
          (i == 0 ? "" : ", ") + program_->inputs()[i].name + " = " + std::to_string(inputs[i]);
    }
    logError() << "expr:" << program_->checks().at(static_cast<std::size_t>(check)).message
               << "-- point" << point << "of the call does not meet the check of the program"
               << (where.empty() ? std::string() : "(" + where + ")");
  }

  private:
  [[nodiscard]] std::string pathOf(GridId grid) const {
    return program_->grids().at(static_cast<std::size_t>(grid)).path;
  }

  const Program* program_;
  GridStore* store_;
  std::vector<std::size_t> storeIds_; // GridId -> GridStore index
};

// --- table seam -------------------------------------------------------------

// Interp.h's TileIo over raw bases, for run(KernelArgs). The permutation and
// every stride still come from the Binding; only the bases move.
template <typename T>
class RawTileIo final : public TileIo<T> {
  public:
  RawTileIo(const Binding& binding, const KernelArgs& args) : binding_(&binding), args_(&args) {}

  void gather(std::size_t first, std::size_t count, T* dst) const override {
    binding_->gatherFrom(args_->inputs, args_->inputCount, first, count, dst);
  }

  void scatter(std::size_t first, std::size_t count, const T* src) override {
    binding_->scatterTo(args_->outputs, args_->outputCount, first, count, src);
  }

  [[nodiscard]] ContractOperands contraction(MatrixId matrix, BlockId block) const override {
    const auto m = static_cast<std::size_t>(matrix);
    const auto b = static_cast<std::size_t>(block);
    return binding_->contraction(
        matrix,
        block,
        m < args_->matrixCount && args_->matrices != nullptr ? args_->matrices[m] : nullptr,
        b < args_->blockCount && args_->blocks != nullptr ? args_->blocks[b] : nullptr);
  }

  [[nodiscard]] const std::size_t* points(std::size_t first) const override {
    return binding_->permutation().empty() ? nullptr : binding_->permutation().data() + first;
  }

  void gatherState(std::size_t first, std::size_t count, T* dst) const override {
    binding_->gatherState(args_->states, args_->stateCount, first, count, dst);
  }

  void scatterState(std::size_t first, std::size_t count, const T* src) override {
    binding_->scatterState(args_->states, args_->stateCount, first, count, src);
  }

  private:
  const Binding* binding_;
  const KernelArgs* args_;
};

// Interp.h's TileIo over a (Binding, DataTable) pair. Two lines, as Interp.h
// predicted: Binding already produces exactly the dst[channel * count + lane]
// layout the interpreter wants, and already resolves the permutation.
template <typename T>
class BoundTileIo final : public TileIo<T> {
  public:
  BoundTileIo(const Binding& binding, const DataTable& table)
      : binding_(&binding), table_(&table) {}

  void gather(std::size_t first, std::size_t count, T* dst) const override {
    binding_->gather(*table_, first, count, dst);
  }

  void scatter(std::size_t first, std::size_t count, const T* src) override {
    binding_->scatter(*table_, first, count, src);
  }

  [[nodiscard]] ContractOperands contraction(MatrixId matrix, BlockId block) const override {
    return binding_->contraction(matrix, block);
  }

  [[nodiscard]] const std::size_t* points(std::size_t first) const override {
    return binding_->permutation().empty() ? nullptr : binding_->permutation().data() + first;
  }

  void gatherState(std::size_t first, std::size_t count, T* dst) const override {
    binding_->gatherState(nullptr, 0, first, count, dst);
  }

  void scatterState(std::size_t first, std::size_t count, const T* src) override {
    binding_->scatterState(nullptr, 0, first, count, src);
  }

  private:
  const Binding* binding_;
  const DataTable* table_;
};

// --- interpreter kernel -----------------------------------------------------

std::vector<PointRange> toPointRanges(const Binding& binding) {
  std::vector<PointRange> ranges;
  ranges.reserve(binding.groupRanges().size());
  for (const auto& range : binding.groupRanges()) {
    ranges.push_back(PointRange{range.begin, range.end});
  }
  return ranges;
}

template <typename T>
class InterpreterKernel final : public Kernel {
  public:
  InterpreterKernel(const Program& program,
                    Binding& binding,
                    LoweredProgram lowered,
                    std::unique_ptr<StoreArrays> grids,
                    const InterpreterOptions& options)
      : binding_(&binding), lowered_(std::move(lowered)), grids_(std::move(grids)),
        interpreter_(program, lowered_, grids_.get(), options), partitions_(toPointRanges(binding)),
        needsPrecompute_(lowered_.hasPrecompute()) {}

  void precompute(const DataTable& table) override {
    if (!needsPrecompute_) {
      return;
    }
    const BoundTileIo<T> io(*binding_, table);
    interpreter_.precompute(io, binding_->numPoints(), persistent(), partitions_);
    precomputed_ = true;
  }

  void run(const KernelArgs& args) override {
    if (!binding_->hostAddressable()) {
      logError() << "expr: this program has a computed column, so it cannot be evaluated from raw "
                    "bases; call run(table) instead.";
      return;
    }
    guardPrecompute();
    RawTileIo<T> io(*binding_, args);
    // One range, exactly what was asked for. The group partitioning is NOT
    // applied here: the caller has already chosen the range, and re-imposing
    // partitions over it would evaluate points outside it.
    const std::vector<PointRange> range{PointRange{args.first, args.first + args.count}};
    interpreter_.run(io, binding_->numPoints(), persistent(), range);
  }

  void run(const DataTable& table) override {
    guardPrecompute();
    BoundTileIo<T> io(*binding_, table);
    interpreter_.run(io, binding_->numPoints(), persistent(), partitions_);
  }

  [[nodiscard]] BackendKind kind() const override { return BackendKind::Interpreter; }

  private:
  void guardPrecompute() const {
    if (needsPrecompute_ && !precomputed_) {
      // Not a warning: the hoisted slots would be read as whatever the
      // allocation zeroed them to, which is a plausible-looking wrong answer
      // rather than a crash.
      logError() << "expr: the kernel has a precompute stage that was never run; call "
                    "Kernel::precompute() from prepare().";
    }
  }

  T* persistent();

  Binding* binding_;
  LoweredProgram lowered_;
  std::unique_ptr<StoreArrays> grids_;
  TileInterpreter<T> interpreter_;
  std::vector<PointRange> partitions_;
  bool needsPrecompute_{false};
  bool precomputed_{false};
};

template <>
double* InterpreterKernel<double>::persistent() {
  return binding_->persistentF64();
}
template <>
float* InterpreterKernel<float>::persistent() {
  return binding_->persistentF32();
}

/// The grids of `store` at `ids` (by GridId) checked against how `program` samples them
/// (gridUses()): an array each, with the dimensions and the components its lookups have. A
/// written-out lookup reads the elements its offsets name, so a mismatch is an error before the
/// first sample, and the configuration's: the script and the file do not agree.
void checkGrids(const Program& program, GridStore& store, const std::vector<std::size_t>& ids) {
  const std::vector<GridUse> uses = gridUses(program);
  for (std::size_t grid = 0; grid < uses.size(); ++grid) {
    const GridUse& use = uses[grid];
    if (use.dimensions == 0) {
      continue;
    }
    const std::string& path = program.grids()[grid].path;
    const reader::datafield::ArrayView* view = store.get(ids.at(grid)).arrayView();
    if (view == nullptr) {
      logError() << "expr: the grid" << path << "keeps no array to sample.";
      continue;
    }
    if (static_cast<std::int32_t>(view->dimensions) != use.dimensions) {
      logError() << "expr: a lookup on the grid" << path << "passes" << use.dimensions
                 << "coordinates, but the grid has" << view->dimensions << "dimensions.";
    }
    if (use.components > static_cast<std::int32_t>(view->components)) {
      logError() << "expr: a lookup on the grid" << path << "requests component"
                 << use.components - 1 << "of" << view->components << ".";
    }
  }
}

/// The same options with the architecture dropped, for handing to a DIFFERENT
/// backend than the one `arch` was meant for. See the note on BackendOptions::arch.
BackendOptions forFallback(const BackendOptions& options) {
  BackendOptions fallback = options;
  fallback.arch.clear();
  return fallback;
}

std::unique_ptr<Kernel> makeInterpreter(const Program& program,
                                        Binding& binding,
                                        GridStore& grids,
                                        const std::vector<std::size_t>& ids,
                                        const BackendOptions& options) {
  // once, at its first use: the reference has to compute as IEEE 754 says
  static const bool Checked = [] {
    const std::string fault = interpreterFault();
    if (!fault.empty()) {
      logError() << "expr: the interpreter of the scripts is built with optimizations that change "
                    "values ("
                 << fault.c_str()
                 << "); build src/Expr without -ffast-math, with -fp-model=precise for icpx.";
    }
    return true;
  }();
  static_cast<void>(Checked);
  LoweredProgram lowered = lower(program, options.lowering);
  if (!options.quiet) {
    logInfo()
        << "expr: interpreter kernel --" << lowered.summary().c_str() << "--"
        << cost(program, lowered, program.computeType()).summary(program.computeType()).c_str();
  }

  binding.allocatePersistent(program, lowered.persistentSlotCount());
  // the arrays of the grids, and who receives what the checks find
  auto arrays = ids.empty() && program.checks().empty()
                    ? nullptr
                    : std::make_unique<StoreArrays>(program, ids, grids);

  InterpreterOptions interpreterOptions;
  interpreterOptions.tileSize = options.tileSize;

  // in the type the data moves in, which the instructions compute from in theirs
  if (lowered.dataType() == ComputeType::F32) {
    return std::make_unique<InterpreterKernel<float>>(
        program, binding, std::move(lowered), std::move(arrays), interpreterOptions);
  }
  return std::make_unique<InterpreterKernel<double>>(
      program, binding, std::move(lowered), std::move(arrays), interpreterOptions);
}

// --- the check of a compiled backend -----------------------------------------

/// What the check of the first call of a compiled backend does where it deviates from the
/// interpreter: SEISSOL_EXPR_CHECK=error (the default) stops the run, `warn` warns, `off` skips the
/// check.
enum class CheckPolicy : std::uint8_t { Off, Warn, Error };

CheckPolicy checkPolicy() {
  const char* env = std::getenv("SEISSOL_EXPR_CHECK");
  const std::string value = env != nullptr ? env : "";
  if (value == "off") {
    return CheckPolicy::Off;
  }
  if (value == "warn") {
    return CheckPolicy::Warn;
  }
  if (!value.empty() && value != "error") {
    logWarning() << "SEISSOL_EXPR_CHECK is" << value
                 << "; error, warn or off decide what a compiled kernel that deviates from the "
                    "interpreter at its first call does.";
  }
  return CheckPolicy::Error;
}

std::vector<std::size_t> spreadPoints(std::size_t first, std::size_t count) {
  constexpr std::size_t Sample = 64;
  std::vector<std::size_t> points;
  if (count <= Sample) {
    for (std::size_t i = 0; i < count; ++i) {
      points.push_back(first + i);
    }
    return points;
  }
  for (std::size_t k = 0; k < Sample; ++k) {
    points.push_back(first + k * (count - 1) / (Sample - 1));
  }
  return points;
}

/// Interp.h's TileIo over the points of a sample: the inputs each point reads, as the binding reads
/// them for the call (from its table or from the bases of its arguments), and the outputs, kept.
template <typename T>
class SampleTileIo final : public TileIo<T> {
  public:
  SampleTileIo(const Binding& binding,
               const DataTable* table,
               const KernelArgs* args,
               const std::vector<std::size_t>& points,
               std::vector<double>& outputs,
               std::size_t outputCount)
      : binding_(&binding), table_(table), args_(args), points_(points), outputs_(&outputs),
        outputCount_(outputCount), inputs_(binding.inputs().size()) {
    // the point each lane is, past the permutation, as a contraction reads it
    for (const std::size_t point : points) {
      permuted_.push_back(binding.permutation().empty() ? point : binding.permutation()[point]);
    }
  }

  void gather(std::size_t first, std::size_t count, T* dst) const override {
    for (std::size_t lane = 0; lane < count; ++lane) {
      const std::size_t point = points_[first + lane];
      if (table_ != nullptr) {
        binding_->gather(*table_, point, 1, inputs_.data());
      } else {
        binding_->gatherFrom(args_->inputs, args_->inputCount, point, 1, inputs_.data());
      }
      for (std::size_t i = 0; i < inputs_.size(); ++i) {
        dst[i * count + lane] = inputs_[i];
      }
    }
  }

  void scatter(std::size_t first, std::size_t count, const T* src) override {
    for (std::size_t o = 0; o < outputCount_; ++o) {
      for (std::size_t lane = 0; lane < count; ++lane) {
        (*outputs_)[o * points_.size() + first + lane] = static_cast<double>(src[o * count + lane]);
      }
    }
  }

  [[nodiscard]] ContractOperands contraction(MatrixId matrix, BlockId block) const override {
    if (args_ == nullptr) {
      return binding_->contraction(matrix, block);
    }
    const auto m = static_cast<std::size_t>(matrix);
    const auto b = static_cast<std::size_t>(block);
    return binding_->contraction(
        matrix,
        block,
        m < args_->matrixCount && args_->matrices != nullptr ? args_->matrices[m] : nullptr,
        b < args_->blockCount && args_->blocks != nullptr ? args_->blocks[b] : nullptr);
  }

  [[nodiscard]] const std::size_t* points(std::size_t first) const override {
    return permuted_.data() + first;
  }

  private:
  const Binding* binding_;
  const DataTable* table_;
  const KernelArgs* args_;
  const std::vector<std::size_t>& points_;
  std::vector<std::size_t> permuted_;
  std::vector<double>* outputs_;
  std::size_t outputCount_;
  mutable std::vector<T> inputs_;
};

/// Interp.h's TileIo over the values of the points of a sample: their inputs, given, and their
/// outputs, kept -- for what a backend on a device read back.
template <typename T>
class ValuesTileIo final : public TileIo<T> {
  public:
  ValuesTileIo(const std::vector<double>& inputs,
               std::size_t inputCount,
               std::size_t points,
               std::vector<double>& outputs,
               std::size_t outputCount)
      : inputs_(&inputs), inputCount_(inputCount), points_(points), outputs_(&outputs),
        outputCount_(outputCount) {}

  void gather(std::size_t first, std::size_t count, T* dst) const override {
    for (std::size_t i = 0; i < inputCount_; ++i) {
      for (std::size_t lane = 0; lane < count; ++lane) {
        dst[i * count + lane] = static_cast<T>((*inputs_)[i * points_ + first + lane]);
      }
    }
  }

  void scatter(std::size_t first, std::size_t count, const T* src) override {
    for (std::size_t o = 0; o < outputCount_; ++o) {
      for (std::size_t lane = 0; lane < count; ++lane) {
        (*outputs_)[o * points_ + first + lane] = static_cast<double>(src[o * count + lane]);
      }
    }
  }

  private:
  const std::vector<double>* inputs_;
  std::size_t inputCount_;
  std::size_t points_;
  std::vector<double>* outputs_;
  std::size_t outputCount_;
};

/// A compiled kernel, checked at its first call (the specification's check of a backend): at a
/// sample of the points of the call, the interpreter evaluates the program as written -- nothing
/// hoisted, so that a hoisted value is checked as well -- from the inputs the call read, and the
/// backend has to give the same bits, as the compiled CPU kernels do. The largest deviation per
/// output goes to the log; a deviation is an error (SEISSOL_EXPR_CHECK).
///
/// Not checked: a program with state, which the evaluation would advance twice; and a call whose
/// outputs cannot be read back (a computed column).
class CheckedKernel final : public Kernel {
  public:
  CheckedKernel(std::unique_ptr<Kernel> inner,
                const Program& program,
                Binding& binding,
                GridStore& grids)
      : inner_(std::move(inner)), program_(program), binding_(&binding), grids_(&grids) {}

  void precompute(const DataTable& table) override { inner_->precompute(table); }

  void run(const DataTable& table) override {
    inner_->run(table);
    if (!checked_) {
      checked_ = true;
      compare(&table, nullptr);
    }
  }

  void run(const KernelArgs& args) override {
    inner_->run(args);
    if (!checked_) {
      checked_ = true;
      compare(nullptr, &args);
    }
  }

  [[nodiscard]] BackendKind kind() const override { return inner_->kind(); }

  private:
  void compare(const DataTable* table, const KernelArgs* args) {
    const auto comparison = compareWithInterpreter(program_, *binding_, *grids_, table, args);
    if (comparison.has_value()) {
      reportComparison(program_, inner_->kind(), *comparison);
    }
  }

  std::unique_ptr<Kernel> inner_;
  Program program_;
  Binding* binding_;
  GridStore* grids_;
  bool checked_{false};
};

/// `kernel`, checked at its first call where SEISSOL_EXPR_CHECK allows.
std::unique_ptr<Kernel> checked(std::unique_ptr<Kernel> kernel,
                                const Program& program,
                                Binding& binding,
                                GridStore& grids,
                                const std::vector<std::size_t>& ids,
                                ComputeType dataType) {
  if (kernel == nullptr || checkPolicy() == CheckPolicy::Off || !binding.hostAddressable()) {
    return kernel;
  }
  static_cast<void>(ids);
  static_cast<void>(dataType);
  return std::make_unique<CheckedKernel>(std::move(kernel), program, binding, grids);
}

} // namespace

namespace {

double ulpOf(double x) {
  const double magnitude = std::fabs(x);
  if (!std::isfinite(magnitude)) {
    return std::numeric_limits<double>::infinity();
  }
  return std::nextafter(magnitude, std::numeric_limits<double>::infinity()) - magnitude;
}

std::string hexOf(double x) {
  std::array<char, 64> buffer{};
  std::snprintf(buffer.data(), buffer.size(), "%a", x);
  return buffer.data();
}

/// The interpreter's outputs at the points of the sample into `outputs`, [output][point], in the
/// type the program moves its data in, through the TileIo `makeIo` makes for them; perturbed and
/// with the signatures of the decisions of the selects as `options` asks.
template <typename T, typename MakeIo>
void interpretAt(const Program& program,
                 const LoweredProgram& lowered,
                 const GridArrays* arrays,
                 std::size_t points,
                 const InterpreterOptions& options,
                 const MakeIo& makeIo,
                 std::vector<double>& outputs) {
  outputs.assign(program.outputs().size() * points, 0.0);
  const auto io = makeIo(outputs);
  TileInterpreter<T> interpreter(program, lowered, arrays, options);
  std::vector<T> persistent(1, T(0));
  interpreter.run(*io, points, persistent.data());
}

/// The reference at the sample: the interpreter's outputs; and, for `perturbed`, R = 8 perturbed
/// evaluations, which give the tolerance of each value and the unstable points.
template <typename T, typename MakeIo>
void interpret(const Program& program,
               const GridArrays* arrays,
               const MakeIo& makeIo,
               bool perturbed,
               BackendComparison& comparison) {
  constexpr int Runs = 8;
  constexpr double Spread = 2.0;
  const LoweredProgram lowered = lower(program, LowerOptions{});
  const std::size_t count = comparison.points.size();
  std::vector<std::uint64_t> decisions(count, 0);
  InterpreterOptions options;
  options.decisions = &decisions;
  interpretAt<T>(program, lowered, arrays, count, options, makeIo, comparison.interpreter);
  if (!perturbed) {
    return;
  }
  comparison.perturbed = true;
  comparison.tolerance.assign(comparison.interpreter.size(), 0.0);
  std::vector<bool> unstable(count, false);
  std::vector<double> values;
  for (int run = 0; run < Runs; ++run) {
    std::vector<std::uint64_t> moved(count, 0);
    InterpreterOptions perturbation;
    perturbation.decisions = &moved;
    perturbation.perturbation.on = true;
    perturbation.perturbation.seed = splitMix64(static_cast<std::uint64_t>(run) + 1);
    interpretAt<T>(program, lowered, arrays, count, perturbation, makeIo, values);
    for (std::size_t v = 0; v < values.size(); ++v) {
      const double deviation = std::fabs(values[v] - comparison.interpreter[v]);
      if (!std::isnan(deviation)) {
        comparison.tolerance[v] = std::max(comparison.tolerance[v], Spread * deviation);
      }
    }
    for (std::size_t k = 0; k < count; ++k) {
      unstable[k] = unstable[k] || moved[k] != decisions[k];
    }
  }
  for (std::size_t k = 0; k < count; ++k) {
    if (unstable[k]) {
      comparison.unstable.push_back(comparison.points[k]);
    }
  }
}

/// Where the backend deviates: in its bits, or beyond the tolerance of the perturbed reference.
void judge(const Program& program, BackendComparison& comparison) {
  const std::size_t sample = comparison.points.size();
  double worstUlps = 0.0;
  for (std::size_t o = 0; o < program.outputs().size(); ++o) {
    for (std::size_t k = 0; k < sample; ++k) {
      const double a = comparison.backend[o * sample + k];
      const double r = comparison.interpreter[o * sample + k];
      if (std::memcmp(&a, &r, sizeof(double)) == 0 || (std::isnan(a) && std::isnan(r))) {
        continue;
      }
      if (comparison.perturbed) {
        // within the tolerance of the value, which the perturbed evaluations set
        const double deviation = std::fabs(a - r);
        const double tolerance = comparison.tolerance[o * sample + k];
        const double ratio =
            tolerance > 0.0 ? deviation / tolerance : std::numeric_limits<double>::infinity();
        comparison.worstRatio = std::max(comparison.worstRatio, ratio);
        if (deviation <= tolerance) {
          continue;
        }
      }
      ++comparison.deviating;
      const double ulps = std::fabs(a - r) / std::max(ulpOf(r), ulpOf(a));
      if (comparison.worst.empty() || ulps > worstUlps || std::isnan(ulps)) {
        worstUlps = ulps;
        comparison.worst = "`" + program.outputs()[o].name + "` at point " +
                           std::to_string(comparison.points[k]) + " is " + hexOf(a) +
                           ", the interpreter's " + hexOf(r);
      }
    }
  }
}

std::string fingerprintOf(const Program& program) {
  std::array<char, 32> buffer{};
  std::snprintf(buffer.data(),
                buffer.size(),
                "%016llx",
                static_cast<unsigned long long>(program.fingerprint()));
  return buffer.data();
}

} // namespace

std::vector<std::size_t> samplePoints(std::size_t first, std::size_t count) {
  return spreadPoints(first, count);
}

bool backendChecks() { return checkPolicy() != CheckPolicy::Off; }

void reportComparison(const Program& program, BackendKind kind, const BackendComparison& found) {
  const std::string name = seissol::expr::name(kind);
  const std::string unstable =
      found.unstable.empty()
          ? std::string()
          : "; " + std::to_string(found.unstable.size()) +
                " of them unstable, where a perturbed evaluation decides a select otherwise";
  if (found.deviating == 0) {
    std::string how = found.perturbed && found.worstRatio > 0.0
                          ? "is within the tolerance of the interpreter, at most " +
                                std::to_string(found.worstRatio) + " of it,"
                          : "gives the bits of the interpreter";
    logInfo() << "expr: the" << name.c_str() << "kernel of program"
              << fingerprintOf(program).c_str() << how.c_str() << "at" << found.points.size()
              << ("points of its first call" + unstable + ".").c_str();
    return;
  }
  const std::string what = "the " + name + " kernel of program " + fingerprintOf(program) +
                           " deviates from the interpreter" +
                           (found.perturbed ? " beyond its tolerance" : "") + " at " +
                           std::to_string(found.deviating) + " values of its first call" +
                           unstable + "; the largest, " + found.worst;
  if (checkPolicy() == CheckPolicy::Error) {
    logError() << "expr:" << what.c_str()
               << ". SEISSOL_EXPR_CHECK=warn lets the run go on with the compiled kernel.";
  } else {
    logWarning() << "expr:" << what.c_str();
  }
}

std::optional<BackendComparison> compareWithInterpreter(const Program& program,
                                                        const Binding& binding,
                                                        GridStore& grids,
                                                        const DataTable* table,
                                                        const KernelArgs* args,
                                                        bool perturbed) {
  const std::size_t first = table != nullptr ? 0 : args->first;
  const std::size_t count = table != nullptr ? binding.numPoints() : args->count;
  if (!program.state().empty() || count == 0) {
    return std::nullopt;
  }
  BackendComparison comparison;
  comparison.points = spreadPoints(first, count);
  const std::size_t outputs = program.outputs().size();
  const std::size_t sample = comparison.points.size();

  // what the backend wrote, and the inputs, for the record
  comparison.backend.assign(outputs * sample, 0.0);
  comparison.inputs.assign(binding.inputs().size() * sample, 0.0);
  std::vector<double> value(std::max(outputs, binding.inputs().size()));
  for (std::size_t k = 0; k < sample; ++k) {
    const bool read =
        table != nullptr
            ? binding.gatherOutputs(*table, comparison.points[k], 1, value.data())
            : binding.gatherOutputsFrom(
                  args->outputs, args->outputCount, comparison.points[k], 1, value.data());
    if (!read) {
      return std::nullopt;
    }
    for (std::size_t o = 0; o < outputs; ++o) {
      comparison.backend[o * sample + k] = value[o];
    }
    if (table != nullptr) {
      binding.gather(*table, comparison.points[k], 1, value.data());
    } else {
      binding.gatherFrom(args->inputs, args->inputCount, comparison.points[k], 1, value.data());
    }
    for (std::size_t i = 0; i < binding.inputs().size(); ++i) {
      comparison.inputs[i * sample + k] = value[i];
    }
  }

  // what the program gives, by the interpreter: the binding reads the inputs, as for the call
  const auto makeIo = [&](auto tag) {
    using T = decltype(tag);
    return [&](std::vector<double>& results) -> std::unique_ptr<TileIo<T>> {
      return std::make_unique<SampleTileIo<T>>(
          binding, table, args, comparison.points, results, outputs);
    };
  };
  const auto makeFloat = makeIo(0.0F);
  const auto makeDouble = makeIo(0.0);
  const auto ids = internGrids(program, grids);
  const StoreArrays arrays(program, ids, grids);
  const GridArrays* sampled = ids.empty() ? nullptr : &arrays;
  if (lower(program, LowerOptions{}).dataType() == ComputeType::F32) {
    interpret<float>(program, sampled, makeFloat, perturbed, comparison);
  } else {
    interpret<double>(program, sampled, makeDouble, perturbed, comparison);
  }
  judge(program, comparison);
  return comparison;
}

std::optional<BackendComparison> compareValues(const Program& program,
                                               GridStore& grids,
                                               const std::vector<std::size_t>& points,
                                               const std::vector<double>& inputs,
                                               const std::vector<double>& backend,
                                               bool perturbed) {
  if (!program.state().empty() || !program.matrices().empty() || !program.blocks().empty() ||
      points.empty()) {
    return std::nullopt;
  }
  BackendComparison comparison;
  comparison.points = points;
  comparison.inputs = inputs;
  comparison.backend = backend;
  const std::size_t inputCount = program.inputs().size();
  const std::size_t outputs = program.outputs().size();
  const auto makeIo = [&](auto tag) {
    using T = decltype(tag);
    return [&](std::vector<double>& results) -> std::unique_ptr<TileIo<T>> {
      return std::make_unique<ValuesTileIo<T>>(
          comparison.inputs, inputCount, points.size(), results, outputs);
    };
  };
  const auto makeFloat = makeIo(0.0F);
  const auto makeDouble = makeIo(0.0);
  const auto ids = internGrids(program, grids);
  const StoreArrays arrays(program, ids, grids);
  const GridArrays* sampled = ids.empty() ? nullptr : &arrays;
  if (lower(program, LowerOptions{}).dataType() == ComputeType::F32) {
    interpret<float>(program, sampled, makeFloat, perturbed, comparison);
  } else {
    interpret<double>(program, sampled, makeDouble, perturbed, comparison);
  }
  judge(program, comparison);
  return comparison;
}

std::vector<std::size_t> internGrids(const Program& program, GridStore& grids) {
  std::vector<std::size_t> ids;
  ids.reserve(program.grids().size());
  for (const auto& desc : program.grids()) {
    ids.push_back(grids.intern(desc));
  }
  if (!ids.empty()) {
    // a no-op for grids that are open and sized already (GridStore::load())
    grids.load();
  }
  return ids;
}

const char* name(BackendKind kind) {
  switch (kind) {
  case BackendKind::Interpreter:
    return "interpreter";
  case BackendKind::RtcCpu:
    return "rtc-cpu";
  case BackendKind::RtcCuda:
    return "rtc-cuda";
  case BackendKind::RtcHip:
    return "rtc-hip";
  case BackendKind::RtcOpenCl:
    return "rtc-opencl";
  case BackendKind::Texture:
    return "texture";
  case BackendKind::Distributed:
    return "distributed";
  }
  return "unknown";
}

std::unique_ptr<Kernel> makeKernel(const Program& program,
                                   Binding& binding,
                                   GridStore& grids,
                                   const BackendOptions& options) {
  // The grids before any kernel: every backend samples them as they lie in the store, the device
  // backends from copies they make of them, and all check them against the program first.
  const auto ids = internGrids(program, grids);
  if (!ids.empty()) {
    checkGrids(program, grids, ids);
  }
  switch (options.preferred) {
  case BackendKind::Interpreter:
    return makeInterpreter(program, binding, grids, ids, options);
  case BackendKind::RtcCpu: {
    LoweredProgram lowered = lower(program, options.lowering);
    const ComputeType dataType = lowered.dataType();
    auto kernel = makeRtcCpuKernel(program, binding, std::move(lowered), options, &grids, ids);
    if (kernel != nullptr) {
      return checked(std::move(kernel), program, binding, grids, ids, dataType);
    }
    if (!options.allowFallback) {
      logError() << "expr: the compiled CPU backend is unusable and fallback was disabled.";
      return nullptr;
    }
    return makeInterpreter(program, binding, grids, ids, options);
  }
  case BackendKind::RtcCuda:
  case BackendKind::RtcHip:
  case BackendKind::RtcOpenCl: {
    LoweredProgram lowered = lower(program, options.lowering);
    auto kernel = makeRtcGpuKernel(program,
                                   binding,
                                   std::move(lowered),
                                   options,
                                   options.preferred == BackendKind::RtcCuda  ? GpuTarget::Cuda
                                   : options.preferred == BackendKind::RtcHip ? GpuTarget::Hip
                                                                              : GpuTarget::OpenCl,
                                   grids,
                                   ids);
    if (kernel != nullptr) {
      return kernel;
    }
    if (!options.allowFallback) {
      logError() << "expr: backend" << name(options.preferred)
                 << "is unusable and fallback was disabled.";
      return nullptr;
    }
    // The device path declines for reasons the model author can act on -- a
    // grid without an array, a computed column, a host pointer -- and
    // makeRtcGpuKernel has already named which. Trying the compiled CPU kernel
    // before the interpreter, because most of those reasons do not apply to it.
    LoweredProgram cpuLowering = lower(program, options.lowering);
    const ComputeType dataType = cpuLowering.dataType();
    auto cpu = makeRtcCpuKernel(
        program, binding, std::move(cpuLowering), forFallback(options), &grids, ids);
    if (cpu != nullptr) {
      return checked(std::move(cpu), program, binding, grids, ids, dataType);
    }
    return makeInterpreter(program, binding, grids, ids, forFallback(options));
  }
  case BackendKind::Texture:
    return makeTextureKernel(program, binding, grids);
  case BackendKind::Distributed:
    return makeDistributedKernel(program, binding, grids);
  }
  logError() << "expr: unknown backend requested.";
  return nullptr;
}

std::unique_ptr<Kernel>
    makeTextureKernel(const Program& /*program*/, Binding& /*binding*/, GridStore& /*grids*/) {
  // See the deferral note in Backend.h: fp32-only with interpolation weights
  // quantised to 8-9 subpixel bits, which breaks CPU/GPU reproducibility.
  logError() << "expr: the texture backend is not implemented.";
  return nullptr;
}

std::unique_ptr<Kernel>
    makeDistributedKernel(const Program& /*program*/, Binding& /*binding*/, GridStore& /*grids*/) {
  // See the deferral note in Backend.h: needs the sub-box/halo analysis on
  // Lookup coordinate arguments before it can be correct.
  logError() << "expr: the distributed backend is not implemented.";
  return nullptr;
}

} // namespace seissol::expr
