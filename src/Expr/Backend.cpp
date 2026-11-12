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
#include "Expr/Record.h"
#include "Expr/Rewrite.h"
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
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <set>
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
  std::unique_ptr<Kernel> kernel;
  if (lowered.dataType() == ComputeType::F32) {
    kernel = std::make_unique<InterpreterKernel<float>>(
        program, binding, std::move(lowered), std::move(arrays), interpreterOptions);
  } else {
    kernel = std::make_unique<InterpreterKernel<double>>(
        program, binding, std::move(lowered), std::move(arrays), interpreterOptions);
  }
  return kernel;
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

bool recording();

/// An interpreter kernel with a sample of its first call taken, for the evaluation record.
std::unique_ptr<Kernel> recorded(std::unique_ptr<Kernel> kernel,
                                 const Program& program,
                                 Binding& binding,
                                 GridStore& grids) {
  if (kernel == nullptr || !recording() || !binding.hostAddressable()) {
    return kernel;
  }
  return std::make_unique<CheckedKernel>(std::move(kernel), program, binding, grids);
}

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

/// The size of the pool of a call, from which the sample of the record is chosen and on which
/// the backend is checked: on the host, and on a device, from which every value is copied apart.
constexpr std::size_t PoolSize = 4096;
constexpr std::size_t DevicePoolSize = 1024;
/// The points of the smallest hash in the sample.
constexpr std::size_t SmallestHashes = 64;

/// The closed pool of a call of `count` points from `first`: its `size` points of the smallest
/// hashes of their identity (all of them, for fewer), in the order of the call, with their
/// identities.
CallPool poolOf(std::size_t first,
                std::size_t count,
                std::size_t size,
                const std::function<std::uint64_t(std::size_t)>& identity) {
  // the `size` smallest hashes, the largest of them on top
  std::priority_queue<std::pair<std::uint64_t, std::size_t>> smallest;
  for (std::size_t point = first; point < first + count; ++point) {
    const std::uint64_t hash = splitMix64(identity(point));
    if (smallest.size() < size) {
      smallest.emplace(hash, point);
    } else if (hash < smallest.top().first) {
      smallest.pop();
      smallest.emplace(hash, point);
    }
  }
  CallPool pool;
  pool.scope.callPoints = count;
  pool.scope.closedAll = count <= size;
  pool.scope.bound = smallest.empty() ? 0 : smallest.top().first;
  while (!smallest.empty()) {
    pool.points.push_back(smallest.top().second);
    smallest.pop();
  }
  std::sort(pool.points.begin(), pool.points.end());
  for (const std::size_t point : pool.points) {
    pool.identities.push_back(identity(point));
  }
  pool.scope.closed = pool.points.size();
  return pool;
}

/// Adds to `pool` the points that the rules over all points of the call choose (the specification's
/// "Evaluation record"): per group the two points of the smallest hash of their identity; per
/// output the points of its smallest and of its largest value, by what the backend wrote, and of
/// two equal values the one of the smaller hash. `readOutputs(first, count, dst)` reads what the
/// backend wrote ([output][point]), false where it cannot; `readInputs(first, count, dst)` the
/// inputs ([input][point], in the order of Program::inputs()). Nothing is added, and the scope
/// stays the pool, where the outputs cannot be read.
template <typename ReadOutputs, typename ReadInputs>
void scanCall(const Program& program,
              std::size_t first,
              std::size_t count,
              const std::function<std::uint64_t(std::size_t)>& identity,
              const ReadOutputs& readOutputs,
              const ReadInputs& readInputs,
              CallPool& pool) {
  constexpr std::size_t Chunk = 1024;
  constexpr std::size_t None = std::numeric_limits<std::size_t>::max();
  struct Candidate {
    std::uint64_t hash{0};
    std::size_t point{None};
    double value{0.0};
  };
  const std::size_t outputs = program.outputs().size();
  const std::size_t inputs = program.inputs().size();
  const auto group = std::find_if(program.inputs().begin(),
                                  program.inputs().end(),
                                  [](const VarSpec& input) { return input.name == "group"; });
  const bool grouped = group != program.inputs().end();
  const auto groupColumn = static_cast<std::size_t>(group - program.inputs().begin());

  std::vector<Candidate> smallest(outputs);
  std::vector<Candidate> largest(outputs);
  // per group, the two of the smallest hash, the smaller first
  std::map<double, std::array<Candidate, 2>> groups;
  std::vector<double> values(outputs * Chunk);
  std::vector<double> read(grouped ? inputs * Chunk : 0);
  for (std::size_t from = first; from < first + count; from += Chunk) {
    const std::size_t n = std::min(Chunk, first + count - from);
    if (!readOutputs(from, n, values.data())) {
      return;
    }
    if (grouped) {
      readInputs(from, n, read.data());
    }
    for (std::size_t lane = 0; lane < n; ++lane) {
      const std::size_t point = from + lane;
      const std::uint64_t hash = splitMix64(identity(point));
      for (std::size_t o = 0; o < outputs; ++o) {
        const double value = values[o * n + lane];
        if (std::isnan(value)) {
          continue;
        }
        auto& low = smallest[o];
        if (low.point == None || value < low.value || (value == low.value && hash < low.hash)) {
          low = Candidate{hash, point, value};
        }
        auto& high = largest[o];
        if (high.point == None || value > high.value || (value == high.value && hash < high.hash)) {
          high = Candidate{hash, point, value};
        }
      }
      if (grouped) {
        const double key = read[groupColumn * n + lane];
        if (std::isnan(key)) {
          continue;
        }
        auto& two = groups[key];
        const Candidate candidate{hash, point, key};
        if (two[0].point == None || hash < two[0].hash) {
          two[1] = two[0];
          two[0] = candidate;
        } else if (two[1].point == None || hash < two[1].hash) {
          two[1] = candidate;
        }
      }
    }
  }

  std::set<std::size_t> chosen(pool.points.begin(), pool.points.end());
  const auto add = [&](const Candidate& candidate) {
    if (candidate.point != None) {
      chosen.insert(candidate.point);
    }
  };
  for (std::size_t o = 0; o < outputs; ++o) {
    add(smallest[o]);
    add(largest[o]);
  }
  for (const auto& [key, two] : groups) {
    static_cast<void>(key);
    add(two[0]);
    add(two[1]);
  }
  pool.points.assign(chosen.begin(), chosen.end());
  pool.identities.clear();
  for (const std::size_t point : pool.points) {
    pool.identities.push_back(identity(point));
  }
  pool.scope.scanned = true;
}

/// The interpreter's outputs at the points of the pool into `outputs`, [output][point], in the
/// type the program moves its data in, through the TileIo `makeIo` makes for them; perturbed, and
/// with the decisions of the selects, as `options` asks.
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

/// What the interpreter gives at the pool: its outputs, the branches of its selects, and, for
/// `perturbed`, from R = 8 perturbed evaluations, the tolerance of each value and the unstable
/// points.
struct PoolReference {
  std::vector<double> outputs;
  std::map<std::size_t, std::vector<std::uint8_t>> branches;
  std::vector<double> tolerance;
  std::vector<bool> unstable;
};

template <typename T, typename MakeIo>
PoolReference interpret(const Program& program,
                        const GridArrays* arrays,
                        std::size_t points,
                        const MakeIo& makeIo,
                        bool perturbed) {
  constexpr int Runs = 8;
  constexpr double Spread = 2.0;
  const LoweredProgram lowered = lower(program, LowerOptions{});
  PoolReference reference;
  std::vector<std::uint64_t> decisions(points, 0);
  InterpreterOptions options;
  options.decisions = &decisions;
  options.branches = &reference.branches;
  interpretAt<T>(program, lowered, arrays, points, options, makeIo, reference.outputs);
  reference.unstable.assign(points, false);
  if (!perturbed) {
    return reference;
  }
  reference.tolerance.assign(reference.outputs.size(), 0.0);
  std::vector<double> values;
  for (int run = 0; run < Runs; ++run) {
    std::vector<std::uint64_t> moved(points, 0);
    InterpreterOptions perturbation;
    perturbation.decisions = &moved;
    perturbation.perturbation.on = true;
    perturbation.perturbation.seed = splitMix64(static_cast<std::uint64_t>(run) + 1);
    interpretAt<T>(program, lowered, arrays, points, perturbation, makeIo, values);
    for (std::size_t v = 0; v < values.size(); ++v) {
      const double deviation = std::fabs(values[v] - reference.outputs[v]);
      if (!std::isnan(deviation)) {
        reference.tolerance[v] = std::max(reference.tolerance[v], Spread * deviation);
      }
    }
    for (std::size_t k = 0; k < points; ++k) {
      reference.unstable[k] = reference.unstable[k] || moved[k] != decisions[k];
    }
  }
  return reference;
}

/// The comparison of the backend with the reference at the pool `callPool` (with the inputs and
/// the backend's outputs there), judged at every point of the pool, and the sample of the record
/// chosen from it (the specification's "Evaluation record"): the 64 points of the smallest hash;
/// in every group the two of the smallest; per output the points of its smallest and its largest
/// value, by what the backend wrote; per select the point of the smallest hash in the closed pool
/// where it takes its first branch and the one where it does not; per output the point where the
/// backend deviates most; and up to 64 of the unstable points.
template <typename MakeFloat, typename MakeDouble>
BackendComparison compareAt(const Program& program,
                            GridStore& grids,
                            const CallPool& callPool,
                            const std::vector<double>& inputs,
                            const std::vector<double>& backend,
                            const MakeFloat& makeFloat,
                            const MakeDouble& makeDouble,
                            bool perturbed) {
  const auto& points = callPool.points;
  const auto& identities = callPool.identities;
  const std::size_t pool = points.size();
  const std::size_t outputs = program.outputs().size();
  const auto ids = internGrids(program, grids);
  const StoreArrays arrays(program, ids, grids);
  const GridArrays* sampled = ids.empty() ? nullptr : &arrays;
  const PoolReference reference =
      lower(program, LowerOptions{}).dataType() == ComputeType::F32
          ? interpret<float>(program, sampled, pool, makeFloat, perturbed)
          : interpret<double>(program, sampled, pool, makeDouble, perturbed);

  BackendComparison comparison;
  comparison.perturbed = perturbed;
  comparison.compared = pool;
  comparison.scope = callPool.scope;
  std::vector<std::uint64_t> hashes(pool);
  for (std::size_t k = 0; k < pool; ++k) {
    hashes[k] = splitMix64(identities[k]);
  }
  // the pool, in the order of the hashes of the identities
  std::vector<std::size_t> byHash(pool);
  for (std::size_t k = 0; k < pool; ++k) {
    byHash[k] = k;
  }
  std::sort(byHash.begin(), byHash.end(), [&](std::size_t a, std::size_t b) {
    return hashes[a] < hashes[b];
  });

  // judged at every point of the pool
  std::vector<bool> chosen(pool, false);
  std::vector<double> worstOfOutput(outputs, 0.0);
  std::vector<std::size_t> worstPoint(outputs, pool);
  double worstUlps = 0.0;
  for (std::size_t o = 0; o < outputs; ++o) {
    for (std::size_t k = 0; k < pool; ++k) {
      const double a = backend[o * pool + k];
      const double r = reference.outputs[o * pool + k];
      if (std::memcmp(&a, &r, sizeof(double)) == 0 || (std::isnan(a) && std::isnan(r))) {
        continue;
      }
      const double deviation = std::fabs(a - r);
      if (!(deviation <= worstOfOutput[o])) {
        worstOfOutput[o] = deviation;
        worstPoint[o] = k;
      }
      if (perturbed) {
        const double tolerance = reference.tolerance[o * pool + k];
        const double ratio =
            tolerance > 0.0 ? deviation / tolerance : std::numeric_limits<double>::infinity();
        comparison.worstRatio = std::max(comparison.worstRatio, ratio);
        if (deviation <= tolerance) {
          continue;
        }
      }
      ++comparison.deviating;
      const double ulps = deviation / std::max(ulpOf(r), ulpOf(a));
      if (comparison.worst.empty() || ulps > worstUlps || std::isnan(ulps)) {
        worstUlps = ulps;
        comparison.worst = "`" + program.outputs()[o].name + "` at point " +
                           std::to_string(points[k]) + " is " + hexOf(a) + ", the interpreter's " +
                           hexOf(r);
      }
    }
  }

  // the sample
  for (std::size_t i = 0; i < std::min(SmallestHashes, pool); ++i) {
    chosen[byHash[i]] = true;
  }
  const auto group = std::find_if(program.inputs().begin(),
                                  program.inputs().end(),
                                  [](const VarSpec& input) { return input.name == "group"; });
  if (group != program.inputs().end()) {
    const auto column = static_cast<std::size_t>(group - program.inputs().begin());
    std::map<double, int> taken;
    for (const std::size_t k : byHash) {
      if (taken[inputs[column * pool + k]]++ < 2) {
        chosen[k] = true;
      }
    }
  }
  // by what the backend wrote, as scanCall() chooses over all points
  for (std::size_t o = 0; o < outputs; ++o) {
    std::size_t smallest = pool;
    std::size_t largest = pool;
    for (const std::size_t k : byHash) {
      const double value = backend[o * pool + k];
      if (std::isnan(value)) {
        continue;
      }
      if (smallest == pool || value < backend[o * pool + smallest]) {
        smallest = k;
      }
      if (largest == pool || value > backend[o * pool + largest]) {
        largest = k;
      }
    }
    for (const std::size_t k : {smallest, largest, worstPoint[o]}) {
      if (k < pool) {
        chosen[k] = true;
      }
    }
  }
  // from the closed pool: the points of the smallest hashes, which come first
  const auto closed = [&](std::size_t k) {
    return callPool.scope.closedAll || hashes[k] <= callPool.scope.bound;
  };
  for (const auto& [select, taken] : reference.branches) {
    static_cast<void>(select);
    for (const std::uint8_t branch : {std::uint8_t{1}, std::uint8_t{0}}) {
      for (const std::size_t k : byHash) {
        if (!closed(k)) {
          break;
        }
        if (k < taken.size() && taken[k] == branch) {
          chosen[k] = true;
          break;
        }
      }
    }
  }
  // per grid, the points of the closed pool beyond its volume, and the two of the smallest hash
  const auto clamps = clampsAt(program, grids, inputs, pool);
  for (std::size_t g = 0; g < clamps.size(); ++g) {
    if (clamps[g].axes.empty()) {
      continue;
    }
    const auto& desc = program.grids()[g];
    ClampCount count;
    count.grid = desc.name.empty() ? desc.path : desc.name;
    count.axes.assign(clamps[g].axes.size(), {0, 0, 0});
    int taken = 0;
    for (const std::size_t k : byHash) {
      if (!closed(k)) {
        break;
      }
      ++count.of;
      bool beyond = false;
      for (std::size_t d = 0; d < clamps[g].axes.size(); ++d) {
        const std::uint8_t sides = clamps[g].axes[d][k];
        count.axes[d][0] += (sides & GridClamps::Below) != 0 ? 1 : 0;
        count.axes[d][1] += (sides & GridClamps::Above) != 0 ? 1 : 0;
        count.axes[d][2] += (sides & GridClamps::Undefined) != 0 ? 1 : 0;
        beyond = beyond || sides != 0;
      }
      if (beyond && taken < 2) {
        chosen[k] = true;
        ++taken;
      }
    }
    comparison.clamps.push_back(std::move(count));
  }
  std::size_t unstable = 0;
  for (const std::size_t k : byHash) {
    if (reference.unstable[k]) {
      comparison.unstable.push_back(points[k]);
      comparison.unstableIdentities.push_back(identities[k]);
      if (unstable++ < SmallestHashes) {
        chosen[k] = true;
      }
    }
  }

  // and their values, in the order of the call
  std::vector<std::size_t> sample;
  for (std::size_t k = 0; k < pool; ++k) {
    if (chosen[k]) {
      sample.push_back(k);
    }
  }
  const std::size_t inputCount = pool == 0 ? 0 : inputs.size() / pool;
  const std::size_t n = sample.size();
  comparison.points.resize(n);
  comparison.identities.resize(n);
  comparison.inputs.assign(inputCount * n, 0.0);
  comparison.interpreter.assign(outputs * n, 0.0);
  comparison.backend.assign(outputs * n, 0.0);
  if (perturbed) {
    comparison.tolerance.assign(outputs * n, 0.0);
  }
  for (std::size_t j = 0; j < n; ++j) {
    const std::size_t k = sample[j];
    comparison.points[j] = points[k];
    comparison.identities[j] = identities[k];
    for (std::size_t i = 0; i < inputCount; ++i) {
      comparison.inputs[i * n + j] = inputs[i * pool + k];
    }
    for (std::size_t o = 0; o < outputs; ++o) {
      comparison.interpreter[o * n + j] = reference.outputs[o * pool + k];
      comparison.backend[o * n + j] = backend[o * pool + k];
      if (perturbed) {
        comparison.tolerance[o * n + j] = reference.tolerance[o * pool + k];
      }
    }
  }
  return comparison;
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

namespace {

struct Records {
  std::mutex mutex;
  std::string directory;
  std::vector<std::string> build;
  int rank{0};
  std::set<std::string> written;
};

Records& records() {
  static Records Instance;
  return Instance;
}

bool recording() {
  auto& state = records();
  const std::lock_guard<std::mutex> lock(state.mutex);
  return !state.directory.empty();
}

/// The evaluation record of `program` on the `kind` backend, once per program and backend.
void writeRecord(const Program& program, BackendKind kind, const BackendComparison& found) {
  auto& state = records();
  const std::lock_guard<std::mutex> lock(state.mutex);
  if (state.directory.empty()) {
    return;
  }
  const std::string fingerprint = fingerprintOf(program);
  const std::string backend = seissol::expr::name(kind);
  if (!state.written.insert(fingerprint + "-" + backend).second) {
    return;
  }
  const std::string text = printRecord(recordOf(program, backend, found, state.build, state.rank));
  const std::filesystem::path path =
      std::filesystem::path(state.directory) /
      (fingerprint + "-" + backend + "-r" + std::to_string(state.rank) + ".txt");
  std::ofstream file(path);
  file << text;
  if (!file) {
    logWarning() << "expr: could not write the evaluation record" << path.string();
  }
}

} // namespace

void setRecordDirectory(const std::string& directory,
                        const std::vector<std::string>& build,
                        int rank) {
  auto& state = records();
  const std::lock_guard<std::mutex> lock(state.mutex);
  state.directory = directory;
  state.build = build;
  state.rank = rank;
  if (!directory.empty()) {
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error) {
      logWarning() << "expr: could not create the directory of the evaluation records" << directory
                   << ":" << error.message();
      state.directory.clear();
    }
  }
}

void reportComparison(const Program& program, BackendKind kind, const BackendComparison& found) {
  writeRecord(program, kind, found);
  if (kind == BackendKind::Interpreter) {
    // compared with itself: only for its record
    return;
  }
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
  // the identities of the points: the call's, where it gives them, else the table's
  const bool given = args != nullptr && args->identity != nullptr;
  const std::function<std::uint64_t(std::size_t)> identity =
      given ? std::function<std::uint64_t(std::size_t)>([args](std::size_t point) {
        return args->identity(args->identityContext, point);
      })
            : std::function<std::uint64_t(std::size_t)>(
                  [&binding](std::size_t point) { return binding.identity(point); });
  CallPool callPool = poolOf(first, count, PoolSize, identity);
  callPool.scope.globalIdentities = given || binding.identified();
  // the groups and the extreme values over all points of the call
  scanCall(
      program,
      first,
      count,
      identity,
      [&](std::size_t from, std::size_t n, double* dst) {
        return table != nullptr
                   ? binding.gatherOutputs(*table, from, n, dst)
                   : binding.gatherOutputsFrom(args->outputs, args->outputCount, from, n, dst);
      },
      [&](std::size_t from, std::size_t n, double* dst) {
        if (table != nullptr) {
          binding.gather(*table, from, n, dst);
        } else {
          binding.gatherFrom(args->inputs, args->inputCount, from, n, dst);
        }
      },
      callPool);
  const auto& points = callPool.points;
  const std::size_t pool = points.size();
  const std::size_t outputs = program.outputs().size();

  // what the backend wrote at the pool, and the inputs there
  std::vector<double> backend(outputs * pool, 0.0);
  std::vector<double> inputs(binding.inputs().size() * pool, 0.0);
  std::vector<double> value(std::max(outputs, binding.inputs().size()));
  for (std::size_t k = 0; k < pool; ++k) {
    const bool read = table != nullptr
                          ? binding.gatherOutputs(*table, points[k], 1, value.data())
                          : binding.gatherOutputsFrom(
                                args->outputs, args->outputCount, points[k], 1, value.data());
    if (!read) {
      return std::nullopt;
    }
    for (std::size_t o = 0; o < outputs; ++o) {
      backend[o * pool + k] = value[o];
    }
    if (table != nullptr) {
      binding.gather(*table, points[k], 1, value.data());
    } else {
      binding.gatherFrom(args->inputs, args->inputCount, points[k], 1, value.data());
    }
    for (std::size_t i = 0; i < binding.inputs().size(); ++i) {
      inputs[i * pool + k] = value[i];
    }
  }

  // the interpreter reads the inputs as the call did, through the binding
  const auto makeIo = [&](auto tag) {
    using T = decltype(tag);
    return [&](std::vector<double>& results) -> std::unique_ptr<TileIo<T>> {
      return std::make_unique<SampleTileIo<T>>(binding, table, args, points, results, outputs);
    };
  };
  const auto makeFloat = makeIo(0.0F);
  const auto makeDouble = makeIo(0.0);
  return compareAt(program, grids, callPool, inputs, backend, makeFloat, makeDouble, perturbed);
}

std::optional<BackendComparison> compareValues(const Program& program,
                                               GridStore& grids,
                                               const CallPool& pool,
                                               const std::vector<double>& inputs,
                                               const std::vector<double>& backend,
                                               bool perturbed) {
  const auto& points = pool.points;
  if (!program.state().empty() || !program.matrices().empty() || !program.blocks().empty() ||
      points.empty()) {
    return std::nullopt;
  }
  const std::size_t inputCount = program.inputs().size();
  const std::size_t outputs = program.outputs().size();
  const auto makeIo = [&](auto tag) {
    using T = decltype(tag);
    return [&](std::vector<double>& results) -> std::unique_ptr<TileIo<T>> {
      return std::make_unique<ValuesTileIo<T>>(inputs, inputCount, points.size(), results, outputs);
    };
  };
  const auto makeFloat = makeIo(0.0F);
  const auto makeDouble = makeIo(0.0);
  return compareAt(program, grids, pool, inputs, backend, makeFloat, makeDouble, perturbed);
}

std::vector<double> interpretValues(const Program& program,
                                    GridStore& grids,
                                    const std::vector<double>& inputs,
                                    std::size_t points,
                                    std::map<std::size_t, std::vector<std::uint8_t>>* branches) {
  const std::size_t outputs = program.outputs().size();
  const auto ids = internGrids(program, grids);
  const StoreArrays arrays(program, ids, grids);
  const GridArrays* sampled = ids.empty() ? nullptr : &arrays;
  const LoweredProgram lowered = lower(program, LowerOptions{});
  InterpreterOptions options;
  options.branches = branches;
  std::vector<double> results;
  const auto run = [&](auto tag) {
    using T = decltype(tag);
    const auto makeIo = [&](std::vector<double>& values) {
      return std::make_unique<ValuesTileIo<T>>(
          inputs, program.inputs().size(), points, values, outputs);
    };
    interpretAt<T>(program, lowered, sampled, points, options, makeIo, results);
  };
  if (lowered.dataType() == ComputeType::F32) {
    run(0.0F);
  } else {
    run(0.0);
  }
  return results;
}

void scanValues(const Program& program,
                std::size_t first,
                std::size_t count,
                const std::function<std::uint64_t(std::size_t)>& identity,
                const std::vector<std::vector<double>>& outputs,
                const std::vector<double>& groups,
                CallPool& pool) {
  const auto group = std::find_if(program.inputs().begin(),
                                  program.inputs().end(),
                                  [](const VarSpec& input) { return input.name == "group"; });
  const auto groupColumn = static_cast<std::size_t>(group - program.inputs().begin());
  if (outputs.size() != program.outputs().size() ||
      (group != program.inputs().end() && groups.size() < count)) {
    return;
  }
  scanCall(
      program,
      first,
      count,
      identity,
      [&](std::size_t from, std::size_t n, double* dst) {
        for (std::size_t o = 0; o < outputs.size(); ++o) {
          if (outputs[o].size() < count) {
            return false;
          }
          std::copy_n(
              outputs[o].begin() + static_cast<std::ptrdiff_t>(from - first), n, dst + o * n);
        }
        return true;
      },
      [&](std::size_t from, std::size_t n, double* dst) {
        if (group != program.inputs().end()) {
          std::copy_n(
              groups.begin() + static_cast<std::ptrdiff_t>(from - first), n, dst + groupColumn * n);
        }
      },
      pool);
}

std::vector<GridClamps> clampsAt(const Program& program,
                                 GridStore& grids,
                                 const std::vector<double>& inputs,
                                 std::size_t points) {
  // the lookups the program evaluates, and whether a coordinate reads a contraction
  const Arena& arena = program.arena();
  std::vector<NodeId> roots = program.roots();
  for (const auto& state : program.state()) {
    roots.push_back(state.root);
  }
  for (const auto& check : program.checks()) {
    roots.push_back(check.root);
  }
  std::vector<char> seen(arena.size(), 0);
  std::vector<NodeId> lookups;
  bool contracts = false;
  std::vector<NodeId> stack(roots.begin(), roots.end());
  std::vector<NodeId> kids;
  while (!stack.empty()) {
    const NodeId id = stack.back();
    stack.pop_back();
    if (id == NoNode || seen[static_cast<std::size_t>(id)] != 0) {
      continue;
    }
    seen[static_cast<std::size_t>(id)] = 1;
    const Node& node = arena[id];
    if (node.kind == Kind::Lookup) {
      lookups.push_back(id);
    }
    contracts = contracts || node.kind == Kind::Contract;
    arena.children(id, kids);
    stack.insert(stack.end(), kids.begin(), kids.end());
  }
  if (lookups.empty() || contracts || points == 0) {
    return {};
  }
  std::sort(lookups.begin(), lookups.end());

  // a program of the coordinates of every lookup, on the inputs of the program
  Program probe;
  probe.setComputeType(program.computeType());
  for (const auto& input : program.inputs()) {
    probe.addInput(input.name, input.type);
  }
  std::vector<NodeId> coordinates;
  for (const NodeId lookup : lookups) {
    const Node& node = arena[lookup];
    const NodeId* args = arena.args(node);
    coordinates.insert(coordinates.end(), args, args + node.argCount);
  }
  const auto imported = importNodes(probe, program, coordinates);
  for (std::size_t c = 0; c < imported.size(); ++c) {
    probe.addOutput("c" + std::to_string(c), reader::scripting::DataType::F64, imported[c]);
  }
  const auto values = interpretValues(probe, grids, inputs, points);

  // their sides, by the geometry of the arrays the kernels sample, as Expand.cpp locates them
  const auto ids = internGrids(program, grids);
  std::vector<GridClamps> clamps(program.grids().size());
  std::size_t c = 0;
  for (const NodeId lookup : lookups) {
    const Node& node = arena[lookup];
    const auto grid = static_cast<std::size_t>(node.grid);
    const reader::datafield::ArrayView* view = grids.get(ids.at(grid)).arrayView();
    auto& axes = clamps[grid].axes;
    axes.resize(std::max<std::size_t>(axes.size(), node.argCount));
    for (std::int32_t d = 0; d < node.argCount; ++d, ++c) {
      auto& sides = axes[static_cast<std::size_t>(d)];
      sides.resize(points, 0);
      if (view == nullptr) {
        continue;
      }
      const auto axis = static_cast<std::size_t>(d);
      const double top = static_cast<double>(view->num[axis]) - 1.0;
      for (std::size_t k = 0; k < points; ++k) {
        const double raw = (values[c * points + k] - view->min[axis]) * view->deltaInv[axis];
        if (std::isnan(raw)) {
          sides[k] |= GridClamps::Undefined;
        } else if (raw < 0.0 && view->before[axis] == 0) {
          sides[k] |= GridClamps::Below;
        } else if (raw > top && view->after[axis] == 0) {
          sides[k] |= GridClamps::Above;
        }
      }
    }
  }
  return clamps;
}

CallPool devicePool(std::size_t first,
                    std::size_t count,
                    const std::function<std::uint64_t(std::size_t)>& identity,
                    bool globalIdentities) {
  CallPool pool = poolOf(first, count, DevicePoolSize, identity);
  pool.scope.globalIdentities = globalIdentities;
  return pool;
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
    return recorded(
        makeInterpreter(program, binding, grids, ids, options), program, binding, grids);
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
    return recorded(
        makeInterpreter(program, binding, grids, ids, options), program, binding, grids);
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
    return recorded(makeInterpreter(program, binding, grids, ids, forFallback(options)),
                    program,
                    binding,
                    grids);
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
