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
#include <cstddef>
#include <cstdint>
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

} // namespace

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
    auto kernel = makeRtcCpuKernel(program, binding, std::move(lowered), options, &grids, ids);
    if (kernel != nullptr) {
      return kernel;
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
    auto cpu = makeRtcCpuKernel(
        program, binding, std::move(cpuLowering), forFallback(options), &grids, ids);
    if (cpu != nullptr) {
      return cpu;
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
