// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Expr/Backend.h"

#include "Expr/Binding.h"
#include "Expr/Interp.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/RtcCpu.h"
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

// Interp.h's GridSampler over a GridStore. Since Package 4 changed
// Grid::sampleBatch to pointer arrays with a null-skip, this is genuinely an
// adapter and not a staging buffer: the coordinate pointers go straight through
// and the one requested component points straight at the destination slot.
//
// The f32 path is the exception, and deliberately so. Grid.h fixes coordinates
// at f64 because index arithmetic in f32 loses exactness long before the value
// does, so a float tile's coordinates are converted here. That is a type
// conversion the design asks for, not a copy the signature forces.
class StoreSampler final : public GridSampler {
  public:
  StoreSampler(std::vector<std::size_t> storeIds, GridStore& store, std::size_t maxComponents)
      : store_(&store), storeIds_(std::move(storeIds)), outF64_(maxComponents, nullptr),
        outF32_(maxComponents, nullptr) {}

  void sampleBatch(GridId grid,
                   std::int32_t component,
                   const double* const* coords,
                   std::int32_t dimension,
                   std::size_t count,
                   double* dst) const override {
    const Grid& target = resolve(grid, component, dimension);
    fill(outF64_, target.components(), component, dst);
    target.sampleBatch(coords, count, outF64_.data());
  }

  void sampleBatch(GridId grid,
                   std::int32_t component,
                   const float* const* coords,
                   std::int32_t dimension,
                   std::size_t count,
                   float* dst) const override {
    const Grid& target = resolve(grid, component, dimension);

    const auto dim = static_cast<std::size_t>(dimension);
    if (coordScratch_.size() < dim * count) {
      coordScratch_.resize(dim * count);
    }
    for (std::size_t d = 0; d < dim; ++d) {
      double* const row = coordScratch_.data() + d * count;
      for (std::size_t lane = 0; lane < count; ++lane) {
        row[lane] = static_cast<double>(coords[d][lane]);
      }
      coordRows_[d] = row;
    }

    fill(outF32_, target.components(), component, dst);
    target.sampleBatch(coordRows_.data(), count, outF32_.data());
  }

  private:
  const Grid& resolve(GridId grid, std::int32_t component, std::int32_t dimension) const {
    const auto id = static_cast<std::size_t>(grid);
    if (grid == NoGrid || id >= storeIds_.size()) {
      logError() << "expr: a lookup references grid" << static_cast<int>(grid)
                 << "which the program does not declare.";
    }
    const Grid& target = store_->get(storeIds_[id]);
    // Checked per call rather than once, because it costs two compares next to
    // a width^d gather and it is the difference between a diagnostic and an
    // out-of-bounds read of the grid array.
    if (static_cast<std::size_t>(dimension) != target.dimensions()) {
      logError() << "expr: a lookup on grid" << static_cast<int>(grid) << "passes" << dimension
                 << "coordinates, but the grid has" << target.dimensions() << "dimensions.";
    }
    if (component < 0 || static_cast<std::size_t>(component) >= target.components()) {
      logError() << "expr: a lookup on grid" << static_cast<int>(grid) << "requests component"
                 << component << "of" << target.components() << ".";
    }
    return target;
  }

  // Every entry null but the one the Lookup names. The nulls are what stop the
  // reduction from running for components nobody reads; see Interpolation.h.
  template <typename T>
  static void fill(std::vector<T*>& out, std::size_t components, std::int32_t component, T* dst) {
    std::fill(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(components), nullptr);
    out[static_cast<std::size_t>(component)] = dst;
  }

  GridStore* store_;
  std::vector<std::size_t> storeIds_; // GridId -> GridStore index
  // Mutable because sampleBatch is const on the GridSampler interface and the
  // interpreter is single-threaded per Kernel (see the note on Kernel::run).
  mutable std::vector<double*> outF64_;
  mutable std::vector<float*> outF32_;
  mutable std::vector<double> coordScratch_;
  mutable std::array<const double*, MaxLookupDimension> coordRows_{};
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
                    std::unique_ptr<StoreSampler> sampler,
                    const InterpreterOptions& options)
      : binding_(&binding), lowered_(std::move(lowered)), sampler_(std::move(sampler)),
        interpreter_(program, lowered_, sampler_.get(), options),
        partitions_(toPointRanges(binding)), needsPrecompute_(lowered_.hasPrecompute()) {}

  void precompute(const DataTable& table) override {
    if (!needsPrecompute_) {
      return;
    }
    const BoundTileIo<T> io(*binding_, table);
    interpreter_.precompute(io, binding_->numPoints(), persistent(), partitions_);
    precomputed_ = true;
  }

  void run(const KernelArgs& args) override {
    if (!binding_->addressable()) {
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
  std::unique_ptr<StoreSampler> sampler_;
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

// --- grid interning ---------------------------------------------------------

std::vector<std::size_t> internGrids(const Program& program, GridStore& store) {
  std::vector<std::size_t> ids;
  ids.reserve(program.grids().size());
  for (const auto& desc : program.grids()) {
    ids.push_back(store.intern(desc));
  }
  if (!ids.empty()) {
    // Idempotent for opening, but the window sizing re-runs, so a caller
    // building several kernels over one store should intern every program
    // before the first makeKernel rather than interleaving the two.
    store.load();
  }
  return ids;
}

std::size_t maxComponents(const std::vector<std::size_t>& ids, GridStore& store) {
  std::size_t most = 1;
  for (const std::size_t id : ids) {
    most = std::max(most, store.get(id).components());
  }
  return most;
}

std::unique_ptr<Kernel> makeInterpreter(const Program& program,
                                        Binding& binding,
                                        GridStore& grids,
                                        const BackendOptions& options) {
  LoweredProgram lowered = lower(program, options.lowering);
  logInfo() << "expr: interpreter kernel --" << lowered.summary().c_str();

  binding.allocatePersistent(program, lowered.persistentSlotCount());

  auto ids = internGrids(program, grids);
  auto sampler =
      ids.empty() ? nullptr : std::make_unique<StoreSampler>(ids, grids, maxComponents(ids, grids));

  InterpreterOptions interpreterOptions;
  interpreterOptions.tileSize = options.tileSize;

  if (program.computeType() == ComputeType::F32) {
    return std::make_unique<InterpreterKernel<float>>(
        program, binding, std::move(lowered), std::move(sampler), interpreterOptions);
  }
  return std::make_unique<InterpreterKernel<double>>(
      program, binding, std::move(lowered), std::move(sampler), interpreterOptions);
}

} // namespace

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
  switch (options.preferred) {
  case BackendKind::Interpreter:
    return makeInterpreter(program, binding, grids, options);
  case BackendKind::RtcCpu: {
    // The grids have to be interned before the kernel exists either way, and a
    // compiled program that turns out to sample one falls through to the
    // interpreter below rather than failing.
    LoweredProgram lowered = lower(program, options.lowering);
    auto kernel = makeRtcCpuKernel(program, binding, std::move(lowered), options);
    if (kernel != nullptr) {
      internGrids(program, grids);
      return kernel;
    }
    if (!options.allowFallback) {
      logError() << "expr: the compiled CPU backend is unusable and fallback was disabled.";
      return nullptr;
    }
    return makeInterpreter(program, binding, grids, options);
  }
  case BackendKind::RtcCuda:
  case BackendKind::RtcHip:
    // Package 5. Not a stub that pretends: with allowFallback the caller gets a
    // working kernel and a warning naming what it did not get, which is the
    // behaviour makeKernel promises ("never returns null").
    if (options.allowFallback) {
      logWarning() << "expr: backend" << name(options.preferred)
                   << "is not implemented yet; using the interpreter instead.";
      return makeInterpreter(program, binding, grids, options);
    }
    logError() << "expr: backend" << name(options.preferred)
               << "is not implemented and fallback was disabled.";
    return nullptr;
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
