// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Interpolation.h"

#include "Grid.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace seissol::reader::datafield {

void interpolationWeights(Interpolation type, double s, double* weights) {
  switch (type) {
  case Interpolation::Nearest:
    weights[0] = 1.0;
    return;
  case Interpolation::Linear:
    weights[0] = 1.0 - s;
    weights[1] = s;
    return;
  case Interpolation::Cubic: {
    // Keys' cubic convolution, a = -1/2. Written in this factored form because
    // easi writes it in this form; a mathematically equivalent regrouping would
    // move the last bits and eat into the cross-path tolerance for nothing.
    const double s2 = s * s;
    const double s3 = s2 * s;
    weights[0] = 0.5 * (-s3 + 2.0 * s2 - s);
    weights[1] = 0.5 * (3.0 * s3 - 5.0 * s2 + 2.0);
    weights[2] = 0.5 * (-3.0 * s3 + 4.0 * s2 + s);
    weights[3] = 0.5 * (s3 - s2);
    return;
  }
  }
}

void linearFallbackWeights(unsigned width, unsigned slotOfBase, double s, double* weights) {
  for (unsigned j = 0; j < width; ++j) {
    weights[j] = 0.0;
  }
  weights[slotOfBase] = 1.0 - s;
  if (slotOfBase + 1 < width) {
    weights[slotOfBase + 1] = s;
  }
}

void SampleScratch::reserve(unsigned dimensions,
                            unsigned width,
                            unsigned components,
                            std::size_t count) {
  std::size_t stencilSize = 1;
  for (unsigned d = 0; d < dimensions; ++d) {
    stencilSize *= width;
  }
  const std::size_t block = stencilSize * components * count;
  if (weights.size() < static_cast<std::size_t>(dimensions) * width * count) {
    weights.resize(static_cast<std::size_t>(dimensions) * width * count);
  }
  if (start.size() < static_cast<std::size_t>(dimensions) * count) {
    start.resize(static_cast<std::size_t>(dimensions) * count);
  }
  if (gather.size() < block) {
    gather.resize(block);
  }
  if (work.size() < block) {
    work.resize(block);
  }
}

namespace {

/// Phase 1 for one axis and one lane. Transcribed from easi's Grid<Derived>
/// locate branches; the branch structure is load-bearing and is documented
/// where it is not obvious.
template <Interpolation Type>
inline bool locateAxis(
    const ArrayView& view, unsigned axis, double x, std::int64_t& windowStart, double* weights) {
  constexpr StencilKernel Kernel = kernelOf(Type);
  constexpr unsigned Width = Kernel.width;

  const std::int64_t num = static_cast<std::int64_t>(view.num[axis]);
  const double top = static_cast<double>(num - 1);

  double raw = (x - view.min[axis]) * view.deltaInv[axis];
  bool clamped = false;
  // Written as `!(raw > 0.0)` rather than `raw < 0.0` so that NaN takes this
  // branch instead of falling through to the index cast.
  if (!(raw > 0.0)) {
    clamped = (raw < 0.0) || std::isnan(raw);
    raw = 0.0;
  } else if (raw > top) {
    clamped = true;
    raw = top;
  }

  if constexpr (Kernel.roundToNearest) {
    windowStart = static_cast<std::int64_t>(std::lround(raw));
    weights[0] = 1.0;
    return clamped;
  } else {
    // The base cell is the last one that has a right neighbour, so a query
    // sitting exactly on the final node uses s = 1 in the previous cell rather
    // than opening a cell that does not exist.
    const std::int64_t maxBase = (num >= 2) ? num - 2 : 0;
    const double lower = std::floor(raw);
    std::int64_t index = static_cast<std::int64_t>(lower);
    double s = 0.0;
    if (index > maxBase) {
      index = maxBase;
      s = (num >= 2) ? 1.0 : 0.0;
    } else {
      s = raw - lower;
    }

    if constexpr (Width <= 2) {
      windowStart = index;
      interpolationWeights(Type, s, weights);
    } else {
      const std::int64_t maxStart = num - static_cast<std::int64_t>(Width);
      windowStart = (num >= static_cast<std::int64_t>(Width)) ? index + Kernel.offset : 0;
      if (windowStart < 0) {
        windowStart = 0;
      } else if (windowStart > maxStart) {
        windowStart = (maxStart > 0) ? maxStart : 0;
      }
      // The full stencil is only valid when the base cell sits in the slot the
      // scheme expects it in. Anywhere else the window has been shifted to stay
      // inside the array, and using the full weights would extrapolate.
      const unsigned slotOfBase = static_cast<unsigned>(index - windowStart);
      if (num >= static_cast<std::int64_t>(Width) &&
          slotOfBase == static_cast<unsigned>(-Kernel.offset)) {
        interpolationWeights(Type, s, weights);
      } else {
        linearFallbackWeights(Width, slotOfBase, s, weights);
      }
    }
    return clamped;
  }
}

template <typename T>
inline double loadAs(const void* data, std::size_t index) {
  return static_cast<double>(static_cast<const T*>(data)[index]);
}

template <Interpolation Type, typename T>
std::size_t sampleBatchTyped(const ArrayView& view,
                             const double* x,
                             std::size_t count,
                             double* out,
                             SampleScratch& scratch) {
  constexpr unsigned Width = kernelOf(Type).width;
  const unsigned dim = view.dimensions;
  const unsigned comps = view.components;

  std::size_t stencilSize = 1;
  for (unsigned d = 0; d < dim; ++d) {
    stencilSize *= Width;
  }
  scratch.reserve(dim, Width, comps, count);

  double* weights = scratch.weights.data();
  std::int64_t* start = scratch.start.data();

  // ---- phase 1: locate and weights, SoA ------------------------------------
  std::size_t clampedLanes = 0;
  for (std::size_t lane = 0; lane < count; ++lane) {
    bool any = false;
    for (unsigned d = 0; d < dim; ++d) {
      double axisWeights[MaxStencilWidth];
      std::int64_t windowStart = 0;
      any |= locateAxis<Type>(
          view, d, x[static_cast<std::size_t>(d) * count + lane], windowStart, axisWeights);
      start[static_cast<std::size_t>(d) * count + lane] = windowStart;
      for (unsigned j = 0; j < Width; ++j) {
        weights[(static_cast<std::size_t>(d) * Width + j) * count + lane] = axisWeights[j];
      }
    }
    clampedLanes += any ? 1 : 0;
  }

  // ---- phase 2: gather -----------------------------------------------------
  // Random access, will not vectorise. Kept separate from the reduction for
  // exactly that reason.
  double* gather = scratch.gather.data();
  unsigned stencilIndex[MaxGridDimensions] = {};
  for (std::size_t f = 0; f < stencilSize; ++f) {
    for (std::size_t lane = 0; lane < count; ++lane) {
      std::size_t offset = 0;
      for (unsigned d = 0; d < dim; ++d) {
        const std::int64_t limit = static_cast<std::int64_t>(view.num[d]) - 1;
        std::int64_t node = start[static_cast<std::size_t>(d) * count + lane] +
                            static_cast<std::int64_t>(stencilIndex[d]);
        // Only the upper end needs clamping here: windowStart is already >= 0.
        if (node > limit) {
          node = limit;
        }
        offset += static_cast<std::size_t>(node) * view.stride[d];
      }
      for (unsigned c = 0; c < comps; ++c) {
        gather[(f * comps + c) * count + lane] =
            loadAs<T>(view.data, offset + static_cast<std::size_t>(c) * view.componentStride);
      }
    }
    for (unsigned d = 0; d < dim; ++d) {
      if (++stencilIndex[d] < Width) {
        break;
      }
      stencilIndex[d] = 0;
    }
  }

  // ---- phase 3: tensor-product reduction, lane loop innermost --------------
  // Ping-pong rather than in place: vectorised over lanes there is no scalar
  // accumulator, so an in-place pass would zero the k == 0 term before reading
  // it. See the note on sampleBatch in Grid.h — the failure is invisible in 1-D.
  const double* source = gather;
  double* target = scratch.work.data();
  std::size_t cnt = stencilSize;
  for (int d = static_cast<int>(dim) - 1; d >= 0; --d) {
    cnt /= Width;
    for (std::size_t p = 0; p < cnt; ++p) {
      for (unsigned c = 0; c < comps; ++c) {
        double* dst = target + (p * comps + c) * count;
        for (std::size_t lane = 0; lane < count; ++lane) {
          dst[lane] = 0.0;
        }
        for (unsigned k = 0; k < Width; ++k) {
          const double* wk = weights + (static_cast<std::size_t>(d) * Width + k) * count;
          const double* sk = source + ((p + k * cnt) * comps + c) * count;
          for (std::size_t lane = 0; lane < count; ++lane) {
            dst[lane] += wk[lane] * sk[lane];
          }
        }
      }
    }
    double* next = const_cast<double*>(source);
    source = target;
    target = next;
  }

  for (unsigned c = 0; c < comps; ++c) {
    const double* src = source + static_cast<std::size_t>(c) * count;
    double* dst = out + static_cast<std::size_t>(c) * count;
    for (std::size_t lane = 0; lane < count; ++lane) {
      dst[lane] = src[lane];
    }
  }
  return clampedLanes;
}

template <Interpolation Type>
std::size_t dispatchType(const ArrayView& view,
                         const double* x,
                         std::size_t count,
                         double* out,
                         SampleScratch& scratch) {
  return (view.type == ElementType::Float)
             ? sampleBatchTyped<Type, float>(view, x, count, out, scratch)
             : sampleBatchTyped<Type, double>(view, x, count, out, scratch);
}

} // namespace

std::size_t sampleBatch(const ArrayView& view,
                        Interpolation interp,
                        const double* x,
                        std::size_t count,
                        double* out,
                        SampleScratch& scratch) {
  if (count == 0 || view.data == nullptr) {
    return 0;
  }
  switch (interp) {
  case Interpolation::Nearest:
    return dispatchType<Interpolation::Nearest>(view, x, count, out, scratch);
  case Interpolation::Linear:
    return dispatchType<Interpolation::Linear>(view, x, count, out, scratch);
  case Interpolation::Cubic:
    return dispatchType<Interpolation::Cubic>(view, x, count, out, scratch);
  }
  return 0;
}

std::size_t sample(const ArrayView& view, Interpolation interp, const double* x, double* out) {
  SampleScratch scratch;
  return sampleBatch(view, interp, x, 1, out, scratch);
}

namespace {

/// Slot range within a window of `count` slices that `kernel` can serve at full
/// accuracy, as offsets from the window start. Empty when the window is too
/// narrow for the scheme at all.
std::pair<std::int64_t, std::int64_t> usableSlots(StencilKernel kernel, std::size_t count) {
  const auto width = static_cast<std::int64_t>(kernel.width);
  const std::int64_t lo = -kernel.offset;
  const std::int64_t hi = static_cast<std::int64_t>(count) - width - kernel.offset;
  return {lo, hi};
}

} // namespace

bool windowServes(const TimeWindow& window,
                  std::int64_t base,
                  StencilKernel kernel,
                  std::size_t numSlices) {
  if (window.count == 0) {
    return false;
  }
  if (window.count >= numSlices) {
    return true; // fully resident: nothing to move
  }
  const auto [lo, hi] = usableSlots(kernel, window.count);
  const std::int64_t local = base - static_cast<std::int64_t>(window.start);
  if (local >= lo && local <= hi) {
    return true;
  }
  // At the file ends the window is already as far as it goes, and the edge
  // fallback is the correct answer rather than a reason to reload.
  const std::int64_t maxStart = static_cast<std::int64_t>(numSlices - window.count);
  if (window.start == 0 && local < lo) {
    return true;
  }
  if (static_cast<std::int64_t>(window.start) == maxStart && local > hi) {
    return true;
  }
  return false;
}

TimeWindow placeTimeWindow(std::int64_t base,
                           StencilKernel kernel,
                           std::size_t resident,
                           std::size_t numSlices) {
  if (numSlices == 0) {
    return {};
  }
  const std::size_t count = (resident < numSlices) ? resident : numSlices;
  if (count >= numSlices) {
    return {0, numSlices};
  }
  const std::int64_t maxStart = static_cast<std::int64_t>(numSlices - count);
  std::int64_t start = base + kernel.offset;
  if (start < 0) {
    start = 0;
  } else if (start > maxStart) {
    start = maxStart;
  }
  return {static_cast<std::size_t>(start), count};
}

} // namespace seissol::reader::datafield
