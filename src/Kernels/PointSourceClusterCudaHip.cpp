// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PointSourceCluster.h"
#include <Memory/MemoryAllocator.h>
#include <cstddef>

#ifdef __HIP__
#include "hip/hip_runtime.h"
#endif

namespace {
constexpr std::size_t SubBlock = 64;
constexpr std::size_t Blocksize = 256;
constexpr auto PerBlock = Blocksize / SubBlock;

using namespace seissol::kernels;

template <typename Cfg>
__launch_bounds__(Blocksize) __global__ void launchKernel(
    std::size_t numElements,
    double from,
    double to,
    sourceterm::CellToPointSourcesMapping* __restrict mappingPtr,
    const seissol::memory::AlignedArray<
        real,
        tensor::mInvJInvPhisAtSources<Cfg>::size()>* __restrict mInvJInvPhisAtSources,
    const std::uint32_t* __restrict simulationIndex,
    const real* __restrict tensor,
    const double* __restrict onsetTime,
    const double* __restrict samplingInterval,
    const std::size_t* __restrict sampleRange,
    const std::size_t* __restrict sampleOffsets,
    const real* __restrict sample) {
  const auto index = threadIdx.y + PerBlock * blockIdx.x;
  if (index < numElements) {
    pointSourceKernelDevice<Cfg, SubBlock>(threadIdx.x,
                                           index,
                                           from,
                                           to,
                                           mappingPtr,
                                           mInvJInvPhisAtSources,
                                           simulationIndex,
                                           tensor,
                                           onsetTime,
                                           samplingInterval,
                                           sampleRange,
                                           sampleOffsets,
                                           sample);
  }
}

} // namespace

namespace seissol::kernels {

template <typename Cfg>
void pointSourceKernel(sourceterm::ClusterMapping& clusterMapping,
                       sourceterm::PointSources<Cfg>& sources,
                       double from,
                       double to,
                       seissol::parallel::runtime::StreamRuntime& runtime) {
  auto& mapping = clusterMapping.cellToSources;
  auto* __restrict mappingPtr = mapping.data();
  if (mapping.size() > 0) {
    const auto* __restrict mInvJInvPhisAtSources = sources.mInvJInvPhisAtSources.data();
    const auto* __restrict tensor = sources.tensor.data();
    const auto* __restrict onsetTime = sources.onsetTime.data();
    const auto* __restrict samplingInterval = sources.samplingInterval.data();
    const auto* __restrict sampleRange = sources.sampleRange.data();
    const auto* __restrict sampleOffsets = sources.sampleOffsets.data();
    const auto* __restrict sample = sources.sample.data();
    const auto* __restrict simulationIndex = sources.simulationIndex.data();

#ifdef __CUDACC__
    using StreamT = cudaStream_t;
#endif
#ifdef __HIP__
    using StreamT = hipStream_t;
#endif
    auto stream = reinterpret_cast<StreamT>(runtime.stream());

    dim3 block(SubBlock, PerBlock);
    dim3 grid((mapping.size() + PerBlock - 1) / PerBlock);

    // special case for a smaller grid
    if (grid.x == 1) {
      block.y = mapping.size();
    }

    launchKernel<Cfg><<<grid, block, 0, stream>>>(mapping.size(),
                                                  from,
                                                  to,
                                                  mappingPtr,
                                                  mInvJInvPhisAtSources,
                                                  simulationIndex,
                                                  tensor,
                                                  onsetTime,
                                                  samplingInterval,
                                                  sampleRange,
                                                  sampleOffsets,
                                                  sample);
  }
}

#define SEISSOL_CONFIGITER(cfg)                                                                    \
  template void pointSourceKernel(sourceterm::ClusterMapping& clusterMapping,                      \
                                  sourceterm::PointSources<cfg>& sources,                          \
                                  double from,                                                     \
                                  double to,                                                       \
                                  seissol::parallel::runtime::StreamRuntime& runtime);
#include "ConfigInclude.h"

} // namespace seissol::kernels
