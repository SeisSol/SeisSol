// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Config.h"
#include "Memory/MemoryAllocator.h"
#include "Parallel/Runtime/Stream.h"
#include "PointSourceCluster.h"
#include "SourceTerm/Typedefs.h"

#include <cstddef>
#include <sycl/sycl.hpp>

namespace {

constexpr std::size_t SubBlock = 64;
constexpr std::size_t Blocksize = 256;
constexpr auto PerBlock = Blocksize / SubBlock;

} // namespace

namespace seissol::kernels {

template <typename Cfg>
void pointSourceKernel(sourceterm::ClusterMapping& clusterMapping,
                       sourceterm::PointSources<Cfg>& sources,
                       double from,
                       double timeStepSize,
                       const double* clock,
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

    auto& queue = *reinterpret_cast<sycl::queue*>(runtime.stream());

    const auto elements = mapping.size();

    // As the CUDA/HIP kernel: a group takes PerBlock cells, with SubBlock work-items each. The
    // work-items of a cell go along the last dimension, which SYCL linearizes fastest and cuts the
    // sub-groups along (CUDA's x), so that neighboring work-items update neighboring degrees of
    // freedom of the same cell.
    const auto groups = (elements + PerBlock - 1) / PerBlock;
    sycl::nd_range<2> rng{{groups * PerBlock, SubBlock}, {PerBlock, SubBlock}};
    queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<2> item) {
        const auto block = item.get_group().get_group_id(0) * PerBlock + item.get_local_id(0);
        const auto thread = item.get_local_id(1);

        if (block < elements) {
          const double start = clock == nullptr ? from : *clock;
          const double end = start + timeStepSize;
          pointSourceKernelDevice<Cfg, SubBlock>(thread,
                                                 block,
                                                 start,
                                                 end,
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
      });
    });
  }
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template void pointSourceKernel<Cfg>(sourceterm::ClusterMapping & clusterMapping,                \
                                       sourceterm::PointSources<Cfg> & sources,                    \
                                       double from,                                                \
                                       double timeStepSize,                                        \
                                       const double* clock,                                        \
                                       seissol::parallel::runtime::StreamRuntime& runtime);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels
