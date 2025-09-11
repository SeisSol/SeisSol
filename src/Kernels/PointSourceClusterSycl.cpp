// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PointSourceCluster.h"

#include <Memory/MemoryAllocator.h>
#include <cstddef>

#include <sycl/sycl.hpp>

namespace {

constexpr std::size_t SubBlock = 64;
constexpr std::size_t Blocksize = 256;
constexpr auto PerBlock = Blocksize / SubBlock;

} // namespace

namespace seissol::kernels {

void pointSourceKernel(sourceterm::ClusterMapping& clusterMapping,
                       sourceterm::PointSources& sources,
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

    auto* queue = reinterpret_cast<sycl::queue*>(runtime.stream());

    const auto elements = mapping.size();

    sycl::nd_range<2> rng{{elements * SubBlock, PerBlock}, {SubBlock, PerBlock}};
    queue->submit([&](sycl::handler& cgh) {
      cgh.parallel_for(rng, [=](sycl::nd_item<2> item) {
        const auto block = item.get_group().get_group_id(0) * PerBlock + item.get_local_id(1);
        const auto thread = item.get_local_id(0);

        if (block < elements) {
          pointSourceKernelDevice<SubBlock>(thread,
                                            block,
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
      });
    });
  }
}

} // namespace seissol::kernels
