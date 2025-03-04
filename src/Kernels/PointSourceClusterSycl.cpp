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
    const auto* __restrict a = sources.A.data();
    const auto* __restrict stiffnessTensor = sources.stiffnessTensor.data();
    const auto* __restrict onsetTime = sources.onsetTime.data();
    const auto* __restrict samplingInterval = sources.samplingInterval.data();
    auto sampleOffsets = memory::AlignedArray<const std::size_t* __restrict, 3>();
    sampleOffsets[0] = sources.sampleOffsets[0].data();
    sampleOffsets[1] = sources.sampleOffsets[1].data();
    sampleOffsets[2] = sources.sampleOffsets[2].data();
    auto sample = memory::AlignedArray<const real* __restrict, 3>();
    sample[0] = sources.sample[0].data();
    sample[1] = sources.sample[1].data();
    sample[2] = sources.sample[2].data();

    auto* queue = reinterpret_cast<sycl::queue*>(runtime.stream());

    sycl::range rng{mapping.size()};
    if (sources.mode == sourceterm::PointSourceMode::Nrf) {
      queue->submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::item<1> id) {
          pointSourceKernelNRF(id[0],
                               from,
                               to,
                               mappingPtr,
                               mInvJInvPhisAtSources,
                               tensor,
                               a,
                               stiffnessTensor,
                               onsetTime,
                               samplingInterval,
                               sampleOffsets,
                               sample);
        });
      });
    } else {
      queue->submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::item<1> id) {
          pointSourceKernelFSRM(id[0],
                                from,
                                to,
                                mappingPtr,
                                mInvJInvPhisAtSources,
                                tensor,
                                a,
                                stiffnessTensor,
                                onsetTime,
                                samplingInterval,
                                sampleOffsets,
                                sample);
        });
      });
    }
  }
}

} // namespace seissol::kernels
