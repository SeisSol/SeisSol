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
constexpr std::size_t Blocksize = 256;

using namespace seissol::kernels;

__launch_bounds__(Blocksize) __global__ void launchKernelNRF(
    std::size_t numElements,
    double from,
    double to,
    sourceterm::CellToPointSourcesMapping* __restrict mappingPtr,
    const seissol::memory::
        AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>* __restrict mInvJInvPhisAtSources,
    const seissol::memory::AlignedArray<real,
                                        sourceterm::PointSources::TensorSize>* __restrict tensor,
    const real* __restrict a,
    const seissol::memory::AlignedArray<real, 81>* __restrict stiffnessTensor,
    const double* __restrict onsetTime,
    const double* __restrict samplingInterval,
    const memory::AlignedArray<const std::size_t* __restrict, 3> sampleOffsets,
    const memory::AlignedArray<const real* __restrict, 3> sample) {
  const auto index = threadIdx.x + blockDim.x * blockIdx.x;
  if (index < numElements) {
    pointSourceKernelNRF(index,
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
  }
}

__launch_bounds__(Blocksize) __global__ void launchKernelFSRM(
    std::size_t numElements,
    double from,
    double to,
    sourceterm::CellToPointSourcesMapping* __restrict mappingPtr,
    const seissol::memory::
        AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>* __restrict mInvJInvPhisAtSources,
    const seissol::memory::AlignedArray<real,
                                        sourceterm::PointSources::TensorSize>* __restrict tensor,
    const real* __restrict a,
    const seissol::memory::AlignedArray<real, 81>* __restrict stiffnessTensor,
    const double* __restrict onsetTime,
    const double* __restrict samplingInterval,
    const memory::AlignedArray<const std::size_t* __restrict, 3> sampleOffsets,
    const memory::AlignedArray<const real* __restrict, 3> sample) {
  const auto index = threadIdx.x + blockDim.x * blockIdx.x;
  if (index < numElements) {
    pointSourceKernelFSRM(index,
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
  }
}

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

#ifdef __CUDACC__
    using StreamT = cudaStream_t;
#endif
#ifdef __HIP__
    using StreamT = hipStream_t;
#endif
    auto stream = reinterpret_cast<StreamT>(runtime.stream());

    dim3 block(Blocksize);
    dim3 grid((mapping.size() + Blocksize - 1) / Blocksize);

    // special case for a smaller grid
    if (grid.x == 1) {
      block.x = mapping.size();
    }

    if (sources.mode == sourceterm::PointSourceMode::Nrf) {
      launchKernelNRF<<<grid, block, 0, stream>>>(mapping.size(),
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
    } else {
      launchKernelFSRM<<<grid, block, 0, stream>>>(mapping.size(),
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
    }
  }
}

} // namespace seissol::kernels
