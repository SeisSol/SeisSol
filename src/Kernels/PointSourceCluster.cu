#include "PointSourceCluster.h"
#include <Initializer/MemoryAllocator.h>
#include <cstddef>

namespace {
constexpr std::size_t Blocksize = 256;

using namespace seissol::kernels;

__launch_bounds__(Blocksize)
__global__ void launchKernelNRF(std::size_t numElements, 
  double from, double to,
  sourceterm::CellToPointSourcesMapping* mappingPtr,
  seissol::memory::AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>* mInvJInvPhisAtSources,
  seissol::memory::AlignedArray<real, sourceterm::PointSources::TensorSize>* tensor,
  real* a,
  seissol::memory::AlignedArray<real, 81>* stiffnessTensor,
  double* onsetTime,
  double* samplingInterval,
  memory::AlignedArray<std::size_t*, 3> sampleOffsets,
  memory::AlignedArray<real*, 3> sample) {
    const auto index = threadIdx.x + blockDim.x * blockIdx.x;
    if (index < numElements) {
        pointSourceKernelNRF(index, from, to, mappingPtr, mInvJInvPhisAtSources, tensor, a, stiffnessTensor, onsetTime, samplingInterval,
                    sampleOffsets, sample);
    }
}

__launch_bounds__(Blocksize)
__global__ void launchKernelFSRM(std::size_t numElements, 
  double from, double to,
  sourceterm::CellToPointSourcesMapping* mappingPtr,
  seissol::memory::AlignedArray<real, tensor::mInvJInvPhisAtSources::size()>* mInvJInvPhisAtSources,
  seissol::memory::AlignedArray<real, sourceterm::PointSources::TensorSize>* tensor,
  real* a,
  seissol::memory::AlignedArray<real, 81>* stiffnessTensor,
  double* onsetTime,
  double* samplingInterval,
  memory::AlignedArray<std::size_t*, 3> sampleOffsets,
  memory::AlignedArray<real*, 3> sample) {
    const auto index = threadIdx.x + blockDim.x * blockIdx.x;
    if (index < numElements) {
        pointSourceKernelFSRM(index, from, to, mappingPtr, mInvJInvPhisAtSources, tensor, a, stiffnessTensor, onsetTime, samplingInterval,
                    sampleOffsets, sample);
    }
}

} // namespace

namespace seissol::kernels {

void pointSourceKernel(sourceterm::ClusterMapping& clusterMapping, sourceterm::PointSources& sources, double from, double to, seissol::parallel::runtime::StreamRuntime& runtime) {
    auto& mapping = clusterMapping.cellToSources;
    auto* mappingPtr = mapping.data();
    if (mapping.size() > 0) {
        auto* mInvJInvPhisAtSources = sources.mInvJInvPhisAtSources.data();
        auto* tensor = sources.tensor.data();
        auto* a = sources.A.data();
        auto* stiffnessTensor = sources.stiffnessTensor.data();
        auto* onsetTime = sources.onsetTime.data();
        auto* samplingInterval = sources.samplingInterval.data();

        auto sampleOffsets = memory::AlignedArray<std::size_t*, 3>();
        sampleOffsets[0] = sources.sampleOffsets[0].data();
        sampleOffsets[1] = sources.sampleOffsets[1].data();
        sampleOffsets[2] = sources.sampleOffsets[2].data();
        auto sample = memory::AlignedArray<real*, 3>();
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
            launchKernelNRF<<<grid, block, 0, stream>>>(mapping.size(), from, to, mappingPtr, mInvJInvPhisAtSources, tensor, a, stiffnessTensor, onsetTime, samplingInterval,
                sampleOffsets, sample);
        } else {
            launchKernelFSRM<<<grid, block, 0, stream>>>(mapping.size(), from, to, mappingPtr, mInvJInvPhisAtSources, tensor, a, stiffnessTensor, onsetTime, samplingInterval,
                sampleOffsets, sample);
        }
    }
}

} // namespace seissol::kernels
