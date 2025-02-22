#include "PointSourceCluster.h"

namespace seissol::kernels {

void pointSourceKernel(sourceterm::ClusterMapping& clusterMapping,
                       sourceterm::PointSources& sources,
                       double from,
                       double to,
                       seissol::parallel::runtime::StreamRuntime& runtime) {
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
