// Copyright (c) 2024 SeisSol Group
// Copyright (c) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "PointSourceClusterOnDevice.h"

#include "generated_code/init.h"
#include "generated_code/tensor.h"

// needs to be loaded after Eigen at the moment, due to SYCL
#include "Parallel/AcceleratorDevice.h"

#include <Kernels/PointSourceCluster.h>
#include <Kernels/precision.hpp>
#include <Parallel/Runtime/Stream.hpp>
#include <SourceTerm/typedefs.hpp>
#include <array>
#include <cstddef>
#include <memory>

#include "Numerical/SyclFunctions.h"

namespace seissol::kernels {

PointSourceClusterOnDevice::PointSourceClusterOnDevice(
    std::shared_ptr<sourceterm::ClusterMapping> mapping,
    std::shared_ptr<sourceterm::PointSources> sources)
    : clusterMapping_(mapping), sources_(sources) {}

unsigned PointSourceClusterOnDevice::size() const { return sources_->numberOfSources; }

void PointSourceClusterOnDevice::addTimeIntegratedPointSources(
    double from, double to, seissol::parallel::runtime::StreamRuntime& runtime) {
  auto& queue = seissol::AcceleratorDevice::getInstance().getSyclDefaultQueue();
  auto& mapping = clusterMapping_->cellToSources;
  if (mapping.size() > 0) {
    runtime.syncToSycl(&queue);

    auto* mappingPtr = mapping.data();
    auto* mInvJInvPhisAtSources = sources_->mInvJInvPhisAtSources.data();
    auto* tensor = sources_->tensor.data();
    auto* a = sources_->A.data();
    auto* stiffnessTensor = sources_->stiffnessTensor.data();
    auto* onsetTime = sources_->onsetTime.data();
    auto* samplingInterval = sources_->samplingInterval.data();
    auto sampleOffsets = std::array<std::size_t*, 3u>{sources_->sampleOffsets[0].data(),
                                                      sources_->sampleOffsets[1].data(),
                                                      sources_->sampleOffsets[2].data()};
    auto sample = std::array<real*, 3u>{
        sources_->sample[0].data(), sources_->sample[1].data(), sources_->sample[2].data()};

    sycl::range rng{mapping.size()};
    if (sources_->mode == sourceterm::PointSources::NRF) {
      queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::item<1> id) {
          const unsigned startSource = mappingPtr[id[0]].pointSourcesOffset;
          const unsigned endSource =
              mappingPtr[id[0]].pointSourcesOffset + mappingPtr[id[0]].numberOfPointSources;
          for (unsigned source = startSource; source < endSource; ++source) {
            std::array<real, 3u> slip;
            for (int i = 0; i < 3; ++i) {
              auto o0 = sampleOffsets[i][source];
              auto o1 = sampleOffsets[i][source + 1];
              slip[i] = computeSampleTimeIntegral<seissol::functions::SyclStdFunctions>(
                  from, to, onsetTime[source], samplingInterval[source], sample[i] + o0, o1 - o0);
            }

            addTimeIntegratedPointSourceNRF(slip,
                                            mInvJInvPhisAtSources[source].data(),
                                            tensor[source].data(),
                                            a[source],
                                            stiffnessTensor[source].data(),
                                            from,
                                            to,
                                            *mappingPtr[id[0]].dofs);
          }
        });
      });
    } else {
      queue.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(rng, [=](sycl::item<1> id) {
          const unsigned startSource = mappingPtr[id[0]].pointSourcesOffset;
          const unsigned endSource =
              mappingPtr[id[0]].pointSourcesOffset + mappingPtr[id[0]].numberOfPointSources;
          for (unsigned source = startSource; source < endSource; ++source) {
            auto o0 = sampleOffsets[0][source];
            auto o1 = sampleOffsets[0][source + 1];
            const real slip = computeSampleTimeIntegral<seissol::functions::SyclStdFunctions>(
                from, to, onsetTime[source], samplingInterval[source], sample[0] + o0, o1 - o0);
            addTimeIntegratedPointSourceFSRM(slip,
                                             mInvJInvPhisAtSources[source].data(),
                                             tensor[source].data(),
                                             from,
                                             to,
                                             *mappingPtr[id[0]].dofs);
          }
        });
      });
    }
    runtime.syncFromSycl(&queue);
  }
}

// workaround for NVHPC (using constexpr arrays directly caused errors in 24.01)
constexpr std::size_t QSpan = init::Q::Stop[0] - init::Q::Start[0];
constexpr std::size_t MomentFsrmSpan = tensor::momentFSRM::Shape[0];
constexpr std::size_t MInvJInvPhisAtSourcesSpan = tensor::mInvJInvPhisAtSources::Shape[0];

void PointSourceClusterOnDevice::addTimeIntegratedPointSourceNRF(const std::array<real, 3>& slip,
                                                                 const real* mInvJInvPhisAtSources,
                                                                 const real* tensor,
                                                                 real a,
                                                                 const real* stiffnessTensor,
                                                                 double from,
                                                                 double to,
                                                                 real dofs[tensor::Q::size()]) {
  real rotatedSlip[3] = {real(0.0)};
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      rotatedSlip[j] += tensor[j + i * 3] * slip[i];
    }
  }

  auto mom = [&](unsigned p, unsigned q) {
    real m = 0.0;
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        m += -a * stiffnessTensor[p + 3 * q + 9 * i + 27 * j] * rotatedSlip[i] * tensor[6 + j];
      }
    }
    return m;
  };

  const real moment[6] = {mom(0, 0), mom(1, 1), mom(2, 2), mom(0, 1), mom(1, 2), mom(0, 2)};
  for (unsigned t = 0; t < 6; ++t) {
    for (unsigned k = 0; k < MInvJInvPhisAtSourcesSpan; ++k) {
      dofs[k + t * QSpan] += mInvJInvPhisAtSources[k] * moment[t];
    }
  }
}

void PointSourceClusterOnDevice::addTimeIntegratedPointSourceFSRM(real slip,
                                                                  const real* mInvJInvPhisAtSources,
                                                                  const real* tensor,
                                                                  double from,
                                                                  double to,
                                                                  real* dofs) {
  for (unsigned p = 0; p < MomentFsrmSpan; ++p) {
    for (unsigned k = 0; k < MInvJInvPhisAtSourcesSpan; ++k) {
      dofs[k + p * QSpan] += slip * mInvJInvPhisAtSources[k] * tensor[p];
    }
  }
}

} // namespace seissol::kernels
