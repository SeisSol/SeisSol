// SPDX-FileCopyrightText: 2015 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PointSourceClusterOnHost.h"

#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"

#include <GeneratedCode/tensor.h>
#include <Kernels/Common.h>
#include <Kernels/PointSourceCluster.h>
#include <Kernels/Precision.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/MultipleSimulations.h>
#include <SourceTerm/Typedefs.h>
#include <array>
#include <cstddef>
#include <memory>
#include <utility>

GENERATE_HAS_MEMBER(sourceToMultSim)

namespace seissol::kernels {

PointSourceClusterOnHost::PointSourceClusterOnHost(
    std::shared_ptr<sourceterm::ClusterMapping> mapping,
    std::shared_ptr<sourceterm::PointSources> sources)
    : clusterMapping_(std::move(mapping)), sources_(std::move(sources)) {}

void PointSourceClusterOnHost::addTimeIntegratedPointSources(
    double from, double to, seissol::parallel::runtime::StreamRuntime& runtime) {
  auto& mapping = clusterMapping_->cellToSources;
  if (mapping.size() > 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t m = 0; m < mapping.size(); ++m) {
      const auto startSource = mapping[m].pointSourcesOffset;
      const auto endSource = mapping[m].pointSourcesOffset + mapping[m].numberOfPointSources;
      for (auto source = startSource; source < endSource; ++source) {
        addTimeIntegratedPointSource(source, from, to, *mapping[m].dofs);
      }
    }
  }
}

std::size_t PointSourceClusterOnHost::size() const { return sources_->numberOfSources; }

void PointSourceClusterOnHost::addTimeIntegratedPointSource(std::size_t source,
                                                            double from,
                                                            double to,
                                                            real dofs[tensor::Q::size()]) {
  std::array<real, Quantities> update{};
  const auto base = sources_->sampleRange[source];
  const auto localSamples = sources_->sampleRange[source + 1] - base;

  const auto* __restrict tensorLocal = sources_->tensor.data() + base * Quantities;

  for (std::size_t i = 0; i < localSamples; ++i) {
    const auto o0 = sources_->sampleOffsets[i + base];
    const auto o1 = sources_->sampleOffsets[i + base + 1];
    const auto slip = computeSampleTimeIntegral(from,
                                                to,
                                                sources_->onsetTime[source],
                                                sources_->samplingInterval[source],
                                                sources_->sample.data() + o0,
                                                o1 - o0);

#pragma omp simd
    for (std::size_t t = 0; t < Quantities; ++t) {
      update[t] += slip * tensorLocal[t + i * Quantities];
    }
  }

  kernel::addPointSource krnl;
  krnl.update = update.data();
  krnl.Q = dofs;
  krnl.mInvJInvPhisAtSources = sources_->mInvJInvPhisAtSources[source].data();

  const auto simulationIndex = sources_->simulationIndex[source];
  std::array<real, seissol::multisim::NumSimulations> sourceToMultSim{};
  sourceToMultSim[simulationIndex] = 1.0;
  set_sourceToMultSim(krnl, sourceToMultSim.data());
  krnl.execute();
}

} // namespace seissol::kernels
