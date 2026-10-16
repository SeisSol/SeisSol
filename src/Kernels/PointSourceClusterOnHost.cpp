// SPDX-FileCopyrightText: 2015 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PointSourceClusterOnHost.h"

#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/PointSourceCluster.h"
#include "Parallel/Runtime/Stream.h"
#include "SourceTerm/Typedefs.h"

#include <array>
#include <cstddef>
#include <memory>
#include <utility>

GENERATE_HAS_MEMBER(oneSimToMultSim)

namespace seissol::kernels {

template <typename Cfg>
PointSourceClusterOnHost<Cfg>::PointSourceClusterOnHost(
    std::shared_ptr<sourceterm::ClusterMapping> mapping,
    std::shared_ptr<sourceterm::PointSources<Cfg>> sources)
    : clusterMapping_(std::move(mapping)), sources_(std::move(sources)) {}

template <typename Cfg>
void PointSourceClusterOnHost<Cfg>::addTimeIntegratedPointSources(
    double from, double to, seissol::parallel::runtime::StreamRuntime& /*runtime*/) {
  auto& mapping = clusterMapping_->cellToSources;
  if (mapping.size() > 0) {

#pragma omp parallel for schedule(static)
    for (std::size_t m = 0; m < mapping.size(); ++m) {
      const auto startSource = mapping[m].pointSourcesOffset;
      const auto endSource = mapping[m].pointSourcesOffset + mapping[m].numberOfPointSources;
      for (auto source = startSource; source < endSource; ++source) {
        addTimeIntegratedPointSource(source, from, to, static_cast<real*>(mapping[m].dofs));
      }
    }
  }
}

template <typename Cfg>
std::size_t PointSourceClusterOnHost<Cfg>::size() const {
  return sources_->numberOfSources;
}

template <typename Cfg>
void PointSourceClusterOnHost<Cfg>::addTimeIntegratedPointSource(
    std::size_t source, double from, double to, real dofs[tensor::Q<Cfg>::size()]) {
  constexpr auto QuantityCount = Quantities<Cfg>;
  std::array<real, QuantityCount> update{};
  const auto base = sources_->sampleRange[source];
  const auto localSamples = sources_->sampleRange[source + 1] - base;

  const auto* __restrict tensorLocal = sources_->tensor.data() + base * QuantityCount;

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
    for (std::size_t t = 0; t < QuantityCount; ++t) {
      update[t] += slip * tensorLocal[t + i * QuantityCount];
    }
  }

  kernel::addPointSource<Cfg> krnl;
  krnl.update = update.data();
  krnl.Q = dofs;
  krnl.mInvJInvPhisAtSources = sources_->mInvJInvPhisAtSources[source].data();

  const auto simulationIndex = sources_->simulationIndex[source];
  std::array<real, Cfg::NumSimulations> sourceToMultSim{};
  sourceToMultSim[simulationIndex] = 1.0;
  set_oneSimToMultSim(krnl, sourceToMultSim.data());
  krnl.execute();
}

#define SEISSOL_INSTANTIATE(Cfg) template class PointSourceClusterOnHost<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels
