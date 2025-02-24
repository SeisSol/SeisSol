// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PointSourceClusterOnHost.h"

#include "generated_code/init.h"
#include "generated_code/kernel.h"

#include <Kernels/PointSourceCluster.h>
#include <Kernels/Precision.h>
#include <Parallel/Runtime/Stream.h>
#include <SourceTerm/Typedefs.h>
#include <memory>
#include <tensor.h>
#include <utility>

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
    for (unsigned m = 0; m < mapping.size(); ++m) {
      const unsigned startSource = mapping[m].pointSourcesOffset;
      const unsigned endSource = mapping[m].pointSourcesOffset + mapping[m].numberOfPointSources;
      if (sources_->mode == sourceterm::PointSourceMode::Nrf) {
        for (unsigned source = startSource; source < endSource; ++source) {
          addTimeIntegratedPointSourceNRF(source, from, to, *mapping[m].dofs);
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
          addTimeIntegratedPointSourceFSRM(source, from, to, *mapping[m].dofs);
        }
      }
    }
  }
}

unsigned PointSourceClusterOnHost::size() const { return sources_->numberOfSources; }

void PointSourceClusterOnHost::addTimeIntegratedPointSourceNRF(unsigned source,
                                                               double from,
                                                               double to,
                                                               real dofs[tensor::Q::size()]) {
  real slip[] = {0.0, 0.0, 0.0};
  for (unsigned i = 0; i < 3; ++i) {
    auto o0 = sources_->sampleOffsets[i][source];
    auto o1 = sources_->sampleOffsets[i][source + 1];
    slip[i] = computeSampleTimeIntegral(from,
                                        to,
                                        sources_->onsetTime[source],
                                        sources_->samplingInterval[source],
                                        sources_->sample[i].data() + o0,
                                        o1 - o0);
  }

  real rotatedSlip[] = {0.0, 0.0, 0.0};
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      rotatedSlip[j] += sources_->tensor[source][j + i * 3] * slip[i];
    }
  }

  kernel::sourceNRF krnl;
  krnl.Q = dofs;
  krnl.mInvJInvPhisAtSources = sources_->mInvJInvPhisAtSources[source].data();
  krnl.stiffnessTensor = sources_->stiffnessTensor[source].data();
  krnl.mSlip = rotatedSlip;
  krnl.mNormal = sources_->tensor[source].data() + 6;
  krnl.mArea = -sources_->A[source];
  krnl.momentToNRF = init::momentToNRF::Values;
#ifdef MULTIPLE_SIMULATIONS
  krnl.oneSimToMultSim = init::oneSimToMultSim::Values;
#endif
  krnl.execute();
}

void PointSourceClusterOnHost::addTimeIntegratedPointSourceFSRM(unsigned source,
                                                                double from,
                                                                double to,
                                                                real dofs[tensor::Q::size()]) {
  auto o0 = sources_->sampleOffsets[0][source];
  auto o1 = sources_->sampleOffsets[0][source + 1];
  auto slip = computeSampleTimeIntegral(from,
                                        to,
                                        sources_->onsetTime[source],
                                        sources_->samplingInterval[source],
                                        sources_->sample[0].data() + o0,
                                        o1 - o0);
  kernel::sourceFSRM krnl;
  krnl.Q = dofs;
  krnl.mInvJInvPhisAtSources = sources_->mInvJInvPhisAtSources[source].data();
  krnl.momentFSRM = sources_->tensor[source].data();
  krnl.stfIntegral = slip;
#ifdef MULTIPLE_SIMULATIONS
  krnl.oneSimToMultSim = init::oneSimToMultSim::Values;
#endif
  krnl.execute();
}

} // namespace seissol::kernels
