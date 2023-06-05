// Copyright (c) 2015-2020 SeisSol Group
// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "PointSourceClusterOnHost.h"

#include <generated_code/kernel.h>
#include <generated_code/init.h>
#include <SourceTerm/PointSource.h>

#include <utility>

namespace seissol::kernels {

PointSourceClusterOnHost::PointSourceClusterOnHost(sourceterm::ClusterMapping mapping,
                                                   sourceterm::PointSources sources)
    : clusterMapping_(std::move(mapping)), sources_(std::move(sources)) {}

void PointSourceClusterOnHost::addTimeIntegratedPointSources(double from, double to) {
  auto& mapping = clusterMapping_.cellToSources;
  if (mapping.size() > 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned m = 0; m < mapping.size(); ++m) {
      unsigned startSource = mapping[m].pointSourcesOffset;
      unsigned endSource = mapping[m].pointSourcesOffset + mapping[m].numberOfPointSources;
      if (sources_.mode == sourceterm::PointSources::NRF) {
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

void PointSourceClusterOnHost::addTimeIntegratedPointSourceNRF(unsigned source,
                                                               double from,
                                                               double to,
                                                               real dofs[tensor::Q::size()]) {
  real slip[] = {0.0, 0.0, 0.0};
  for (unsigned i = 0; i < 3; ++i) {
    if (sources_.slipRates[source][i].slopes.size() > 0) {
      slip[i] = sourceterm::computePwLFTimeIntegral(sources_.slipRates[source][i], from, to);
    }
  }

  real rotatedSlip[] = {0.0, 0.0, 0.0};
  for (unsigned i = 0; i < 3; ++i) {
    for (unsigned j = 0; j < 3; ++j) {
      rotatedSlip[j] += sources_.tensor[source][j + i * 3] * slip[i];
    }
  }

  kernel::sourceNRF krnl;
  krnl.Q = dofs;
  krnl.mInvJInvPhisAtSources = sources_.mInvJInvPhisAtSources[source].data();
  krnl.stiffnessTensor = sources_.stiffnessTensor[source].data();
  krnl.mSlip = rotatedSlip;
  krnl.mNormal = sources_.tensor[source].data() + 6;
  krnl.mArea = -sources_.A[source];
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
  kernel::sourceFSRM krnl;
  krnl.Q = dofs;
  krnl.mInvJInvPhisAtSources = sources_.mInvJInvPhisAtSources[source].data();
  krnl.momentFSRM = sources_.tensor[source].data();
  krnl.stfIntegral = sourceterm::computePwLFTimeIntegral(sources_.slipRates[source][0], from, to);
#ifdef MULTIPLE_SIMULATIONS
  krnl.oneSimToMultSim = init::oneSimToMultSim::Values;
#endif
  krnl.execute();
}

} // namespace seissol::kernels
