// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "PointSourceClusterOnDevice.h"

#include <generated_code/tensor.h>
#include <generated_code/init.h>
#include <Parallel/AcceleratorDevice.h>
#include <SourceTerm/PointSource.h>

#include <utility>
#include <sycl/sycl.hpp>

namespace seissol::kernels {

PointSourceClusterOnDevice::PointSourceClusterOnDevice(sourceterm::ClusterMapping mapping,
                                                       sourceterm::PointSources sources)
    : clusterMapping_(std::move(mapping)), sources_(std::move(sources)) {}

void PointSourceClusterOnDevice::addTimeIntegratedPointSources(double from, double to) {
  auto& queue = seissol::AcceleratorDevice::getInstance().getSyclDefaultQueue();
  auto& mapping = clusterMapping_.cellToSources;
  if (mapping.size() > 0) {
    auto* mapping_ptr = mapping.data();
    auto* slipRates0 = sources_.slipRates[0].data();
    auto* slipRates1 = sources_.slipRates[1].data();
    auto* slipRates2 = sources_.slipRates[2].data();
    auto* mInvJInvPhisAtSources = sources_.mInvJInvPhisAtSources.data();
    auto* tensor = sources_.tensor.data();
    auto* A = sources_.A.data();
    auto* stiffnessTensor = sources_.stiffnessTensor.data();
    if (sources_.mode == sourceterm::PointSources::NRF) {
      queue
          .parallel_for({mapping.size()},
                        [=](sycl::item<1> id) {
                          unsigned startSource = mapping_ptr[id[0]].pointSourcesOffset;
                          unsigned endSource = mapping_ptr[id[0]].pointSourcesOffset +
                                               mapping_ptr[id[0]].numberOfPointSources;
                          for (unsigned source = startSource; source < endSource; ++source) {
                            addTimeIntegratedPointSourceNRF(
                                {&slipRates0[source], &slipRates1[source], &slipRates2[source]},
                                mInvJInvPhisAtSources[source].data(),
                                tensor[source].data(),
                                A[source],
                                stiffnessTensor[source].data(),
                                from,
                                to,
                                *mapping_ptr[id[0]].dofs);
                          }
                        })
          .wait();
    } else {
      queue
          .parallel_for({mapping.size()},
                        [=](sycl::item<1> id) {
                          unsigned startSource = mapping_ptr[id[0]].pointSourcesOffset;
                          unsigned endSource = mapping_ptr[id[0]].pointSourcesOffset +
                                               mapping_ptr[id[0]].numberOfPointSources;
                          for (unsigned source = startSource; source < endSource; ++source) {
                            addTimeIntegratedPointSourceFSRM(&slipRates0[source],
                                                             mInvJInvPhisAtSources[source].data(),
                                                             tensor[source].data(),
                                                             from,
                                                             to,
                                                             *mapping_ptr[id[0]].dofs);
                          }
                        })
          .wait();
    }
  }
}

void PointSourceClusterOnDevice::addTimeIntegratedPointSourceNRF(
    std::array<sourceterm::PiecewiseLinearFunction1D const*, 3> slipRates,
    real* mInvJInvPhisAtSources,
    real* tensor,
    real A,
    real* stiffnessTensor,
    double from,
    double to,
    real dofs[tensor::Q::size()]) {
  real slip[3] = {real(0.0)};
  for (unsigned i = 0; i < 3; ++i) {
    if (slipRates[i]->slopes.size() > 0) {
      slip[i] = sourceterm::computePwLFTimeIntegral(*slipRates[i], from, to);
    }
  }

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
        m += -A * stiffnessTensor[p + 3 * q + 9 * i + 27 * j] * rotatedSlip[i] * tensor[6 + j];
      }
    }
    return m;
  };

  real moment[6] = {mom(0, 0), mom(1, 1), mom(2, 2), mom(0, 1), mom(1, 2), mom(0, 2)};
  for (unsigned t = 0; t < 6; ++t) {
    for (unsigned k = 0; k < tensor::mInvJInvPhisAtSources::Shape[0]; ++k) {
      dofs[k + t * (init::Q::Stop[0] - init::Q::Start[0])] += mInvJInvPhisAtSources[k] * moment[t];
    }
  }
}

void PointSourceClusterOnDevice::addTimeIntegratedPointSourceFSRM(
    sourceterm::PiecewiseLinearFunction1D const* slipRate0,
    real* mInvJInvPhisAtSources,
    real* tensor,
    double from,
    double to,
    real* dofs) {
  auto stfIntegral = sourceterm::computePwLFTimeIntegral(*slipRate0, from, to);
  for (unsigned p = 0; p < tensor::momentFSRM::Shape[0]; ++p) {
    for (unsigned k = 0; k < tensor::mInvJInvPhisAtSources::Shape[0]; ++k) {
      dofs[k + p * (init::Q::Stop[0] - init::Q::Start[0])] +=
          stfIntegral * mInvJInvPhisAtSources[k] * tensor[p];
    }
  }
}

} // namespace seissol::kernels
