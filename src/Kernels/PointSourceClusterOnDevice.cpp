// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PointSourceClusterOnDevice.h"

#include "generated_code/init.h"
#include "generated_code/tensor.h"

// needs to be loaded after Eigen at the moment, due to SYCL
#include "Parallel/AcceleratorDevice.h"

#include <Kernels/PointSourceCluster.h>
#include <Kernels/Precision.h>
#include <Parallel/Runtime/Stream.h>
#include <SourceTerm/Typedefs.h>
#include <array>
#include <cstddef>
#include <memory>

namespace seissol::kernels {

PointSourceClusterOnDevice::PointSourceClusterOnDevice(
    std::shared_ptr<sourceterm::ClusterMapping> mapping,
    std::shared_ptr<sourceterm::PointSources> sources)
    : clusterMapping_(mapping), sources_(sources) {}

unsigned PointSourceClusterOnDevice::size() const { return sources_->numberOfSources; }

void PointSourceClusterOnDevice::addTimeIntegratedPointSources(
    double from, double to, seissol::parallel::runtime::StreamRuntime& runtime) {
  pointSourceKernel(*clusterMapping_, *sources_, from, to, runtime);
}

} // namespace seissol::kernels
