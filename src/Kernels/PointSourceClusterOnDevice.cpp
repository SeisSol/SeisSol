// SPDX-FileCopyrightText: 2023 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PointSourceClusterOnDevice.h"

#include <Kernels/PointSourceCluster.h>
#include <Kernels/Precision.h>
#include <Parallel/Runtime/Stream.h>
#include <SourceTerm/Typedefs.h>
#include <memory>
#include <utility>

namespace seissol::kernels {

template <typename Cfg>
PointSourceClusterOnDevice<Cfg>::PointSourceClusterOnDevice(
    std::shared_ptr<sourceterm::ClusterMapping> mapping,
    std::shared_ptr<sourceterm::PointSources<Cfg>> sources)
    : clusterMapping_(std::move(mapping)), sources_(std::move(sources)) {}

template <typename Cfg>
std::size_t PointSourceClusterOnDevice<Cfg>::size() const {
  return sources_->numberOfSources;
}

template <typename Cfg>
void PointSourceClusterOnDevice<Cfg>::addTimeIntegratedPointSources(
    double from, double to, seissol::parallel::runtime::StreamRuntime& runtime) {
  pointSourceKernel(*clusterMapping_, *sources_, from, to, runtime);
}

#define _H_(cfg) template class PointSourceClusterOnDevice<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::kernels
