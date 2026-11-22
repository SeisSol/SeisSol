// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLCLUSTERFACTORY_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLCLUSTERFACTORY_H_

#include "Core/Memory/Descriptor/LTS.h"
#include "Core/TimeStepping/Compute/CellCluster.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"

#include <memory>

namespace seissol::solver {

/**
 * Makes the clusters of the cells, with the numerics that compute on them (e.g. those of the wave
 * propagation): one per layer of the cells that are not ghosts.
 */
class CellClusterFactory {
  public:
  CellClusterFactory() = default;
  virtual ~CellClusterFactory() = default;

  CellClusterFactory(const CellClusterFactory&) = delete;
  CellClusterFactory(CellClusterFactory&&) = delete;
  auto operator=(const CellClusterFactory&) = delete;
  auto operator=(CellClusterFactory&&) = delete;

  /// the cluster of the cells of `layer`, with the time steps of its time cluster in
  /// `clusterLayout`; `printProgress` makes it report its steps
  virtual std::unique_ptr<CellCluster> makeCluster(LTS::Layer& layer,
                                                   const initializer::ClusterLayout& clusterLayout,
                                                   bool printProgress,
                                                   const ClusterMonitoring& monitoring) = 0;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLCLUSTERFACTORY_H_
