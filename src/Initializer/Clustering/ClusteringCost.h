// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERINGCOST_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERINGCOST_H_

#include <cstdint>
#include <mpi.h>
#include <vector>

namespace seissol::initializer {

double computeLocalCostOfClustering(const std::vector<int>& clusterIds,
                                    const std::vector<int>& cellCosts,
                                    const std::vector<uint64_t>& rate,
                                    double wiggleFactor,
                                    double minimalTimestep);

double computeGlobalCostOfClustering(const std::vector<int>& clusterIds,
                                     const std::vector<int>& cellCosts,
                                     const std::vector<uint64_t>& rate,
                                     double wiggleFactor,
                                     double minimalTimestep,
                                     MPI_Comm comm);

/// Local cost of every clustering obtained by capping `clusterIds` at 0, 1, ... `maxClusterId`.
///
/// Entry `m` equals `computeLocalCostOfClustering(enforceMaxClusterId(clusterIds, m), ...)`,
/// bit for bit: each entry accumulates over the cells in the same order. Computing them
/// together replaces one full clustering copy and one reduction per candidate with a single
/// pass and a single reduction.
std::vector<double> computeLocalCostsOfCappedClusterings(const std::vector<int>& clusterIds,
                                                         const std::vector<int>& cellCosts,
                                                         const std::vector<uint64_t>& rate,
                                                         double wiggleFactor,
                                                         double minimalTimestep,
                                                         int maxClusterId);

std::vector<int> enforceMaxClusterId(const std::vector<int>& clusterIds, int maxClusterId);

int computeMaxClusterIdAfterAutoMerge(const std::vector<int>& clusterIds,
                                      const std::vector<int>& cellCosts,
                                      const std::vector<uint64_t>& rate,
                                      double maximalAdmissibleCost,
                                      double wiggleFactor,
                                      double minimalTimestep);

std::uint64_t ratepow(const std::vector<std::uint64_t>& rate, std::uint64_t a, std::uint64_t b);

/// Bins a cell timestep into a cluster id of the ladder described by `rate`.
///
/// Cluster k covers `[ratepow(rate,0,k), ratepow(rate,0,k+1))` in units of
/// `wiggleFactor * globalMinTimestep`; a ratio of 1 at index >= 1 terminates the ladder, and
/// a rate vector shorter than the ladder is extended with its last entry.
std::uint64_t getCluster(double timestep,
                         double globalMinTimestep,
                         double wiggleFactor,
                         const std::vector<std::uint64_t>& rate);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERINGCOST_H_
