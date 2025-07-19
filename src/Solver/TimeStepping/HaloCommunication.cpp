// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HaloCommunication.h"
#include <Config.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/Typedefs.h>
#include <cstddef>

namespace seissol::solver {
HaloCommunication getHaloCommunication(const initializer::ClusterLayout& layout,
                                       const MeshStructure* structure) {
  HaloCommunication communication;
  communication.resize(layout.globalClusterCount);
  for (std::size_t i = 0; i < layout.globalClusterCount; ++i) {
    communication[i].resize(layout.globalClusterCount);
  }

  for (std::size_t i = 0; i < layout.globalClusterCount; ++i) {
    const auto& clusterStructure = structure[i];
    const std::size_t localClusterIndex = i;
    for (std::size_t j = 0; j < clusterStructure.numberOfRegions; ++j) {
      const std::size_t remoteClusterIndex = clusterStructure.neighboringClusters[j][1];

      communication.at(localClusterIndex)
          .at(remoteClusterIndex)
          .emplace_back(RemoteClusterPair{
              RemoteCluster{clusterStructure.copyRegions[j],
                            clusterStructure.copyRegionSizes[j],
                            Config::Precision,
                            clusterStructure.neighboringClusters[j][0],
                            DataTagOffset + clusterStructure.sendIdentifiers[j]},
              RemoteCluster{clusterStructure.ghostRegions[j],
                            clusterStructure.ghostRegionSizes[j],
                            Config::Precision,
                            clusterStructure.neighboringClusters[j][0],
                            DataTagOffset + clusterStructure.receiveIdentifiers[j]}});
    }
  }
  return communication;
}
} // namespace seissol::solver
