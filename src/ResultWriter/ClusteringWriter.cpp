// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ClusteringWriter.h"

#include "Common/Filesystem.h"
#include <Memory/Tree/Layer.h>
#include <cstddef>
#include <fstream>
#include <ios>
#include <string>
#include <utils/logger.h>

#include "Parallel/MPI.h"
namespace seissol::writer {

ClusteringWriter::ClusteringWriter(const std::string& outputPrefix) : outputPrefix(outputPrefix) {}

void ClusteringWriter::addCluster(unsigned profilingId,
                                  unsigned localClusterId,
                                  LayerType layerType,
                                  unsigned size,
                                  unsigned dynRupSize) {
  clusteringInformation.profilingIds.push_back(profilingId);
  clusteringInformation.localClusterIds.push_back(localClusterId);
  clusteringInformation.layerTypes.push_back(layerType);
  clusteringInformation.sizes.push_back(size);
  clusteringInformation.dynamicRuptureSizes.push_back(dynRupSize);
}

void ClusteringWriter::write() const {
  using namespace seissol::filesystem;
  const auto& mpi = MPI::mpi;

  const auto localRanks = mpi.collect(mpi.sharedMemMpiRank());
  const auto profilingIds = mpi.collectContainer(clusteringInformation.profilingIds);
  const auto localClusterIds = mpi.collectContainer(clusteringInformation.localClusterIds);
  const auto layerTypes = mpi.collectContainer(clusteringInformation.layerTypes);
  const auto sizes = mpi.collectContainer(clusteringInformation.sizes);
  const auto dynamicRuptureSizes = mpi.collectContainer(clusteringInformation.dynamicRuptureSizes);

  if (mpi.rank() == 0) {

    auto filepath = path(outputPrefix);
    filepath += path("-clustering.csv");

    auto fileStream = std::ofstream(filepath, std::ios::out);

    fileStream << "profilingId,localId,layerType,size,dynamicRuptureSize,rank,localRank\n";

    for (int rank = 0; rank < mpi.size(); ++rank) {
      const auto localRank = localRanks[rank];
      const auto& curProfilingIds = profilingIds[rank];
      const auto& curLocalClusterIds = localClusterIds[rank];
      const auto& curLayerTypes = layerTypes[rank];
      const auto& curSizes = sizes[rank];
      const auto& curDynamicRuptureSizes = dynamicRuptureSizes[rank];

      for (std::size_t i = 0; i < curProfilingIds.size(); ++i) {
        const auto layerType = static_cast<LayerType>(curLayerTypes[i]);
        if (layerType != LayerType::Interior && layerType != LayerType::Copy) {
          logError() << "Encountered illegal layer type in ClusteringWriter.";
        }
        const auto* layerTypeStr = layerType == Interior ? "Interior" : "Copy";
        fileStream << curProfilingIds[i] << "," << curLocalClusterIds[i] << "," << layerTypeStr
                   << "," << curSizes[i] << "," << curDynamicRuptureSizes[i] << "," << rank << ","
                   << localRank << "\n";
      }
    }

    fileStream.close();
  }
}

} // namespace seissol::writer
