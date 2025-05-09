// SPDX-FileCopyrightText: 2023 SeisSol Group
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
                                  ClusterType clusterType,
                                  LayerType layerType,
                                  unsigned size) {
  clusteringInformation.profilingIds.push_back(profilingId);
  clusteringInformation.localClusterIds.push_back(localClusterId);
  clusteringInformation.clusterTypes.push_back(static_cast<int>(clusterType));
  clusteringInformation.layerTypes.push_back(layerType);
  clusteringInformation.sizes.push_back(size);
}

void ClusteringWriter::write() const {
  using namespace seissol::filesystem;
  const auto& mpi = MPI::mpi;

  const auto localRanks = mpi.collect(mpi.sharedMemMpiRank());
  const auto profilingIds = mpi.collectContainer(clusteringInformation.profilingIds);
  const auto localClusterIds = mpi.collectContainer(clusteringInformation.localClusterIds);
  const auto clusterTypes = mpi.collectContainer(clusteringInformation.clusterTypes);
  const auto layerTypes = mpi.collectContainer(clusteringInformation.layerTypes);
  const auto sizes = mpi.collectContainer(clusteringInformation.sizes);

  if (mpi.rank() == 0) {

    auto filepath = path(outputPrefix);
    filepath += path("-clustering.csv");

    auto fileStream = std::ofstream(filepath, std::ios::out);

    fileStream << "profilingId,localId,clusterType,layerType,size,rank,localRank\n";

    for (int rank = 0; rank < mpi.size(); ++rank) {
      const auto localRank = localRanks[rank];
      const auto& curProfilingIds = profilingIds[rank];
      const auto& curLocalClusterIds = localClusterIds[rank];
      const auto& curClusterTypes = clusterTypes[rank];
      const auto& curLayerTypes = layerTypes[rank];
      const auto& curSizes = sizes[rank];

      for (std::size_t i = 0; i < curProfilingIds.size(); ++i) {
        const auto layerType = static_cast<LayerType>(curLayerTypes[i]);
        const auto clusterType = static_cast<ClusterType>(curClusterTypes[i]);
        if (layerType != LayerType::Interior && layerType != LayerType::Copy) {
          logError() << "Encountered illegal layer type in ClusteringWriter.";
        }
        const auto* layerTypeStr = layerType == Interior ? "Interior" : "Copy";
        const auto* clusterTypeStr = clusterType == ClusterType::Cell ? "Cell" : "Face";
        fileStream << curProfilingIds[i] << "," << curLocalClusterIds[i] << "," << layerTypeStr
                   << "," << clusterTypeStr << "," << curSizes[i] << "," << rank << "," << localRank
                   << "\n";
      }
    }

    fileStream.close();
  }
}

} // namespace seissol::writer
