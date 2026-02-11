// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ClusteringWriter.h"

#include "Common/Filesystem.h"
#include "Initializer/BasicTypedefs.h"
#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"

#include <cstddef>
#include <fstream>
#include <ios>
#include <string>
#include <type_traits>
#include <utils/logger.h>
#include <vector>
namespace seissol::writer {

ClusteringWriter::ClusteringWriter(const std::string& outputPrefix) : outputPrefix(outputPrefix) {}

void ClusteringWriter::addCluster(unsigned profilingId,
                                  unsigned localClusterId,
                                  HaloType layerType,
                                  std::size_t size,
                                  std::size_t dynRupSize) {
  clusteringInformation.profilingIds.push_back(profilingId);
  clusteringInformation.localClusterIds.push_back(localClusterId);
  clusteringInformation.layerTypes.push_back(
      static_cast<std::underlying_type_t<HaloType>>(layerType));
  clusteringInformation.sizes.push_back(size);
  clusteringInformation.dynamicRuptureSizes.push_back(dynRupSize);
}

void ClusteringWriter::write() const {
  using namespace seissol::filesystem;
  const auto& mpi = Mpi::mpi;

  const auto localRanks = mpi.collect(mpi.sharedMemMpiRank());
  const auto profilingIds = mpi.collectContainer(clusteringInformation.profilingIds);
  const auto localClusterIds = mpi.collectContainer(clusteringInformation.localClusterIds);
  const auto layerTypes = mpi.collectContainer(clusteringInformation.layerTypes);
  const auto sizes = mpi.collectContainer(clusteringInformation.sizes);
  const auto dynamicRuptureSizes = mpi.collectContainer(clusteringInformation.dynamicRuptureSizes);

  if (mpi.rank() == 0) {
    logInfo() << "Cluster statistics:";
    for (std::size_t i = 0; i < clusteringInformation.profilingIds.size(); ++i) {
      std::vector<double> sizestat(mpi.size());
      for (std::size_t j = 0; j < sizestat.size(); ++j) {
        sizestat[j] = sizes[j][i];
      }
      const auto sizeSummary = statistics::Summary(sizestat);
      const auto layerType = static_cast<HaloType>(clusteringInformation.layerTypes[i]);
      const std::string layerTypeStr = layerType == HaloType::Interior ? "interior" : "copy";
      logInfo() << "cell" << layerTypeStr.c_str() << localClusterIds[0][i] << ":" << sizeSummary.sum
                << "(per rank:" << sizeSummary.mean << "±" << sizeSummary.std << "; range: ["
                << sizeSummary.min << ";" << sizeSummary.max << "])";
    }
    for (std::size_t i = 0; i < clusteringInformation.profilingIds.size(); ++i) {
      std::vector<double> sizestat(mpi.size());
      for (std::size_t j = 0; j < sizestat.size(); ++j) {
        sizestat[j] = dynamicRuptureSizes[j][i];
      }
      const auto sizeSummary = statistics::Summary(sizestat);
      const auto layerType = static_cast<HaloType>(clusteringInformation.layerTypes[i]);
      const std::string layerTypeStr = layerType == HaloType::Interior ? "interior" : "copy";
      logInfo() << "DR" << layerTypeStr.c_str() << localClusterIds[0][i] << ":" << sizeSummary.sum
                << "(per rank:" << sizeSummary.mean << "±" << sizeSummary.std << "; range: ["
                << sizeSummary.min << ";" << sizeSummary.max << "])";
    }

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
        const auto layerType = static_cast<HaloType>(curLayerTypes[i]);
        if (layerType != HaloType::Interior && layerType != HaloType::Copy) {
          logError() << "Encountered illegal layer type in ClusteringWriter.";
        }
        const auto* layerTypeStr = layerType == HaloType::Interior ? "Interior" : "Copy";
        fileStream << curProfilingIds[i] << "," << curLocalClusterIds[i] << "," << layerTypeStr
                   << "," << curSizes[i] << "," << curDynamicRuptureSizes[i] << "," << rank << ","
                   << localRank << "\n";
      }
    }

    fileStream.close();
  }
}

} // namespace seissol::writer
