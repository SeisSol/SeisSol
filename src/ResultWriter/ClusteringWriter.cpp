// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ClusteringWriter.h"

#include "Common/Filesystem.h"
#include "IO/Instance/Point/Csv.h"
#include "Initializer/BasicTypedefs.h"
#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"

#include <cstddef>
#include <string>
#include <type_traits>
#include <utils/logger.h>
#include <vector>
namespace seissol::writer {

ClusteringWriter::ClusteringWriter(const std::string& outputPrefix) : outputPrefix_(outputPrefix) {}

void ClusteringWriter::addCluster(unsigned profilingId,
                                  unsigned localClusterId,
                                  HaloType layerType,
                                  std::size_t size,
                                  std::size_t dynRupSize) {
  clusteringInformation_.profilingIds.push_back(profilingId);
  clusteringInformation_.localClusterIds.push_back(localClusterId);
  clusteringInformation_.layerTypes.push_back(
      static_cast<std::underlying_type_t<HaloType>>(layerType));
  clusteringInformation_.sizes.push_back(size);
  clusteringInformation_.dynamicRuptureSizes.push_back(dynRupSize);
}

void ClusteringWriter::write() const {
  using namespace seissol::filesystem;
  const auto& mpi = Mpi::mpi;

  const auto localRanks = mpi.collect(mpi.sharedMemMpiRank());
  const auto profilingIds = mpi.collectContainer(clusteringInformation_.profilingIds);
  const auto localClusterIds = mpi.collectContainer(clusteringInformation_.localClusterIds);
  const auto layerTypes = mpi.collectContainer(clusteringInformation_.layerTypes);
  const auto sizes = mpi.collectContainer(clusteringInformation_.sizes);
  const auto dynamicRuptureSizes = mpi.collectContainer(clusteringInformation_.dynamicRuptureSizes);

  if (mpi.rank() == 0) {
    logInfo() << "Cluster statistics:";
    for (std::size_t i = 0; i < clusteringInformation_.profilingIds.size(); ++i) {
      std::vector<double> sizestat(mpi.size());
      for (std::size_t j = 0; j < sizestat.size(); ++j) {
        sizestat[j] = sizes[j][i];
      }
      const auto sizeSummary = statistics::Summary(sizestat);
      const auto layerType = static_cast<HaloType>(clusteringInformation_.layerTypes[i]);
      const std::string layerTypeStr = layerType == HaloType::Interior ? "interior" : "copy";
      logInfo() << "cell" << layerTypeStr.c_str() << localClusterIds[0][i] << ":" << sizeSummary.sum
                << "(per rank:" << sizeSummary.mean << "±" << sizeSummary.std << "; range: ["
                << sizeSummary.min << ";" << sizeSummary.max << "])";
    }
    for (std::size_t i = 0; i < clusteringInformation_.profilingIds.size(); ++i) {
      std::vector<double> sizestat(mpi.size());
      for (std::size_t j = 0; j < sizestat.size(); ++j) {
        sizestat[j] = dynamicRuptureSizes[j][i];
      }
      const auto sizeSummary = statistics::Summary(sizestat);
      const auto layerType = static_cast<HaloType>(clusteringInformation_.layerTypes[i]);
      const std::string layerTypeStr = layerType == HaloType::Interior ? "interior" : "copy";
      logInfo() << "DR" << layerTypeStr.c_str() << localClusterIds[0][i] << ":" << sizeSummary.sum
                << "(per rank:" << sizeSummary.mean << "±" << sizeSummary.std << "; range: ["
                << sizeSummary.min << ";" << sizeSummary.max << "])";
    }

    seissol::io::instance::point::Csv table("clustering");
    table.addColumn<int>("profilingId");
    table.addColumn<int>("localId");
    // "Interior" or "Copy", so the widest of the two
    table.addTextColumn("layerType", 8);
    table.addColumn<std::size_t>("size");
    table.addColumn<std::size_t>("dynamicRuptureSize");
    table.addColumn<int>("rank");
    table.addColumn<int>("localRank");

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
        table.addCell<int>(curProfilingIds[i]);
        table.addCell<int>(curLocalClusterIds[i]);
        table.addText(layerType == HaloType::Interior ? "Interior" : "Copy");
        table.addCell<std::size_t>(curSizes[i]);
        table.addCell<std::size_t>(curDynamicRuptureSizes[i]);
        table.addCell<int>(rank);
        table.addCell<int>(localRank);
      }
    }

    auto filepath = path(outputPrefix_);
    filepath += path("-clustering.csv");
    table.writeFile(filepath.string());
  }
}

} // namespace seissol::writer
