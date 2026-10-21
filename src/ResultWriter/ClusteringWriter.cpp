// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ClusteringWriter.h"

#include "Common/Filesystem.h"
#include "Common/RunId.h"
#include "IO/Instance/Point/Csv.h"
#include "Initializer/BasicTypedefs.h"
#include "Numerical/Statistics.h"
#include "Parallel/MPI.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utils/logger.h>
#include <vector>
namespace seissol::writer {

ClusteringWriter::ClusteringWriter(const std::string& outputPrefix, const std::string& runId)
    : outputPrefix_(outputPrefix), runId_(runId) {}

void ClusteringWriter::addCluster(unsigned profilingId,
                                  unsigned localClusterId,
                                  HaloType layerType,
                                  std::size_t size,
                                  std::size_t dynRupSize,
                                  double timestep,
                                  long updateFactor) {
  clusteringInformation_.profilingIds.push_back(profilingId);
  clusteringInformation_.localClusterIds.push_back(localClusterId);
  clusteringInformation_.layerTypes.push_back(
      static_cast<std::underlying_type_t<HaloType>>(layerType));
  clusteringInformation_.sizes.push_back(size);
  clusteringInformation_.dynamicRuptureSizes.push_back(dynRupSize);
  clusteringInformation_.timesteps.push_back(timestep);
  clusteringInformation_.updateFactors.push_back(updateFactor);
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
  const auto timesteps = mpi.collectContainer(clusteringInformation_.timesteps);
  const auto updateFactors = mpi.collectContainer(clusteringInformation_.updateFactors);

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
    table.addColumn<double>("timestep");
    table.addColumn<std::int64_t>("updateFactor");
    table.addTextColumn("runId", RunIdLength);

    for (int rank = 0; rank < mpi.size(); ++rank) {
      const auto localRank = localRanks[rank];
      const auto& curProfilingIds = profilingIds[rank];
      const auto& curLocalClusterIds = localClusterIds[rank];
      const auto& curLayerTypes = layerTypes[rank];
      const auto& curSizes = sizes[rank];
      const auto& curDynamicRuptureSizes = dynamicRuptureSizes[rank];
      const auto& curTimesteps = timesteps[rank];
      const auto& curUpdateFactors = updateFactors[rank];

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
        table.addCell<double>(curTimesteps[i]);
        table.addCell<std::int64_t>(curUpdateFactors[i]);
        table.addText(runId_);
      }
    }

    auto filepath = path(outputPrefix_);
    filepath += path("-clustering.csv");
    table.writeFile(filepath.string());
  }
}

} // namespace seissol::writer
