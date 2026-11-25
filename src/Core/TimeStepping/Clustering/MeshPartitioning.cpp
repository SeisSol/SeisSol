// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MeshPartitioning.h"

#include "Core/Geometry/PUMLReader.h"
#include "Core/TimeStepping/Clustering/Clustering.h"
#include "Core/TimeStepping/Clustering/VertexWeights/VertexWeightModel.h"

#include <string>
#include <utils/logger.h>

namespace seissol::initializer {

const ClusteringResult& clusterAndPartition(geometry::PUMLReader& reader,
                                            Clustering& clustering,
                                            VertexWeightModel& weightModel,
                                            const std::string& partitioningLib,
                                            double tpwgt) {
  logInfo() << "Compute clustering.";
  const auto& result = clustering.compute(reader.topologyMesh(),
                                          reader.topologyArrays(),
                                          reader.geometryMesh(),
                                          reader.geometryArrays(),
                                          reader.vertexOrders());
  weightModel.build(result);
  reader.partition(
      geometry::PUMLReader::Partitioning{.vertexWeights = weightModel.vertexWeights(),
                                         .weightsPerVertex = weightModel.nWeightsPerVertex(),
                                         .imbalance = weightModel.imbalances()[0] - 1.0,
                                         .clusterIds = result.clusterIds.data(),
                                         .timesteps = result.timesteps.cellTimeStepWidths.data()},
      partitioningLib,
      tpwgt);
  return result;
}

} // namespace seissol::initializer
