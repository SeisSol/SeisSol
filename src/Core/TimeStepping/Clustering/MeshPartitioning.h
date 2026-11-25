// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_CLUSTERING_MESHPARTITIONING_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_CLUSTERING_MESHPARTITIONING_H_

#include "Core/Geometry/PUMLReader.h"
#include "Core/TimeStepping/Clustering/Clustering.h"
#include "Core/TimeStepping/Clustering/VertexWeights/VertexWeightModel.h"

#include <string>

namespace seissol::initializer {

/**
 * Clusters the cells of the mesh `reader` has read (see Clustering::compute()), weighs them with
 * `weightModel`, and partitions the mesh by those weights (see geometry::PUMLReader::partition()),
 * with the node weight `tpwgt` of this rank and the library `partitioningLib`; their clusters and
 * time steps go along with the cells. Collective. Returns what the clustering decided.
 */
const ClusteringResult& clusterAndPartition(geometry::PUMLReader& reader,
                                            Clustering& clustering,
                                            VertexWeightModel& weightModel,
                                            const std::string& partitioningLib,
                                            double tpwgt = 1.0);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_CLUSTERING_MESHPARTITIONING_H_
