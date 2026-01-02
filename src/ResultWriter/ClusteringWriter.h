// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_CLUSTERINGWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_CLUSTERINGWRITER_H_

#include "Memory/Tree/Layer.h"

#include <string>
#include <type_traits>
#include <vector>

namespace seissol::writer {

class ClusteringWriter {
  public:
  explicit ClusteringWriter(const std::string& outputPrefix);
  void addCluster(unsigned profilingId,
                  unsigned localClusterId,
                  HaloType layerType,
                  std::size_t size,
                  std::size_t dynRupSize);
  void write() const;

  // SoA that contains info about clusters
  struct ClusteringInformation {
    std::vector<int> ranks;
    std::vector<int> localRanks;
    std::vector<int> profilingIds;
    std::vector<int> localClusterIds;
    std::vector<std::underlying_type_t<HaloType>> layerTypes;
    std::vector<std::size_t> sizes;
    std::vector<std::size_t> dynamicRuptureSizes;
  };

  private:
  std::string outputPrefix;
  ClusteringInformation clusteringInformation;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_CLUSTERINGWRITER_H_
