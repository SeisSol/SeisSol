#ifndef SEISSOL_CLUSTERINGWRITER_H
#define SEISSOL_CLUSTERINGWRITER_H

#include <string>
#include <vector>
#include <type_traits>
#include "Initializer/tree/Layer.hpp"

namespace seissol::writer {

class ClusteringWriter {
public:
  ClusteringWriter(const std::string& outputPrefix);
    void addCluster(unsigned profilingId,
                    unsigned localClusterId,
                    LayerType layerType,
                    unsigned size,
                    unsigned dynRupSize);
    void write() const;

  // SoA that contains info about clusters
  struct ClusteringInformation {
    std::vector<int> ranks;
    std::vector<int> localRanks;
    std::vector<int> profilingIds;
    std::vector<int> localClusterIds;
    std::vector<std::underlying_type<LayerType>::type> layerTypes;
    std::vector<unsigned> sizes;
    std::vector<unsigned> dynamicRuptureSizes;
  };

private:
  std::string outputPrefix;
  ClusteringInformation clusteringInformation;
};

} // namespace seissol::writer

#endif // SEISSOL_CLUSTERINGWRITER_H
