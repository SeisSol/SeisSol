#ifndef SEISSOL_NEIGHBORCLUSTER_H
#define SEISSOL_NEIGHBORCLUSTER_H

#include <memory>
#include "PostOffice.h"

namespace seissol {
namespace time_stepping {
struct NeighborCluster {
  std::shared_ptr<PostOffice> messagesFromNeighbor;
  std::shared_ptr<PostOffice> messagesToNeighbor;
  double predictedTime;
  double maxTimeStepSize;
};

}
}
#endif //SEISSOL_NEIGHBORCLUSTER_H
