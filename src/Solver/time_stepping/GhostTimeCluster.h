#ifndef SEISSOL_GHOSTTIMECLUSTER_H
#define SEISSOL_GHOSTTIMECLUSTER_H

#include "AbstractTimeCluster.h"

namespace seissol {
namespace time_stepping {

class GhostTimeCluster : public AbstractTimeCluster {
public:
  bool act();
};

}
}



#endif //SEISSOL_GHOSTTIMECLUSTER_H
