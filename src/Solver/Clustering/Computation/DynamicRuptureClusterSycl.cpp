#include "DynamicRuptureCluster.hpp"
#include <DynamicRupture/FrictionLaws/FrictionSolver.h>
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>
#include <Parallel/Helper.hpp>

namespace seissol::time_stepping {
void DynamicRuptureCluster::initFrictionSolverDevice() {
#ifdef ACL_DEVICE
  if (concurrentClusters()) {
    frictionSolverDeviceStorage = this->frictionSolverDevice->clone();
    this->frictionSolverDevice = frictionSolverDeviceStorage.get();

    if (auto* impl = dynamic_cast<dr::friction_law::gpu::FrictionSolverInterface*>(
            this->frictionSolverDevice)) {
      impl->initSyclQueue();

      impl->setMaxClusterSize(layer->getNumberOfCells());

      impl->allocateAuxiliaryMemory();

      // TODO(David): remove duplicate constant data
      impl->copyStaticDataToDevice();
    }
  }
#endif
}
} // namespace seissol::time_stepping
