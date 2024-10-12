// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DynamicRuptureCluster.h"
#include <DynamicRupture/FrictionLaws/FrictionSolver.h>
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>
#include <Parallel/Helper.h>

namespace seissol::solver::clustering::computation {
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
} // namespace seissol::solver::clustering::computation
