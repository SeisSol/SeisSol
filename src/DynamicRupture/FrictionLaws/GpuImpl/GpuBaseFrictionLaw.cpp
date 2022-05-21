#include "DynamicRupture/FrictionLaws/GpuImpl/GpuBaseFrictionLaw.h"
#include "utils/logger.h"
#include <omp.h>
#include <algorithm>
#include <sstream>

namespace seissol::dr::friction_law::gpu {
GpuBaseFrictionLaw::GpuBaseFrictionLaw(dr::DRParameters& drParameters)
    : FrictionSolver(drParameters) {
  checkOffloading();
}

GpuBaseFrictionLaw::~GpuBaseFrictionLaw() {
  //#pragma omp target exit data map(release:faultStresses[0:maxClusterSize])
  delete[] faultStresses;
  delete[] tractionResults;
  delete[] stateVariableBuffer;
  delete[] strengthBuffer;
}

void GpuBaseFrictionLaw::allocateAuxiliaryMemory(seissol::initializers::LTSTree* drTree,
                                                 seissol::initializers::DynamicRupture* drDescr,
                                                 int currDiviceId) {
  for (auto it = drTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != drTree->endLeaf();
       ++it) {
    size_t currClusterSize = static_cast<size_t>(it->getNumberOfCells());
    maxClusterSize = std::max(currClusterSize, maxClusterSize);
  }

  diviceId = currDiviceId;
  faultStresses = new FaultStresses[maxClusterSize];
  tractionResults = new TractionResults[maxClusterSize];

  using ArrayType = std::remove_pointer<decltype(stateVariableBuffer)>::type;
  stateVariableBuffer = new ArrayType[maxClusterSize];
  strengthBuffer = new ArrayType[maxClusterSize];
  //#pragma omp target enter data map(alloc:faultStresses[0:maxClusterSize])
}

void GpuBaseFrictionLaw::checkOffloading() {
  bool canOffload = false;
  #pragma omp target map(tofrom : canOffload)
  {
    if (!omp_is_initial_device()) {
      canOffload = true;
    }
  }
  std::ostringstream info;
  info << "Device offloading: " << std::boolalpha << canOffload;
  logInfo() << info.str();
}
} // namespace seissol::dr::friction_law::gpu
