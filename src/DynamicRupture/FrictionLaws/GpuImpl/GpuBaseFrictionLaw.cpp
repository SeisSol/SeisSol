#include "DynamicRupture/FrictionLaws/GpuImpl/GpuBaseFrictionLaw.h"
#include "utils/logger.h"
#include <omp.h>
#include <algorithm>
#include <sstream>

#pragma omp declare target(                                                                        \
    seissol::dr::friction_law::FrictionSolver::precomputeStressFromQInterpolated)
#pragma omp declare target(                                                                        \
    seissol::dr::friction_law::FrictionSolver::postcomputeImposedStateFromNewStress)

namespace seissol::dr::friction_law::gpu {
GpuBaseFrictionLaw::GpuBaseFrictionLaw(dr::DRParameters& drParameters)
    : FrictionSolver(drParameters){};

GpuBaseFrictionLaw::~GpuBaseFrictionLaw() {
  //#pragma omp target exit data map(release:faultStresses[0:maxClusterSize])
  delete[] faultStresses;
  delete[] tractionResults;
}

void GpuBaseFrictionLaw::evaluate(seissol::initializers::Layer& layerData,
                                  seissol::initializers::DynamicRupture* dynRup,
                                  real fullUpdateTime,
                                  double timeWeights[CONVERGENCE_ORDER]) {
  FrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  this->copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);

#pragma omp target data map(to : this)
  {
    auto layerSize = layerData.getNumberOfCells();

// clang-format off
    #pragma omp target teams loop map(from: faultStresses [0:layerSize]) \
    is_device_ptr(qInterpolatedPlus, qInterpolatedMinus, impAndEta) device(diviceId)
    // clang-format on
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
      precomputeStressFromQInterpolated(faultStresses[ltsFace], ltsFace);
    }

// loop over all dynamic rupture faces, in this LTS layer
#pragma omp parallel for schedule(static)
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {

      // define some temporary variables
      std::array<real, misc::numPaddedPoints> stateVariableBuffer{0};
      std::array<real, misc::numPaddedPoints> strengthBuffer{0};

      this->preHook(stateVariableBuffer, ltsFace);

      // loop over sub time steps (i.e. quadrature points in time)
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
        this->updateFrictionAndSlip(faultStresses[ltsFace],
                                    tractionResults[ltsFace],
                                    stateVariableBuffer,
                                    strengthBuffer,
                                    ltsFace,
                                    timeIndex);
      }

      this->postHook(stateVariableBuffer, ltsFace);

      // output rupture front
      this->saveRuptureFrontOutput(ltsFace);

      // output time when shear stress is equal to the dynamic stress after rupture arrived
      this->saveDynamicStressOutput(ltsFace);

      // output peak slip rate
      this->savePeakSlipRateOutput(ltsFace);

      // output average slip
      // TODO: What about outputSlip
      // this->saveAverageSlipOutput(outputSlip, ltsFace);
    }

    // clang-format off
    #pragma omp target teams loop map(to: faultStresses [0:layerSize], tractionResults[0:layerSize], timeWeights[0:CONVERGENCE_ORDER]) \
    is_device_ptr(imposedStatePlus, imposedStateMinus, qInterpolatedPlus, qInterpolatedMinus, impAndEta) \
    device(diviceId)
    // clang-format on
    for (unsigned ltsFace = 0; ltsFace < layerSize; ++ltsFace) {
      this->postcomputeImposedStateFromNewStress(
          faultStresses[ltsFace], tractionResults[ltsFace], timeWeights, ltsFace);
    }
  }
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
  //#pragma omp target enter data map(alloc:faultStresses[0:maxClusterSize])
}
} // namespace seissol::dr::friction_law::gpu
