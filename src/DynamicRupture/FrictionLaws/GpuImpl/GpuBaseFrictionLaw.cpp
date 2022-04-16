#include "DynamicRupture/FrictionLaws/GpuImpl/GpuBaseFrictionLaw.h"
#include "utils/logger.h"
#include <omp.h>
#include <sstream>


namespace seissol::dr::friction_law::gpu {
void GpuBaseFrictionLaw::evaluate(seissol::initializers::Layer& layerData,
                                  seissol::initializers::DynamicRupture* dynRup,
                                  real fullUpdateTime,
                                  double timeWeights[CONVERGENCE_ORDER]) {
  FrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  this->copySpecificLtsDataTreeToLocal(layerData, dynRup, fullUpdateTime);

  // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
    FaultStresses faultStresses = this->precomputeStressFromQInterpolated(ltsFace);

    // define some temporary variables
    std::array<real, misc::numPaddedPoints> stateVariableBuffer{0};
    std::array<real, misc::numPaddedPoints> strengthBuffer{0};

    this->preHook(stateVariableBuffer, ltsFace);

    TractionResults tractionResults = {};

    // loop over sub time steps (i.e. quadrature points in time)
    for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
      this->updateFrictionAndSlip(
          faultStresses, tractionResults, stateVariableBuffer, strengthBuffer, ltsFace, timeIndex);
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

    // compute output
    this->postcomputeImposedStateFromNewStress(
        faultStresses, tractionResults, timeWeights, ltsFace);
  }
}

void GpuBaseFrictionLaw::checkOffloading() {
  bool canOffload = false;
  #pragma omp target map(tofrom: canOffload)
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
