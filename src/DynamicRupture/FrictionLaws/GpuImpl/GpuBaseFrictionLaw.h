#ifndef SEISSOL_GPUBASEFRICTIONLAW_H
#define SEISSOL_GPUBASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
class GpuBaseFrictionLaw : public FrictionSolver {
  public:
  GpuBaseFrictionLaw(dr::DRParameters& drParameters) : FrictionSolver(drParameters) {
    checkOffloading();
  };

  void checkOffloading();
  void allocateAuxiliaryMemory(seissol::initializers::LTSTree* drTree,
                               seissol::initializers::DynamicRupture* drDescr);

  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture* dynRup,
                real fullUpdateTime,
                double timeWeights[CONVERGENCE_ORDER]) override;

  virtual void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                              seissol::initializers::DynamicRupture* dynRup,
                                              real fullUpdateTime) = 0;

  virtual void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                       unsigned ltsFace) = 0;
  virtual void updateFrictionAndSlip(FaultStresses& faultStresses,
                                     TractionResults& tractionResults,
                                     std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                                     std::array<real, misc::numPaddedPoints>& strengthBuffer,
                                     unsigned int ltsFace,
                                     unsigned int timeIndex) = 0;

  virtual void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                        unsigned int ltsFace) = 0;

  virtual void calcStrengthHook(FaultStresses& faultStresses,
                                std::array<real, misc::numPaddedPoints>& strength,
                                unsigned int timeIndex,
                                unsigned int ltsFace) = 0;

  virtual void calcStateVariableHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                                     unsigned int timeIndex,
                                     unsigned int ltsFace) = 0;

  virtual void saveDynamicStressOutput(unsigned int ltsFace) = 0;

  protected:
  size_t maxClusterSize{};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPUBASEFRICTIONLAW_H
