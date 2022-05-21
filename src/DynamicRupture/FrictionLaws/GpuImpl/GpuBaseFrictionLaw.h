#ifndef SEISSOL_GPUBASEFRICTIONLAW_H
#define SEISSOL_GPUBASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
class GpuBaseFrictionLaw : public FrictionSolver {
  public:
  GpuBaseFrictionLaw(dr::DRParameters& drParameters);
  ~GpuBaseFrictionLaw();

  void allocateAuxiliaryMemory(seissol::initializers::LTSTree* drTree,
                               seissol::initializers::DynamicRupture* drDescr,
                               int diviceId);

  virtual void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                              seissol::initializers::DynamicRupture* dynRup,
                                              real fullUpdateTime) = 0;

  protected:
  void checkOffloading();

  size_t maxClusterSize{};
  size_t currLayerSize{};
  int diviceId{};

  FaultStresses* faultStresses{nullptr};
  TractionResults* tractionResults{nullptr};
  real (*stateVariableBuffer)[misc::numPaddedPoints]{nullptr};
  real (*strengthBuffer)[misc::numPaddedPoints]{nullptr};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPUBASEFRICTIONLAW_H
