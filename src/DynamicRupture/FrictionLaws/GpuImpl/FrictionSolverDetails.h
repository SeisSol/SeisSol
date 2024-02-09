#ifndef SEISSOL_FRICTION_SOLVER_DETAILS_H
#define SEISSOL_FRICTION_SOLVER_DETAILS_H

#include <yaml-cpp/yaml.h>
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"

namespace seissol::dr::friction_law::gpu {
class FrictionSolverDetails : public FrictionSolverInterface {
  public:
  explicit FrictionSolverDetails(seissol::initializer::parameters::DRParameters* drParameters);
  ~FrictionSolverDetails() override;

  void initSyclQueue() override;
  void allocateAuxiliaryMemory() override;
  void copyStaticDataToDevice() override;

  virtual void
      copySpecificLtsDataTreeToLocal(seissol::initializer::Layer& layerData,
                                     seissol::initializer::DynamicRupture const* const dynRup,
                                     real fullUpdateTime) = 0;

  protected:
  size_t currLayerSize{};

  FaultStresses* faultStresses{nullptr};
  TractionResults* tractionResults{nullptr};
  real (*stateVariableBuffer)[misc::numPaddedPoints]{nullptr};
  real (*strengthBuffer)[misc::numPaddedPoints]{nullptr};
  real* resampleMatrix{nullptr};
  double* devTimeWeights{nullptr};
  real* devSpaceWeights{nullptr};
  int* queue{nullptr};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_FRICTION_SOLVER_DETAILS_H
