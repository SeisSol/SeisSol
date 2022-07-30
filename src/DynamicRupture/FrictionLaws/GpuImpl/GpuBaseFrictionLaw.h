#ifndef SEISSOL_GPUBASEFRICTIONLAW_H
#define SEISSOL_GPUBASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
class GpuBaseFrictionLaw : public FrictionSolver {
  public:
  GpuBaseFrictionLaw(dr::DRParameters* drParameters);
  ~GpuBaseFrictionLaw();

  void setDeviceId(int currDeviceId);
  void setMaxClusterSize(size_t size) { maxClusterSize = size; }
  void allocateAuxiliaryMemory();
  void copyStaticDataToDevice();

  virtual void copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                              seissol::initializers::DynamicRupture* dynRup,
                                              real fullUpdateTime) = 0;

  protected:
  void checkOffloading();

  size_t maxClusterSize{};
  size_t currLayerSize{};
  int hostId{};
  int deviceId{};

  FaultStresses* faultStresses{nullptr};
  TractionResults* tractionResults{nullptr};
  real (*stateVariableBuffer)[misc::numPaddedPoints]{nullptr};
  real (*strengthBuffer)[misc::numPaddedPoints]{nullptr};
  real* resampleMatrix{nullptr};
  double* devTimeWeights{nullptr};
  real* devDeltaT{nullptr};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPUBASEFRICTIONLAW_H
