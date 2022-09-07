#ifndef SEISSOL_GPUBASEFRICTIONLAW_H
#define SEISSOL_GPUBASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include <CL/sycl.hpp>

namespace sycl = cl::sycl;

namespace seissol::dr::friction_law::gpu {
class GpuBaseFrictionLaw : public FrictionSolver {
  public:
  GpuBaseFrictionLaw(dr::DRParameters* drParameters);
  ~GpuBaseFrictionLaw();

  void initSyclQueue();
  void setMaxClusterSize(size_t size) { maxClusterSize = size; }
  void allocateAuxiliaryMemory();
  void copyStaticDataToDevice();

  virtual void
      copySpecificLtsDataTreeToLocal(seissol::initializers::Layer& layerData,
                                     seissol::initializers::DynamicRupture const* const dynRup,
                                     real fullUpdateTime) = 0;

  protected:
  size_t maxClusterSize{};
  size_t currLayerSize{};

  FaultStresses* faultStresses{nullptr};
  TractionResults* tractionResults{nullptr};
  real (*stateVariableBuffer)[misc::numPaddedPoints]{nullptr};
  real (*strengthBuffer)[misc::numPaddedPoints]{nullptr};
  real* resampleMatrix{nullptr};
  double* devTimeWeights{nullptr};
  real* devDeltaT{nullptr};

  sycl::device device;
  sycl::queue queue;
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPUBASEFRICTIONLAW_H
