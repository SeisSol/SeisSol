#ifndef SEISSOL_FRICTION_SOLVER_DETAILS_H
#define SEISSOL_FRICTION_SOLVER_DETAILS_H

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"
#include <yaml-cpp/yaml.h>

// #define CURRCHUNK std::min(this->currLayerSize,this->chunksize*chunk):std::min(this->currLayerSize,this->chunksize*(chunk+1))

// #define CCHUNK(var) var[std::min(this->currLayerSize,this->chunksize*chunk):std::min(this->currLayerSize,this->chunksize*(chunk+1))-std::min(this->currLayerSize,this->chunksize*chunk)]
// #define CCHUNK(var) (var + std::min(this->currLayerSize,this->chunksize*chunk))[0:std::min(this->currLayerSize,this->chunksize*(chunk+1))]

#define CCHUNK(var) var[std::min(layerSize,chunksize*chunk):std::min(layerSize,chunksize*(chunk+1))-std::min(layerSize,chunksize*chunk)]
#define CCHUNKLOOP(var) for (int var = std::min(layerSize,chunksize*chunk); var < std::min(layerSize,chunksize*(chunk+1)); ++var)

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
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) = 0;

  protected:
  size_t currLayerSize{};

  FaultStresses* faultStresses{nullptr};
  TractionResults* tractionResults{nullptr};
  real (*stateVariableBuffer)[misc::NumPaddedPoints]{nullptr};
  real (*strengthBuffer)[misc::NumPaddedPoints]{nullptr};
  real* resampleMatrix{nullptr};
  double* devTimeWeights{nullptr};
  real* devSpaceWeights{nullptr};
  int* queue{nullptr};
  std::size_t chunksize{0};
  std::size_t chunkcount{1};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_FRICTION_SOLVER_DETAILS_H
