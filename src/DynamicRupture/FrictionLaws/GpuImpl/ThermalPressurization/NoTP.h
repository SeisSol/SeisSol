#ifndef SEISSOL_GPU_NOTP_H
#define SEISSOL_GPU_NOTP_H

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
class NoTP {
  public:
  NoTP(DRParameters* drParameters){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {}

  #pragma omp declare target
  void calcFluidPressure(real (*normalStress)[misc::numPaddedPoints],
                         real (*mu)[misc::numPaddedPoints],
                         real (*slipRateMagnitude)[misc::numPaddedPoints],
                         real deltaT,
                         bool saveTmpInTP) {}

  struct Details {};
  Details getCurrentLayerDetails() { return Details{}; }
  static real getFluidPressure(Details, unsigned, unsigned) { return static_cast<real>(0.0); };
  #pragma omp end declare target
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_NOTP_H
