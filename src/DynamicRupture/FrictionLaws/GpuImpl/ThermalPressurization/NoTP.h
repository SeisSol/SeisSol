#ifndef SEISSOL_GPU_NOTP_H
#define SEISSOL_GPU_NOTP_H

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
class NoTP {
  public:
  NoTP(seissol::initializer::parameters::DRParameters* drParameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

  void calcFluidPressure(real (*normalStress)[misc::NumPaddedPoints],
                         real (*mu)[misc::NumPaddedPoints],
                         real (*slipRateMagnitude)[misc::NumPaddedPoints],
                         real deltaT,
                         bool saveTmpInTP) {}

  struct Details {};
  Details getCurrentLayerDetails() { return Details{}; }
  static real getFluidPressure(Details, unsigned, unsigned) { return static_cast<real>(0.0); };
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_NOTP_H
