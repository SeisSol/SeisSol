#ifndef SEISSOL_GPU_NOTP_H
#define SEISSOL_GPU_NOTP_H

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
template <typename Config>
class NoTP {
  public:
  using RealT = typename Config::RealT;
  NoTP(DRParameters* drParameters){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime) {}

  void calcFluidPressure(RealT (*normalStress)[misc::numPaddedPoints<Config>],
                         RealT (*mu)[misc::numPaddedPoints<Config>],
                         RealT (*slipRateMagnitude)[misc::numPaddedPoints<Config>],
                         RealT deltaT,
                         bool saveTmpInTP) {}

  struct Details {};
  Details getCurrentLayerDetails() { return Details{}; }
  static RealT getFluidPressure(Details, unsigned, unsigned) { return static_cast<RealT>(0.0); };
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_GPU_NOTP_H
