#ifndef SEISSOL_NOTP_H
#define SEISSOL_NOTP_H

namespace seissol::dr::friction_law {
template <typename Config>
class NoTP {
  public:
  using RealT = typename Config::RealT;
  NoTP(DRParameters* drParameters){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture<Config> const* const dynRup,
                          RealT fullUpdateTime) {}

  void calcFluidPressure(std::array<RealT, misc::numPaddedPoints<Config>>& normalStress,
                         RealT (*mu)[misc::numPaddedPoints<Config>],
                         std::array<RealT, misc::numPaddedPoints<Config>>& slipRateMagnitude,
                         RealT deltaT,
                         bool saveTmpInTP,
                         unsigned int timeIndex,
                         unsigned int ltsFace) {}

  RealT getFluidPressure(unsigned, unsigned) const { return 0; };
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_NOTP_H
