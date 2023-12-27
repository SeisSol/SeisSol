#ifndef SEISSOL_NOTP_H
#define SEISSOL_NOTP_H

namespace seissol::dr::friction_law {
class NoTP {
  public:
  NoTP(seissol::initializers::parameters::DRParameters* drParameters){};

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture const* const dynRup,
                          real fullUpdateTime) {}

  void calcFluidPressure(std::array<real, misc::numPaddedPoints>& normalStress,
                         real (*mu)[misc::numPaddedPoints],
                         std::array<real, misc::numPaddedPoints>& slipRateMagnitude,
                         real deltaT,
                         bool saveTmpInTP,
                         unsigned int timeIndex,
                         unsigned int ltsFace) {}

  real getFluidPressure(unsigned, unsigned) const { return 0; };
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_NOTP_H
