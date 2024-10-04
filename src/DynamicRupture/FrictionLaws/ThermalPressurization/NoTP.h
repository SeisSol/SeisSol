#ifndef SEISSOL_NOTP_H
#define SEISSOL_NOTP_H

namespace seissol::dr::friction_law {
class NoTP {
  public:
  NoTP(seissol::initializer::parameters::DRParameters* drParameters) {};

  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

  void calcFluidPressure(std::array<real, misc::NumPaddedPoints>& normalStress,
                         real (*mu)[misc::NumPaddedPoints],
                         std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
                         real deltaT,
                         bool saveTmpInTP,
                         unsigned int timeIndex,
                         unsigned int ltsFace) {}

  [[nodiscard]] real getFluidPressure(unsigned /*unused*/, unsigned /*unused*/) { return 0; };
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_NOTP_H
