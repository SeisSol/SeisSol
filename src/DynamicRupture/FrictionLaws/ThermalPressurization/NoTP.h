#ifndef SEISSOL_NOTP_H
#define SEISSOL_NOTP_H

namespace seissol::dr::friction_law {
class NoTP {
  public:
  NoTP(DRParameters& drParameters){};

  void setInitialFluidPressure(unsigned ltsFace) {}

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {}

  void calcFluidPressure(FaultStresses const& faultStresses,
                         real (*initialStressInFaultCS)[misc::numPaddedPoints][6],
                         real (*mu)[misc::numPaddedPoints],
                         std::array<real, misc::numPaddedPoints>&,
                         real deltaT,
                         bool saveTmpInTP,
                         unsigned int timeIndex,
                         unsigned int ltsFace) {}

  real fluidPressure(unsigned, unsigned) const { return 0; };
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_NOTP_H
