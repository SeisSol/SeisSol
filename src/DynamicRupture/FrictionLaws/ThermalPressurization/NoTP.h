#ifndef SEISSOL_NOTP_H
#define SEISSOL_NOTP_H

namespace seissol::dr::friction_law {
class NoTP {
  public:
  NoTP(DRParameters& drParameters){};

  /**
   * Initialize local attributes (used in initializer class respectively)
   */
  void initializeTP(seissol::Interoperability& eInteroperability){};

  void setInitialFluidPressure(unsigned ltsFace) {}

  /**
   * Copy all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {}

  /**
   * Set initial value of thermal pressure
   */
  void setInitialFluidPressureHook(std::array<real, misc::numPaddedPoints>& fluidPressure,
                                   unsigned int ltsFace) {}

  /**
   * Compute thermal pressure according to Noda and Lapusta 2010
   * bool saveTmpInTP is used to save final thermal pressure values for theta and sigma
   */
  void calcFluidPressure(FaultStresses const& faultStresses,
                         real (*initialStressInFaultCS)[misc::numPaddedPoints][6],
                         real (*mu)[misc::numPaddedPoints],
                         real (*slipRateMagnitude)[misc::numPaddedPoints],
                         real deltaT,
                         bool saveTmpInTP,
                         unsigned int timeIndex,
                         unsigned int ltsFace) {}

  /**
   * compute thermal pressure according to Noda and Lapusta 2010
   */
  void updateTemperatureAndPressure(unsigned int pointIndex,
                                    unsigned int timeIndex,
                                    unsigned int ltsFace) {}

  real heatSource(real tmp, real alpha, unsigned int tpGridPointIndex, unsigned int timeIndex) {
    return 0;
  }

  real fluidPressure(unsigned pointIndex) const { return 0; };
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_NOTP_H
